
//###################################################################################
//                     Meteor Trajectory Iterative Solution
//###################################################################################
//
// This meteor trajectory module solves for the atmospheric track of a meteor given
// video measurements from multiple cameras. It uses a multi-parameter fit approach
// that assumes that for each camera, the timing information for measurements are
// precisely known. That is, the measurement angles are as extracted from the focal
// plane meteor positions and can be noisy, but the timing is provided by the user,
// whom most often assumes a constant video rate. So a given measurement's time, 
// relative to others for a given camera, is well known. But if the precise time
// of each measurement is obtainable, then that can be used instead and it does 
// not need to be uniformly spaced.
//
// The multi-parameter fit algorithm and implemented software can handle two or more 
// cameras and also deal with multiple cameras from the same site. It uses a 
// particle swarm optimzation (PSO) algorithm to fit the best propagating linear-path 
// motion-model to the measurements provided. The parameters solved for are the 
// radiant direction, begin point velocity, deceleration terms, timing offsets  
// between the unsynchronized camera measurements, begin and end point LLA positions, 
// and LLA position, range, and velocity at each measurement and motion model point.
//
// It solves the problem in several steps via a bootstrapping approach to help ensure
// the non-linear iterative solver finds the global minimum of the cost function. The
// PSO is used to find the minimum of the residual angles between the measurement rays 
// and points along the 3D linear-path motion-model.
//
// For video systems, the time of the measurements are typically assumed to be  
// uniformly spaced, but can be non-uniformly spaced and/or have drop-outs. Also each 
// camera can be unsynchronized with respect to the others, and the solver will
// estimate the sub-frame timing offsets between them. In a post-solution step, 
// a Monte Carlo loop provides an estimate of the standard deviations for the 
// principle parameters by adding 2D Gaussian noise (normally distributed) to the  
// measurements and re-solving the trajectory. Note the final trajectory returned is
// for the no added noise case and is NOT the mean of the Monte Carlo solutions.
//
// IMPORTANT NOTES:
//
// Results are for geocentric radiants and velocities and have not been corrected 
// for zenith attraction. The trajectory estimator also does not correct for the
// relative motion between the camera sites during passage of the meteor, nor for
// stellar aberration.
//
//###################################################################################
//
// Software Interface: The trajectory solver relies on a number of functions that are  
//                     self contained within this file. In addition, the processing
//                     approach relies heavily on the particle swarm optimization
//                     module (via the include file ParticleSwarmFunctions.h). The 
//                     following lines of pseudo-code gives a calling sequence
//                     example. For each trajectory to be solved, the user should
//                     RESET the trajectory structure, then feed measurements to the
//                     trajectory_info structure per camera via INFILL. After all the
//                     measurements have been populated into the trajectory_info
//                     structure, then the function MeteorTrajectory should be 
//                     invoked to perform a simultaneous solve for all parameters
//                     via a low fidelity to high fidelity bootstrapping technique.
//                     --------------------------------------------------------------
//
//  #include  TrajectorySolution.h
//
//  #include  ParticleSwarmFunctions.h
//
//
//  struct trajectory_info  traj;
//
//
//  InitTrajectoryStructure( maxcameras_expected, &traj );  //... Called at program startup
//
//  ...
//
//  //-------- For every solution to be computed, always RESET the trajectory structure providing:
//                  reference julian date/time, max timing offset, velocity model mode, 
//                  #monte carlo trials, measurement type mode, verbose flag
//
//	ResetTrajectoryStructure( jdt_ref, max_toffset, velmodel, nummonte, meastype, verbose, &traj );
//
//  //--------- Next loop over each camera track (set of measurements) to INFILL the trajectory structure
//  //              where the loop count will set the camera count in the structure. User to provide the
//  //              number of measurements for each camera, the two anglar measurement vectors (radians),
//  //              time relative to the jdt reference time (seconds), and noise std dev estimate per
//  //              measurement (radians). Also the camera's site latitude, longitude, and height 
//  //              in standard GPS coordinates (radians and kilometers).
//
//  Loop over cameras (each multi-frame measurement sequence) to infill the trajectory structure
//
//       InfillTrajectoryStructure( #measurements, ra, dec, deltatime, noise, latitude, longitude, height, &traj );
//
//  End camera loop
//
//
//  MeteorTrajectory( &traj );   //... Solve for best trajectory fit parameters
//
//
//  //--------- Solution complete, start next trajectory solution with a RESET call ...
//  ...
//  ...
//  ...
//  
//
//  FreeTrajectoryStructure( &traj );  //... Called at program completion to free memory
//
//
//###################################################################################
//
// Bootstrapping Algorithm:
//
//          Step 1 uses the intersecting planes method between all pairs of camera
//                 measurements, obtaining the position and radiant for the pair
//                 with the largest convergence angle.
//          Step 2 is based upon Jiri Borivicka's LMS approach from the paper "The
//                 Comparison of Two Methods of Determining Meteor Trajectories 
//                 from Photographs" where the position and radiant direction are 
//                 refined from step 1 relying solely on minimizing the sum of
//                 angles (vice distances in the paper) between all camera 
//                 measurements and a 3D line, ignoring any temporal information.
//          Step 3 uses a least mean square solution to the velocity by projecting
//                 all measurement rays to the radiant line and fitting the CPA
//                 points nearest the line for ALL camera measurement sets
//                 simultaneously assuming a common velocity. 
//          Step 4 uses the PSO on a subset of parameters to iteratively solve
//                 for the timing offsets between the unsynchronized cameras 
//                 keeping all other parameters fixed.
//          Step 5 uses the PSO to refine both the velocity and timing offsets 
//                 assuming a constant velocity model.
//          Step 6 uses the PSO to refine the velocity and make estimates of the 
//                 deceleration parameters keeping all other parameters fixed.
//          Step 7 uses the PSO to iteratively solve for ALL the variable parameters 
//                 simultaneously given that a good estimate is now available 
//                 for each of the unknowns. Produces the final "no additive 
//                 noise" solution.
//          Step 8 loops over Monte Carlo trials using PSO to quickly solve for
//                 a solution having the initial PSO guess start around the no  
//                 noise solution. Each Monte Carlo trial adds limited extent 
//                 2D Gaussian noise to the measurements, and obtain a set of 
//                 new solutions which feeds a statistical standard deviation 
//                 for a subset of the critical parameters (Radiant Rah/Dec, 
//                 Vbegin, Decel1, Decel2).
//
//###################################################################################
//
// Date         Ver   Author        Comment
// -----------  ---   ------------  -------------------------------------------------
// Feb 06 2011  1.0   Pete Gural    Initial C implementation with ameoba-simplex
// Oct 15 2016  2.0   Pete Gural    Revised with PSO, modularized, new I/O interface, 
//                                  weighted LMS velocity, weighted cost functions 
//
//###################################################################################



//###################################################################################
//                        Includes and Constants
//###################################################################################

#include "ParticleSwarmFunctions.h"

#define   X   0
#define   Y   1
#define   Z   2
#define   T   3

#define   RADEC         1   // "meastype" is RA and Dec
#define   NAZEL         2   // Azim +east of north, altitude angle
#define   SAZZA         3   // Azim +west of south, zenith angle

#define   POSITION      0   // "propagation_state" 
#define   VELOCITY      1

#define   CONSTANT      0   // "velmodel" velocity motion model
#define   LINEAR        1
#define   QUADRATIC     2
#define   EXPONENT      3

#define   MEAS2LINE     0   // cost function selection
#define   MEAS2MODEL    1

#define   ACCURATE_FIT  0   // "pso_accuracy_config"  setting
#define   QUICK_FIT     1

#define   LLA_BEG       0   // "LLA_position"  extreme begin or end point
#define   LLA_END       1


//###################################################################################
//                      Trajectory Structure Definition
//###################################################################################

struct trajectory_info   
{
	//======= PARAMETERS REFLECTING INPUT SETTINGS and CONTROL ==========================================

	//----------------------------- Max memory handling, intermediate reporting
	int       maxcameras;        // Maximum number of cameras expected (to initially allocate arrays)
	int       verbose;           // Flag to turn on intermediate product reporting during function call
                                 //     0 = no writes to the console
                                 //     1 = all step's intermediate values displayed on console
                                 //     2 = only final step solution parameters displayed on console
                                 //     3 = TBD measurements and model positions sent to console

	//----------------------------- Modeling parameters
	int       nummonte;          // Number of Monte Carlo trials for standard deviation calculation
	int       velmodel;          // Velocity propagation model
                                 //     0 = constant   v(t) = vinf
                                 //     1 = linear     v(t) = vinf - |acc1| * t
                                 //     2 = quadratic  v(t) = vinf - |acc1| * t + acc2 * t^2
                                 //     3 = exponent   v(t) = vinf - |acc1| * |acc2| * exp( |acc2| * t )
	int       meastype;          // Flag to indicate the type of angle measurements the user is providing
	                             //    for meas1 and meas2 below. The following are all in radians:
	                             //        1 = Right Ascension for meas1, Declination for meas2.
	                             //        2 = Azimuth +east of due north for meas1, Elevation angle 
	                             //            above the horizon for meas2
	                             //        3 = Azimuth +west of south for meas1, Zenith angle for meas2

	//----------------------------- Reference timing and offset constraint 
	double    jdt_ref;           // Reference Julian date/time that the measurements times are provided
                                 //     relative to. This is user selectable and can be the time of the 
	                             //     first camera, or the first measurement, or some average time for 
	                             //     the meteor, but should be close to the time of passage of 
	                             //     the meteor. This same reference date/time will be used on all
                                 //     camera measurements for the purposes of computing local sidereal
	                             //     time and making  geocentric coordinate transformations.	
	double    max_toffset;       // Maximum allowed time offset between cameras in seconds


	//======= SOLUTION STARTUP and DERIVED PRODUCTS ==========================================================

	//----------------------------- Memory allocation handling 
	int      *malloced;          // Vectors and arrays memory allocation flag per camera (0=freed,1=allocated)

	//----------------------------- Camera/site information 
	int       numcameras;        // Number of cameras having measurements passed in and is determined by
	                             //    the number of calls to InfillTrajectoryStructure. Measurements are
	                             //    typically from multiple sites but can include multiple cameras from
	                             //    the same site. There must be one pair of well separated site cameras
	                             //    for triangulation to proceed successfully.
    double   *camera_lat;        // Vector of GEODETIC latitudes of the cameras in radians (GPS latitude)
	double   *camera_lon;        // Vector of longitudes (positive east) of the cameras in radians
	double   *camera_hkm;        // Vector of camera heights in km above a WGS84 ellipsoid (GPS height)
	double   *camera_LST;        // Vector of local sidereal times for each camera (based on jdt_ref)
	double  **rcamera_ECI;       // Camera site ECI radius vector from Earth center (#cameras x XYZ)

	//----------------------------- Measurement information 
	int      *nummeas;           // Vector containing the number of measurements per camera

	                             // ------The following are dimensioned #cameras x #measurements(for that camera)
	double  **meas1;             // Array of 1st measurement type (see meastype), typically RA, typically in radians
	double  **meas2;             // Array of 2nd measurement type (see meastype), typically Dec, typically in radians
	double  **dtime;             // Array of time measurements in seconds relative to jdt_ref
	double  **noise;             // Array of measurement standard deviations in radians (per measurement)

	double ***meashat_ECI;       // Measurement ray unit vectors (#cameras x #measurements x XYZT)


	//----------------------------- Trajectory fitting parameters to feed to the particle swarm optimizer
	int       numparams;         // Number of total parameters for a full fit
	double   *params;            // Working vector of parameter values of length numparams
	double   *limits;            // Each parameter's value search limit of length numparams
	double   *xguess;            // Initial starting guess for the parameters of length numparams
	double   *xshift;            // Plus/minus limit around the guess values of length numparams (if set
	                             //    to zero, this parameter will remain fixed during the PSO)


	//======= SOLUTION OUTPUT PRODUCTS ====================================================================

    //------------------ Note: The output RA/Dec equatorial coordinates epoch is the same as that  
    //                         of the input measurements and represent Geocentric values, but they   
    //                         are NOT corrected for zenith attraction.
    //                         The velocity at the begin point "vbegin" could be considered equivalent
    //                         to the velocity at the top of the atmosphere Vinfinity, and is also NOT
    //                         corrected for zenith attraction.

	//----------------------------- Best solution vector of the parameter values of length numparams
	double   *solution;          // { Rx, Ry, Rz, Vx, Vy, Vz, Decel1, Decel2, tref_offsets[*] }
	                             // Note that R and V are in ECI (ECEF)

	//----------------------------- Primary output products and their standard deviations (sigma) 

    double    ra_radiant;        // Radiant right ascension in radians (multi-parameter fit)
    double    dec_radiant;       // Radiant declination in radians (multi-parameter fit)
    double    vbegin;            // Meteor solution velocity at the begin point in km/sec
    double    decel1;            // Deceleration term 1 defined by the given velocity model
    double    decel2;            // Deceleration term 2 defined by the given velocity model

    double    ra_sigma;          // Standard deviation of radiant right ascension in radians
    double    dec_sigma;         // Standard deviation of radiant declination in radians
    double    vbegin_sigma;      // Standard deviation of vbegin in km/sec
    double    decel1_sigma;      // Standard deviation of decceleration term 1
    double    decel2_sigma;      // Standard deviation of decceleration term 2

	//----------------------------- Intermediate bootstrapping solutions

	                             // Intersecting planes solution for best convergence angle pair
	double    max_convergence;   // Max convergence angle between all camera pairs in radians
    double    ra_radiant_IP;     // Radiant right ascension in radians
    double    dec_radiant_IP;    // Radiant declination in radians

	                             // Intersecting planes solution for weighted multiple tracks
    double    ra_radiant_IPW;    // Radiant right ascension in radians
    double    dec_radiant_IPW;   // Radiant declination in radians

	                             // Borovicka's least mean squares solution for the radiant
    double    ra_radiant_LMS;    // Radiant right ascension in radians
    double    dec_radiant_LMS;   // Radiant declination in radians

	//----------------------------- Relative timing output
    double    begtime;           // Begin position time relative to jdt_ref in seconds
    double    endtime;           // End position time relative to jdt_ref in seconds
    double   *tref_offsets;      // Vector of timing offsets in seconds re to jdt_ref per camera

	//----------------------------- Measurement and model LLA, range and velocity arrays 
	//                                  with dimension #cameras x #measurements(camera)
    double  **meas_lat;          // Array of geodetic latitudes closest to trail for each measurement
    double  **meas_lon;          // Array of +east longitudes closest to trail for each measurement
    double  **meas_hkm;          // Array of heights re WGS84 closest to trail for each measurement
    double  **meas_range;        // Array of ranges from site along measurement to the CPA of the trail
    double  **meas_vel;          // Array of velocity along the trail for each measurement

    double  **model_lat;         // Array of geodetic latitudes for the model positions
    double  **model_lon;         // Array of +east longitudes for the model positions
    double  **model_hkm;         // Array of heights re WGS84 for the model positions
    double  **model_range;       // Array of ranges from site to the model positions
    double  **model_vel;         // Array of velocity on the trail at each model position

	//----------------------------- Model fit vectors which follow the same "meastype" on output
	//                                  with dimension #cameras x #measurements(camera)
    double  **model_fit1;        // Array of 1st data sequence containing the model fit in meastype format
    double  **model_fit2;        // Array of 2nd data sequence containing the model fit in meastype format
    double  **model_time;        // Array of model time which includes offsets relative to the reference time

	//----------------------------- BEGIN position and standard deviation in LLA
    double    rbeg_lat;          // Position on radiant line as GEODETIC latitude in radians
    double    rbeg_lon;          // Position on radiant line as +EAST longitude in radians
    double    rbeg_hkm;          // Position on radiant line as height in km relative WGS84

    double    rbeg_lat_sigma;    // Standard deviation of rbeg_lat in radians
    double    rbeg_lon_sigma;    // Standard deviation of rbeg_lon in radians
    double    rbed_hkm_sigma;    // Standard deviation of rbeg_hkm in km

    //----------------------------- END position and standard deviation in LLA
    double    rend_lat;          // Position on radiant line as GEODETIC latitude in radians
    double    rend_lon;          // Position on radiant line as +EAST longitude in radians
    double    rend_hkm;          // Position on radiant line as height in km relative WGS84

    double    rend_lat_sigma;    // Standard deviation of rend_lat in radians
    double    rend_lon_sigma;    // Standard deviation of rend_lon in radians
    double    rend_hkm_sigma;    // Standard deviation of rend_hkm in km


};  //... end trajectory_info structure definition


//###################################################################################
//                          Prototype Definitions
//###################################################################################

int     MeteorTrajectory( struct trajectory_info *traj );

void    InitTrajectoryStructure( int maxcameras, struct trajectory_info *traj );

void    FreeTrajectoryStructure( struct trajectory_info *traj );

void    ResetTrajectoryStructure( double jdt_ref, double max_toffset, 
	                              int velmodel, int nummonte, int meastype, int verbose, 
                                  struct trajectory_info *traj );

void    InfillTrajectoryStructure( int nummeas, double *meas1, double *meas2, double *dtime, double *noise,
	                               double site_latitude, double site_longitude, double site_height,
                                   struct trajectory_info *traj );

void    AllocateTrajectoryMemory4Infill( int kcamera, int nummeas, struct trajectory_info *traj );

void    IntermediateConsoleOutput( char *ctitle, double *params, int ncams );

void    EnsureNonZeroStdDev( struct trajectory_info *traj );

void    Angles2SiteMeasurements( double   camera_lat,  double  camera_lon,  double  camera_hkm,  double camera_LST,
		                         int      nummeas,     int     meastype,
								 double  *rah,         double *dec,         double *tsec,
								 int      noiseflag,   double *noise,
								 double **meashat_ECI, double *rcamera_ECI  );

void    IntersectingPlanes_MultiTrack( struct trajectory_info *traj, double *radiant_hat );

void    IntersectingPlanes_BestConvergence( struct trajectory_info *traj, 
	                                        int    *kcamera_bestconv, double *max_convergence, 
										    double *radiant_bestconv, double *r_bestconv  );

void    Normal2MeteorMeasurements( int nmeas, double **meashat_ECI, double *planenormal );

void    TwoLineCPA( double *position1, double *vector1, double *r1, 
				    double *position2, double *vector2, double *r2, double *d21 );

double  Propagation( int propmodel, int velmodel, double t, double vinf, double acc1, double acc2 );

double  VelocityFit_Differencing( struct trajectory_info *traj );

double  VelocityFit_LMS( struct trajectory_info *traj );

void    ConstantVelocity_MultiTrackFit( int ncameras, int *nmeas_per_camera, double *pos, double *tim, double *noi,
	                                     double *velocityLMS, double *xo, double *rmserror, double *err );

void    ParameterRefinementViaPSO( struct trajectory_info *traj, int cost_function, int velocity_model, int pso_accuracy_config );

double  Anglesum_Measurements2Line( double *ro_ECI, double *radiant_ECI, struct trajectory_info *traj );

double  Anglesum_Measurements2Model( int velmodel,  double *ro_ECI, double *radiant_ECI, 
				                     double decel1,  double decel2, double *dt,
				                     struct trajectory_info *traj );

void    ReportFill_LLAVT_Meas_Model( struct trajectory_info *traj );

void    ReportFill_LLA_Beg_End( int LLA_position, struct trajectory_info *traj );

void    MonteCarlo_ErrorEstimate( struct trajectory_info *traj );

double  VectorDotProduct( double *a, double *b );

void    VectorCrossProduct( double *a, double *b, double *c );

void    VectorNorm( double *a, double *b, double *c );

double  RandomGauss( double sigma, double maxvalue );

void    ECI2RahDec( double *eci, double *rah, double *dec );

void    RahDec2ECI( double rah, double dec, double *eci );

void    AzimuthElevation2RADec( double  azim, 
                                double  elev,
                                double  geodetic_latitude,
                                double  LST,
				                double *RA,
					            double *DEC );

void    RADec2AzimuthElevation( double  RA,
                                double  DEC,
                                double  geodetic_latitude, 
                                double  LST,
					            double *Azim,
					            double *Elev  );

double  LocalSiderealTimeE( double jdt, double longitude_east );

double  LST2LongitudeEast( double jdt, double LST );

void    LatLonAlt2ECEF( double lat, double lon, double alt_km, double *ecef_km  );

void    ECEF2LatLonAlt( double *ecef_km, double *lat, double *lon, double *alt_km );


//############################################################################################
//                                   Functions
//############################################################################################

int   MeteorTrajectory( struct trajectory_info *traj )
{
int     kcamera, kcamera_bestconv, k, kbeg, noise_flag;
double  radiant_hat[3], radiant_bestconv[3];
double  r_bestconv[3], rdummy[3], rbeg[3], dist, magdummy;
double  max_convergence, vbegin;


	//======== Zero all the parameter's initial values and set their search limits

	traj->numparams = 8 + traj->numcameras;

	for( k=0; k<traj->numparams; k++ )  traj->params[k] = 0.0;  //... All parameters initialized to zero

	for( k=0; k<3;               k++ )  traj->limits[k] = 5.0;                //... ECI begin position search limit in km
	for( k=3; k<6;               k++ )  traj->limits[k] = 5.0;                //... ECI begin velocity search limit in km/sec
	for( k=6; k<8;               k++ )  traj->limits[k] = 1.0;                //... Deceleration coefficients search limit
	for( k=8; k<traj->numparams; k++ )  traj->limits[k] = traj->max_toffset;  //... Timing offset search limit in seconds


	//======== Since we now use an inverse of the noise variance per measurement to weight
	//         the minimzation cost function, look for any zero valued standard deviations
	//         input by the user and set them to the minimum sigma found.

	EnsureNonZeroStdDev( traj );


    //======== Site positions and measured coordinates are converted into measurement unit vectors
	//         that have a common coordinate system (Earth Centered Inertial - ECI) where the
	//         conversion depends on the user's choice for input measurement type.
	//         Zero measurement noise is added to the measurements at this stage of processing.

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {
		
		noise_flag = 0;; //... No noise added

	    Angles2SiteMeasurements( traj->camera_lat[kcamera], 
								 traj->camera_lon[kcamera], 
								 traj->camera_hkm[kcamera],
								 traj->camera_LST[kcamera],
								 traj->nummeas[kcamera], 
								 traj->meastype,
								 traj->meas1[kcamera], 
								 traj->meas2[kcamera], 
								 traj->dtime[kcamera],
								 noise_flag,
								 traj->noise[kcamera],
								 traj->meashat_ECI[kcamera], 
								 traj->rcamera_ECI[kcamera]  );

	} //... end of camera loop


    //======== Get a weighted multi-track intersecting planes solution (currently not used)

	IntersectingPlanes_MultiTrack( traj, radiant_hat );

	ECI2RahDec( radiant_hat, &traj->ra_radiant_IPW, &traj->dec_radiant_IPW );

	if( traj->verbose == 1 )  printf(" Multi-Track Weighted IP:  RA = %lf   Dec = %lf\n", traj->ra_radiant_IPW * 57.296, traj->dec_radiant_IPW * 57.296 );


    //======== Get an initial guess for the radiant by pairing all track combinations 
	//         and selecting the intersecting planes solution with the largest 
	//         convergence angle.

    IntersectingPlanes_BestConvergence( traj, &kcamera_bestconv, &max_convergence, radiant_bestconv, r_bestconv  );

	ECI2RahDec( radiant_bestconv, &traj->ra_radiant_IP, &traj->dec_radiant_IP );

	if( traj->verbose == 1 )  printf(" Best Convergence Ang IP:  RA = %lf   Dec = %lf\n", traj->ra_radiant_IP * 57.296, traj->dec_radiant_IP * 57.296 );


	traj->max_convergence = max_convergence;

	if( traj->verbose == 1 )  printf(" Best camera used for Ro is %i with convergence angle %lf\n\n", kcamera_bestconv, traj->max_convergence * 57.296 );


	//======== Use the first measurement of the best convergence angle camera,
	//         to find a starting position in 3D space (ECI coordinates). The 
	//         position vector is based on the measurement CPA to the estimated 
	//         radiant line and lies along the measurement unit vector.

	kbeg = 0;  //... 1st measurement index

    TwoLineCPA( traj->rcamera_ECI[kcamera_bestconv], traj->meashat_ECI[kcamera_bestconv][kbeg], rbeg,
		        r_bestconv, radiant_bestconv, rdummy, &dist );

	traj->params[0] = rbeg[X];   //... Initial 3D position state vector in ECI spatial coordinates
	traj->params[1] = rbeg[Y];
	traj->params[2] = rbeg[Z];


	//======== Use the best radiant from the intersecting planes solution as an initial guess for the 
	//         unity velocity state vector (radiant direction only) in ECI coordinates

	traj->params[3] = radiant_bestconv[X];  //... Initial 3D velocity for the IP's best convergence angle
	traj->params[4] = radiant_bestconv[Y];
	traj->params[5] = radiant_bestconv[Z];

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "Intersecting Planes Solution", traj->params, traj->numcameras );
	

	//======== First pass PSO to refine position and velocity (radiant drection) via Jiri Borivicka's LMS method

	for( k=0; k<traj->numparams; k++ )  traj->xguess[k] = traj->params[k];
	for( k=0; k<traj->numparams; k++ )  traj->xshift[k] = 0.0;

    for( k=0; k<3; k++ )  traj->xshift[k] = traj->limits[k];  //... Vary starting position
    for( k=3; k<6; k++ )  traj->xshift[k] = 0.02;             //... Vary the velocity (radiant unit vector at this stage)

	ParameterRefinementViaPSO( traj, MEAS2LINE, CONSTANT, ACCURATE_FIT );

	VectorNorm( &traj->params[3], &traj->params[3], &magdummy );  //... Ensure radiant direction is a unit vector (set |V|=1 for now)

	ECI2RahDec( &traj->params[3], &traj->ra_radiant_LMS, &traj->dec_radiant_LMS );

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "Position and Radiant LMS Solution", traj->params, traj->numcameras );


	//-------- Use the measurements of all the camera tracks to compute the initial
	//         velocity estimate working from absolute positions (CPA) of the rays
	//         to the current 3D radiant line. Least mean squares estimate is made
	//         of the velocity with outlier removal.

	////vbegin = VelocityFit_Differencing( traj ); 

	vbegin = VelocityFit_LMS( traj ); 
	
	traj->params[3] *= vbegin;  //... Velocity state vector is scaled up from a radiant unit vector
	traj->params[4] *= vbegin; 
	traj->params[5] *= vbegin; 
    

	//======== Second pass PSO to estimate timing offsets only (all other parameters fixed)

	for( k=0; k<traj->numparams; k++ )  traj->xguess[k] = traj->params[k];
	for( k=0; k<traj->numparams; k++ )  traj->xshift[k] = 0.0;

    for( k=0; k<traj->numcameras; k++ )  traj->xshift[8+k] = traj->limits[8+k];  //... Vary timing offsets

	ParameterRefinementViaPSO( traj, MEAS2MODEL, CONSTANT, ACCURATE_FIT );

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "PSO: Timing Offsets Estimation", traj->params, traj->numcameras );


	//======== Third pass PSO to estimate velocity and the timing offsets (all other parameters fixed)

	for( k=0; k<traj->numparams;  k++ )  traj->xguess[k] = traj->params[k];
	for( k=0; k<traj->numparams;  k++ )  traj->xshift[k] = 0.0;

                                         traj->xshift[3]   = traj->limits[3];    //... Vary velocity vector
                                         traj->xshift[4]   = traj->limits[4];  
                                         traj->xshift[5]   = traj->limits[5];  
    for( k=0; k<traj->numcameras; k++ )  traj->xshift[8+k] = traj->limits[8+k];  //... Vary timing offsets

	ParameterRefinementViaPSO( traj, MEAS2MODEL, CONSTANT, ACCURATE_FIT );

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "PSO: Timing Offsets and Vinf Estimation", traj->params, traj->numcameras );


	//======== Fourth pass PSO to estimate velocity and deceleration terms only (all other parameters fixed)

	for( k=0; k<traj->numparams; k++ )  traj->xguess[k] = traj->params[k];
	for( k=0; k<traj->numparams; k++ )  traj->xshift[k] = 0.0;

                               traj->xshift[3] = traj->limits[3];  //... Vary velocity vector
                               traj->xshift[4] = traj->limits[4];  
                               traj->xshift[5] = traj->limits[5];  
    if( traj->velmodel >= 1 )  traj->xshift[6] = traj->limits[6];  //... Vary decel1
    if( traj->velmodel >= 2 )  traj->xshift[7] = traj->limits[7];  //... Vary decel2

	ParameterRefinementViaPSO( traj, MEAS2MODEL, traj->velmodel, ACCURATE_FIT );

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "PSO: Velocity and Deceleration Estimation", traj->params, traj->numcameras );


	//======== Fifth pass PSO to refine all parameters simultaneously to get FINAL SOLUTION

	for( k=0; k<traj->numparams; k++ )  traj->xguess[k] = traj->params[k];
	for( k=0; k<traj->numparams; k++ )  traj->xshift[k] = traj->limits[k];

	if( traj->velmodel <= 0 )  traj->xshift[6] = 0.0;  //... Constrain decel1 due to velocity model    
	if( traj->velmodel <= 1 )  traj->xshift[7] = 0.0;  //... Constrain decel2 due to velocity model  

	ParameterRefinementViaPSO( traj, MEAS2MODEL, traj->velmodel, ACCURATE_FIT );

	if( traj->verbose == 1 )  IntermediateConsoleOutput( "PSO: Solution on ALL parameters", traj->params, traj->numcameras );


	//======== Assign return parameters for final solution, radiant, velocity, deceleration, and time offsets

    for( k=0; k<traj->numparams; k++ )  traj->solution[k] = traj->params[k];

	ECI2RahDec( &traj->solution[3], &traj->ra_radiant, &traj->dec_radiant );

	VectorNorm( &traj->solution[3], radiant_hat, &traj->vbegin );

	traj->decel1 = fabs( traj->solution[6] );
	traj->decel2 = fabs( traj->solution[7] );
    
 	for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) traj->tref_offsets[kcamera] = traj->solution[8+kcamera];
   
	////for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) tref_offsets[kcamera] = params[8+kcamera] - params[8];


	//======== Assign return LLA, velocity and time parameters for the measurements 
	//         at their closest point of approach (CPA) to the radiant line for
	//         both measurements and model positions.

	ReportFill_LLAVT_Meas_Model( traj );


	//======== Assign LLA parameters for the BEGIN and END positions and associated
	//         time offsets relative to the reference Julian date/time.

	ReportFill_LLA_Beg_End( LLA_BEG, traj );

	ReportFill_LLA_Beg_End( LLA_END, traj );


	//======== Now run a noisy measurement Monte Carlo set of trials to estimate
	//         the error in each parameter. The no noise "solution" parameters
	//         will be used as the starting point for each minimization.

	MonteCarlo_ErrorEstimate( traj );


	//======== Diagnostic standard deviation output to console

	if( traj->verbose != 0 )  {
		printf("     Rah Dec sigma  %lf  %lf\n",        traj->ra_sigma * 57.29577951,     traj->dec_sigma * 57.29577951 );
		printf("     Vel Acc sigma  %lf  %lf  %lf\n",   traj->vbegin_sigma,               traj->decel1_sigma,               traj->decel2_sigma );
		printf("     LLA beg sigma  %lf  %lf  %lf\n",   traj->rbeg_lon_sigma * 57.29577951, traj->rbeg_lat_sigma * 57.29577951, traj->rbed_hkm_sigma );
		printf("     LLA end sigma  %lf  %lf  %lf\n\n", traj->rend_lon_sigma * 57.29577951, traj->rend_lat_sigma * 57.29577951, traj->rend_hkm_sigma );
	}


    return(0);
}


//##################################################################################
//
//======== Function to report on critical parameters during intermediate steps using
//         the "verbose" mode.
//==================================================================================

void  IntermediateConsoleOutput( char *ctitle, double *params, int ncams )
{
int     kcam;
double  rah, dec, radiant_hat[3], Vbegin;

	ECI2RahDec( &params[3], &rah, &dec );

	VectorNorm( &params[3], radiant_hat, &Vbegin );

    printf(" %s\n", ctitle );
    printf(" ------------------------------------------------------\n" );
	printf("     Radiant at %lf %lf \n", rah * 57.29577951, dec * 57.29577951 );
    printf("     Ro         %lf %lf %lf \n", params[0], params[1], params[2] );
    printf("     Vel Acc    %lf %lf %lf \n", Vbegin, fabs(params[6]), fabs(params[7]) );
    for( kcam=0; kcam<ncams; kcam++ ) printf("     Timing Offsets %lf \n", params[8+kcam] );
	printf("\n");

}


//##################################################################################
//
//======== Function to reset any standard deviations that were input as zero to the 
//         minimum non-zero value found amongst the full data set.
//==================================================================================

void    EnsureNonZeroStdDev( struct trajectory_info *traj )
{
int     kmeas, kcamera;
double  sigma, minsigma, maxsigma;


    //======== Find the smallest (nonzero) and largest standard deviation

    minsigma = +1.0e+30;
	maxsigma = -1.0e+30;

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {

		for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ )  {

			sigma = fabs( traj->noise[kcamera][kmeas] );

			if( sigma < minsigma  &&  sigma > 0.0 )  minsigma = sigma;

			if( sigma > maxsigma )  maxsigma = sigma;

		}  //... end of measurement loop per camera

	} //... end of camera loop


	//======== Set the smallest standard deviation and infill any zeroes

	if( maxsigma == 0.0 )  sigma = 1.0;       //... All input sigmas were zero
	else                   sigma = minsigma;  //... Smallest sigma greater than zero


	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {

		for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ )  {

			if( traj->noise[kcamera][kmeas] == 0.0 )  traj->noise[kcamera][kmeas] = sigma;

		}  //... end of measurement loop per camera

	} //... end of camera loop


}


//##################################################################################
//
//======== Function to convert angle-angle-time to measurement rays and obtain the
//         site position, both in ECI coordinates.
//==================================================================================

void   Angles2SiteMeasurements( double    camera_lat,  double   camera_lon,  double  camera_hkm,  double camera_LST,
		                        int       nummeas,     int      meastype,
								double   *meas_angle1, double  *meas_angle2, double *tsec,  
								int       noise_flag,  double  *noise,
								double  **mhat_camera, double  *rs_camera )
{
int     kmeas;
double  ecef[3], rhat[3], uhat[3], vhat[3], zhat[3], mvec[3], mag;
double  Rcenter, camera_lat_geocentric, sigma, azim, elev, ra, dec;


    //-------- Site position vectors in ECI coords relative to center of Earth
    
    LatLonAlt2ECEF( camera_lat, camera_lon, camera_hkm, ecef );
    
    Rcenter = sqrt( ecef[X]*ecef[X] + ecef[Y]*ecef[Y] + ecef[Z]*ecef[Z] );
    
    camera_lat_geocentric = atan( ecef[Z] / sqrt( ecef[X]*ecef[X] + ecef[Y]*ecef[Y] ) );
    
    rs_camera[X] = Rcenter * cos( camera_lat_geocentric ) * cos( camera_LST );
    rs_camera[Y] = Rcenter * cos( camera_lat_geocentric ) * sin( camera_LST );
    rs_camera[Z] = Rcenter * sin( camera_lat_geocentric );


    //-------- Measurement unit vectors in ECI (ECEF) coords

	for( kmeas=0; kmeas<nummeas; kmeas++ )  {

		if( meastype == RADEC )  {  //... no conversion

			ra  = meas_angle1[kmeas];

			dec = meas_angle2[kmeas];

		}
		else if( meastype == NAZEL )  {  //... convert from azimuth +east of north and altitude angle

			azim = meas_angle1[kmeas];

			elev = meas_angle2[kmeas];

			AzimuthElevation2RADec( azim, elev, camera_lat, camera_LST, &ra, &dec );

		}
		else if( meastype == SAZZA )  {  //... convert from azimuth +west of south and zenith angle

			azim = meas_angle1[kmeas] + 3.141592654;

			if( azim >= 2.0 * 3.141592654 )  azim -= 2.0 * 3.141592654;

			elev = 3.141592654 / 2.0 - meas_angle2[kmeas];

			AzimuthElevation2RADec( azim, elev, camera_lat, camera_LST, &ra, &dec );

		}
		else  {
		    printf(" ====> ERROR in Angles2SiteMeasurements: meastype %i not implemented\n", meastype );
		    exit(1);
	    }


		//-------- Convert equatorial to ECI coordinates

		RahDec2ECI( ra, dec, rhat );  //...Returns a unit vector

		if( rhat[Z] < 0.0 )  {  //... Southern Hemisphere
			zhat[X] =  0.0;
			zhat[Y] =  0.0;
			zhat[Z] = +1.0;
			VectorCrossProduct( rhat, zhat, uhat );
			VectorNorm( uhat, uhat, &mag );
			VectorCrossProduct( uhat, rhat, vhat );
			VectorNorm( vhat, vhat, &mag );
		}
		else  {                 //... Northern Hemisphere
			zhat[X] =  0.0;
			zhat[Y] =  0.0;
			zhat[Z] = -1.0;
			VectorCrossProduct( zhat, rhat, uhat );
			VectorNorm( uhat, uhat, &mag );
			VectorCrossProduct( uhat, rhat, vhat );
			VectorNorm( vhat, vhat, &mag );
		}

		//-------- Noise added assumes a small angle approximation that tangent(sigma in radians) ~ sigma

		if( noise_flag == 0 )  sigma = 0.0;
		else                   sigma = noise[kmeas] / sqrt(2.0);  //...sqrt(2)/2 * noisesigma in each orthogonal dimension

        mvec[X] = rhat[X] + RandomGauss( sigma, 3.0*sigma ) * uhat[X] + RandomGauss( sigma, 3.0*sigma ) * vhat[X];
        mvec[Y] = rhat[Y] + RandomGauss( sigma, 3.0*sigma ) * uhat[Y] + RandomGauss( sigma, 3.0*sigma ) * vhat[Y];
        mvec[Z] = rhat[Z] + RandomGauss( sigma, 3.0*sigma ) * uhat[Z] + RandomGauss( sigma, 3.0*sigma ) * vhat[Z];

        VectorNorm( mvec, mhat_camera[kmeas], &mag );  //... Normalize to a unit vector

	    mhat_camera[kmeas][T] = tsec[kmeas];   //... seconds relative to jdt reference time
			 
	}


}


//##################################################################################
//
//======== Function to compute the intersecting planes solution for two measurement
//         sets given number of measurements, ECI measurement rays X, Y, Z, and T,
//         and the site coords in ECI. Returns the radiant direction, convergence
//         angle, and CPA range vector to the 3D line from each site.
//==================================================================================

void   IntersectingPlanes( int     nmeas_camera1, double **meashat_camera1, double *rcamera1,
						   int     nmeas_camera2, double **meashat_camera2, double *rcamera2,
						   double *radiant_hat, double  *convergence_angle, double *rcpa_camera1, double *rcpa_camera2  )
{
int     kbeg, kend;
double  n1_hat[3], n2_hat[3], radiant[3], w1[3], w2[3], rcamera_diff[3];
double  dotw1w2, range_cpa1, range_cpa2, mag, cosangle;


    //======== This is a multi-ray solution to obtain the normal to each measurement plane.

	Normal2MeteorMeasurements( nmeas_camera1, meashat_camera1, n1_hat ); 

	Normal2MeteorMeasurements( nmeas_camera2, meashat_camera2, n2_hat );


	//======== Get the radiant unit vector via cross product of the normals

	VectorCrossProduct( n1_hat, n2_hat, radiant );

	VectorNorm( radiant, radiant_hat, &mag );

	//======== If closer to the anti-radiant, then reverse sign

    if( nmeas_camera1 >= 4 )  {
		kbeg = 1;
		kend = nmeas_camera1-2;
	}
	else  {
		kbeg = 0;
		kend = nmeas_camera1-1;
	}

	if( VectorDotProduct( meashat_camera1[kbeg], radiant_hat ) 
	  < VectorDotProduct( meashat_camera1[kend], radiant_hat ) )  {
	    radiant_hat[X] = -radiant_hat[X];
	    radiant_hat[Y] = -radiant_hat[Y];
	    radiant_hat[Z] = -radiant_hat[Z];
	}


	//======== Compute the convergence angle

	cosangle = VectorDotProduct( n1_hat, n2_hat );

	if( cosangle > +1.0 )  cosangle = +1.0;
	if( cosangle < -1.0 )  cosangle = -1.0;

	*convergence_angle = acos( fabs( cosangle ) );


	//======== Compute the cpa distance from the two cameras to the radiant line
	//         (i.e the range from a camera site to the line along the normal to
	//         the line) and return as a vector with that distance as magnitude.

	VectorCrossProduct( radiant_hat, n1_hat, w1 );

	VectorNorm( w1, w1, &mag );

	if( VectorDotProduct( w1, meashat_camera1[kbeg] ) < 0.0 )  {
		w1[X] = -w1[X];
		w1[Y] = -w1[Y];
		w1[Z] = -w1[Z];
	}


	VectorCrossProduct( radiant_hat, n2_hat, w2 );

	VectorNorm( w2, w2, &mag );

	if( VectorDotProduct( w2, meashat_camera2[kbeg] ) < 0.0 )  {
		w2[X] = -w2[X];
		w2[Y] = -w2[Y];
		w2[Z] = -w2[Z];
	}


	rcamera_diff[X] = rcamera1[X] - rcamera2[X];
	rcamera_diff[Y] = rcamera1[Y] - rcamera2[Y];
	rcamera_diff[Z] = rcamera1[Z] - rcamera2[Z];

	dotw1w2 = VectorDotProduct( w1, w2 );

    range_cpa1 = ( dotw1w2 * VectorDotProduct( rcamera_diff, w2 ) - VectorDotProduct( rcamera_diff, w1 ) )
		        / ( 1.0 - dotw1w2 * dotw1w2 );

	rcpa_camera1[X] = range_cpa1 * w1[X];
	rcpa_camera1[Y] = range_cpa1 * w1[Y];
	rcpa_camera1[Z] = range_cpa1 * w1[Z];


    range_cpa2 = ( VectorDotProduct( rcamera_diff, w2 ) - dotw1w2 * VectorDotProduct( rcamera_diff, w1 ) )
		        / ( 1.0 - dotw1w2 * dotw1w2 );

	rcpa_camera2[X] = range_cpa2 * w2[X];
	rcpa_camera2[Y] = range_cpa2 * w2[Y];
	rcpa_camera2[Z] = range_cpa2 * w2[Z];

}

//##################################################################################
//    
//======== Pair all track combinations and select the intersecting planes solution 
//         corresponding to the pair with the largest convergence angle.
//==================================================================================

void   IntersectingPlanes_BestConvergence( struct trajectory_info *traj, 
	                                       int    *kcamera_bestconv, double *max_convergence, 
										   double *radiant_bestconv, double *r_bestconv  )
{
int     kcamera1, kcamera2;
double  radiant_hat[3], convergence_angle, rcpa_camera1[3], rcpa_camera2[3];


	*max_convergence = 0.0;

	*kcamera_bestconv = 0;

	for( kcamera1=0; kcamera1<traj->numcameras; kcamera1++ ) {

		for( kcamera2=kcamera1+1; kcamera2<traj->numcameras; kcamera2++ ) {

	         IntersectingPlanes( traj->nummeas[kcamera1], traj->meashat_ECI[kcamera1], traj->rcamera_ECI[kcamera1],
					             traj->nummeas[kcamera2], traj->meashat_ECI[kcamera2], traj->rcamera_ECI[kcamera2],
					             radiant_hat, &convergence_angle, rcpa_camera1, rcpa_camera2 );

		     if( convergence_angle > *max_convergence )  {

                 *max_convergence = convergence_angle;

				 radiant_bestconv[X] = radiant_hat[X];
				 radiant_bestconv[Y] = radiant_hat[Y];
				 radiant_bestconv[Z] = radiant_hat[Z];

				 //... CPA position on the line w.r.t the camera in ECI (used later for starting position)

				 r_bestconv[X] = traj->rcamera_ECI[kcamera1][X] + rcpa_camera1[X];
				 r_bestconv[Y] = traj->rcamera_ECI[kcamera1][Y] + rcpa_camera1[Y];
				 r_bestconv[Z] = traj->rcamera_ECI[kcamera1][Z] + rcpa_camera1[Z];

				 *kcamera_bestconv = kcamera1;

			 }

		}  //... end kcamera2 loop

	}  //... end kcamera1 loop

}


//##################################################################################
//    
//======== Uniquely pair all track combinations and perform a weighted combination
//         of intersecting planes solutions. 
//==================================================================================

void   IntersectingPlanes_MultiTrack( struct trajectory_info *traj, double *radiant_hat  )
{
int     kcamera1, kcamera2, lastmeas;
double  radiant_sum[3], convergence_angle, rcpa_camera1[3], rcpa_camera2[3];
double  obsangle1, obsangle2, weight, wsum, rmagnitude;


    wsum = 0.0;

	radiant_sum[X] = 0.0;
	radiant_sum[Y] = 0.0;
	radiant_sum[Z] = 0.0;

	for( kcamera1=0; kcamera1<traj->numcameras; kcamera1++ ) {

		lastmeas = traj->nummeas[kcamera1] - 1;

		obsangle1 = acos( VectorDotProduct( traj->meashat_ECI[kcamera1][0], traj->meashat_ECI[kcamera1][lastmeas] ) );

		for( kcamera2=kcamera1+1; kcamera2<traj->numcameras; kcamera2++ ) {

		    lastmeas = traj->nummeas[kcamera2] - 1;

		    obsangle2 = acos( VectorDotProduct( traj->meashat_ECI[kcamera2][0], traj->meashat_ECI[kcamera2][lastmeas] ) );

	        IntersectingPlanes( traj->nummeas[kcamera1], traj->meashat_ECI[kcamera1], traj->rcamera_ECI[kcamera1],
						        traj->nummeas[kcamera2], traj->meashat_ECI[kcamera2], traj->rcamera_ECI[kcamera2],
						        radiant_hat, &convergence_angle, rcpa_camera1, rcpa_camera2 );

			weight = obsangle1 * obsangle2 * sin( convergence_angle ) * sin( convergence_angle );

			wsum += weight;

			radiant_sum[X] += weight * radiant_hat[X];
			radiant_sum[Y] += weight * radiant_hat[Y];
			radiant_sum[Z] += weight * radiant_hat[Z];

		}  //... end kcamera2 loop

	}  //... end kcamera1 loop

	radiant_sum[X] /= wsum;
	radiant_sum[Y] /= wsum;
	radiant_sum[Z] /= wsum;

	VectorNorm( radiant_sum, radiant_hat, &rmagnitude ); 

}


//##################################################################################
//
//======== Function to compute the normal to the plane defined by a set of 
//         measurement rays from one camera that followed the meteor track. The 
//         normal is in the same coords as the measurement unit vectors (typ. ECI).
//==================================================================================

void  Normal2MeteorMeasurements( int nmeas, double **meashat_ECI, double *planenormal )
{
int     k, kbeg, kend;
double  xdotx, xdoty, xdotz, ydoty, ydotz;
double  sx, sy, smag, denom;


    //======== The plane is really undefined for a single measurement ray

    if( nmeas == 1 )  {  
		smag = sqrt( meashat_ECI[0][X] * meashat_ECI[0][X] + meashat_ECI[0][Y] * meashat_ECI[0][Y] );
		planenormal[X] = -meashat_ECI[0][Y] / smag;
		planenormal[Y] = +meashat_ECI[0][X] / smag;
		planenormal[Z] =  0.0;
		return;
	}


    //======== Compute the plane using the second through second-to-last measurements
    //         unless there are fewer than four measurements, then first to last.

    if( nmeas >= 4 )  {
		kbeg = 1;
		kend = nmeas-2;
	}
	else  {
		kbeg = 0;
		kend = nmeas-1;
	}


    //======== Compute running sums of dot products

    xdotx = 0.0;
    xdoty = 0.0;
    xdotz = 0.0;
    ydoty = 0.0;
    ydotz = 0.0;

	for( k=kbeg; k<=kend; k++ )  {

		 xdotx += meashat_ECI[k][X] * meashat_ECI[k][X];
		 xdoty += meashat_ECI[k][X] * meashat_ECI[k][Y];
		 xdotz += meashat_ECI[k][X] * meashat_ECI[k][Z];
		 ydoty += meashat_ECI[k][Y] * meashat_ECI[k][Y];
		 ydotz += meashat_ECI[k][Y] * meashat_ECI[k][Z];
    }


	//======== Solve for the unit vector normal to the meteor plane

	denom = xdotx * ydoty - xdoty * xdoty;

	sx = ( ydoty * xdotz - xdoty * ydotz ) / denom;

	sy = ( xdotx * ydotz - xdoty * xdotz ) / denom;

	smag = sqrt( 1.0 + sx*sx + sy*sy );

	planenormal[X] =   sx / smag;
	planenormal[Y] =   sy / smag;
	planenormal[Z] = -1.0 / smag;

}

//##################################################################################
//
//======== Function to compute the points on two lines that are the closest point
//         of approach between the two lines. Each line is defined by a 3D point 
//         and a 3D unit vector along the line. The points on the two lines that 
//         are at CPA are "r1" and "r2" with the distance between them "d21".
//==================================================================================

void   TwoLineCPA( double *position1, double *vector1, double *r1, 
				   double *position2, double *vector2, double *r2, double *d21 )
{
double  d1, d2, mag_dummy, unitvector1[3], unitvector2[3], p21[3], rxm[3], pxm[3], mdotm;


    VectorNorm( vector1, unitvector1, &mag_dummy );
    VectorNorm( vector2, unitvector2, &mag_dummy );

    p21[0] = position1[0] - position2[0];
    p21[1] = position1[1] - position2[1];
    p21[2] = position1[2] - position2[2];

	VectorCrossProduct( unitvector1, unitvector2, rxm );

	VectorCrossProduct( p21, rxm, pxm );

	mdotm = VectorDotProduct( rxm, rxm );

	d1 = VectorDotProduct( pxm, unitvector2 ) / mdotm;

	r1[0] = position1[0] + d1 * unitvector1[0];
	r1[1] = position1[1] + d1 * unitvector1[1];
	r1[2] = position1[2] + d1 * unitvector1[2];

	d2 = VectorDotProduct( pxm, unitvector1 ) / mdotm;

	r2[0] = position2[0] + d2 * unitvector2[0];
	r2[1] = position2[1] + d2 * unitvector2[1];
	r2[2] = position2[2] + d2 * unitvector2[2];

	*d21 = fabs( VectorDotProduct( p21, rxm ) ) / sqrt( mdotm );

}


//##################################################################################
//
//======== Function to compute a motion model for either position or velocity
//         propagation of the meteor. It uses either a constant velocity motion model,
//         a linear deceleration model, a quadratic deceleration model, or Jacchia's
//         exponential motion model.
//===================================================================================


double  Propagation( int propagation_state, int velmodel, double t, double vbegin, double decel1, double decel2 )
{
    
	if( propagation_state == POSITION )  {

		if(      velmodel == CONSTANT  ) return( fabs(vbegin) * t );
		else if( velmodel == LINEAR    ) return( fabs(vbegin) * t - fabs(decel1) * t * t / 2.0 );
		else if( velmodel == QUADRATIC ) return( fabs(vbegin) * t - fabs(decel1) * t * t / 2.0 + decel2 * t * t * t / 3.0 );
		else if( velmodel == EXPONENT  ) return( fabs(vbegin) * t - fabs(decel1) * exp( fabs(decel2) * t ) );

		else printf("ERROR--> In Propagation of MeteorTrajectory: velocity model %d not implemented \n", velmodel );
	}


	else if( propagation_state == VELOCITY )  {

		if(      velmodel == CONSTANT  ) return( fabs(vbegin) );
		else if( velmodel == LINEAR    ) return( fabs(vbegin) - fabs(decel1) * t );
		else if( velmodel == QUADRATIC ) return( fabs(vbegin) - fabs(decel1) * t + decel2 * t * t );
		else if( velmodel == EXPONENT  ) return( fabs(vbegin) - fabs(decel1 * decel2) * exp( fabs(decel2) * t ) );

        else printf("ERROR--> In Propagation of MeteorTrajectory: velocity model %d not implemented \n", velmodel );
	}

	else  printf("ERROR--> In Propagation of MeteorTrajectory: propagation_state %d not implemented \n", propagation_state );

	return(-999.0);

}

//##################################################################################
//
//======== Function to find the velocity from a set of positions and times along
//         the radiant line by differencing adjacent-in-time positions. The set
//         of velocities are averaged and their standard deviation computed so
//         an iterative outlier removal can be used to obtain a more robust mean
//         velocity estimate.
//==================================================================================

double   VelocityFit_Differencing( struct trajectory_info *traj )
{
int       kcamera, kmeas, ntotalmeas, nmeas, k, kbeg, kend, kloop;
double   *vel, velsum, velssq, velave, velstd;
double    time_difference, rbeg[3], rend[3], rdummy[3], distance;


    //======== Allocate memory for velocity measurements

	ntotalmeas = 0;

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  ntotalmeas += traj->nummeas[kcamera];

    vel = (double*) malloc( ntotalmeas * sizeof(double) );


	//======== Compute the positional and time differencs to get the velocity mean and std dev
	            
	velsum = 0.0;
	velssq = 0.0;
	kmeas = 0;
    
	for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) {
        
		//-------- Avoid first and last measurement, if possible

		kbeg = 1;                        
		kend = traj->nummeas[kcamera]-2;
        
 		if( traj->nummeas[kcamera] == 3 )  {
			kbeg = 0;
			kend = 2;
		}

		if( traj->nummeas[kcamera] == 2 )  {
			kbeg = 0;
			kend = 1;
		}

		if( traj->nummeas[kcamera] <= 1 )  continue;

        //-------- Since this is a vbegin estimate - use no more than first
        //            60 interleaved frames (~1 second duration)

        if( traj->nummeas[kcamera]-2 > 60 )  kend = 60;

		//-------- Compute velocity for each temporally adjacent pair of measurements

		for( k=kbeg; k<kend; k++ )  {

			//........ Get an absolute ECI position along the measurement ray that has a CPA to the radiant line

            TwoLineCPA( traj->rcamera_ECI[kcamera], traj->meashat_ECI[kcamera][k], rbeg, 
		                &traj->params[0], &traj->params[3], rdummy, &distance );

            TwoLineCPA( traj->rcamera_ECI[kcamera], traj->meashat_ECI[kcamera][k+1], rend, 
		                &traj->params[0], &traj->params[3], rdummy, &distance );

			distance = sqrt( ( rend[X] - rbeg[X] ) * ( rend[X] - rbeg[X] ) + 
                             ( rend[Y] - rbeg[Y] ) * ( rend[Y] - rbeg[Y] ) +
		                     ( rend[Z] - rbeg[Z] ) * ( rend[Z] - rbeg[Z] )   ); 

	        time_difference = traj->meashat_ECI[kcamera][k+1][T] - traj->meashat_ECI[kcamera][k][T];

	        vel[kmeas] = distance / time_difference;   // Velocity estimate

			velsum += vel[kmeas];
			velssq += vel[kmeas] * vel[kmeas];

			kmeas++;

		}
	}

	nmeas = kmeas;
	velave = velsum / (double)kmeas;
	velstd = sqrt( ( velssq - (double)kmeas * velave * velave ) / (double)(kmeas-1) );


	//........ Iterate to remove outliers from mean velocity calculation

	for( kloop=0; kloop<4; kloop++ )  {

	    velsum = 0.0;
        velssq = 0.0;
	    kmeas  = 0;

	    for( k=0; k<nmeas; k++ )  {
		    if( fabs( vel[k] - velave ) < 2.0 * velstd )  {
			    velsum += vel[k];
			    velssq += vel[k] * vel[k];
			    kmeas++;
		    }
	    }

		if( kmeas < 3 )  break;

	    velave = velsum / (double)kmeas;
	    velstd = sqrt( ( velssq - (double)kmeas * velave * velave ) / (double)(kmeas-1) );

	}

    free( vel );

	return( velave );

}

//##################################################################################
//
//======== Function to find the velocity from a set of positions and times along
//         the radiant line by using an LMS fit to multiple measurement sequences. 
//         The RMS error is used on a second pass to remove outliers in the final 
//         LMS fit.
//==================================================================================

double   VelocityFit_LMS( struct trajectory_info *traj )
{
int       kcamera, kmeas, ntotalmeas, k, kbeg, kend, kmeas_kept, nmeas_kept;
int      *nmeas_per_camera;
double   *pos, *tim, *noi, *xo, velocityLMS, rmserror, *err;
double    rbeg[3], rend[3], rdummy[3], distance;



    //======== Allocate memory for the positional measurements and time stamps actually used, 
    //         plus the measurement count and starting position per camera 

	ntotalmeas = 0;

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  ntotalmeas += traj->nummeas[kcamera];

    pos           = (double*) malloc( ntotalmeas       * sizeof(double) );
    tim           = (double*) malloc( ntotalmeas       * sizeof(double) );
    noi           = (double*) malloc( ntotalmeas       * sizeof(double) );
    err           = (double*) malloc( ntotalmeas       * sizeof(double) );
    xo            = (double*) malloc( traj->numcameras * sizeof(double) );
	nmeas_per_camera = (int*) malloc( traj->numcameras * sizeof(int)    );


	//======== Infill the positional and temporal measurements making a single vector
	//         of concatenated measurements for all the cameras.
	            
	kmeas = 0;
    
	for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) {
        
		//-------- Avoid first and last measurement, unless 3 or less measurements

		kbeg = 1;                        
		kend = traj->nummeas[kcamera] - 2;
        
 		if( traj->nummeas[kcamera] <= 3 )  {
			kbeg = 0;
			kend = traj->nummeas[kcamera] - 1;
		}

        if( kend > 30 )  kend = 30; //... vbegin calculation, so restrict to early frames


		nmeas_per_camera[kcamera] = kend - kbeg + 1;


		//-------- Set reference position to estimated begin point on the radiant line

		rbeg[X] = traj->params[0];
		rbeg[Y] = traj->params[1];
		rbeg[Z] = traj->params[2];


		//-------- Compute each positional and temporal measurement

		for( k=kbeg; k<=kend; k++ )  {

			//........ Get an absolute ECI position along the measurement ray at its CPA to the radiant line

            TwoLineCPA( traj->rcamera_ECI[kcamera], traj->meashat_ECI[kcamera][k], rend, 
		                &traj->params[0], &traj->params[3], rdummy, &distance );

			//........ Compute the distance from the reference point along the line to the measurement CPA

			pos[kmeas] = sqrt( ( rend[X] - rbeg[X] ) * ( rend[X] - rbeg[X] ) + 
                               ( rend[Y] - rbeg[Y] ) * ( rend[Y] - rbeg[Y] ) +
		                       ( rend[Z] - rbeg[Z] ) * ( rend[Z] - rbeg[Z] )   ); 

			noi[kmeas] = traj->noise[kcamera][k];

			//........ Time of the measurement can be relative for the LMS

	        tim[kmeas] = traj->meashat_ECI[kcamera][k][T];

			//printf(" x t  %lf  %lf\n", pos[kmeas], tim[kmeas] );

			kmeas++;

		}  // end of measurement loop per camera

	}  //... end of camera loop


	//======== LMS solution for velocity and starting positions per camera

	ConstantVelocity_MultiTrackFit( traj->numcameras, nmeas_per_camera, pos, tim, noi, &velocityLMS, xo, &rmserror, err );


	//======== Perform outlier removal set at 2x the rmserror

	kmeas      = 0;
	kmeas_kept = 0;

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) {

		nmeas_kept = 0;
        
		for( k=0; k<nmeas_per_camera[kcamera]; k++ )  {

			if( fabs( err[kmeas] ) <= 2.0 * rmserror )  {

				pos[kmeas_kept] = pos[kmeas];
				tim[kmeas_kept] = tim[kmeas];
				noi[kmeas_kept] = noi[kmeas];
				kmeas_kept++;
				nmeas_kept++;

			}

			kmeas++;

		}  // end of measurement loop per camera

		nmeas_per_camera[kcamera] = nmeas_kept;

	}  //... end of camera loop


	//======== LMS solution for velocity with outliers removed

	ConstantVelocity_MultiTrackFit( traj->numcameras, nmeas_per_camera, pos, tim, noi, &velocityLMS, xo, &rmserror, err );


	//======== Free memory and return LMS velocity estimate

	free( pos );
    free( tim );
    free( noi );
    free( err );
	free( xo  );
	free( nmeas_per_camera );


	return ( velocityLMS );

}

//############################################################################################

//======== This function computes a weighted LMS fit to velocity and starting positions of 
//         sequences of position and time measurements obtained from several independent 
//         tracks. Assumes all the tracks have the same constant motion velocity "velocityLMS"
//         but that each track is not synchronized in any way to another track sequence in 
//         time or starting position "xo". The input data is comprised of a set of concatenated 
//         track measurements in vectors for position "pos", time "tim", and standard deviation
//         "noi". The weight is the inverse of the variance per measurement.
//
//         Thus the user inputs the number of sequences or tracks "ncameras" and the number of
//         measurements for each camera "nmeas_per_camera" which does NOT have to be the same
//         measurement count per camera.
//============================================================================================

void  ConstantVelocity_MultiTrackFit( int ncameras, int *nmeas_per_camera, double *pos, double *tim, double *noi, double *velocityLMS, double *xo, double *rmserror, double *err )
{
int     kcamera, kmeas, k;
double  sumw, sumwx, sumwt, sumwxt, sumwtt, sumwtsumwtsumw, sumwxsumwtsumw, sumesq, weight;


    //======== Compute sums for the LMS solution
	            
    sumwxt         = 0.0;
    sumwtt         = 0.0;
	sumwtsumwtsumw = 0.0;
	sumwxsumwtsumw = 0.0;

	kmeas = 0;
    
	for( kcamera=0; kcamera<ncameras; kcamera++ ) {
        
		sumw  = 0.0;
        sumwx = 0.0;
        sumwt = 0.0;

		if( nmeas_per_camera[kcamera] > 0 )  {

		    for( k=0; k<nmeas_per_camera[kcamera]; k++ )  {

				if( noi[kmeas] != 0.0 )  weight = 1.0 / ( noi[kmeas] * noi[kmeas] );
				else                     weight = 1.0;

				sumw   += weight;
                sumwx  += weight * pos[kmeas];        
                sumwt  += weight * tim[kmeas];
        
                sumwxt += weight * pos[kmeas] * tim[kmeas];
                sumwtt += weight * tim[kmeas] * tim[kmeas];

			    kmeas++;

		    }  // end of measurement loop per camera

            sumwtsumwtsumw += sumwt * sumwt / sumw;
            sumwxsumwtsumw += sumwx * sumwt / sumw;

		}

	}  //... end of camera loop

        
    //======== LMS solve for the velocity

    *velocityLMS = ( sumwxt - sumwxsumwtsumw ) / ( sumwtt - sumwtsumwtsumw );


    //======== LMS solve for the starting position
	
	kmeas = 0;
    
	for( kcamera=0; kcamera<ncameras; kcamera++ ) {
        
		sumw  = 0.0;
        sumwx = 0.0;
        sumwt = 0.0;

		if( nmeas_per_camera[kcamera] > 0 )  {

			for( k=0; k<nmeas_per_camera[kcamera]; k++ )  {
        
				if( noi[kmeas] != 0.0 )  weight = 1.0 / ( noi[kmeas] * noi[kmeas] );
				else                     weight = 1.0;

				sumw  += weight;
                sumwx += weight * pos[kmeas];
                sumwt += weight * tim[kmeas];

			    kmeas++;

		    }  // end of measurement loop per camera

            xo[kcamera] = ( sumwx - *velocityLMS * sumwt ) / sumw;

		}

		else  xo[kcamera] = 0.0;

	}  //... end of camera loop

	
	//======== Determine the root-mean-square error

	sumesq = 0.0;

	kmeas  = 0;

	for( kcamera=0; kcamera<ncameras; kcamera++ ) {
        
		for( k=0; k<nmeas_per_camera[kcamera]; k++ )  {

			err[kmeas] = pos[kmeas] - xo[kcamera] - *velocityLMS * tim[kmeas];

			sumesq += err[kmeas] * err[kmeas];

			kmeas++;

		}

	}

	*rmserror = sqrt( sumesq / (double)kmeas );

	
}


//##################################################################################
//
//======== Function used to refine estimation parameters by calling the PSO
//         minimization module. Two levels of accuracy are used based on need:
//         a high quality "noise-free" solution and a multiple Monte Carlo low
//         quality soutions but done very quickly.
//=================================================================================


void   ParameterRefinementViaPSO( struct trajectory_info *traj, int cost_function, int velocity_model, int pso_accuracy_config )
{
int     k, number_particles, maximum_iterations, boundary_flag, limits_flag, particle_distribution_flag;
double  epsilon_convergence, weight_inertia, weight_stubborness, weight_grouppressure;
double  cost_func_value;

struct particleswarming  pso;


    //======== Particle Swarm Optimizer (PSO) default settings

    if( pso_accuracy_config == ACCURATE_FIT )  {
        number_particles           = 100;
	    maximum_iterations         = 1000;
	    boundary_flag              = BOUNDARY_REFLECTIVE;
	    limits_flag                = LIMITS_ARE_LOOSE;
	    particle_distribution_flag = PARTICLEDISTRO_RANDOM;
	    epsilon_convergence        = 1.0e-10;
	    weight_inertia             = 0.8;
	    weight_stubborness         = 1.0;
	    weight_grouppressure       = 2.0;
	}
	else if( pso_accuracy_config == QUICK_FIT )  {
        number_particles           = 20;
	    maximum_iterations         = 300;
	    boundary_flag              = BOUNDARY_REFLECTIVE;
	    limits_flag                = LIMITS_ARE_STRICT;
	    particle_distribution_flag = PARTICLEDISTRO_GAUSS;
	    epsilon_convergence        = 1.0e-7;
	    weight_inertia             = 0.8;
	    weight_stubborness         = 1.0;
	    weight_grouppressure       = 2.0;
	}
	else  {
		printf("ERROR--> Option not implemented for pso_accuracy_config = %d in ParameterRefinementViaPSO\n");
		exit(1);
	}


	//======== Pre-populate the PSO and initialize with guess and shifts (zero shifts for fixed parameters)

	ParticleSwarm_PrePopulate( number_particles, traj->numparams, maximum_iterations, 
	                           boundary_flag, limits_flag, particle_distribution_flag, epsilon_convergence,
	                           weight_inertia, weight_stubborness, weight_grouppressure,  								   
						       &pso );

	ParticleSwarm_Initialize( traj->xguess, traj->xshift, &pso );


	//======== PSO minimization loop

	while( pso.processing == CONTINUE_PROCESSING )  {

		if(      cost_function == MEAS2LINE  )  cost_func_value = Anglesum_Measurements2Line( &pso.xtest[0], 
			                                                                                  &pso.xtest[3], 
																						      traj );

		else if( cost_function == MEAS2MODEL )  cost_func_value = Anglesum_Measurements2Model( velocity_model, 
			                                                                                   &pso.xtest[0], 
																						       &pso.xtest[3], 
																						       pso.xtest[6], 
																						       pso.xtest[7], 
																						       &pso.xtest[8], 
																						       traj );

		else  {
			printf("ERROR--> cost_function %d not implemented in ParameterRefinementViaPSO\n");
			exit(1);
		}


		ParticleSwarm_Update( cost_func_value, &pso );

	}

	////printf(" Halt Condition %d,  Number of iterations %d\n", pso.processing, pso.niteration );


	//======== Assign refined parameters, cleanup memory, and return

	for( k=0; k<traj->numparams; k++ )  traj->params[k] = pso.gbest[k];

	ParticleSwarm_PostCleanup( &pso );

}


//##################################################################################
//
//======== Function to sum the angles between the radiant line and the measurement
//         rays to provide a cost function value to minimize against. Used for
//         Borovicka's least squares radiant solution but replaces the distance
//         sum for an angle sum for the cost function that is being minimized.
//         The cost function is weighted by the inverse of the noise variance.
//==================================================================================


double  Anglesum_Measurements2Line( double *ro_ECI, double *vo_ECI, struct trajectory_info *traj )
{
int     kcamera, kmeas;
double  fsum, cosangle, r[3], rmeas[3], radiant_hat_ECI[3], rmagnitude, cpa_distance, wsum, weight;



    fsum = 0.0;
	wsum = 0.0;

	VectorNorm( vo_ECI, radiant_hat_ECI, &rmagnitude );  //... Unit vector for the radiant direction

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {

		for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ )  {

			//-------- Find the ECI vector to point "r" on the radiant line closest to the the measurement ray

            TwoLineCPA( traj->rcamera_ECI[kcamera], traj->meashat_ECI[kcamera][kmeas], rmeas, 
				        ro_ECI, radiant_hat_ECI, r,   
					    &cpa_distance );

			//-------- The weighting of 1/variance is based on the per measurment noise standard deviation

			weight = 1.0 / ( traj->noise[kcamera][kmeas] * traj->noise[kcamera][kmeas] + 1.0e-10 );

			wsum += weight;

			//-------- CPA ray's unit vector to the point on the line from the camera position

            r[X] = r[X] - traj->rcamera_ECI[kcamera][X];
            r[Y] = r[Y] - traj->rcamera_ECI[kcamera][Y];
            r[Z] = r[Z] - traj->rcamera_ECI[kcamera][Z];

	        VectorNorm( r, r, &rmagnitude );

			//-------- Angle between measurement ray and the CPA ray as seen from the camera

			cosangle = VectorDotProduct( traj->meashat_ECI[kcamera][kmeas], r );

	        if( cosangle > +1.0 )  cosangle = +1.0;
	        if( cosangle < -1.0 )  cosangle = -1.0;

			fsum += weight * acos( cosangle );
		}

	}

	fsum /= wsum;

	return( fsum );

}

//##################################################################################
//
//======== Function to sum the angles between the motion model positions and the 
//         measurement rays to provide a cost function value to minimize against.
//         The cost function is weighted by the inverse of the noise variance.
//==================================================================================

double  Anglesum_Measurements2Model( int velmodel,  double *ro_ECI, double *vo_ECI, 
				                     double decel1,  double decel2, double *dt,
				                     struct trajectory_info *traj )
{
int     kcamera, kmeas;
double  fsum, tt, length_km, r[3], radiant_hat_ECI[3];
double  vbegin, rmagnitude, cosangle, weight, wsum;



    fsum = 0.0;
	wsum = 0.0;

	VectorNorm( vo_ECI, radiant_hat_ECI, &vbegin );

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {

		for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ )  {

			//-------- Compute the model time from the measurement time plsu the camera time offset

			tt = traj->meashat_ECI[kcamera][kmeas][T] + dt[kcamera];

			//-------- The weighting of 1/variance is based on the per measurment noise standard deviation 

			weight = 1.0 / ( traj->noise[kcamera][kmeas] * traj->noise[kcamera][kmeas] + 1.0e-10 );

			wsum += weight;

			//-------- Propagation distance at that time

		    length_km = Propagation( POSITION, velmodel, tt, vbegin, decel1, decel2 );

			//-------- Model ray's unit vector from the camera position

            r[X] = ro_ECI[X] - traj->rcamera_ECI[kcamera][X] - radiant_hat_ECI[X] * length_km;
            r[Y] = ro_ECI[Y] - traj->rcamera_ECI[kcamera][Y] - radiant_hat_ECI[Y] * length_km;
            r[Z] = ro_ECI[Z] - traj->rcamera_ECI[kcamera][Z] - radiant_hat_ECI[Z] * length_km;

			VectorNorm( r, r, &rmagnitude );

			//-------- Angle between the measurement ray and model ray as seen from the camera

			cosangle = VectorDotProduct( traj->meashat_ECI[kcamera][kmeas], r );

	        if( cosangle > +1.0 )  cosangle = +1.0;
	        if( cosangle < -1.0 )  cosangle = -1.0;

			fsum += weight * acos( cosangle );

		}

	}

	fsum /= wsum;

	return( fsum );

}


//##################################################################################

//======== Assign return LLA parameters for measurement directions at their closest
//         point of approach (CPA) to the radiant line, plus the associated range.
//         Also get the model position along the radiant line and the associated
//         modeled velocity computed ON the radiant line at the TIMES input and 
//         adjusted for the timing offsets and NOT directly associated with the
//         measurements. Note that LLA is latitude, longitude, altitude in 
//         geodetic WGS84 coordinates for the reference Julian date/time.
//==================================================================================

void   ReportFill_LLAVT_Meas_Model( struct trajectory_info *traj )
{
int     kcamera, kmeas;
double  radiant_hat[3], rdummy[3], r[3], rfit[3], r_lat, r_lon, r_LST, r_hkm;
double  dist, vbegin, tt, length_km, rafit, decfit;
double  azim, elev, sazim, zenangle;


	VectorNorm( &traj->solution[3], radiant_hat, &vbegin );


	//======== Loop over all cameras and their associated measurements

	for( kcamera=0; kcamera<traj->numcameras; kcamera++ ) {

		for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ ) {

			 //-------- Get measurement ray positional info at the closest point of approach to the radiant line

			 TwoLineCPA( traj->rcamera_ECI[kcamera], traj->meashat_ECI[kcamera][kmeas], r, 
		                 &traj->solution[0], radiant_hat, rdummy, &dist );

	         ECEF2LatLonAlt( r, &r_lat, &r_LST, &r_hkm );

	         r_lon = LST2LongitudeEast( traj->jdt_ref, r_LST );

			 traj->meas_lat[kcamera][kmeas] = r_lat;
			 traj->meas_lon[kcamera][kmeas] = r_lon;
			 traj->meas_hkm[kcamera][kmeas] = r_hkm;

			 traj->meas_range[kcamera][kmeas] = sqrt( (r[X] - traj->rcamera_ECI[kcamera][X]) * (r[X] - traj->rcamera_ECI[kcamera][X])
				                                    + (r[Y] - traj->rcamera_ECI[kcamera][Y]) * (r[Y] - traj->rcamera_ECI[kcamera][Y])
				                                    + (r[Z] - traj->rcamera_ECI[kcamera][Z]) * (r[Z] - traj->rcamera_ECI[kcamera][Z]) );


			 //-------- Compute the measurement time with timing offsets

			 tt = traj->meashat_ECI[kcamera][kmeas][T] + traj->solution[8+kcamera];

			 traj->model_time[kcamera][kmeas] = tt;


			 //-------- Get model positional info at each time position along the radiant line

		     length_km = Propagation( POSITION, traj->velmodel, tt, traj->vbegin, traj->decel1, traj->decel2 );

             r[X] = traj->solution[0] - radiant_hat[X] * length_km;
             r[Y] = traj->solution[1] - radiant_hat[Y] * length_km;
             r[Z] = traj->solution[2] - radiant_hat[Z] * length_km;

 	         ECEF2LatLonAlt( r, &r_lat, &r_LST, &r_hkm );

	         r_lon = LST2LongitudeEast( traj->jdt_ref, r_LST );

			 traj->model_lat[kcamera][kmeas] = r_lat;
			 traj->model_lon[kcamera][kmeas] = r_lon;
			 traj->model_hkm[kcamera][kmeas] = r_hkm;

			 traj->model_range[kcamera][kmeas] = sqrt( (r[X] - traj->rcamera_ECI[kcamera][X]) * (r[X] - traj->rcamera_ECI[kcamera][X])
				                                     + (r[Y] - traj->rcamera_ECI[kcamera][Y]) * (r[Y] - traj->rcamera_ECI[kcamera][Y])
				                                     + (r[Z] - traj->rcamera_ECI[kcamera][Z]) * (r[Z] - traj->rcamera_ECI[kcamera][Z]) );
			 

			 //-------- Compute the model velocity at eash time (duplicate this value for the measurement)

			 traj->model_vel[kcamera][kmeas] = Propagation( VELOCITY, traj->velmodel, tt, traj->vbegin, traj->decel1, traj->decel2 );

			 traj->meas_vel[kcamera][kmeas] = traj->model_vel[kcamera][kmeas];


			 //-------- Compute the model ray unit vector "fit" in ECI and convert to "meastype" for output reporting

			 rfit[X] = r[X] - traj->rcamera_ECI[kcamera][X];
             rfit[Y] = r[Y] - traj->rcamera_ECI[kcamera][Y];
             rfit[Z] = r[Z] - traj->rcamera_ECI[kcamera][Z];

			 ECI2RahDec( rfit, &rafit, &decfit );

			 if( traj->meastype == RADEC )  {
                 traj->model_fit1[kcamera][kmeas] = rafit;
			     traj->model_fit2[kcamera][kmeas] = decfit;
			 }

	         else if( traj->meastype == NAZEL )  {
				 RADec2AzimuthElevation( rafit, decfit, traj->camera_lat[kcamera], traj->camera_LST[kcamera], &azim, &elev );
                 traj->model_fit1[kcamera][kmeas] = azim;
			     traj->model_fit2[kcamera][kmeas] = elev;
		     }

	         else if( traj->meastype == SAZZA )  {  
 				 RADec2AzimuthElevation( rafit, decfit, traj->camera_lat[kcamera], traj->camera_LST[kcamera], &azim, &elev );
				 sazim = azim + 3.141592654;
				 if( sazim >= 2.0 * 3.141592654 )  sazim -= 2.0 * 3.141592654;
				 zenangle = 3.141592654 / 2.0 - elev;
                 traj->model_fit1[kcamera][kmeas] = sazim;
			     traj->model_fit2[kcamera][kmeas] = zenangle;
			 }

	         else  {
		         printf(" ====> ERROR in ReportFill_LLAVT_Meas_Model: meastype %i not implemented for meas_fit*\n", traj->meastype );
		         exit(1);
	         }


		}  //... end of measurement loop per camera
	
	}  //... end of camera loop


}

//##################################################################################

//======== Assign return LLA parameters for the BEGIN or END position and the time 
//         offset relative to the reference Julian date/time by finding the 
//         earliest or latest measurement respectively. This point is on the radiant
//         line as specified by the adjusted time - thus is NOT tied to the
//         measurement unit vector which may be of poor quality because it has
//         low illumination level and inaccurate centroid.
//==================================================================================

void    ReportFill_LLA_Beg_End( int LLA_position, struct trajectory_info *traj )
{
int     kcamera, kcamera1, kmeas;
double  tt, length_km, vbegin;
double  r_lat, r_lon, r_LST, r_hkm, r[3], radiant_hat[3];


    //======== Get the earliest first measurement time

    if( LLA_position == LLA_BEG )  {

	    tt = +1.0e+20;

	    for( kcamera1=0; kcamera1<traj->numcameras; kcamera1++ ) {

			kmeas = 0;

		    if( tt > traj->meashat_ECI[kcamera1][kmeas][T] + traj->tref_offsets[kcamera1] )  {
			    tt = traj->meashat_ECI[kcamera1][kmeas][T] + traj->tref_offsets[kcamera1];
			    kcamera = kcamera1;
		    }
	    }

	    traj->begtime = tt;

		kmeas = 0;
	}


    //======== Get the latest last measurement time

    if( LLA_position == LLA_END )  {

	    tt = -1.0e+20;

	    for( kcamera1=0; kcamera1<traj->numcameras; kcamera1++ ) {

	        kmeas = traj->nummeas[kcamera1] - 1;

		    if( tt < traj->meashat_ECI[kcamera1][kmeas][T] + traj->tref_offsets[kcamera1] )  {
			    tt = traj->meashat_ECI[kcamera1][kmeas][T] + traj->tref_offsets[kcamera1];
			    kcamera = kcamera1;
		    }
	    }

	    traj->endtime = tt;

		kmeas = traj->nummeas[kcamera] - 1;
	}


	//======== Get the time for the desired measurement point, compute its ECI coords, convert to LLA

	length_km = Propagation( POSITION, traj->velmodel, tt, traj->vbegin, traj->decel1, traj->decel2 );

	VectorNorm( &traj->solution[3], radiant_hat, &vbegin );

	r[X] = traj->solution[0] - radiant_hat[X] * length_km; 
	r[Y] = traj->solution[1] - radiant_hat[Y] * length_km;  
	r[Z] = traj->solution[2] - radiant_hat[Z] * length_km;

	ECEF2LatLonAlt( r, &r_lat, &r_LST, &r_hkm );

	r_lon = LST2LongitudeEast( traj->jdt_ref, r_LST );


	//======== Assign to output report structure elements

	if( LLA_position == LLA_BEG )  {
	    traj->rbeg_lat = r_lat;
	    traj->rbeg_lon = r_lon;
	    traj->rbeg_hkm = r_hkm;
	}

	if( LLA_position == LLA_END )  {
	    traj->rend_lat = r_lat;
	    traj->rend_lon = r_lon;
	    traj->rend_hkm = r_hkm;
	}


}

//##################################################################################

//======== Function to run a noisy measurement Monte Carlo set of trials to estimate
//         the error in each parameter. Adds noise to the original measurements and
//         starts the solver with the no noise "solution" information to be used as
//         initialization point for each minimization. Returns the standard 
//         deviations of various critical parameters.
//===================================================================================

void  MonteCarlo_ErrorEstimate( struct trajectory_info *traj )
{
int     kmonte, kcamera, k, noise_flag;
double  ra, dec, vbegin, length_km, kount;
double  radiant_hat[3], r[3], r_lat, r_lon, r_hkm, r_LST;


    //======== Initialize the standard deviation (squared accumulation) to zero

	traj->ra_sigma       = 0.0;
	traj->dec_sigma      = 0.0;
	traj->vbegin_sigma   = 0.0;
	traj->decel1_sigma   = 0.0;
	traj->decel2_sigma   = 0.0;

	traj->rbed_hkm_sigma = 0.0;
	traj->rbeg_lat_sigma = 0.0;
	traj->rbeg_lon_sigma = 0.0;

	traj->rend_hkm_sigma = 0.0;
	traj->rend_lat_sigma = 0.0;
	traj->rend_lon_sigma = 0.0;


    //======== Monte Carlo loop that adds measurement noise to find error estimate

	for( kmonte=0; kmonte<traj->nummonte; kmonte++ )  {

		//-------- Retrieve measurement unit vectors again but this time add noise

	    for( kcamera=0; kcamera<traj->numcameras; kcamera++ )  {

			noise_flag = 1;
		
	        Angles2SiteMeasurements( traj->camera_lat[kcamera], 
									 traj->camera_lon[kcamera], 
									 traj->camera_hkm[kcamera],
									 traj->camera_LST[kcamera],
									 traj->nummeas[kcamera],
									 traj->meastype,
									 traj->meas1[kcamera], 
									 traj->meas2[kcamera], 
									 traj->dtime[kcamera],
									 noise_flag,
									 traj->noise[kcamera],
									 traj->meashat_ECI[kcamera], 
									 traj->rcamera_ECI[kcamera]  );

		}  //... end of camera loop to add noise to the measurements


		//-------- Find noise added trajectory with the no noise solution as the starting guess

	    for( k=0; k<traj->numparams; k++ )  traj->xguess[k] = traj->solution[k];
	    for( k=0; k<traj->numparams; k++ )  traj->xshift[k] = traj->limits[k];

	    if( traj->velmodel <= 0 )  traj->xshift[6] = 0.0;  //... Constrain decel1 due to velocity model    
	    if( traj->velmodel <= 1 )  traj->xshift[7] = 0.0;  //... Constrain decel2 due to velocity model  

	    ParameterRefinementViaPSO( traj, MEAS2MODEL, traj->velmodel, QUICK_FIT );  // --> traj-params[*]


        //-------- Form statistics on radiant, velocity, deceleration estimates

	    ECI2RahDec( &traj->params[3], &ra, &dec );

	    VectorNorm( &traj->params[3], radiant_hat, &vbegin );

		if( ra - traj->ra_radiant > +3.14159265359 )  ra -= 2.0 * 3.14159265359;
		if( ra - traj->ra_radiant < -3.14159265359 )  ra += 2.0 * 3.14159265359;

		traj->ra_sigma   += (ra  - traj->ra_radiant ) * (ra  - traj->ra_radiant );
		traj->dec_sigma  += (dec - traj->dec_radiant) * (dec - traj->dec_radiant);

		traj->vbegin_sigma += (               vbegin - traj->vbegin) * (               vbegin - traj->vbegin);
		traj->decel1_sigma += (fabs(traj->params[6]) - traj->decel1) * (fabs(traj->params[6]) - traj->decel1);
		traj->decel2_sigma += (fabs(traj->params[7]) - traj->decel2) * (fabs(traj->params[7]) - traj->decel2);


		//-------- Form statistics on begin LLA

		length_km = Propagation( POSITION, traj->velmodel, traj->begtime, vbegin, fabs(traj->params[6]), fabs(traj->params[7]) );

	    r[X] = traj->params[0] - radiant_hat[X] * length_km; 
	    r[Y] = traj->params[1] - radiant_hat[Y] * length_km;  
	    r[Z] = traj->params[2] - radiant_hat[Z] * length_km;

	    ECEF2LatLonAlt( r, &r_lat, &r_LST, &r_hkm );

	    r_lon = LST2LongitudeEast( traj->jdt_ref, r_LST );

		if( r_lon - traj->rbeg_lon > +3.14159265359 )  r_lon -= 2.0 * 3.14159265359;
		if( r_lon - traj->rbeg_lon < -3.14159265359 )  r_lon += 2.0 * 3.14159265359;

		traj->rbeg_lon_sigma += (r_lon - traj->rbeg_lon) * (r_lon - traj->rbeg_lon);
		traj->rbeg_lat_sigma += (r_lat - traj->rbeg_lat) * (r_lat - traj->rbeg_lat);

		traj->rbed_hkm_sigma += (r_hkm - traj->rbeg_hkm) * (r_hkm - traj->rbeg_hkm);


		//..... Form statistics on end LLA

		length_km = Propagation( POSITION, traj->velmodel, traj->endtime, vbegin, fabs(traj->params[6]), fabs(traj->params[7]) );

	    r[X] = traj->params[0] - radiant_hat[X] * length_km; 
	    r[Y] = traj->params[1] - radiant_hat[Y] * length_km;  
	    r[Z] = traj->params[2] - radiant_hat[Z] * length_km;

	    ECEF2LatLonAlt( r, &r_lat, &r_LST, &r_hkm );

	    r_lon = LST2LongitudeEast( traj->jdt_ref, r_LST );

		if( r_lon - traj->rend_lon > +3.14159265359 )  r_lon -= 2.0 * 3.14159265359;
		if( r_lon - traj->rend_lon < -3.14159265359 )  r_lon += 2.0 * 3.14159265359;

		traj->rend_lon_sigma += (r_lon - traj->rend_lon) * (r_lon - traj->rend_lon);
		traj->rend_lat_sigma += (r_lat - traj->rend_lat) * (r_lat - traj->rend_lat);

		traj->rend_hkm_sigma += (r_hkm - traj->rend_hkm) * (r_hkm - traj->rend_hkm);


	}  //... end of Monte Carlo loop


	//======== Compute standard deviations

	if( traj->nummonte > 1 ) kount = (double)(traj->nummonte - 1);
	else                     kount = 1.0;

    traj->ra_sigma       = sqrt( traj->ra_sigma       / kount );
    traj->dec_sigma      = sqrt( traj->dec_sigma      / kount );
    traj->vbegin_sigma   = sqrt( traj->vbegin_sigma   / kount );
    traj->decel1_sigma   = sqrt( traj->decel1_sigma   / kount );
    traj->decel2_sigma   = sqrt( traj->decel2_sigma   / kount );

	traj->rbeg_lon_sigma = sqrt( traj->rbeg_lon_sigma / kount );
	traj->rend_lon_sigma = sqrt( traj->rend_lon_sigma / kount );

	traj->rbeg_lat_sigma = sqrt( traj->rbeg_lat_sigma / kount );
	traj->rend_lat_sigma = sqrt( traj->rend_lat_sigma / kount );

	traj->rbed_hkm_sigma = sqrt( traj->rbed_hkm_sigma / kount );
	traj->rend_hkm_sigma = sqrt( traj->rend_hkm_sigma / kount );

}


//##################################################################################

double  VectorDotProduct( double *a, double *b )
{
   double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   return(c);
}

//==============================================================

void  VectorCrossProduct( double *a, double *b, double *c )
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

//==============================================================

void  VectorNorm( double *a, double *b, double *c )
{
   *c = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
   if( *c != 0 ) {
       b[0] = a[0] / *c;
       b[1] = a[1] / *c;
       b[2] = a[2] / *c;
   }
}

//##################################################################################

double  RandomGauss( double sigma, double maxvalue )
{
static short   kthcall = 0;
static double  rangauss1, rangauss2;
double         ran1, ran2, v1, v2, fac, vsq;


     if( kthcall == 0 )  {
		 do {

             do {
                  ran1 = (double)rand() / (double)RAND_MAX;
                  ran2 = (double)rand() / (double)RAND_MAX;
                  v1   = 2.0 * ran1 - 1.0;
                  v2   = 2.0 * ran2 - 1.0;
                  vsq  = v1*v1 + v2*v2;
             }
             while( vsq >= 1.0 );
         
             fac = sqrt( -2.0*log(vsq)/vsq );
         
             rangauss1 = sigma * v1 * fac;
             rangauss2 = sigma * v2 * fac;
		 }
		 while( rangauss1 > maxvalue || rangauss2 > maxvalue );

         kthcall = 1;
         return( rangauss1 );
     }
     else {
         kthcall = 0;
         return( rangauss2 );
     }
     
} 

//##################################################################################

//=======================================================================================
// Conversion from geodetic lat, lon, alt to ECEF assuming WGS84 height in km and radians
//=======================================================================================

void  LatLonAlt2ECEF( double lat, double lon, double alt_km, double *ecef_km  )
{
double  a, e, N;

  a = 6378.137;   //... WGS84 constants
  e = 0.081819190842621;
    
  N = a / sqrt( 1.0 - e*e * sin(lat)*sin(lat) );
  
  ecef_km[0] = (N + alt_km) * cos(lat) * cos(lon);
  ecef_km[1] = (N + alt_km) * cos(lat) * sin(lon);
  ecef_km[2] =  ((1-e*e) * N + alt_km) * sin(lat);
   
}
  
//============================================================================
// Conversion from ECEF to geodetic lat, lon, alt using WGS84 (radians and km)
//============================================================================

void  ECEF2LatLonAlt( double *ecef_km, double *lat, double *lon, double *alt_km )
{
double a, b, e, ep, N, p, theta;


  a = 6378.137;   //... WGS84 constants
  e = 0.081819190842621;
  
  b = sqrt( a*a * (1.0 - e*e) );
  ep = sqrt( (a*a - b*b) / (b*b) );
  
  *lon = atan2( ecef_km[1], ecef_km[0] );
  
  p = sqrt( ecef_km[0] * ecef_km[0]  +  ecef_km[1] * ecef_km[1] );
  theta = atan2( ecef_km[2] * a, p * b );
  
  *lat = atan2( ecef_km[2] + ep*ep*b*sin(theta)*sin(theta)*sin(theta), p - e*e*a*cos(theta)*cos(theta)*cos(theta) );
              
  N = a / sqrt( 1.0 - e*e * sin(*lat)*sin(*lat) );
  
  *alt_km = p / cos(*lat) - N;
  
  if( fabs(ecef_km[0]) < 0.001  &&  fabs(ecef_km[1]) < 0.001 )  *alt_km = fabs(ecef_km[2]) - b;  
   
}

//============================================================================
// Conversion from ECI to Right Ascension and Declination
//============================================================================

void  ECI2RahDec( double *eci, double *rah, double *dec )
{
double  pi, mag_dummy, rad[3];


    pi = 3.14159265359;

	VectorNorm( eci, rad, &mag_dummy );

	*dec = asin( rad[Z] );

    *rah = atan2( rad[Y], rad[X] );

	if( *rah < 0.0 )  *rah += 2.0 * pi;

}

//============================================================================
// Conversion from Right Ascension and Declination to ECI Coordinates
//============================================================================

void  RahDec2ECI( double rah, double dec, double *eci )
{
    eci[X] = cos(dec) * cos(rah);
    eci[Y] = cos(dec) * sin(rah);
    eci[Z] = sin(dec);

}

//============================================================================
// Conversion from Azimuth (east of north) and Elevation (altitude angle) to
//     Right Ascension and Declination. All angles in radians.
//============================================================================

void   AzimuthElevation2RADec( double  azim, 
                               double  elev,
                               double  geodetic_latitude,
                               double  LST,
				               double *RA,
					           double *DEC    )
{
double pi, hourangle, sinlat, coslat;


       pi = 3.14159265359;

       //... Note formulae signs assume azim measured positive east of north

       sinlat = sin( geodetic_latitude );
       coslat = cos( geodetic_latitude );

	   *DEC = asin( sinlat * sin( elev )
                  + coslat * cos( elev ) * cos( azim ) );
                        
       hourangle = atan2( -sin( azim ), 
                           tan( elev ) * coslat - cos( azim ) * sinlat );                
       
       *RA = LST - hourangle;
       
       while( *RA < 0.0      )  *RA += 2.0 * pi;
       while( *RA > 2.0 * pi )  *RA -= 2.0 * pi;

} 

//============================================================================
// Conversion from Right Ascension and Declination to Azimuth (east of north) 
//     and Elevation (altitude angle). All angles in radians.
//============================================================================

void     RADec2AzimuthElevation( double  RA,
                                 double  DEC,
                                 double  geodetic_latitude, 
                                 double  LST,
						         double *Azim,
						         double *Elev  )
{
double pi, hourangle, sinelev, sinlat, coslat;


       pi = 3.14159265359;

       sinlat = sin( geodetic_latitude );
       coslat = cos( geodetic_latitude );

       hourangle = LST - RA;

       while( hourangle < -pi )  hourangle += 2.0 * pi;
       while( hourangle > +pi )  hourangle -= 2.0 * pi;
       
       *Azim = pi + atan2( sin( hourangle ),
                           cos( hourangle ) * sinlat 
                         - tan( DEC ) * coslat );
                                     
       sinelev = sinlat * sin( DEC )
               + coslat * cos( DEC ) * cos( hourangle );
               
       if( sinelev > +1.0 )  sinelev = +1.0;
       if( sinelev < -1.0 )  sinelev = -1.0;

       *Elev = asin( sinelev );
       
} 

//============================================================================
// Computation of Local Sidereal Time - Based on Meeus
//============================================================================

double  LocalSiderealTimeE( double jdt, double longitude_east )
{
double  pi, tim, stg, lst;  // longitude assumed +east and in radians

        //.................. sidereal time at Greenwich

        pi = 3.14159265359;
        
        tim = ( jdt - 2451545.0 ) / 36525.0;
                
        stg = 280.46061837 
            + 360.98564736629 * ( jdt - 2451545.0 )
            + tim * tim * 0.000387933
            - tim * tim * tim / 38710000.0;
                    
        //.................. local sidereal time
            
        lst = stg * pi / 180.0  +  longitude_east;
        
        //.................. set value between 0 and 2pi
        
        while( lst >= 2.0 * pi )  lst -= 2.0 * pi;
        while( lst <       0.0 )  lst += 2.0 * pi;
       
        return( lst ); 
}

//============================================================================
// Conversion from LST to Longitude using Julian date/time
//============================================================================

double  LST2LongitudeEast( double jdt, double LST )
{
double  pi, tim, stg, eastlongitude;  // longitude assumed +east

        //.................. sidereal time at Greenwich
        
        pi = 3.14159265359;
        
        tim = ( jdt - 2451545.0 ) / 36525.0;
                
        stg = 280.46061837 
            + 360.98564736629 * ( jdt - 2451545.0 )
            + tim * tim * 0.000387933
            - tim * tim * tim / 38710000.0;
                    
        //.................. local sidereal time
            
		eastlongitude = LST - stg * pi / 180.0;
        
        //.................. set value between 0 and 360 degrees
        
        while( eastlongitude >= +pi )  eastlongitude -= 2.0 * pi;
        while( eastlongitude <  -pi )  eastlongitude += 2.0 * pi;
       
        return( eastlongitude ); 

}


//##################################################################################
//
//======== Function to initially allocate memory for 1D and 2D arrays based on the 
//         maximum number of camera measurement sets expected. This function should
//         be called only once at program start-up.
//==================================================================================

void    InitTrajectoryStructure( int maxcameras, struct trajectory_info *traj )
{
int  kcamera, maxparams;


    maxparams = 9 + maxcameras;   //... Multi-parameter fitting 
	                              //    position vector, radiant unit vector, velocity, decel1, decel2, dtime(#cameras) 

	traj->maxcameras = maxcameras;
	traj->numcameras = 0;

    traj->camera_lat     =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->camera_lon     =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->camera_hkm     =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->camera_LST     =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->tref_offsets   =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->nummeas        =      (int*) malloc( maxcameras * sizeof(    int  ) );
    traj->malloced       =      (int*) malloc( maxcameras * sizeof(    int  ) );

    traj->params         =   (double*) malloc( maxparams  * sizeof( double  ) );
    traj->solution       =   (double*) malloc( maxparams  * sizeof( double  ) );
    traj->limits         =   (double*) malloc( maxparams  * sizeof( double  ) );
    traj->xguess         =   (double*) malloc( maxcameras * sizeof( double  ) );
    traj->xshift         =   (double*) malloc( maxcameras * sizeof( double  ) );

    traj->meas1          =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->meas2          =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->dtime          =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->noise          =  (double**) malloc( maxcameras * sizeof( double* ) );

    traj->meas_lat       =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->meas_lon       =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->meas_hkm       =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->meas_range     =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->meas_vel       =  (double**) malloc( maxcameras * sizeof( double* ) );

	traj->model_lat      =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_lon      =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_hkm      =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_range    =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_vel      =  (double**) malloc( maxcameras * sizeof( double* ) );

	traj->model_fit1     =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_fit2     =  (double**) malloc( maxcameras * sizeof( double* ) );
    traj->model_time     =  (double**) malloc( maxcameras * sizeof( double* ) );

    traj->meashat_ECI    = (double***) malloc( maxcameras * sizeof( double**) );
    traj->rcamera_ECI    =  (double**) malloc( maxcameras * sizeof( double* ) );



	if( traj->camera_lat     == NULL  ||
		traj->camera_lon     == NULL  ||
		traj->camera_hkm     == NULL  ||
		traj->camera_LST     == NULL  ||
		traj->tref_offsets   == NULL  ||
		traj->nummeas        == NULL  ||
		traj->params         == NULL  ||
		traj->solution       == NULL  ||
		traj->limits         == NULL  ||
		traj->xguess         == NULL  ||
		traj->xshift         == NULL  ||
		traj->malloced       == NULL  ||
		traj->meas1          == NULL  ||
		traj->meas2          == NULL  ||
		traj->dtime          == NULL  ||
		traj->noise          == NULL  ||
		traj->meas_lat       == NULL  ||
		traj->meas_lon       == NULL  ||
		traj->meas_hkm       == NULL  ||
		traj->meas_range     == NULL  ||
		traj->meas_vel       == NULL  ||
		traj->model_lat      == NULL  ||
		traj->model_lon      == NULL  ||
		traj->model_hkm      == NULL  ||
		traj->model_range    == NULL  ||
		traj->model_vel      == NULL  ||
		traj->model_fit1     == NULL  ||
		traj->model_fit2     == NULL  ||
		traj->model_time     == NULL  ||
		traj->rcamera_ECI    == NULL  ||
		traj->meashat_ECI    == NULL      )  {

		printf("ERROR--> Memory not allocated for vectors and arrays in InitTrajectoryStructure\n");
		exit(1);

	}


	//======== Set the memory allocated flag to "freed"

	for( kcamera=0; kcamera<traj->maxcameras; kcamera++ )  traj->malloced[kcamera] = 0; 


	//======== The second dimension of each 2D array will be allocated as the measurments are ingested
	//              except for "rcamera_ECI" where the second dimension size is known for the three 
	//              components XYZ.

	for( kcamera=0; kcamera<traj->maxcameras; kcamera++ )  { //... Allocate up to MAX cameras

		traj->rcamera_ECI[kcamera] = (double*) malloc( 3 * sizeof(double) ); 

		if( traj->rcamera_ECI[kcamera] == NULL )  {
			printf("ERROR--> Memory not allocated for rcamera_ECI in InitTrajectoryStructure\n");
		    exit(1);
		}

	} //... end of "kcamera" loop


}


//##################################################################################
//
//======== Function to free up ALL allocated memory for 1D, 2D and 3D arrays. This  
//         should be called only once at program completion
//==================================================================================

void    FreeTrajectoryStructure( struct trajectory_info *traj )
{
int  kcamera;

	//... First free up the column dimensions of the 2D arrays

	for( kcamera=0; kcamera<traj->maxcameras; kcamera++ )  free( traj->rcamera_ECI[kcamera] );

	ResetTrajectoryStructure( 0.0, 0.0, 0, 0, 0, 0, traj ); 

    free( traj->meas1          );     //... Now free up the row dimension of the 2D arrays
    free( traj->meas2          );
    free( traj->dtime          );
    free( traj->noise          );

    free( traj->meas_lat       );
    free( traj->meas_lon       );
    free( traj->meas_hkm       );
    free( traj->meas_range     );
    free( traj->meas_vel       );

	free( traj->model_lat      );
    free( traj->model_lon      );
    free( traj->model_hkm      );
    free( traj->model_range    );
    free( traj->model_vel      );

	free( traj->model_fit1     );
    free( traj->model_fit2     );
    free( traj->model_time     );

	free( traj->rcamera_ECI    );
    free( traj->meashat_ECI    );

	free( traj->malloced       );


	free( traj->camera_lat     );     //... Free up the vectors
    free( traj->camera_lon     );
    free( traj->camera_hkm     );
    free( traj->camera_LST     );
    free( traj->tref_offsets   );
    free( traj->nummeas        );
    free( traj->params         );
    free( traj->solution       );
    free( traj->limits         );
    free( traj->xguess         );
    free( traj->xshift         );


}

//##################################################################################
//
//======== Function to reset the trajectory structure by freeing up the column 
//         dimension (measurement dimension) of the 2D arrays.
//==================================================================================

void    ResetTrajectoryStructure( double jdt_ref, double max_toffset, 
	                              int velmodel, int nummonte, int meastype, int verbose, 
                                  struct trajectory_info *traj )
{
int  kcamera, kmeas;


    //======== Set up some initial parameters for this trajectory solution

	traj->jdt_ref     = jdt_ref;
	traj->max_toffset = max_toffset;
	traj->velmodel    = velmodel;
	traj->nummonte    = nummonte;
	traj->meastype    = meastype;
	traj->verbose     = verbose;

	traj->numcameras = 0;

    for( kcamera=0; kcamera<traj->maxcameras; kcamera++ )  traj->tref_offsets[kcamera] = 0.0;


	//======== Free memory through all possible camera arrays of 2D or 3D dimensions

    for( kcamera=0; kcamera<traj->maxcameras; kcamera++ )  {

		if( traj->malloced[kcamera] == 1 )  {

            free( traj->meas1[kcamera]          );
            free( traj->meas2[kcamera]          );
            free( traj->dtime[kcamera]          );
            free( traj->noise[kcamera]          );

            free( traj->meas_lat[kcamera]       );
            free( traj->meas_lon[kcamera]       );
            free( traj->meas_hkm[kcamera]       );
            free( traj->meas_range[kcamera]     );
            free( traj->meas_vel[kcamera]       );

	        free( traj->model_lat[kcamera]      );
            free( traj->model_lon[kcamera]      );
            free( traj->model_hkm[kcamera]      );
            free( traj->model_range[kcamera]    );
            free( traj->model_vel[kcamera]      );

	        free( traj->model_fit1[kcamera]     );
            free( traj->model_fit2[kcamera]     );
            free( traj->model_time[kcamera]     );

		    //... free first the XYZT 3rd dimension of meashat_ECI then the 2nd measurement dimension

		    for( kmeas=0; kmeas<traj->nummeas[kcamera]; kmeas++ )  free( traj->meashat_ECI[kcamera][kmeas] );

			free( traj->meashat_ECI[kcamera] );

			traj->malloced[kcamera] = 0;  //... memory freed

		} //... end of IF preiously malloc'ed

	} //... end of "kcamera" loop


}

//##################################################################################
//
//======== Function to infill the trajectory structure with a new set of measurements 
//         from a single camera. Multiple calls to this function builds up the full 
//         measurement set for multiple sites and multiple cameras.
//==================================================================================

void    InfillTrajectoryStructure( int nummeas, double *meas1, double *meas2, double *dtime, double *noise,
	                               double site_latitude, double site_longitude, double site_height,
                                   struct trajectory_info *traj )
{
int  kcamera, kmeas;


    kcamera = traj->numcameras;

	AllocateTrajectoryMemory4Infill( kcamera, nummeas, traj );

	traj->nummeas[kcamera] = nummeas;


	//======== Assign the measurements to the working arrays

	for( kmeas=0; kmeas<nummeas; kmeas++ )  {

		traj->meas1[kcamera][kmeas] = meas1[kmeas];
		traj->meas2[kcamera][kmeas] = meas2[kmeas];
		traj->dtime[kcamera][kmeas] = dtime[kmeas];
		traj->noise[kcamera][kmeas] = noise[kmeas];

	}

	traj->camera_lat[kcamera]   = site_latitude;
	traj->camera_lon[kcamera]   = site_longitude;
	traj->camera_hkm[kcamera]   = site_height;
	traj->camera_LST[kcamera]   = LocalSiderealTimeE( traj->jdt_ref, site_longitude );

	traj->tref_offsets[kcamera] = 0.0;
	

    traj->numcameras = kcamera + 1;  //... Increment the active camera counter for next camera's measurments


}


//##################################################################################
//
//======== Function to allocate memory for the column dimension of the 2D arrays
//==================================================================================

void    AllocateTrajectoryMemory4Infill( int kcamera, int nummeas, struct trajectory_info *traj )
{
int  kmeas;

	//======== Check to make sure memory was previously freed for each array column

	if( traj->malloced[kcamera] == 1 )  {

		printf("ERROR--> You must call ResetTrajectoryStructure prior to the first InfillTrajectoryStructure call\n");
		Sleep(10000);
	    exit(1);

    }


	//======== Allocate memory through all the working arrays of this camera index

    traj->meas1[kcamera]          = (double*) malloc( nummeas * sizeof(double) );
    traj->meas2[kcamera]          = (double*) malloc( nummeas * sizeof(double) );
    traj->dtime[kcamera]          = (double*) malloc( nummeas * sizeof(double) );
    traj->noise[kcamera]          = (double*) malloc( nummeas * sizeof(double) );

    traj->meas_lat[kcamera]       = (double*) malloc( nummeas * sizeof(double) );
    traj->meas_lon[kcamera]       = (double*) malloc( nummeas * sizeof(double) );
    traj->meas_hkm[kcamera]       = (double*) malloc( nummeas * sizeof(double) );
    traj->meas_range[kcamera]     = (double*) malloc( nummeas * sizeof(double) );
    traj->meas_vel[kcamera]       = (double*) malloc( nummeas * sizeof(double) );
		
    traj->model_lat[kcamera]      = (double*) malloc( nummeas * sizeof(double) );
    traj->model_lon[kcamera]      = (double*) malloc( nummeas * sizeof(double) );
    traj->model_hkm[kcamera]      = (double*) malloc( nummeas * sizeof(double) );
    traj->model_range[kcamera]    = (double*) malloc( nummeas * sizeof(double) );
    traj->model_vel[kcamera]      = (double*) malloc( nummeas * sizeof(double) );

	traj->model_fit1[kcamera]     = (double*) malloc( nummeas * sizeof(double) );
    traj->model_fit2[kcamera]     = (double*) malloc( nummeas * sizeof(double) );
    traj->model_time[kcamera]     = (double*) malloc( nummeas * sizeof(double) );

    traj->meashat_ECI[kcamera]        = (double**) malloc( nummeas * sizeof(double*) );


	//======== Check the 2D memory was allocated

	if( traj->meas1[kcamera]          == NULL  ||
		traj->meas2[kcamera]          == NULL  ||
		traj->dtime[kcamera]          == NULL  ||
		traj->noise[kcamera]          == NULL  ||
		traj->meas_lat[kcamera]       == NULL  ||
		traj->meas_lon[kcamera]       == NULL  ||
		traj->meas_hkm[kcamera]       == NULL  ||
		traj->meas_range[kcamera]     == NULL  ||
		traj->meas_vel[kcamera]       == NULL  ||
		traj->model_lat[kcamera]      == NULL  ||
		traj->model_lon[kcamera]      == NULL  ||
		traj->model_hkm[kcamera]      == NULL  ||
		traj->model_range[kcamera]    == NULL  ||
		traj->model_vel[kcamera]      == NULL  ||
		traj->model_fit1[kcamera]     == NULL  ||
		traj->model_fit2[kcamera]     == NULL  ||
		traj->model_time[kcamera]     == NULL  ||
		traj->meashat_ECI[kcamera]        == NULL      )  {

		printf("ERROR--> Memory not allocated for 2D array columns in AllocateTrajectoryMemory\n");
		exit(1);

	}

	//======== Allocate the 3rd dimension for components XYZT of the measurement unit vectors
	
	for( kmeas=0; kmeas<nummeas; kmeas++ )  {

		traj->meashat_ECI[kcamera][kmeas] = (double*) malloc( 4 * sizeof(double) ); 

		if( traj->meashat_ECI[kcamera][kmeas] == NULL )  {
			printf("ERROR--> Memory not allocated for meashat_ECI in AllocateTrajectoryMemory\n");
		    exit(1);
		}

	}


	//======== Set the memory allocated flag for this camera

	traj->malloced[kcamera] = 1;


}

//################################################################################