# Bolide Research Repository

Contains the code for my dissertation research and the prototype for the
data management system for the SkySentinel project (goskysentinel.com).

## Structure

This has been constructed as a standalone Django web app to make use of the
django database engine and model system. This has an advantage over the use
of pure python with sqlalchemy in that it allowed me to create a management
portal through which additional students can import the data from nodes after
events are detected. Eventually, this lends itself to creating a fully
automated data pipeline whereby the nodes can automatically import their
event data into the models.

## Database

This uses a db.sqlite3 instance locally for now. Eventually, this will be
deployed to something like Heroku or AWS and make use of a native RDS instance
of some kind (bonus: allows remote data importing and cloud-based backups).

