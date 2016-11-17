from django.db import models


class Node(models.Model):
    name = models.TextField(max_length=255, unique=True)
    latitude = models.FloatField()
    longitude = models.FloatField()

    class Meta:
        verbose_name = "Camera Node"
        verbose_name_plural = "Camera Nodes"

    def __unicode__(self):
        return self.name


class Event(models.Model):
    event_date = models.DateField()
    node = models.ForeignKey(Node)

    class Meta:
        verbose_name = "Bolide Event"
        verbose_name_plural = "Bolide Events"

    def __unicode__(self):
        name = self.node.name + ': ' + str(self.event_date)
        return name


class EventData(models.Model):
    event = models.ForeignKey(Event)
    date_time = models.DateTimeField()
    pixels_above_threshold = models.IntegerField()
    sum_above_threshold = models.IntegerField()
    centroid_x = models.FloatField()
    centroid_y = models.FloatField()
    computed_az = models.FloatField()
    computed_el = models.FloatField()

    class Meta:
        verbose_name = "Data Point"
        verbose_name_plural = "Data Points"
