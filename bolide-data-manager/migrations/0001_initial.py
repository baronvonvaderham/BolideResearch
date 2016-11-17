# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Event',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('event_date', models.DateField()),
            ],
            options={
                'verbose_name': 'Bolide Event',
                'verbose_name_plural': 'Bolide Events',
            },
        ),
        migrations.CreateModel(
            name='EventData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('date_time', models.DateTimeField()),
                ('pixels_above_threshold', models.IntegerField()),
                ('sum_above_threshold', models.IntegerField()),
                ('centroid_x', models.FloatField()),
                ('centroid_y', models.FloatField()),
                ('computed_az', models.FloatField()),
                ('computed_el', models.FloatField()),
                ('event', models.ForeignKey(to='bolide-data-manager.Event')),
            ],
            options={
                'verbose_name': 'Data Point',
                'verbose_name_plural': 'Data Points',
            },
        ),
        migrations.CreateModel(
            name='Node',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.TextField(unique=True, max_length=255)),
                ('latitude', models.FloatField()),
                ('longitude', models.FloatField()),
            ],
            options={
                'verbose_name': 'Camera Node',
                'verbose_name_plural': 'Camera Nodes',
            },
        ),
        migrations.AddField(
            model_name='event',
            name='node',
            field=models.ForeignKey(to='bolide-data-manager.Node'),
        ),
    ]
