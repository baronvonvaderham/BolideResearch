from django import forms
from models import Node


class NodeForm(forms.Form):
    name = forms.CharField(required=True)
    latitude = forms.FloatField(required=True)
    longitude = forms.FloatField(required=True)


class EventForm(forms.Form):
    event_date = forms.DateField(required=True)
    node = forms.Select(Node)
    data = forms.FilePathField(path=None, required=True)
