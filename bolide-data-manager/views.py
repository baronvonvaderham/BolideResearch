from models import Node, EventData, Event
from django.views.generic import FormView, ListView,\
    DetailView, DeleteView, UpdateView
from forms import NodeForm, EventForm


class CreateNodeView(FormView):
    """
        Generic FormView for the creation of new Nodes
    """
    template_name = 'node_create.html'
    form_class = NodeForm
    success_url = '/node_list/'

    def form_valid(self, form):
        new_node = Node.objects.create(
            name=form.name,
            latitude=form.latitude,
            longitude=form.longitude,
        )
        new_node.save()
        return super(CreateNodeView, self).form_valid(form)


class EditNodeView(UpdateView):
    """
        Generic UpdateView for editing a Node
    """
    form_class = NodeForm
    template_name = 'node_edit.html'
    success_url = '/node_list/'

    def get_object(self, **kwargs):
        return Node.objects.get(name=self.kwargs['node_name'])


class DeleteNodeView(DeleteView):
    """
        Generic DeleteView for deleting a Node
    """
    model = Node
    template_name = 'confirm_node_delete.html'
    success_url = '/node_list/'

    def get_object(self, **kwargs):
        return Node.objects.get(name=self.kwargs['node_name'])


class NodeListView(ListView):
    """
        Generic ListView to display the list of all existing Nodes
    """
    model = Node
    template_name = 'node_list.html'
    context_object_name = 'node_list'

    def get_queryset(self):
        return Node.objects.all().order_by('name')


class NodeDetailView(DetailView):
    """
        Generic DetailView to display information about a node, access
        its edit function (link in template), and see a list of all existing
        Events associated with this Node in reverse chronological order.

        @todo - Add thumbnail display of images?
    """
    model = Node
    template_name = 'node_detail.html'

    def get_object(self, queryset=None):
        return Node.objects.get(name=self.kwargs['node_name'])

    def get_context_data(self, **kwargs):
        context = super(NodeDetailView, self).get_context_data(**kwargs)
        context['event_list'] = Event.objects.filter(self.get_object()).\
            order_by('-event_date')
        return context


class CreateEventView(FormView):
    """
        Generic FormView for the creation of new Events
    """
    template_name = 'event_create.html'
    form_class = EventForm
    success_url = '/event_list/'

    def form_valid(self, form):
        new_event = Event.objects.create(
            event_date=form.event_date,
            node=form.node,
        )
        new_event.save()
        """
            @todo - Insert processing of csv to create EventData objects here!
        """
        return super(CreateEventView, self).form_valid(form)


class EditEventView(UpdateView):
    """
        Generic UpdateView for editing an Event
    """
    form_class = EventForm
    template_name = 'event_edit.html'
    success_url = '/event_list/'

    def get_object(self, **kwargs):
        return Node.objects.get(id=self.kwargs['event_id'])


class DeleteEventView(DeleteView):
    """
        Generic DeleteView for deleting an Event
        and its associated EventData objects
    """
    model = Event
    template_name = 'confirm_event_delete.html'
    success_url = '/event_list/'

    def get_object(self, **kwargs):
        return Node.objects.get(id=self.kwargs['event_id'])

class EventListView(ListView):
    """
        Generic ListView to display the list of all existing Events
    """
    model = Event
    template_name = 'event_list.html'
    context_object_name = 'event_list'

    def get_queryset(self):
        return Event.objects.all().order_by('event_date')
