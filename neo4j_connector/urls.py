from django.urls import path
from .views import test_endpoint,get_structs, get_all_outgoing_struct 


urlpatterns = [
    path('test/', test_endpoint),
    path('get_structs/', get_structs),
    path('get_all_struct/', get_all_outgoing_struct),
]

app_name = 'neo4j_connector'