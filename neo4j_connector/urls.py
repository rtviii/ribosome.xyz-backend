from django.urls import path
from .views import test_endpoint,get_struct, get_all_outgoing_struct, custom_cypher


urlpatterns = [
    path('test/', test_endpoint),
    path('get_struct/', get_struct),
    path('get_all_struct/', get_all_outgoing_struct),
    path('cypher/', custom_cypher)
]

app_name = 'neo4j_connector'