from django.urls import path
from .views import get_struct, custom_cypher, anything, get_homologs
from .pdb_connector import FileDownloadListAPIView


urlpatterns = [
    path('get_struct/', get_struct),
    path('get_homologs/', get_homologs),
    path('cypher/', custom_cypher),
    path('anything/', anything),
    path('get_pdbsubchain/', FileDownloadListAPIView.as_view())
]

app_name = 'neo4j_connector'