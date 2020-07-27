from django.contrib import admin
from django.urls import path, include





urlpatterns = [
    path('admin/',      admin.site.urls),
    path('neo4j/', include('neo4j_connector.urls','neo4j_connector') )
]
