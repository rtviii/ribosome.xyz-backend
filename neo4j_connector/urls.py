from django.urls import path
from .views import *


urlpatterns = [

    path('get_struct/',                 get_struct),
    path('get_homologs/',               get_homologs),
    path('cypher/',                     custom_cypher),
    path('anything/',                   anything),

    path('list_nom_classes/',           list_nom_classes),

    path('get_all_structs/',            get_all_structs) ,
    path('gmo_nom_class/',              gmo_nom_class),
    path('get_all_ligands/',            get_all_ligands),
    path('get_all_rnas/',               get_all_rnas),
    path('get_individual_ligand/',               get_individual_ligand),
    path('get_rnas_by_struct/' ,        get_rnas_by_struct) ,
    path('get_ligands_by_struct/',      get_ligands_by_struct),
    path('match_structs/',              match_structs)          ,
    path('get_surface_ratios/',         get_surface_ratios)          ,
    path('TEMP_classification_sample/',         TEMP_classification_sample)          ,

]

app_name = 'neo4j_connector'