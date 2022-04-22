from django.urls import path
from .views import *


urlpatterns = [
      path('get_struct/'              , get_struct              ),
      path('get_full_struct/'         , get_full_structure      ),
    # path('get_homologs/'            , get_homologs            ) ,
      path('cypher/'                  , custom_cypher           ) ,
      path('anything/'                , anything                ) ,
      path('list_nom_classes/'        , list_nom_classes        ) ,
      path('get_all_structs/'         , get_all_structs         ) ,
      path('gmo_nom_class/'           , gmo_nom_class           ) ,
      path('get_banclasses_metadata/' , get_banclasses_metadata ) ,
      path('nomclass_visualize/'      , nomclass_visualize      ) ,
      path('banclass_annotation/'     , banclass_annotation     ) ,
      path('get_rna_class/'           , get_rna_class           ) ,
      path('get_individual_ligand/'   , get_individual_ligand   ) ,
      path('get_rnas_by_struct/'      , get_rnas_by_struct      ) ,
      path('get_ligands_by_struct/'   , get_ligands_by_struct   ) ,
      path('get_all_ligands/'         , get_all_ligands         ) ,
      path('get_all_ligandlike/'      , get_all_ligandlike      ) ,
      path('match_structs/'           , match_structs           ) ,
      path('proteins_number/'         , proteins_number         ),
      path('tax_ids/'                 , tax_ids                 ),
      path('nomenclature/'            , nomenclature            )
]

app_name = 'neo4j_connector'