from django.urls import path
from .views import * 


urlpatterns = [
    path('get_ligand_nbhd/'     , get_ligand_nbhd     ) ,
    path('download_ligand_nbhd/', download_ligand_nbhd) ,
    path('cif_chain/'           , cif_chain           ) ,
    path('ligand_prediction/'   , ligand_prediction   ) ,
    path('download_structure/'  , download_structure  ) ,
    path('cif_chain_by_class/'  , cif_chain_by_class  ) ,
    path('tunnel/'              , tunnel              ) ,
    path('align_3d/'            , align_3d            ) ,
    path('downloadArchive/'     , downloadArchive     ) ,
]


app_name = 'static_files'