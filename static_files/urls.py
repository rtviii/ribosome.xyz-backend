from django.urls import path
from .views import * 


urlpatterns = [
    path('get_ligand_nbhd/',            get_ligand_nbhd) ,
    path('download_ligand_nbhd/',       download_ligand_nbhd) ,
    path('cif_chain/',                  cif_chain),
    path('tunnel/',                     tunnel),
]

app_name = 'static_files'