from django.contrib import admin
from .models import Molecule

# Register your models here.


class MoleculeAdmin(admin.ModelAdmin):
    list_display=('pdbid','description','reviewed')
    
admin.site.register(Molecule,MoleculeAdmin)
