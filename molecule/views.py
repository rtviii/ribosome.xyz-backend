from django.shortcuts import render
from rest_framework import viewsets          # add this
from .serializers import     MoleculeSerializer
from .models import Molecule     

# Create your views here.




class MoleculeView(viewsets.ModelViewSet):
    serializer_class = MoleculeSerializer
    queryset = Molecule.objects.all()
