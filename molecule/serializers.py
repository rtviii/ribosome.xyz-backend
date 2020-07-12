from rest_framework import serializers
from .models import Molecule

class MoleculeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Molecule
        fields = ('pdbid', 'description', 'reviewed')