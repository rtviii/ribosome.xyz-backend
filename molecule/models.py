from django.db import models

# Create your models here.



class Molecule(models.Model):
    pdbid = models.CharField(max_length=120)
    description = models.TextField()
    reviewed = models.BooleanField(default=False)


    def _str_(self):
        return self.pdbid