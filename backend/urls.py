from django.contrib import admin
from django.urls import path, include
from rest_framework import routers
from molecule import views



router = routers.DefaultRouter()
router.register(r'molecules',views.MoleculeView, 'molecule')

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/', include(router.urls)) 
]
