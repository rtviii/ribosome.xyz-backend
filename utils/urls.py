from django.urls import path
from .views import * 


urlpatterns = [
    path('number_of_structures/', number_of_structures)
]
app_name = 'utils'