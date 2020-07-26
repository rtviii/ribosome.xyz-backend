from django.contrib import admin
from django.urls import path, include
from molecule import views





urlpatterns = [
    path('admin/', admin.site.urls),
    path('test/', views.test_endp),
    path('get_struct/', views.get_struct    )
]
