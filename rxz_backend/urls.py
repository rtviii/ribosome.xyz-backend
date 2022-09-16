from django.contrib import admin
from django.urls import path, include
from rest_framework.schemas import get_schema_view

from django.contrib import admin
from django.contrib.staticfiles.views import serve
from django.views.generic import RedirectView


from rest_framework import permissions
from drf_yasg.views import get_schema_view
from drf_yasg import openapi
from django.shortcuts import redirect
def view_404(request, exception=None):
    return redirect('/docs') # or redirect('name-of-index-url')
handler404  = 'rxz_backend.urls.view_404'
schema_view = get_schema_view(
   openapi.Info(
      title="Ribosome XYZ",
      default_version='v1',
      description="""The API is freely available for public requests, but the development is still ongoing. For any sizeable requests or long-term integrations we do recommend getting in touch at kdd@math.ubc.ca. At this stage we are happy to provide the data statically and consider implementing custom endpoints.
      Cite us at: [PULBICATION IN REVIEW] """,
      contact=openapi.Contact(email="kdd@math.ubc.ca"),
      license=openapi.License(name="BSD License"),
   ),
   public=True,
   permission_classes=[permissions.AllowAny],
)

urlpatterns = [

   #  path('admin/'        , admin  .site.urls                                ),
    path('neo4j/'        , include('neo4j_connector.urls','neo4j_connector')),
    path('static_files/' , include('static_files.urls', 'static_files')     ),
    path('utils/'        , include('utils.urls', 'utils')                   ),
    
   #  path('/docs', schema_view.without_ui(cache_timeout=0), name='schema-json'),

   #  path('docs/', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    path('docs/', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
   #  path('redoc/', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),
]
