import os,sys
import json
from dotenv import load_dotenv


load_dotenv(dotenv_path='/home/rxz/dev/ribosome.xyz-backend/rxz_backend/.env')
print("yup", os.getenv('SECRET_KEY'))
SECRET_KEY = os.getenv('SECRET_KEY')
DEBUG      = os.getenv('DEBUG')

from pathlib import Path
BASE_DIR               =  Path(__file__).resolve().parent.parent

ALLOWED_HOSTS          =  ["*"]
CORS_ORIGIN_ALLOW_ALL  =  True

INSTALLED_APPS = [
    'corsheaders',
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'rest_framework',
    'neo4j_connector',
    'static_files'
]

MIDDLEWARE = [
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'rxz_backend.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'rxz_backend.wsgi.application'


DATABASES = {
}
AUTH_PASSWORD_VALIDATORS = [
]


LANGUAGE_CODE  =  'en-us'
TIME_ZONE      =  'UTC'
USE_I18N       =  True
USE_L10N       =  True
USE_TZ         =  True


STATIC_ROOT  =  os.path.join(BASE_DIR, 'static')
STATIC_URL   =  '/static/'
