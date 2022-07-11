import os,sys
from dotenv import dotenv_values, load_dotenv, find_dotenv

PROJECT_PATH = os.path.abspath(os.path.dirname(__name__))
DOTENV_PATH  = os.path.join(str( PROJECT_PATH ), 'rxz_backend','.env')
load_dotenv(DOTENV_PATH)

STATIC_ROOT              = os.getenv( "STATIC_ROOT" )
DEBUG                    = os.getenv( "DEBUG"       )
SECRET_KEY               = os.getenv('SECRET_KEY')
# DEBUG                    = os.getenv( 'DEBUG'                    )
# NEO4J_USER               = os.getenv( 'NEO4J_USER'               )
# NEO4J_URI                = os.getenv( 'NEO4J_URI'                )
# NEO4J_CURRENTDB          = os.getenv( 'NEO4J_CURRENTDB'          )
# NEO4J_PASSWORD           = os.getenv( 'NEO4J_PASSWORD'           )
# PROTEIN_ALIGNMENT_SCRIPT = os.getenv( 'PROTEIN_ALIGNMENT_SCRIPT' )
# LIGAND_PREDICTION_SCRIPT = os.getenv( 'LIGAND_PREDICTION_SCRIPT' )
# TEMP_CHAIN               = os.getenv( 'TEMP_CHAIN'               )
# STATIC_ROOT              = os.getenv( 'STATIC_ROOT'              )

from pathlib import Path
BASE_DIR               =  Path(__file__).resolve().parent.parent
# print("BASE DIR IS ", BASE_DIR)
sys.path.append(str(BASE_DIR))
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
    'static_files',
    'utils'
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

STATIC_URL   =  '/static/'
