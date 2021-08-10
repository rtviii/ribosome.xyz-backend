import os,sys
import json

with open('./rxz_backend/.env.json', 'rb') as envf:env=json.load(envf)

SECRET_KEY = env[ 'SECRET_KEY' ]
DEBUG      = env[ 'DEBUG' ]



os.environ["SECRET_KEY"      ]           = env[ 'SECRET_KEY'      ]
os.environ["DEBUG"           ]           = env[ 'DEBUG'           ]
os.environ["NEO4J_USER"      ]           = env[ 'NEO4J_USER'      ]
os.environ["NEO4J_URI"       ]           = env[ 'NEO4J_URI'       ]
os.environ["NEO4J_CURRENTDB" ]           = env[ 'NEO4J_CURRENTDB' ]
os.environ["NEO4J_PASSWORD"  ]           = env[ 'NEO4J_PASSWORD'  ]
os.environ["PROTEIN_ALIGNMENT_SCRIPT"  ] = env[ 'PROTEIN_ALIGNMENT_SCRIPT'  ]
os.environ["TEMP_CHAIN"  ]               = env[ 'TEMP_CHAIN'  ]
os.environ['STATIC_ROOT']                = env['STATIC_ROOT']
STATIC_ROOT  =  env['STATIC_ROOT']

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


# STATIC_ROOT  =  os.path.join(BASE_DIR, 'static')
STATIC_URL   =  '/static/'
