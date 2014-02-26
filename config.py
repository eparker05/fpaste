from os import path

basedir = path.abspath(path.dirname(__file__))

CSRF_ENABLED = True
SECRET_KEY = "twentythousandorangeponies"

OPENID_PROVIDERS = [{'name':'Google', 'url':'https://www.google.com/accounts/o8/id'},
                    {'name':'Yahoo', 'url':'https://me.yahoo.com'},
                    {'name':'AOL', 'url':'http://openid.aol.com/<username>'},
                    {'name':'Flickr', 'url':'http://www.flickr.com/<username>'},
                    {'name':'MyOpenID', 'url':'https://www.myopenid.com'}
                   ]
                   
SQLALCHEMY_DATABASE_URI = 'sqlite:///' + path.join(basedir, 'app.db')
SQLALCHEMY_MIGRATE_REPO = path.join(basedir, 'db_repository')
