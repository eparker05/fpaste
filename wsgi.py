activate_this = '/home/wwwadmin/fastadb/python27/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import sys
sys.path.insert(0, '/home/wwwadmin/fastadb')
from fpaste import app as application

if __name__ == '__main__':
    application.run()
