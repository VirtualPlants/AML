
#scons parameters file
#use this file to pass custom parameter to SConstruct script

import sys,os

if('win32' in sys.platform):
    compiler='mingw'

elif('posix' in os.name):
    QTDIR="/usr"
    QT_VERSION=4
