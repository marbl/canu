#!/usr/local/packages/python-2.7.3/bin/python2.7

import sys
import os
import getopt
from distutils import sysconfig

print sysconfig.get_python_inc()

#        flags = ['-I' + ,
#                 '-I' + sysconfig.get_python_inc(plat_specific=True)]
