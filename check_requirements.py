import pkg_resources
import sys

# check Python version
print('Checking Python version')
if not (sys.version_info.major==3 and sys.version_info.minor >= 5):
    print('ERROR: Python version >=3.5 required, current Python version is %i.%i.%i' % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
    sys.exit(1)

# if required package cannot be found, then a DistributionNotFound or VersionConflict error
print('Checking for required Python packages')
pkg_resources.require(['numpy>=1.15', 'scipy>=0.19'])
try:
    pkg_resources.require(['numpy>=1.15', 'scipy>=0.19'])
except pkg_resources.VersionConflict:
    sys.exit(1)
