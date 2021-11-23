import pkg_resources
import sys

# check Python version
print('Checking Python version')
if not (sys.version_info.major==3 and sys.version_info.minor >= 5):
    sys.exit('ERROR: Python version >=3.5 required, current Python version is %i.%i.%i' % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro))

# if required package cannot be found, then a DistributionNotFound or VersionConflict error
print('Checking for required Python packages')
pkg_resources.require(['numpy>=1.18', 'scipy>=1.4'])
