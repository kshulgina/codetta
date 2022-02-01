import pkg_resources
import sys

# check Python version
print('Checking Python version')
if not (sys.version_info.major==3 and sys.version_info.minor >= 5):
    print('ERROR: Python version >=3.5 required, current Python version is %i.%i.%i' % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
    sys.exit(1)

# if required package cannot be found, then a DistributionNotFound or VersionConflict error
print('Checking for required Python packages')

try:
    pkg_resources.require(['scipy>=0.19'])
except pkg_resources.DistributionNotFound:
    print("ERROR: scipy>=0.19 required, none found.\nYou can install scipy using \
        \n \tconda install 'scipy>=0.19' \n or \n\tpip install 'scipy>=0.19'")
    sys.exit(1)
except pkg_resources.VersionConflict:
    print("ERROR: scipy>=0.19 required, scipy %s found.\nYou can update scipy using \
        \n \tconda install 'scipy>=0.19'\n or \n\tpip install -U 'scipy>=0.19'" % pkg_resources.get_distribution("scipy").version )
    sys.exit(1)

try:
    pkg_resources.require(['numpy>=1.15'])
except pkg_resources.DistributionNotFound:
    print("ERROR: numpy>=1.15 required, none found.\nYou can install numpy using \
        \n \tconda install 'numpy>=1.15' \n or \n\tpip install 'numpy>=1.15'")
    sys.exit(1)
except pkg_resources.VersionConflict:
    print("ERROR: numpy>=1.15 required, numpy %s found.\nYou can update numpy using \
        \n \tconda install 'numpy>=1.15' \n or \n\tpip install -U 'numpy>=1.15'" % pkg_resources.get_distribution("numpy").version)
    sys.exit(1)

