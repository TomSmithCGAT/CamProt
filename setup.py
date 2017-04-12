import sys
import os
import glob

from ez_setup import use_setuptools
use_setuptools("10.0")
import setuptools

from setuptools import setup, Extension, find_packages

########################################################################
########################################################################
# collect camprot version
def getVersion():
    return "0.0.13"

version = getVersion()

###############################################################
###############################################################
# Define dependencies 
# Perform a camprot Installation

major, minor1, minor2, s, tmp = sys.version_info

if (major == 2 and minor1 < 7) or major < 2:
    raise SystemExit("""camprot requires Python 2.7 or later.""")

camprot_packages = find_packages()
camprot_package_dirs = {'camprot': 'camprot'}

install_requires = [
    "urllib3>=1.0",
    "pandas>=0.12.0",
    "requests>=2.0"]

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
Operating System :: Microsoft
"""

setup(
    # package information
    name='camprot',
    version=version,
    description='camprot: Tools for Computational proteomics',
    author='Tom Smith',
    author_email='tss38@cam.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational proteomics",
    long_description='camprot: Scripts and modules for computational proteomics',
    classifiers=list(filter(None, classifiers.split("\n"))),
    url="https://github.com/TomSmithCGAT/CamProt",
    download_url="https://github.com/TomSmithCGAT/CamProt/tarball/%s" % version,
    # package contents
    packages=camprot_packages,
    package_dir=camprot_package_dirs,
    include_package_data=True,
    # dependencies
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['camprot = camprot.camprot:main']
    },
    # other options
    zip_safe=False,
)
