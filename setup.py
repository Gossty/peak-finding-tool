import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'peakFinding/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )

setup(
    include_package_data=True,
    name='peakFinding',
    version=VERSION,
    description="find peaks in a provided tag directory. It takes as input a tag directory containing tag files, and it optionally takes a control directory for peak finding as well as other parameters. The code performs various filtering and analysis steps to identify peaks and outputs the results in a BED file format and a txt file with statistical information.",
    author="Stephen Golzari, Yashwin Madakamutil, Yang Han",
    author_email="sgolzari@ucsd.edu, ymadakamutil@ucsd.edu, yah015@ucsd.edu",
    packages=find_packages(),
    install_requires=['pandas', 'numpy', 'argparse', 'scipy'],
    entry_points= {
        "console_scripts": [
            "peakFinding=peakFinding:main"
        ],
    },
)