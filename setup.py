# -*- coding: utf-8 -*-

import sys
import subprocess

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install
import epigeneticgraphanalyzer._version as version

# Function to parse requirements.txt
def parse_requirements(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        return [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]

# Reading dependencies
requirements = parse_requirements('requirements.txt')


def get_version():
    return version.__version__


class sdist(_sdist):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

# Install class to check for external dependencies from OS environment


class install(_install):

    def run(self):
        # update_version_py()
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            # the config file is installed in:
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))

setup(
    name='Epigenetic Graph Analyzer',
    version=get_version(),
    author='Joachim Wolff',
    author_email='wolff.joachim@mh-hannover.de',
    long_description_content_type='text/markdown',
    url='https://github.com/joachimwolff/Epigenetic-Graph-Analyzer',
    packages=find_packages(),
    install_requires=requirements,
    scripts=[
        'bin/egaCreateGraph',
        'bin/egaMergeGraphs',
        'bin/egaClusterBaits',
        'bin/egaPlotlyGraph',
        'bin/egaComputeRLM',
        'bin/egaMergeRLM',
        'bin/egaValidateWithSPRITE',
        'bin/egaPrepareInteractions',
        'bin/egaFetchPromoters',
        'bin/egaReconstructCutSites',
        'bin/egaTADStats',
        'bin/egaCheckDownloads',
    ],
    include_package_data=True,
    license='LICENSE',
    description='Set of programs to integrate the analysis of capture Hi-C data RNA-seq data and epigenetic data',
    long_description=open('README.md', encoding='utf-8').read(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    zip_safe=False,
    cmdclass={'sdist': sdist, 'install': install}
)