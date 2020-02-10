from setuptools import setup
import re
import os
import sys

ver_info = sys.version_info
if ver_info < (3,6,0):
    raise RuntimeError("signatureanalyzer requires at least python 3.6.0")

with open(os.path.join(os.path.dirname(__file__), 'signatureanalyzer', '__init__.py')) as r:
    version = re.search(r'__version__ = \'(\d+\.\d+\.\d+[-_a-zA-Z0-9]*)\'', r.read()).group(1)

setup(
    name = 'signatureanalyzer',
    version = version,
    author = 'Shankara Anand & Justin Cha - Broad Institute - Cancer Genome Computational Analysis',
    author_email = 'sanand@broadinstitute.org',
    url = 'https://github.com/broadinstitute/getzlab-SignatureAnalyzer',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description = 'Bayesian NMF methods for mutational signature analysis & transcriptomic profiling on GPUs (Getz Lab).',
    packages = [
        'signatureanalyzer',
        'signatureanalyzer.plotting',
        'signatureanalyzer.pathways',
        'signatureanalyzer.signatureanalyzer_gpu'
    ],
    install_requires = [
        "twobitreader>=3.1.7",
        "torch>=1.2.0",
        "seaborn>=0.9.0",
        "pandas>=0.25.0",
        "pyarrow>=0.14.1",
        "scikit-learn>=0.21.3",
        "scikit-image>=0.15.0",
        "tqdm>=4.33.0",
        "h5py>=2.9.0",
        "tables>=3.6.1",
        "missingpy",
        "gprofiler",
        "numpy",
        "matplotlib",
        "scipy"
    ],
    package_data = {
    "":[
        "ref/cosmic_v2/sa*",
        "ref/cosmic_v3/sa*"
    ]
    },
    entry_points = {
        'console_scripts': [
            'signatureanalyzer = signatureanalyzer.__main__:main'
        ]
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Interface Engine/Protocol Translator",
    ],
    license="MIT"
)
