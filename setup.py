from pathlib import Path
import re
import sys

from setuptools import find_packages, setup

ROOT = Path(__file__).parent

ver_info = sys.version_info
if ver_info < (3, 9, 0):
    raise RuntimeError("signatureanalyzer requires Python 3.9 or newer")

init_path = ROOT / "signatureanalyzer" / "__init__.py"
version_match = re.search(r"__version__\s*=\s*'([^']+)'", init_path.read_text(encoding="utf-8"))
if version_match is None:
    raise RuntimeError("Unable to determine package version")
version = version_match.group(1)

readme = (ROOT / "README.md").read_text(encoding="utf-8")

setup(
    name="signatureanalyzer",
    version=version,
    author="Shankara Anand & Justin Cha - Broad Institute - Cancer Genome Computational Analysis",
    author_email="sanand@broadinstitute.org",
    url="https://github.com/broadinstitute/getzlab-SignatureAnalyzer",
    long_description=readme,
    long_description_content_type="text/markdown",
    description="Bayesian NMF methods for mutational signature analysis & transcriptomic profiling on GPUs (Getz Lab).",
    python_requires=">=3.9",
    packages=find_packages(exclude=("tests", "tests.*")),
    include_package_data=True,
    install_requires=[
        "gprofiler>=0.1.0",
        "h5py>=3.10.0",
        "matplotlib>=3.8.0",
        "numpy>=2.0.0",
        "pandas>=2.2.0",
        "pyarrow>=16.0.0",
        "scikit-image>=0.24.0",
        "scikit-learn>=1.4.0",
        "scipy>=1.13.0",
        "seaborn>=0.13.0",
        "tables>=3.6.1",
        "tqdm>=4.33.0",
        "twobitreader>=3.1.7",
    ],
    extras_require={
        "gpu": ["torch>=2.3.0"],
    },
    package_data={
        "signatureanalyzer": [
            "ref/cosmic_v2/sa*",
            "ref/cosmic_v3/sa*",
        ]
    },
    entry_points={
        "console_scripts": [
            "signatureanalyzer = signatureanalyzer.__main__:main"
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Interface Engine/Protocol Translator",
    ],
    license="MIT",
)
