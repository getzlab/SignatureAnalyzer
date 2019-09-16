from setuptools import setup

setup(
    name='siganalyzer',
    author='Shankara Anand',
    author_email='sanand@broadinstitute.org',
    version='0.0.0',
    description='Bayesian NMF methods for mutational signature analysis on GPUs (Getz Lab).',
    install_requires=[
        "twobitreader>=3.1.7"
    ]
)
