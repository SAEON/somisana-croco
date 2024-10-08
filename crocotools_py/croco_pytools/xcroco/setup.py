import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="xcroco",

    description="Python post-processing for Croco, based on xgcm",

    packages=find_packages(exclude=[]),

    long_description=read('README.md'),
)
