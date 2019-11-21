try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='DataGen',
    packages=find_packages(),
    version='1.0.0',
    description='A python API for conformation generation, optimization, calculation, and data prepartion',
    long_description=long_description,
    author='Jianing Lu ',
    license='MIT')
