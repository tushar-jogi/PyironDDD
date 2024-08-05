"""
Setuptools based setup module
"""
from setuptools import setup, find_packages


setup(
    name='PyironDDD',
    url='https://github.com/tushar-jogi/PyironDDD/',
    license='BSD',

    packages=find_packages(exclude=["*tests*"]),
    install_requires=[
        'pyiron_base==0.9.10',
    ],
)
