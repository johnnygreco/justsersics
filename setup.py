#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as required_file:
    required = required_file.read()

requirements = required.split('\n')[:-1]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Johnny Greco",
    author_email='jgreco.astro@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Just Sersic fitting.",
    entry_points={
        'console_scripts': [
            'justsersics=justsersics.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='justsersics',
    name='justsersics',
    packages=find_packages(include=['justsersics', 'justsersics.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/johnnygreco/justsersics',
    version='0.1.0',
    zip_safe=False,
)
