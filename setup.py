# file: setup.py
# content: setup file for capycity package
# created: 2022 November 19 Monday
# author: roch schanen

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "capycity",
    version = "0.0.1",
    author = "Roch Schanen",
    author_email = "r.schanen@lancaster.ac.uk",
    description = "Capacitance Calculator",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/RochSchanen/capycity_dev",
    packages = ['capycity'],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['numpy', 'matplotlib'],
    python_requires = '>=3.10'
)
