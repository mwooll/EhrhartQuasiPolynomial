## -*- encoding: utf-8 -*-
from setuptools import setup

# Get information from separate files (README)
def readfile(filename):
    with open(filename,  encoding="utf-8") as f:
        return f.read()

setup(
    name="ehrhart_polynomial",
    version="0.1.0",
    description="Package to calculate the Ehrhart (quasi)polynomial of convex polytopes",
    long_description=readfile("README.md"), # get the long description from the README

    url="https://github.com/mwooll/Ehrhart",
    author="Mark Woolley",
    author_email="mark.woolley@uzh.ch",

    license="GPLv2+",
    classifiers=[
      "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
      "Programming Language :: Python :: 3"
    ], 
    keywords = "SageMath polytope ehrhart quasi-polynomial",
    packages = ["ehrhart_polynomial"],
)
