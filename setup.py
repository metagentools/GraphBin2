#!/usr/bin/env python3

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

packages = setuptools.find_packages()
package_data = {"src": ["src/*", "src/bidirectionalmap/*", "src/support/*"]}

data_files = [(".", ["LICENSE", "README.md"])]

setuptools.setup(
    name="graphbin2",
    version="1.2.0",
    zip_safe=True,
    author="Vijini Mallawaarachchi, Anuradha Wickramarachchi and Yu Lin",
    author_email="viji.mallawaarachchi@gmail.com",
    description="GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs Using Assembly Graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/metagentools/GraphBin2",
    license="BSD-3",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["graphbin2"],
    entry_points={
        "console_scripts": [
            "gfa2fasta=src.support.gfa2fasta:main",
            "prep_result=src.support.prepResult:main",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        "biopython",
        "python-igraph",
        "scipy",
        "tqdm",
        "click",
    ],
    python_requires=">=3.7",
)
