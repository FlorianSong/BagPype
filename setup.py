#! /usr/bin/env python3

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bagpype", # Replace with your own username
    version="1.0",
    author="Ben Amor, Francesca Vianello, Florian Song",
    author_email="florian.song@gmail.com",
    description="Biomolecular, atomistic graph construction software in Python for proteins, etc.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FlorianSong/BagPype",
    # project_urls={
    #     "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    # },
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    install_requires=(
        'numpy',
        'scipy',
        'networkx',
        'pandas',
        'PyCifRW',
        'requests'
    ),
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    package_data={"bagpype": ["dependencies/*", "dependencies/mmcif/*"]},
    python_requires=">=3.6",
)
