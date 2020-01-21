import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="windel",
    version="0.0.1",
    author="Jesse Wallace",
    author_email="jesse.wallace@anu.edu.au",
    description="Sliding-window indel correction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/J-Wall/windel",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=["fire==0.2", "numpy==1.17", "pysam==0.15", "scipy==1.4"],
)
