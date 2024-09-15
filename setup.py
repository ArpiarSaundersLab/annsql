import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='AnnSQL',
    version='0.1',
    author="Kenny Pavan",
    author_email="pavan@ohsu.edu",
    description="A SQL wrapper for querying AnnData objects.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kennypavan/AnnSQL",
    packages=setuptools.find_packages(where='src/AnnSQL'),  
    package_dir={'': 'src/AnnSQL'},  
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scanpy>=1',
    ],
)
