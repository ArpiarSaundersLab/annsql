import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='AnnSQL',
    version='v0.9.2',
    author="Kenny Pavan",
    author_email="pavan@ohsu.edu",
    description="A Python SQL tool for converting Anndata objects to a relational DuckDb database. Methods are included for querying and basic single-cell preprocessing (experimental). ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kennypavan/AnnSQL",
    packages=setuptools.find_packages(where='src'),  
    package_dir={'': 'src'},  
	python_requires='>=3.9',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scanpy>=1.10.3',
		'duckdb>=1.1.2',
    ],
)
