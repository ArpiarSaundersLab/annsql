import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name='AnnSQL',
	version='v0.9.9.1',
	author="Kenny Pavan",
	author_email="pavan@ohsu.edu",
	description="A Python SQL tool for converting Anndata objects to a relational DuckDb database. Methods are included for querying and basic single-cell preprocessing (experimental). ",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/ArpiarSaundersLab/annsql",
	packages=setuptools.find_packages(where='src'),  
	package_dir={'': 'src'},  
	python_requires='>=3.12',
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	install_requires=[
		'scanpy>=1.10.3',
		'duckdb>=1.1.2',
		'memory-profiler>=0.61.0',
		'psutil>=6.0.0',
		'pyarrow>=17.0.0',
	],
)
