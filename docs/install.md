# Installation

## PyPi
AnnSQL requires a Python 3.12 or greater environment. **Note:** Higher memory consumption using Apple M-Series is expected when building AnnSQL databases. Please use the *make_buffer_file* parameter when using the *MakeDb* class if you experience memory related issues. 
```
pip install annsql
```

## Conda
Currently, there isn't a conda specific install, however, you may of course create a fresh conda environment, then install the AnnSql package using pip inside the conda environment.
```
#create the env
conda create -n annsql python=3.12

#activate the environment
conda activate annsql

#install in the conda env using pip
pip install annsql
```