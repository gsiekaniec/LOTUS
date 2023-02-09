# LOTUS installation

## First way

The easiest way to install it is through conda:
```
conda create -p lotus_env
conda activate lotus_env
conda install -c gsiekaniec -c conda-forge -c bioconda lotus
```

It is then possible to execute LOTUS with:
```
lotus -h
```

---

## Second way

A second possibility is simply to clone the github directory: 
```
git clone git@github.com:gsiekaniec/LOTUS.git
``` 

Then install python 3 with the necessary packages for LOTUS to work properly (see the list below in section [Python packages]()). 
To simplify the installation, it is possible to create a conda environment containing python 3 and the necessary packages via the file ``environment.yml`` with the command:
```
conda env create -f {PATH_TO_LOTUS}/conda/environment.yml
conda activate lotus
```

> __Note__
The ```environment.yml``` file can be found [here](https://github.com/gsiekaniec/LOTUS/blob/main/conda/environment.yml).

Using this method, it is then possible to call LOTUS with:
```
python3 {PATH_TO_LOTUS}/lotus.py -h
```

Unlike the first method, this way of doing things also allows you to run unit tests with:
```
python3 -m py.test {PATH_TO_LOTUS}/tests
```

## Python packages

:file_folder: Packages from the [Python Standard Library](https://docs.python.org/3/library/) used:

  - [argparse](https://docs.python.org/3/library/argparse.html)
  - [collections](https://docs.python.org/3/library/collections.html)
  - [copy](https://docs.python.org/3/library/copy.html)
  - [gzip](https://docs.python.org/3/library/gzip.html)
  - [itertools](https://docs.python.org/3/library/itertools.html)
  - [json](https://docs.python.org/3/library/json.html)
  - [logging](https://docs.python.org/3/library/logging.html)
  - [os](https://docs.python.org/3/library/os.html)
  - [pathlib](https://docs.python.org/3/library/pathlib.html)
  - [pickle](https://docs.python.org/3/library/pickle.html)
  - [re](https://docs.python.org/3/library/re.html)
  - [sys](https://docs.python.org/3/library/sys.html)
  - [uuid](https://docs.python.org/3/library/uuid.html)
  - [warnings](https://docs.python.org/3/library/warnings.html)
  
:file_folder: Required python packages to run LOTUS:
  
  - [matplotlib](https://matplotlib.org/)
  - [more_itertools](https://more-itertools.readthedocs.io/en/stable/)
  - [numpy](https://numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [pyfastx](https://pyfastx.readthedocs.io/en/latest/)
  - [pytest](https://docs.pytest.org/en/7.2.x/)
  - [requests](https://requests.readthedocs.io/en/latest/)
  - [tqdm](https://tqdm.github.io/)
  - [uspsetplot](https://upsetplot.readthedocs.io/en/stable/)

<sub>:warning: These packages must be installed before you can use LOTUS.</sub>
