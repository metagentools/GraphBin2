# Installing GraphBin2

## Using Conda

You can install GraphBin2 using the [bioconda](https://anaconda.org/bioconda/graphbin2) distribution. You can download 
[Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains `conda`.

```shell
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment and install
conda create -n graphbin2 graphbin2

# activate conda environment
conda activate graphbin2
```

## Using pip

You can install GraphBin2 using `pip` from the [PyPI](https://pypi.org/project/graphbin2/) distribution.

```shell
# install graphbin2
pip install graphbin2
```

## Test the setup

After installing, run the following command to ensure that GraphBin2 is working.

```shell
graphbin2 --help
```

Now you are ready to run GraphBin2.