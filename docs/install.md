# Setting up GraphBin2

## Downloading GraphBin2

You can download the latest release of GraphBin2 from [Releases](https://github.com/Vini2/GraphBin2/releases) or clone the GraphBin2 repository to your machine.

```bash
git clone https://github.com/Vini2/GraphBin2.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```bash
unzip [file_name].zip
```

Now go in to the GraphBin2 folder using the command

```bash
cd GraphBin2/
```

## Setting up the environment

We recommend that you use [Conda](https://docs.conda.io/en/latest/) to run GraphBin2. You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, make sure you are in the GraphBin2 folder. Now run the following commands to create a Conda environment and activate it to run GraphBin2.

```bash
conda env create -f environment.yml
conda activate graphbin2
```

Now install GraphBin2 using the following command.

```bash
flit install
```

## Test the setup

After installing, run the following command to ensure that GraphBin2 is working.

```bash
graphbin2 -h
```

Now you are ready to run GraphBin2.