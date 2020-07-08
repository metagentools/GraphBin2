<p align="center">
  <img src="GraphBin2_Logo.png" width="450" title="GraphBin2 Logo" alt="GraphBin2 Logo">
</p>

# GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs Using Assembly Graphs

![GitHub](https://img.shields.io/github/license/Vini2/GraphBin2) 
![GitHub top language](https://img.shields.io/github/languages/top/Vini2/GraphBin2)

GraphBin2 is an extension of [GraphBin](https://github.com/Vini2/GraphBin) which refines the binning results obtained from existing tools and, more importantly, is able to assign contigs to multiple bins. GraphBin2 uses the connectivity and coverage information from assembly graphs to adjust existing binning results on contigs and to infer contigs shared by multiple species.

## Getting Started

### Dependencies
You will need the following python packages installed.
* [Biopython](https://biopython.org/)
* [python-igraph](https://igraph.org/python/)

### Downloading GraphBin2
You can download the latest release of GraphBin from [Releases](https://github.com/Vini2/GraphBin/releases) or clone the GraphBin repository to your machine.

```
git clone https://github.com/Vini2/GraphBin2.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```
unzip [file_name].zip
```

Now go in to the GraphBin folder using the command

```
cd GraphBin2/src/
```

## Preprocessing

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes or SGA.

### metaSPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble reads into contigs. A sample command is given below.

```
spades --meta -1 Reads_1.fastq -2 Reads_2.fastq -o /path/output_folder -t 16
```

### SGA
[**SGA**](https://github.com/jts/sga) (String Graph Assembler) is an assembler based on the overlap-layout-consensus (more recently string graph) approach. Use SGA software to assemble reads into contigs. Sample commands are given below. You may change the parameters to suit your datasets.

```
sga preprocess -o reads.fastq --pe-mode 1 Reads_1.fastq Reads_2.fastq
sga index -a ropebwt -t 16 --no-reverse reads.fastq
sga correct -k 41 --learn -t 16 -o reads.k41.fastq reads.fastq
sga index -a ropebwt -t 16 reads.k41.fastq
sga filter -x 2 -t 16 reads.k41.fastq
sga fm-merge -m 45 -t 16  reads.k41.filter.pass.fa
sga index -t 16 reads.k41.filter.pass.merged.fa
sga overlap -m 55 -t 16 reads.k41.filter.pass.merged.fa
sga assemble -m 95 reads.k41.filter.pass.merged.asqg.gz
```
Next, you have to bin the resulting contigs using an existing contig-binning tool. We have used the following tools with their commands for the experiments.

#### [MaxBin2](https://sourceforge.net/projects/maxbin2/)

```
perl MaxBin-2.2.5/run_MaxBin.pl -contig contigs.fasta -abund abundance.abund -thread 8 -out /path/output_folder
```

#### [SolidBin](https://github.com/sufforest/SolidBin)

```
python scripts/gen_kmer.py /path/to/data/contig.fasta 1000 4 
sh gen_cov.sh 
python SolidBin.py --contig_file /path/to/contigs.fasta --composition_profiles /path/to/kmer_4.csv --coverage_profiles /path/to/cov_inputtableR.tsv --output /output/result.tsv --log /output/log.txt --use_sfs
```

## Using GraphBin2
You can see the usage options of GraphBin by typing ```python graphbin2.py -h``` on the command line. For example,

```
python graphbin2.py -h
usage: graphbin2.py [-h] --assembler ASSEMBLER --graph GRAPH --contigs CONTIGS
                    [--paths PATHS] [--abundance ABUNDANCE] --binned BINNED
                    --output OUTPUT [--prefix PREFIX] [--depth DEPTH]
                    [--threshold THRESHOLD] [--nthreads NTHREADS]

GraphBin2 Help. GraphBin2 is a tool which refines the binning results obtained
from existing tools and, more importantly, is able to assign contigs to
multiple bins. GraphBin2 uses the connectivity and coverage information from
assembly graphs to adjust existing binning results on contigs and to infer
contigs shared by multiple species.

optional arguments:
  -h, --help            show this help message and exit
  --assembler ASSEMBLER
                        name of the assembler used (SPAdes or SGA)
  --graph GRAPH         path to the assembly graph file
  --contigs CONTIGS     path to the contigs file
  --paths PATHS         path to the contigs.paths file
  --abundance ABUNDANCE
                        path to the abundance file
  --binned BINNED       path to the .csv file with the initial binning output
                        from an existing tool
  --output OUTPUT       path to the output folder
  --prefix PREFIX       prefix for the output file
  --depth DEPTH         maximum depth for the breadth-first-search. [default:
                        5]
  --threshold THRESHOLD
                        threshold for determining inconsistent vertices.
                        [default: 1.5]
  --nthreads NTHREADS   number of threads to use. [default: 8]
```

## Input Format

For the SPAdes version of `graphbin2.py` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Assembly graph file (in `.gfa` format)
* Paths of contigs (in `.paths` format)
* Binning output from an existing tool (in `.csv` format)

For the SGA version of `graphbin2.py` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Abundance file (tab separated file with contig ID and coverage in each line)
* Assembly graph file (in `.asqg` format)
* Binning output from an existing tool (in `.csv` format)

**Note:** The binning output file should have comma separated values ```(contig_number, bin_number)``` for each contig. The contents of the binning output file should look similar to the example given below. Contigs are named according to their number starting from 0 and the numbering of bins starts from 1.

Example binned input
```
0,1
1,2
2,1
3,1
4,2
...
```

## Example Usage

```
python graphbin2.py --assembler spades --contigs /path/to/contigs.fasta --graph /path/to/graph_file.gfa --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
python graphbin2.py --assembler sga --contigs /path/to/contigs.fa --abundance /path/to/abundance.tsv --graph /path/to/graph_file.asqg --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Visualization of the metaSPAdes Assembly Graph of the Sim-5G Dataset

### Initial Binning Result
<p align="center">
  <img src="images/initial_binning_result.svg" width="500" title="Initial assembly graph" alt="Initial binning result">
</p>

### Assembly Graph with Refined Labels
<p align="center">
  <img src="images/label_refined.svg" width="500" title="Initial assembly graph" alt="Labels refined">
</p>

### Assembly Graph after Label Propagation
<p align="center">
  <img src="images/label_propagated.svg" width="500" title="Initial assembly graph" alt="Labels propagated">
</p>

### Assembly Graph with Multi-labelled Vertices
<p align="center">
  <img src="images/multiple_marked.svg" width="500" title="Initial assembly graph" alt="Multi-labelled">
</p>

## References

[1] Barnum, T.P., _et al._: Genome-resolved metagenomics identifies genetic mobility, metabolic interactions, and unexpected diversity in perchlorate-reducing communities. The ISME Journal 12, 1568-1581 (2018)

[2] Mallawaarachchi, V., Wickramarachchi, A., Lin, Y.: GraphBin: Refined binning of metagenomic contigs using assembly graphs. Bioinformatics, btaa180 (2020)

[3] Nurk, S., _et al._: metaSPAdes: a new versatile metagenomic assembler. Genome Researcg 5, 824-834 (2017)

[4] Simpson, J. T. and Durbin, R.: Efficient de novo assembly of large genomes using compressed data structures. Genome Research, 22(3), 549–556 (2012).

[5] Wang, Z., _et al._:  SolidBin: improving metagenome binning withsemi-supervised normalized cut. Bioinformatics 35(21), 4229–4238 (2019).

[6] Wu, Y.W., _et al._: MaxBin: an automated binning method to recover individual genomes from metagenomes using an expectation-maximization algorithm. Microbiome 2(1), 26 (2014)

[7] Wu, Y.W., _et al._: MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32(4), 605–607 (2016)
