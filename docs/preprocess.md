# Preprocessing

## Assemble datasets

Firstly, you will have to assemble your set of reads into contigs. For this purpose, you can use metaSPAdes, SGA or metaFlye.

### metaSPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is a short-read assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble short reads into contigs. A sample command is given below.

```
spades --meta -1 Reads_1.fastq -2 Reads_2.fastq -o /path/output_folder -t 16
```

### SGA
[**SGA**](https://github.com/jts/sga) (String Graph Assembler) is a short-read assembler based on the overlap-layout-consensus (more recently string graph) approach. Use SGA software to assemble short reads into contigs. Sample commands are given below. You may change the parameters to suit your datasets.

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

### metaFlye
[**Flye**](https://github.com/fenderglass/Flye) is a long-read assembler based on the de Bruijn graph approach. [**metaFlye**](https://www.nature.com/articles/s41592-020-00971-x) is the dedicated metagenomic assembler of Flye. Use metaFlye (Flye in metagenomics mode) software to assemble long reads into contigs. A sample command is given below.

```
flye --meta --pacbio-raw reads.fasta --genome-size estimated_metagenome_size --out-dir /path/output_folder --threads 16
```

## Bin contigs

Next, you have to bin the resulting contigs using an existing contig-binning tool. We have used the following tools with their commands for the experiments.

### [MaxBin2](https://sourceforge.net/projects/maxbin2/)

```
perl MaxBin-2.2.5/run_MaxBin.pl -contig contigs.fasta -abund abundance.abund -thread 8 -out /path/output_folder
```

### [SolidBin](https://github.com/sufforest/SolidBin)

```
python scripts/gen_kmer.py /path/to/data/contig.fasta 1000 4 
sh gen_cov.sh 
python SolidBin.py --contig_file /path/to/contigs.fasta --composition_profiles /path/to/kmer_4.csv --coverage_profiles /path/to/cov_inputtableR.tsv --output /output/result.tsv --log /output/log.txt --use_sfs
```

## Prepare binning results

The binning output file should have delimiter separated (e.g., comma separated) values ```(contig_identifier, bin_number)``` for each contig. The contents of the binning output file should look similar to the example given below. Contigs are named according to their original identifier and the numbering of bins starts from 1. You can use the [`prepResult`](https://github.com/Vini2/GraphBin2/blob/master/src/support/prepResult.py) command to format an initial binning result in to the .csv format with contig identifiers and bin ID. Further details can be found [here](https://github.com/Vini2/GraphBin2/blob/master/src/graphbin2/support/README.md#prepresultpy) and in the next page.

### Example binned inputs

Example metaSPAdes binned input
```
NODE_1_length_507141_cov_16.465306,1
NODE_2_length_487410_cov_94.354557,1
NODE_3_length_483145_cov_59.410818,1
NODE_4_length_468490_cov_20.967912,2
NODE_5_length_459607_cov_59.128379,2
...
```
Example SGA binned input
```
contig-0,1
contig-1,2
contig-2,1
contig-3,1
contig-4,2
...
```
Example Flye binned input
```
edge_1,1
edge_2,2
edge_3,1
edge_4,1
edge_5,2
...
```

### Before using Flye assemblies for binning

Before using Flye assemblies for binning, please use [`gfa2fasta`](https://github.com/Vini2/GraphBin2/blob/master/src/graphbin2/support/gfa2fasta.py) command to get the edge sequences. Further details can be found [here](https://github.com/Vini2/GraphBin2/blob/master/src/graphbin2/support/README.md#gfa2fastapy). More details can be found in the next page.

## Obtain the coverage of contigs (`abundance.tsv`)

You can use [CoverM](https://github.com/wwood/CoverM) to get the coverage of contigs. You can run the following commands to get the `abundance.tsv` file. Please make sure that there are **no headers** in the `abundance.tsv` file.

```
coverm contig -1 reads_1.fastq -2 reads_2.fastq -r contigs.fasta -o abundance.tsv -t 8
sed -i '1d' abundance.tsv   # remove the header of the file
```

The resulting `abundance.tsv` file can be directly used in GraphBin2.

Once you have obtained the assembly output, binning results and the coverage information file, you can run GraphBin2.