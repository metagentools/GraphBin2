# Using GraphBin2

You can see the usage options of GraphBin2 by typing `./graphbin2 -h` on the command line. For example,

```shell
Usage: graphbin2 [OPTIONS]

  GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs Using
  Assembly Graphs GraphBin2 is a tool which refines the binning results
  obtained from existing tools and, is able to  assign contigs to multiple
  bins. GraphBin2 uses the connectivity and coverage information from
  assembly graphs to adjust existing binning results on contigs and to infer
  contigs shared by multiple species.

Options:
  --assembler [spades|megahit|sga|flye]
                                  name of the assembler used. (Supports
                                  SPAdes, SGA and Flye)  [required]
  --graph PATH                    path to the assembly graph file  [required]
  --contigs PATH                  path to the contigs file  [required]
  --paths PATH                    path to the contigs.paths (metaSPAdes) or
                                  assembly.info (metaFlye) file
  --abundance PATH                path to the abundance file  [required]
  --binned PATH                   path to the .csv file with the initial
                                  binning output from an existing toole
                                  [required]
  --output PATH                   path to the output folder  [required]
  --prefix TEXT                   prefix for the output file
  --depth INTEGER                 maximum depth for the breadth-first-search.
                                  [default: 5]
  --threshold FLOAT               threshold for determining inconsistent
                                  vertices.  [default: 1.5]
  --delimiter [,|;|$'\t'|" "]     delimiter for output results. Supports a
                                  comma (,), a semicolon (;), a tab ($'\t'), a
                                  space (" ") and a pipe (|) .  [default: ,]
  --nthreads INTEGER              number of threads to use.  [default: 8]
  -v, --version                   Show the version and exit.
  --help                          Show this message and exit.
```

# Input Format

The SPAdes version of `graphbin2` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Assembly graph file (in `.gfa` format)
* Paths of contigs (in `.paths` format)
* Binning output from an existing tool (in `.csv` format)

The SGA version of `graphbin2` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Abundance file (tab separated file with contig ID and coverage in each line)
* Assembly graph file (in `.asqg` format)
* Binning output from an existing tool (in `.csv` format)

The MEGAHIT version of `graphbin2` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Abundance file (tab separated file with contig ID and coverage in each line)
* Assembly graph file (in `.gfa` format)
* Binning output from an existing tool (in `.csv` format)

The Flye version of `graphbin2` takes in 4 files as inputs (required).
* Contigs file (in `.fasta` format)
* Abundance file (tab separated file with contig ID and coverage in each line)
* Assembly graph file (in `.gfa` format)
* Binning output from an existing tool (in `.csv` format)

**Note:** The abundance file (e.g., `abundance.abund`) is a tab separated file with contig ID and the coverage for each contig in the assembly. metaSPAdes provides the coverage of each contig in the contig identifier of the final assembly. We can directly extract these values to create the abundance.abund file. However, no such information is provided for contigs produced by SGA. Hence, reads should be mapped back to the assembled contigs in order to determine the coverage of SGA contigs.

**Note:** Make sure that the initial binning result consists of contigs belonging to only one bin. GraphBin2 is designed to handle initial contigs which belong to only one bin.

**Note:** You can specify the delimiter for the initial binning result file and the final output file using the `delimiter` paramter. Enter the following values for different delimiters; `,` for a comma, `;` for a semicolon, `$'\t'` for a tab, `" "` for a space and `|` for a pipe.


## Example Usage

```shell
# metaSPAdes assembly
graphbin2 --assembler spades --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --graph /path/to/graph_file.gfa  --binned /path/to/binning_result.csv --abundance /path/to/abundance.tsv --output /path/to/output_folder
```
```shell
# SGA assembly
graphbin2 --assembler sga --contigs /path/to/contigs.fa --graph /path/to/graph_file.asqg --binned /path/to/binning_result.csv --abundance /path/to/abundance.tsv --output /path/to/output_folder
```
```shell
# MEGAHIT version
graphbin2 --assembler megahit --graph /path/to/final.gfa --contigs /path/to/final.contigs.fa --binned /path/to/binning_result.csv --abundance /path/to/abundance.tsv --output /path/to/output_folder
```
```shell
# metaFlye assembly
graphbin2 --assembler flye --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --graph /path/to/graph_file.gfa --binned /path/to/binning_result.csv --abundance /path/to/abundance.tsv --output /path/to/output_folder
```