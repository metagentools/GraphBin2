# Support scripts for GraphBin2

## prepResult.py

`prepResult.py` is a support script that allows you to format an initial binning result in to the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered starting from 1. You can run `prepResult.py` as follows.

```
python prepResult.py --binned /path/to/folder_with_binning_result --output /path/to/output_folder
```
You can see the usage options of `prepResult.py` by typing `python prepResult.py -h` on the command line.

```
usage: prepResult.py [-h] --binned BINNED --output OUTPUT
                     [--delimiter DELIMITER] [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --binned BINNED       path to the folder containing the initial binning
                        result from an existing tool
  --output OUTPUT       path to the output folder
  --delimiter DELIMITER
                        delimiter for results. Supports a comma (,), a
                        semicolon (;), a tab ($'\t'), a space (" ") and a pipe
                        (|) [default: , (comma)]
  --prefix PREFIX       prefix for the output file
```

Formatted binning result will be stored in a file named `initial_contig_bins.csv` in the output folder provided. Bin IDs and corresponding fasta files for each bin will be recorded in a file named `bin_ids.csv` in the output folder provided.

You can also specify the delimiter for the initial binning result file using the `delimiter` paramter. Enter the following values for different delimiters; 
* `,` for a comma
* `;` for a semicolon
* `$'\t'` for a tab
* `" "` for a space 
* `|` for a pipe.

## gfa2fasta.py

Please note that, if you are using Flye/Miniasm assemblies, you should provide the edge sequences for the initial binning tool (not the contigs output from Flye/Miniasm). To get the edge sequences from the GFA file, you can use the script [`gfa2fasta.py`](https://github.com/Vini2/GraphBin2/blob/master/support/gfa2fasta.py) as the assembly graph consists of these edge sequences and not contigs.

You can see the usage options of `gfa2fasta.py` by typing `python gfa2fasta.py -h` on the command line.

```
usage: gfa2fasta.py [-h] --graph GRAPH --assembler ASSEMBLER --output OUTPUT
                    [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --graph GRAPH         path to the assembly graph file
  --assembler ASSEMBLER
                        type of the assembler (Flye or Miniasm)
  --output OUTPUT       path to the output folder
  --prefix PREFIX       prefix for the output file
```
