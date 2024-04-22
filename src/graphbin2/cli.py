#!/usr/bin/env python3

"""graphbin2: Refined and overlapped binning of metagenomic contigs using assembly graphs."""

import logging
import sys

import click

from graphbin2 import graphbin2_Flye, graphbin2_SGA, graphbin2_SPAdes

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2020, GraphBin2 Project"
__license__ = "BSD"
__version__ = "1.2.3"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Stable Release"


class ArgsObj:
    def __init__(
        self,
        assembler,
        graph,
        contigs,
        paths,
        abundance,
        binned,
        output,
        prefix,
        depth,
        threshold,
        delimiter,
        nthreads,
    ):
        self.assembler = assembler
        self.graph = graph
        self.contigs = contigs
        self.paths = paths
        self.abundance = abundance
        self.binned = binned
        self.output = output
        self.prefix = prefix
        self.depth = depth
        self.threshold = threshold
        self.delimiter = delimiter
        self.nthreads = nthreads


@click.command()
@click.option(
    "--assembler",
    help="name of the assembler used. (Supports SPAdes, SGA and Flye)",
    type=click.Choice(["spades", "sga", "flye"], case_sensitive=False),
    required=True,
)
@click.option(
    "--graph",
    help="path to the assembly graph file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--contigs",
    help="path to the contigs file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--paths",
    help="path to the contigs.paths (metaSPAdes) or assembly.info (metaFlye) file",
    type=click.Path(exists=True),
    required=False,
)
@click.option(
    "--abundance",
    help="path to the abundance file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--binned",
    help="path to the .csv file with the initial binning output from an existing toole",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--output",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)
@click.option(
    "--prefix",
    help="prefix for the output file",
    type=str,
    required=False,
)
@click.option(
    "--depth",
    help="maximum depth for the breadth-first-search.",
    type=int,
    default=5,
    show_default=True,
    required=False,
)
@click.option(
    "--threshold",
    help="threshold for determining inconsistent vertices.",
    type=float,
    default=1.5,
    show_default=True,
    required=False,
)
@click.option(
    "--delimiter",
    help="delimiter for output results. Supports a comma (,), a semicolon (;), a tab ($'\\t'), a space (\" \") and a pipe (|) .",
    type=click.Choice([",", ";", "$'\\t'", '" "'], case_sensitive=False),
    default=",",
    show_default=True,
    required=False,
)
@click.option(
    "--nthreads",
    help="number of threads to use.",
    type=int,
    default=8,
    show_default=True,
    required=False,
)
@click.version_option(__version__, "-v", "--version", is_flag=True)
def main(
    assembler,
    graph,
    contigs,
    paths,
    abundance,
    binned,
    output,
    prefix,
    depth,
    threshold,
    delimiter,
    nthreads,
):
    """
    GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs Using Assembly Graphs
    GraphBin2 is a tool which refines the binning results obtained from existing tools and, is able to
    assign contigs to multiple bins. GraphBin2 uses the connectivity and coverage information from
    assembly graphs to adjust existing binning results on contigs and to infer contigs shared by multiple species.
    """

    assembly_graph_file = graph
    contigs_file = contigs
    contig_paths = paths
    abundance_file = abundance
    contig_bins_file = binned
    output_path = output

    # Validate prefix
    if prefix != None:
        if not prefix.endswith("_"):
            prefix = prefix + "_"
    else:
        prefix = ""

    # Setup logger
    # -----------------------
    logger = logging.getLogger("GraphBin2 1.2.0")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    # ---------------------------------------------------

    fileHandler = logging.FileHandler(f"{output_path}/{prefix}graphbin2.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    # Validation of inputs
    # ---------------------------------------------------

    # Check if paths file is provided when the assembler type is SPAdes
    if assembler.lower() == "spades" and contig_paths is None:
        logger.error("Please make sure to provide the path to the contigs.paths file.")
        logger.info("Exiting GraphBin2...\nBye...!\n")
        sys.exit(1)

    # Check if paths file is provided when the assembler type is SPAdes
    if assembler.lower() == "flye" and contig_paths is None:
        logger.error(
            "Please make sure to provide the path to the assembly_info.txt file."
        )
        logger.info("Exiting GraphBin2...\nBye...!\n")
        sys.exit(1)

    # Validate paths
    if assembler.lower() == "sga":
        contig_paths = "None"

    # Validate depth
    if depth < 1:
        logger.error("Please enter a valid number for depth")
        logger.info("Exiting GraphBin2...\nBye...!\n")
        sys.exit(1)

    # Validate threshold
    if threshold < 1.0:
        logger.error("Please enter a valid number for threshold")
        logger.info("Exiting GraphBin2...\nBye...!\n")
        sys.exit(1)

    # Validate number of threads
    if nthreads <= 0:
        logger.error("Please enter a valid number for the number of threads")
        logger.info("Exiting GraphBin2...\nBye...!\n")
        sys.exit(1)

    # Start GraphBin2
    # ---------------------------------------------------

    logger.info(
        "Welcome to GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs using Assembly Graphs."
    )

    if assembler.lower() == "spades":
        logger.info(
            "This version of GraphBin2 makes use of the assembly graph produced by SPAdes which is based on the de Bruijn graph approach."
        )
    elif assembler.lower() == "sga":
        logger.info(
            "This version of GraphBin2 makes use of the assembly graph produced by SGA which is based on the string graph approach."
        )
    elif assembler.lower() == "flye":
        logger.info(
            "This version of GraphBin2 makes use of the assembly graph produced by metaFlye which is a long reads assembler based on the de Bruijn graph approach."
        )

    logger.info("Input arguments:")
    logger.info(f"Contigs file: {contigs_file}")
    logger.info(f"Assembly graph file: {assembly_graph_file}")
    logger.info(f"Contig paths file: {contig_paths}")
    logger.info(f"Existing binning output file: {contig_bins_file}")
    logger.info(f"Final binning output file: {output_path}")
    logger.info(f"Depth: {depth}")
    logger.info(f"Threshold: {threshold}")
    logger.info(f"Number of threads: {nthreads}")

    logger.info("GraphBin2 started")

    # Make args object
    args = ArgsObj(
        assembler,
        graph,
        contigs,
        paths,
        abundance,
        binned,
        output,
        prefix,
        depth,
        threshold,
        delimiter,
        nthreads,
    )

    # Run GraphBin2
    # ---------------------------------------------------
    if assembler.lower() == "spades":
        graphbin2_SPAdes.main(args)

    elif assembler.lower() == "sga":
        graphbin2_SGA.main(args)

    elif assembler.lower() == "flye":
        graphbin2_Flye.main(args)


if __name__ == "__main__":
    main()
