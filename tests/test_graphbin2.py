import subprocess
from pathlib import Path

import pytest

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2020, GraphBin2 Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD"
__version__ = "1.2.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


TEST_ROOTDIR = Path(__file__).parent


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_graphbin_version():
    """test graphbin2 help message"""
    cmd = "graphbin2 --help"
    exec_command(cmd)


def test_graphbin_on_spades_5g_dataset(tmp_dir):
    """test graphbin on spades assembly"""
    dir_name = TEST_ROOTDIR / "data" / "Sim-5G+metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    abundance = dir_name / "abundance.abund"
    binned = dir_name / "initial_contig_bins.csv"
    cmd = f"graphbin2 --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --abundance {abundance} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_spades_10g_dataset(tmp_dir):
    """test graphbin on spades assembly"""
    dir_name = TEST_ROOTDIR / "data" / "Sim-10G+metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    abundance = dir_name / "abundance.abund"
    binned = dir_name / "initial_contig_bins.csv"
    cmd = f"graphbin2 --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --abundance {abundance} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


# def test_graphbin_on_sga_5g_dataset(tmp_dir):
#     """test graphbin on spades assembly"""
#     dir_name = TEST_ROOTDIR / "data" / "Sim-5G+SGA"
#     graph = dir_name / "default-graph.asqg"
#     contigs = dir_name / "default-contigs.fa"
#     abundance = dir_name / "abundance.abund"
#     binned = dir_name / "initial_contig_bins.csv"
#     cmd = f"graphbin2 --assembler sga --graph {graph} --contigs {contigs} --abundance {abundance} --binned {binned} --output {tmp_dir}"
#     exec_command(cmd)
