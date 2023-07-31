import re

from Bio import SeqIO
from .bidirectionalmap.bidirectionalmap import BidirectionalMap

def get_contig_lengths_spades(contigs_file):
    # Get length and coverage of contigs
    contig_lengths = {}
    coverages = {}

    my_map = BidirectionalMap()

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        start = 'NODE_'
        end = '_length'
        contig_num = int(re.search('%s(.*)%s' % (start, end), record.id).group(1))
        
        start = '_length_'
        end = '_cov'
        length = int(re.search('%s(.*)%s' % (start, end), record.id).group(1))
        
        start = '_cov_'
        end = ''
        coverage = int(float(re.search('%s(.*)%s' % (start, end), record.id).group(1)))
        
        contig_lengths[contig_num] = length
        coverages[contig_num] = coverage

    return contig_lengths, coverages


def get_contig_paths_spades(contig_paths):

    paths = {}
    segment_contigs = {}
    node_count = 0

    contig_names = {}

    my_map = BidirectionalMap()

    current_contig_num = ""

    with open(contig_paths) as file:
        name = file.readline()
        path = file.readline()
        
        while name != "" and path != "":
                
            while ";" in path:
                path = path[:-2]+","+file.readline()
            
            start = 'NODE_'
            end = '_length_'
            contig_num = str(int(re.search('%s(.*)%s' % (start, end), name).group(1)))
            
            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                my_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1
            
            if contig_num not in paths:
                paths[contig_num] = [segments[0], segments[-1]]
            
            for segment in segments:
                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)
            
            name = file.readline()
            path = file.readline()

    return my_map, 