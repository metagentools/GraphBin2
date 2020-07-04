#!/usr/bin/env python3

import sys
import csv
import time
import argparse
import re
import heapq
import itertools as it

from multiprocessing import Pool
from Bio import SeqIO
from igraph import *
from collections import defaultdict


# Setup argument parser
#---------------------------------------------------

ap = argparse.ArgumentParser(description="""GraphBin2 Help. GraphBin2 is a tool which refines the binning results obtained from existing tools and, 
more importantly, is able to assign contigs to multiple bins. GraphBin2 uses the connectivity and coverage information from assembly graphs to 
adjust existing binning results on contigs and to infer contigs shared by multiple species.""")

ap.add_argument("--contigs", required=True, help="path to the contigs file")
ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--paths", required=True, help="path to the contigs.paths file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='', help="prefix for the output file")
ap.add_argument("--depth", required=False, type=int, default=5, help="maximum depth for the breadth-first-search. [default: 5]")
ap.add_argument("--threshold", required=False, type=float, default=1.5, help="threshold for determining inconsistent vertices. [default: 1.5]")
ap.add_argument("--nthreads", required=False, type=int, default=8, help="number of threads to use. [default: 8]")

args = vars(ap.parse_args())

contigs_file = args["contigs"]
assembly_graph_file = args["graph"]
contig_paths = args["paths"]
contig_bins_file = args["binned"]
output_path = args["output"]
prefix = args["prefix"]
depth = args["depth"]
threshold = args["threshold"]
nthreads = args["nthreads"]

n_bins = 0

print("\nWelcome to GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs using Assembly Graphs.")
print("This version of GraphBin2 makes use of the assembly graph produced by SPAdes which is based on the de Bruijn graph approach.")

print("\nInput arguments:", contigs_file)
print("Assembly graph file:", assembly_graph_file)
print("Contig paths file:", contig_paths)
print("Existing binning output file:", contig_bins_file)
print("Final binning output file:", output_path)
print("Depth:", depth)
print("Threshold:", threshold)
print("Number of threads:", nthreads)

print("\nGraphBin2 started\n-------------------")

start_time = time.time()

# Get length and coverage of contigs
#--------------------------------------------------------

contig_lengths = {}
coverages = {}

for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
    start = 'NODE_'
    end = '_length'
    contig_num = int(re.search('%s(.*)%s' % (start, end), record.id).group(1))-1
    
    start = '_length_'
    end = '_cov'
    length = int(re.search('%s(.*)%s' % (start, end), record.id).group(1))
    
    start = '_cov_'
    end = ''
    coverage = int(float(re.search('%s(.*)%s' % (start, end), record.id).group(1)))
    
    contig_lengths[contig_num] = length
    coverages[contig_num] = coverage



# Build the assembly graph
#--------------------------------------------------------

paths = {}
segment_contigs = {}
node_count = 0

# Get contig paths from contigs.paths
with open(contig_paths) as file:
    name = file.readline()
    path = file.readline()
    
    while name != "" and path != "":
            
        while ";" in path:
            path = path[:-2]+","+file.readline()
        
        start = 'NODE_'
        end = '_length_'
        contig_num = str(int(re.search('%s(.*)%s' % (start, end), name).group(1))-1)
        
        segments = path.rstrip().split(",")
        
        if contig_num not in paths:
            node_count += 1
            paths[contig_num] = [segments[0], segments[-1]]
        
        for segment in segments:
            if segment not in segment_contigs:
                segment_contigs[segment] = set([contig_num])
            else:
                segment_contigs[segment].add(contig_num)
        
        name = file.readline()
        path = file.readline()

print("\nTotal number of contigs available:", node_count)

links = []
links_map = defaultdict(set)

# Get contig paths from contigs.paths
with open(assembly_graph_file) as file:
    line = file.readline()
    
    while line != "":
        
        # Identify lines with link information
        if "L" in line:
            strings = line.split("\t")
            f1, f2 = strings[1]+strings[2], strings[3]+strings[4]
            links_map[f1].add(f2)
            links_map[f2].add(f1)
            links.append(strings[1]+strings[2]+" "+strings[3]+strings[4])
        line = file.readline()
        

# Create graph
assembly_graph = Graph()

# Add vertices
assembly_graph.add_vertices(node_count)

# Get edges
edges = []

for i in range(len(assembly_graph.vs)):
    assembly_graph.vs[i]["id"]= i
    assembly_graph.vs[i]["label"]= "NODE_"+str(i)+"\nCov: "+str(coverages[i])+"\nLen: "+str(contig_lengths[i])
    
for i in range(len(paths)):
    segments = paths[str(i)]
    
    start = segments[0]
    start_rev = ""
    if start.endswith("+"):
        start_rev = start[:-1]+"-"
    else:
        start_rev = start[:-1]+"+"
        
    end = segments[1]
    end_rev = ""
    if end.endswith("+"):
        end_rev = end[:-1]+"-"
    else:
        end_rev = end[:-1]+"+"
    
    new_links = []
    
    if start in links_map:
        new_links.extend(list(links_map[start]))
    if start_rev in links_map:
        new_links.extend(list(links_map[start_rev]))
    if end in links_map:
        new_links.extend(list(links_map[end]))
    if end_rev in links_map:
        new_links.extend(list(links_map[end_rev]))
    
    for new_link in new_links:
        if new_link in segment_contigs:
            for contig in segment_contigs[new_link]:
                if i!=int(contig):
                    edges.append((i,int(contig)))
                    
assembly_graph.add_edges(edges)
assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)


# Get the number of bins from the initial binning result
#--------------------------------------------------------

try:
    all_bins_list = []

    with open(contig_bins_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            all_bins_list.append(row[1])
            
    bins_list = list(set(all_bins_list))
    bins_list.sort()

    n_bins = len(bins_list)
    print("Number of bins available in binning result:", n_bins)
except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin... Bye...!")
    sys.exit(1)


# Get initial binning result
#----------------------------

bins = [[] for x in range(n_bins)]

try:
    with open(contig_bins_file) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            bin_num = int(row[1])-1
            contig_num = int(row[0])
            bins[bin_num].append(contig_num)

except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin... Bye...!")
    sys.exit(1)


# Get binned and unbinned contigs
#-----------------------------------------------------

binned_contigs = []

for n in range(n_bins):
    binned_contigs = sorted(binned_contigs+bins[n])
    
unbinned_contigs = []

for i in range(node_count):
    if i not in binned_contigs:
        unbinned_contigs.append(i)

binned_contigs.sort()
unbinned_contigs.sort()

print("No. of binned contigs:", len(binned_contigs))
print("No. of unbinned contigs:", len(unbinned_contigs))


# Get isolated vertices and components without labels
#-----------------------------------------------------

isolated=[]

for i in range(node_count):
    
    neighbours = assembly_graph.neighbors(i, mode=ALL)
    
    if len(neighbours)==0:
        isolated.append(i)


non_isolated = []

for i in range(node_count):
    
    if i not in non_isolated and i in binned_contigs:

        component = []
        component.append(i)
        length = len(component)
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        for neighbor in neighbours:
            if neighbor not in component:
                component.append(neighbor)

        component = list(set(component))

        while length!= len(component):

            length = len(component)

            for j in component:

                neighbours = assembly_graph.neighbors(j, mode=ALL)

                for neighbor in neighbours:
                    if neighbor not in component:
                        component.append(neighbor)

        labelled = False
        for j in component:
            if j in binned_contigs:
                labelled = True
                break

        if labelled:
            for j in component:
                if j not in non_isolated:
                    non_isolated.append(j)

print("\nNumber of non-isolated contigs:", len(non_isolated))


# The BFS function to search labelled nodes
#-----------------------------------------------------

def runBFS(node, threhold=depth):
    queue = []
    visited = set()
    queue.append(node)
    depth = {}
    
    depth[node] = 0
    
    labelled_nodes = set()
    
    while (len(queue) > 0):
        active_node = queue.pop(0)
        visited.add(active_node)
        
        if active_node in binned_contigs and len(visited) > 1:
            
            contig_bin = -1

            # Get the bin of the current contig
            for n in range(n_bins):
                if active_node in bins[n]:
                    contig_bin = n
                    break
            
            labelled_nodes.add((node, active_node, contig_bin, depth[active_node], abs(coverages[node]-coverages[active_node])))
            
        else:
            for neighbour in assembly_graph.neighbors(active_node, mode=ALL):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
                    
    return labelled_nodes



# Remove labels of unsupported vertices
#-----------------------------------------------------

print("\nRemoving labels of unsupported vertices...")

iter_num = 1

while True:
    
    print("Iteration:", iter_num)
    
    remove_labels = {}

    for my_node in binned_contigs:

        if my_node in non_isolated:

            my_contig_bin = -1

            # Get the bin of the current contig
            for n in range(n_bins):
                if my_node in bins[n]:
                    my_contig_bin = n
                    break    

            BFS_labelled_nodes = list(runBFS(my_node))
            
            if len(BFS_labelled_nodes)>0:

                # Get the count of nodes in the closest_neighbours that belongs to each bin
                BFS_labelled_bin_counts = [0 for x in range(n_bins)]

                for i in range(len(BFS_labelled_nodes)):
                    BFS_labelled_bin_counts[BFS_labelled_nodes[i][2]] += 1
                
                zero_bin_count = 0

                # Count the number of bins which have no BFS_labelled_contigs
                for j in BFS_labelled_bin_counts:
                    if j == 0:
                        zero_bin_count += 1

                # Get the bin number which contains the maximum number of BFS_labelled_contigs
                max_index = BFS_labelled_bin_counts.index(max(BFS_labelled_bin_counts))
                    
                # If there are no BFS nodes of same label as contig, remove label
                if my_contig_bin!=-1 and BFS_labelled_bin_counts[my_contig_bin]==0:
                    remove_labels[my_node] = my_contig_bin
                
                # Check if all the BFS_labelled_contigs are in one bin
                elif zero_bin_count == (len(BFS_labelled_bin_counts)-1):

                    # If contig is not in the bin with maximum number of BFS_labelled_contigs 
                    if max_index!=my_contig_bin and BFS_labelled_bin_counts[max_index] > 1 and contig_lengths[my_node]<10000:
                        remove_labels[my_node] = my_contig_bin
            
    
    if len(remove_labels)==0:
        break
    else:
        
        for contig in remove_labels:
            bins[remove_labels[contig]].remove(contig)
            binned_contigs.remove(contig)
            unbinned_contigs.append(contig)
    
    iter_num += 1



# Refine labels of inconsistent vertices
#-----------------------------------------------------

print("\nRefining labels of inconsistent vertices...")

iter_num = 1

once_moved = []

while True:
    
    print("Iteration:", iter_num)
    
    contigs_to_correct = {}

    for my_node in binned_contigs:

        if my_node in non_isolated and my_node not in once_moved:

            my_contig_bin = -1

            # Get the bin of the current contig
            for n in range(n_bins):
                if my_node in bins[n]:
                    my_contig_bin = n
                    break    

            BFS_labelled_nodes = list(runBFS(my_node))

            # Get the count of nodes in the closest_neighbours that belongs to each bin
            BFS_labelled_bin_counts = [0 for x in range(n_bins)]

            for i in range(len(BFS_labelled_nodes)):
                BFS_labelled_bin_counts[BFS_labelled_nodes[i][2]] += 1

            zero_bin_count = 0

            # Count the number of bins which have no BFS_labelled_contigs
            for j in BFS_labelled_bin_counts:
                if j == 0:
                    zero_bin_count += 1
            
            # Get the bin number which contains the maximum number of BFS_labelled_contigs
            max_index = BFS_labelled_bin_counts.index(max(BFS_labelled_bin_counts))

            weighted_bin_count = [0 for x in range(n_bins)]

            for i in range(len(BFS_labelled_nodes)):
                
                path_length = BFS_labelled_nodes[i][3]
                weighted_bin_count[BFS_labelled_nodes[i][2]] += 1/(2**path_length)
            
            should_move = False
            
            max_weight = -1
            max_weight_bin = -1
            
            for i in range(len(weighted_bin_count)):
                if len(BFS_labelled_nodes)>0 and my_contig_bin!=-1 and i!=my_contig_bin and weighted_bin_count[i]>0 and weighted_bin_count[i] > weighted_bin_count[my_contig_bin]*threshold:
                    should_move = True
                    if max_weight < weighted_bin_count[i]:
                        max_weight = weighted_bin_count[i]
                        max_weight_bin = i
            
            if should_move and max_weight_bin!=-1:
                contigs_to_correct[my_node] = (my_contig_bin, max_weight_bin)
                once_moved.append(my_node)

    
    if len(contigs_to_correct)==0:
        break
    else:
        for contig in contigs_to_correct:
            old_bin = contigs_to_correct[contig][0]
            new_bin = contigs_to_correct[contig][1]
            bins[old_bin].remove(contig)
            bins[new_bin].append(contig)
            bins[new_bin].sort()
    
    iter_num += 1


# Propagate labels to unlabelled vertices
#-----------------------------------------------------

print("\nPropagating labels to unlabelled vertices...")

class DataWrap:
    def __init__(self, data):
        self.data = data
        
    def __lt__(self, other):
        return (self.data[3], self.data[-1])  < (other.data[3], other.data[-1]) 
    
contigs_to_bin = set()

for contig in binned_contigs:
    if contig in non_isolated:
        closest_neighbours = filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(contig, mode=ALL))
        contigs_to_bin.update(closest_neighbours)


sorted_node_list = []
sorted_node_list_ = [list(runBFS(x, threhold=depth)) for x in contigs_to_bin]
sorted_node_list_ = [item for sublist in sorted_node_list_ for item in sublist]

for data in sorted_node_list_:
    heapObj = DataWrap(data)
    heapq.heappush(sorted_node_list, heapObj)


while sorted_node_list:
    best_choice = heapq.heappop(sorted_node_list)    
    to_bin, binned, bin_, dist, cov_diff = best_choice.data
    
    
    if to_bin in unbinned_contigs:
        bins[bin_].append(to_bin)
        binned_contigs.append(to_bin)
        unbinned_contigs.remove(to_bin)
        
        # Discover to_bin's neighbours
        unbinned_neighbours = set(filter(lambda x: x not in binned_contigs, assembly_graph.neighbors(to_bin, mode=ALL)))
        sorted_node_list = list(filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list))
        heapq.heapify(sorted_node_list)
    
        for n in unbinned_neighbours:
            candidates = list(runBFS(n, threhold=depth))
            for c in candidates:
                heapq.heappush(sorted_node_list, DataWrap(c))



# Determine contigs belonging to multiple bins
#-----------------------------------------------------

print("\nDetermining multi-binned contigs...")

bin_cov_sum = [0 for x in range(n_bins)]
bin_contig_len_total = [0 for x in range(n_bins)]

for i in range(n_bins):
    for j in range(len(bins[i])):
        if bins[i][j] in non_isolated:
            bin_cov_sum[i] += coverages[bins[i][j]]*contig_lengths[bins[i][j]]
            bin_contig_len_total[i] += contig_lengths[bins[i][j]]

def is_multi(contig):
    if contig in non_isolated and contig in binned_contigs:
        
        contig_bin = -1

        # Get the bin of the current contig
        for n in range(n_bins):
            if contig in bins[n]:
                contig_bin = n
                break

        # Get average coverage of each connected component representing a bin excluding the contig
        bin_coverages = list(bin_cov_sum)
        bin_contig_lengths = list(bin_contig_len_total)

        bin_coverages[contig_bin] = bin_coverages[contig_bin] - (coverages[contig]*contig_lengths[contig])
        bin_contig_lengths[contig_bin] = bin_contig_lengths[contig_bin] - contig_lengths[contig]

        for i in range(n_bins):
            if bin_contig_lengths[i] != 0:
                bin_coverages[i] = bin_coverages[i]/bin_contig_lengths[i]

        # Get coverages of neighbours
        neighbour_bins = [[] for x in range(n_bins)]

        neighbour_bin_coverages = [[] for x in range(n_bins)]

        neighbours = assembly_graph.neighbors(contig, mode=ALL)

        for neighbour in neighbours:

            for n in range(n_bins):
                if neighbour in bins[n]:
                    neighbour_bins[n].append(neighbour)
                    neighbour_bin_coverages[n].append(coverages[neighbour])
                    break

        zero_bin_count = 0

        non_zero_bins = []

        # Count the number of bins which have no labelled neighbouring contigs
        for i in range(len(neighbour_bins)):
            if len(neighbour_bins[i]) == 0:
                zero_bin_count += 1
            else:
                non_zero_bins.append(i)

        if zero_bin_count < n_bins-1:

            bin_combinations = []

            for i in range(len(non_zero_bins)):
                bin_combinations += list(it.combinations(non_zero_bins, i+1))

            min_diff = sys.maxsize
            min_diff_combination = -1

            for combination in bin_combinations:

                comb_cov_total = 0

                for i in range(len(combination)):
                    comb_cov_total += bin_coverages[combination[i]]

                cov_diff = abs(comb_cov_total-coverages[contig])

                if cov_diff < min_diff:
                    min_diff = cov_diff
                    min_diff_combination = combination

            if min_diff_combination!=-1 and len(min_diff_combination) > 1 and contig_lengths[contig]>1000:
                # return True
                return contig, min_diff_combination

    return None

# Threads and multi-processing
p = Pool(nthreads)
mapped = p.map(is_multi, list(range(node_count)))
p.close()
multi_bins = list(filter(lambda x: x is not None, mapped))

# Add contigs to multiplt bins
for contig, min_diff_combination in multi_bins:
    print("Contig", contig, "belongs to bins", min_diff_combination)
    for mybin in min_diff_combination:
        if contig not in bins[mybin]:
            bins[mybin].append(contig)


# Determine elapsed time
elapsed_time = time.time() - start_time

# Print elapsed time for the process
print("\nElapsed time: ", elapsed_time, " seconds")

# Sort contigs in bins
for i in range(n_bins):
    bins[i].sort()


# Write result to output file
#-----------------------------------

output_bins = []

for i in range(node_count):
    for k in range(n_bins):
        if i in bins[k]:
            line = []
            line.append("NODE_"+str(i+1))
            line.append(k+1)
            output_bins.append(line)

output_file = output_path + prefix + 'graphbin2_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)

print("\nFinal binning results can be found at", output_file.name)


# Exit program
#-----------------------------------

print("\nThank you for using GraphBin2!\n")