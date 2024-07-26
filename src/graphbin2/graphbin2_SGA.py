#!/usr/bin/env python3

import concurrent.futures
import csv
import heapq
import itertools as it
import logging
import re
import sys
import time

from Bio import SeqIO
from igraph import *
from tqdm import tqdm

from .bidirectionalmap.bidirectionalmap import BidirectionalMap


def run(args):
    # Get arguments
    # ---------------------------------------------------

    contigs_file = args.contigs
    assembly_graph_file = args.graph
    abundance_file = args.abundance
    contig_bins_file = args.binned
    output_path = args.output
    prefix = args.prefix
    depth = args.depth
    threshold = args.threshold
    delimiter = args.delimiter
    nthreads = args.nthreads

    n_bins = 0

    # Setup logger
    # -----------------------
    logger = logging.getLogger("GraphBin2")
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    # ---------------------------------------------------

    fileHandler = logging.FileHandler(output_path + "/" + prefix + "graphbin2.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to GraphBin2: Refined and Overlapped Binning of Metagenomic Contigs using Assembly Graphs."
    )
    logger.info(
        "This version of GraphBin2 makes use of the assembly graph produced by SGA which is based on the string graph approach."
    )

    logger.info("Input arguments:")
    logger.info("Contigs file: " + contigs_file)
    logger.info("Assembly graph file: " + assembly_graph_file)
    logger.info("Existing binning output file: " + contig_bins_file)
    logger.info("Final binning output file: " + output_path)
    logger.info("Depth: " + str(depth))
    logger.info("Threshold: " + str(threshold))
    logger.info("Number of threads: " + str(nthreads))

    logger.info("GraphBin2 started")

    start_time = time.time()

    # Get length and coverage of contigs
    # --------------------------------------------------------

    contig_lengths = {}

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        start_n = "contig-"
        end_n = ""

        contig_num = int(re.search("%s(.*)%s" % (start_n, end_n), record.id).group(1))

        length = len(record.seq)
        contig_lengths[contig_num] = length

    coverages = {}

    with open(abundance_file, "r") as my_file:
        line = my_file.readline()

        while line != "":
            strings = line.split("\t")

            start_n = "contig-"
            end_n = ""

            contig_num = int(
                re.search("%s(.*)%s" % (start_n, end_n), strings[0]).group(1)
            )

            coverages[contig_num] = int(strings[1])
            line = my_file.readline()

    # Get the links from the .asqg file
    # -----------------------------------

    links = []

    contig_names = {}

    my_map = BidirectionalMap()

    node_count = 0

    try:
        # Get contig connections from .asqg file
        with open(assembly_graph_file) as file:
            line = file.readline()

            while line != "":
                # Count the number of contigs
                if "VT" in line:
                    start = "contig-"
                    end = ""
                    contig_num = int(
                        re.search(
                            "%s(.*)%s" % (start, end), str(line.split()[1])
                        ).group(1)
                    )
                    my_map[node_count] = contig_num
                    contig_names[node_count] = line.split()[1].strip()
                    node_count += 1

                # Identify lines with link information
                elif "ED" in line:
                    link = []
                    strings = line.split("\t")[1].split()
                    link.append(int(strings[0][7:]))
                    link.append(int(strings[1][7:]))
                    links.append(link)
                line = file.readline()

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin2... Bye...!")
        sys.exit(1)

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    logger.info("Total number of contigs available: " + str(node_count))

    ## Construct the assembly graph
    # -------------------------------

    try:
        # Create the graph
        assembly_graph = Graph()

        # Create list of edges
        edge_list = []

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Name vertices
        for i in range(len(assembly_graph.vs)):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(i)

        # Iterate links
        for link in links:
            # Remove self loops
            if link[0] != link[1]:
                # Add edge to list of edges
                edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

        # Add edges to the graph
        assembly_graph.add_edges(edge_list)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin2... Bye...!")
        sys.exit(1)

    logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

    # Get the number of bins from the initial binning result
    # --------------------------------------------------------

    try:
        all_bins_list = []

        with open(contig_bins_file) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=delimiter)
            for row in readCSV:
                all_bins_list.append(row[1])

        bins_list = list(set(all_bins_list))
        bins_list.sort()

        n_bins = len(bins_list)
        logger.info("Number of bins available in binning result: " + str(n_bins))
    
    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the binning result file is provided and it is having the correct format"
        )
        logger.info("Exiting GraphBin2... Bye...!")
        sys.exit(1)

    # Get initial binning result
    # ----------------------------

    bins = [[] for x in range(n_bins)]

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                start = "contig-"
                end = ""
                contig_num = contigs_map_rev[
                    int(re.search("%s(.*)%s" % (start, end), row[0]).group(1))
                ]

                bin_num = int(row[1]) - 1
                bins[bin_num].append(contig_num)

        for i in range(n_bins):
            bins[i].sort()

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format."
        )
        logger.info("Exiting GraphBin2... Bye...!")
        sys.exit(1)

    # Get binned and unbinned contigs
    # -----------------------------------------------------

    binned_contigs = []

    for n in range(n_bins):
        binned_contigs = sorted(binned_contigs + bins[n])

    unbinned_contigs = []

    for i in range(node_count):
        if i not in binned_contigs:
            unbinned_contigs.append(i)

    binned_contigs.sort()
    unbinned_contigs.sort()

    logger.info("No. of binned contigs: " + str(len(binned_contigs)))
    logger.info("No. of unbinned contigs: " + str(len(unbinned_contigs)))

    # Get isolated vertices
    # -----------------------------------------------------

    isolated = []

    for i in range(node_count):
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        if len(neighbours) == 0:
            isolated.append(i)

    # The BFS function to search labelled nodes
    # -----------------------------------------------------

    def runBFS(node, threhold=depth):
        queue = []
        visited = set()
        queue.append(node)
        depth = {}

        depth[node] = 0

        labelled_nodes = set()

        while len(queue) > 0:
            active_node = queue.pop(0)
            visited.add(active_node)

            if active_node in binned_contigs and len(visited) > 1:
                contig_bin = -1

                # Get the bin of the current contig
                for n in range(n_bins):
                    if active_node in bins[n]:
                        contig_bin = n
                        break

                labelled_nodes.add(
                    (
                        node,
                        active_node,
                        contig_bin,
                        depth[active_node],
                        abs(
                            coverages[contigs_map[node]]
                            - coverages[contigs_map[active_node]]
                        ),
                    )
                )

            else:
                for neighbour in assembly_graph.neighbors(active_node, mode=ALL):
                    if neighbour not in visited:
                        depth[neighbour] = depth[active_node] + 1
                        if depth[neighbour] > threhold:
                            continue
                        queue.append(neighbour)

        return labelled_nodes

    # Remove labels of unsupported vertices
    # -----------------------------------------------------

    logger.info("Removing labels of unsupported vertices")

    iter_num = 1

    while True:
        logger.debug("Iteration: " + str(iter_num))

        remove_labels = {}

        # Initialise progress bar
        pbar = tqdm(total=len(binned_contigs))

        for my_node in binned_contigs:
            if my_node not in isolated:
                my_contig_bin = -1

                # Get the bin of the current contig
                for n in range(n_bins):
                    if my_node in bins[n]:
                        my_contig_bin = n
                        break

                BFS_labelled_nodes = list(runBFS(my_node))

                if len(BFS_labelled_nodes) > 0:
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
                    max_index = BFS_labelled_bin_counts.index(
                        max(BFS_labelled_bin_counts)
                    )

                    # If there are no BFS nodes of same label as contig, remove label
                    if (
                        my_contig_bin != -1
                        and BFS_labelled_bin_counts[my_contig_bin] == 0
                    ):
                        remove_labels[my_node] = my_contig_bin

                    # Check if all the BFS_labelled_contigs are in one bin
                    elif zero_bin_count == (len(BFS_labelled_bin_counts) - 1):
                        # If contig is not in the bin with maximum number of BFS_labelled_contigs
                        if (
                            max_index != my_contig_bin
                            and BFS_labelled_bin_counts[max_index] > 1
                            and contig_lengths[contigs_map[my_node]] < 10000
                        ):
                            remove_labels[my_node] = my_contig_bin

            # Update progress bar
            pbar.update(1)

        # Close progress bar
        pbar.close()

        if len(remove_labels) == 0:
            break
        else:
            for contig in remove_labels:
                bins[remove_labels[contig]].remove(contig)
                binned_contigs.remove(contig)
                unbinned_contigs.append(contig)

        iter_num += 1

    # Refine labels of inconsistent vertices
    # -----------------------------------------------------

    logger.info("Refining labels of inconsistent vertices")

    iter_num = 1

    once_moved = []

    while True:
        logger.debug("Iteration: " + str(iter_num))

        contigs_to_correct = {}

        # Initialise progress bar
        pbar = tqdm(total=len(binned_contigs))

        for my_node in binned_contigs:
            if my_node not in isolated and my_node not in once_moved:
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
                    weighted_bin_count[BFS_labelled_nodes[i][2]] += 1 / (2**path_length)

                should_move = False

                max_weight = -1
                max_weight_bin = -1

                for i in range(len(weighted_bin_count)):
                    if (
                        len(BFS_labelled_nodes) > 0
                        and my_contig_bin != -1
                        and i != my_contig_bin
                        and weighted_bin_count[i] > 0
                        and weighted_bin_count[i]
                        > weighted_bin_count[my_contig_bin] * threshold
                    ):
                        should_move = True
                        if max_weight < weighted_bin_count[i]:
                            max_weight = weighted_bin_count[i]
                            max_weight_bin = i

                if should_move and max_weight_bin != -1:
                    contigs_to_correct[my_node] = (my_contig_bin, max_weight_bin)
                    once_moved.append(my_node)

            # Update progress bar
            pbar.update(1)

        # Close progress bar
        pbar.close()

        if len(contigs_to_correct) == 0:
            break
        else:
            for contig in contigs_to_correct:
                old_bin = contigs_to_correct[contig][0]
                new_bin = contigs_to_correct[contig][1]
                bins[old_bin].remove(contig)
                bins[new_bin].append(contig)
                bins[new_bin].sort()

        iter_num += 1

    # Get non isolated contigs

    logger.info("Obtaining non isolated contigs")

    # Initialise progress bar
    pbar = tqdm(total=node_count)

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

            while length != len(component):
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

        # Update progress bar
        pbar.update(1)

    # Close progress bar
    pbar.close()

    logger.info("Number of non-isolated contigs: " + str(len(non_isolated)))

    non_isolated_unbinned = list(set(non_isolated).intersection(set(unbinned_contigs)))

    logger.info(
        "Number of non-isolated unbinned contigs: " + str(len(non_isolated_unbinned))
    )

    # Propagate labels to unlabelled vertices
    # -----------------------------------------------------

    logger.info("Propagating labels to unlabelled vertices")

    # Initialise progress bar
    pbar = tqdm(total=len(non_isolated_unbinned))

    class DataWrap:
        def __init__(self, data):
            self.data = data

        def __lt__(self, other):
            return (self.data[3], self.data[-1]) < (other.data[3], other.data[-1])

    contigs_to_bin = set()

    for contig in binned_contigs:
        if contig in non_isolated:
            closest_neighbours = filter(
                lambda x: x not in binned_contigs,
                assembly_graph.neighbors(contig, mode=ALL),
            )
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

        if to_bin in non_isolated_unbinned:
            bins[bin_].append(to_bin)
            binned_contigs.append(to_bin)
            non_isolated_unbinned.remove(to_bin)
            unbinned_contigs.remove(to_bin)

            # Update progress bar
            pbar.update(1)

            # Discover to_bin's neighbours
            unbinned_neighbours = set(
                filter(
                    lambda x: x not in binned_contigs,
                    assembly_graph.neighbors(to_bin, mode=ALL),
                )
            )
            sorted_node_list = list(
                filter(lambda x: x.data[0] not in unbinned_neighbours, sorted_node_list)
            )
            heapq.heapify(sorted_node_list)

            for n in unbinned_neighbours:
                candidates = list(runBFS(n, threhold=depth))
                for c in candidates:
                    heapq.heappush(sorted_node_list, DataWrap(c))

    # Close progress bar
    pbar.close()

    # Determine contigs belonging to multiple bins
    # -----------------------------------------------------

    logger.info("Determining multi-binned contigs")

    bin_cov_sum = [0 for x in range(n_bins)]
    bin_contig_len_total = [0 for x in range(n_bins)]

    for i in range(n_bins):
        for j in range(len(bins[i])):
            if bins[i][j] in non_isolated:
                bin_cov_sum[i] += (
                    coverages[contigs_map[bins[i][j]]]
                    * contig_lengths[contigs_map[bins[i][j]]]
                )
                bin_contig_len_total[i] += contig_lengths[contigs_map[bins[i][j]]]

    mapped = [None for itr in range(node_count)]

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=nthreads)

    # Thread function
    def thread_function(
        n,
        non_isolated,
        binned_contigs,
        n_bins,
        bins,
        bin_cov_sum,
        bin_contig_len_total,
        coverages,
        contigs_map,
        contig_lengths,
        assembly_graph,
    ):
        bin_result = is_multi(
            contig=n,
            non_isolated=non_isolated,
            binned_contigs=binned_contigs,
            n_bins=n_bins,
            bins=bins,
            bin_cov_sum=bin_cov_sum,
            bin_contig_len_total=bin_contig_len_total,
            coverages=coverages,
            contigs_map=contigs_map,
            contig_lengths=contig_lengths,
            assembly_graph=assembly_graph,
        )
        mapped[n] = bin_result

    # Set up execution args for thread function
    exec_args = []

    for n in range(node_count):
        exec_args.append(
            (
                n,
                non_isolated,
                binned_contigs,
                n_bins,
                bins,
                bin_cov_sum,
                bin_contig_len_total,
                coverages,
                contigs_map,
                contig_lengths,
                assembly_graph,
            )
        )

    # Thread executor
    for itr in tqdm(
        executor.map(lambda p: thread_function(*p), exec_args), total=node_count
    ):
        pass

    executor.shutdown(wait=True)

    multi_bins = list(filter(lambda x: x is not None, mapped))

    if len(multi_bins) == 0:
        logger.info("No multi-labelled contigs were found")
    else:
        logger.info("Found " + str(len(multi_bins)) + " multi-labelled contigs ==>")

    # Add contigs to multiplt bins
    for contig, min_diff_combination in multi_bins:
        logger.info(
            contig_names[contig]
            + " belongs to bins "
            + ", ".join(str(s + 1) for s in min_diff_combination)
        )
        for mybin in min_diff_combination:
            if contig not in bins[mybin]:
                bins[mybin].append(contig)

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Show elapsed time for the process
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    # Sort contigs in bins
    for i in range(n_bins):
        bins[i].sort()

    # Write result to output file
    # -----------------------------------

    output_bins = []

    for i in range(node_count):
        for k in range(n_bins):
            if i in bins[k]:
                line = []
                line.append(contig_names[i])
                line.append(k + 1)
                output_bins.append(line)

    output_file = output_path + prefix + "graphbin2_output.csv"

    with open(output_file, mode="w") as output_file:
        output_writer = csv.writer(
            output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
        )

        for row in output_bins:
            output_writer.writerow(row)

    logger.info("Final binning results can be found at " + str(output_file.name))

    # Exit program
    # -----------------------------------

    logger.info("Thank you for using GraphBin2!")


def is_multi(
    contig,
    non_isolated,
    binned_contigs,
    n_bins,
    bins,
    bin_cov_sum,
    bin_contig_len_total,
    coverages,
    contigs_map,
    contig_lengths,
    assembly_graph,
):
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

        bin_coverages[contig_bin] = bin_coverages[contig_bin] - (
            coverages[contigs_map[contig]] * contig_lengths[contigs_map[contig]]
        )
        bin_contig_lengths[contig_bin] = (
            bin_contig_lengths[contig_bin] - contig_lengths[contigs_map[contig]]
        )

        for i in range(n_bins):
            if bin_contig_lengths[i] != 0:
                bin_coverages[i] = bin_coverages[i] / bin_contig_lengths[i]

        # Get coverages of neighbours
        neighbour_bins = [[] for x in range(n_bins)]

        neighbour_bin_coverages = [[] for x in range(n_bins)]

        neighbours = assembly_graph.neighbors(contig, mode=ALL)

        for neighbour in neighbours:
            for n in range(n_bins):
                if neighbour in bins[n]:
                    neighbour_bins[n].append(neighbour)
                    neighbour_bin_coverages[n].append(coverages[contigs_map[neighbour]])
                    break

        zero_bin_count = 0

        non_zero_bins = []

        # Count the number of bins which have no labelled neighbouring contigs
        for i in range(len(neighbour_bins)):
            if len(neighbour_bins[i]) == 0:
                zero_bin_count += 1
            else:
                non_zero_bins.append(i)

        if zero_bin_count < n_bins - 1:
            bin_combinations = []

            for i in range(len(non_zero_bins)):
                bin_combinations += list(it.combinations(non_zero_bins, i + 1))

            min_diff = sys.maxsize
            min_diff_combination = -1

            for combination in bin_combinations:
                comb_cov_total = 0

                for i in range(len(combination)):
                    comb_cov_total += bin_coverages[combination[i]]

                cov_diff = abs(comb_cov_total - coverages[contigs_map[contig]])

                if cov_diff < min_diff:
                    min_diff = cov_diff
                    min_diff_combination = combination

            if (
                min_diff_combination != -1
                and len(min_diff_combination) > 1
                and contig_lengths[contigs_map[contig]] > 1000
            ):
                # return True
                return contig, min_diff_combination

    return None


def main(args):
    run(args)


if __name__ == "__main__":
    main()
