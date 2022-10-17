"""Perform assembly based on debruijn graph."""

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

# ===============
# FROM * IMPORT *
# ===============
# [O]
from operator import itemgetter
# [R]
from random import randint

# ===============
# IMPORT *
# ===============
# STANDARD IMPORT
# [P]
import pickle
import os
# [R]
import random
# [S]
import statistics
import sys

# OTHER IMPORT
# [A]
import argparse
# [N]
import networkx as nx
# [M]
import matplotlib
import matplotlib.pyplot as plt

random.seed(9001)
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", dest="fastq_file", type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument("-k", dest="kmer_size", type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument("-o", dest="output_file", type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument("-f", dest="graphimg_file", type=str,
                        help="Save graph as image (png)")
    return vars(parser.parse_args())


def read_fastq(fastq_file):
    """Read a fastq file.
    """
    with open(fastq_file, "r", encoding="utf-8") as file:
        for line in file:
            if line[0] in ["A", "T", "C", "G"]:
                yield line.strip()


def cut_kmer(read, kmer_size):
    """Cut a sequence in "k" mers.
    """
    length_read = len(read)

    for i in range(length_read):
        if i + kmer_size <= length_read:
            yield read[i:i + kmer_size]
        else:
            break


def build_kmer_dict(fastq_file, kmer_size):
    """Input a fastq file, output a kmer_dic.
    """
    kmer_dic = {}
    whole_seq = ""
    kmer_list = []

    fastq_seq = read_fastq(fastq_file)

    for seq in fastq_seq:
        whole_seq += seq + " "
        kmer_list += [cut_kmer(seq, kmer_size)]

    whole_seq = whole_seq.strip()

    for kmer in kmer_list:
        for kmer_unit in kmer:
            if kmer_unit not in kmer_dic:
                kmer_dic[kmer_unit] = 0

            kmer_dic[kmer_unit] += 1

    return kmer_dic


def build_graph(kmer_dict):
    """Build a `networkx` graph.
    """
    diagram = nx.DiGraph()
    kmer_keys = list(kmer_dict.keys())

    for i, key in enumerate(kmer_keys):
        if i + 1 > len(kmer_keys):
            break
        else:
            left = key[:-1]
            right = key[1:]

            diagram.add_edge(left, right, weight=kmer_dict[key])

    return diagram


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a given path from the graph.
    """
    left = 1 - int(delete_entry_node)

    for path in path_list:
        if delete_sink_node:
            graph.remove_nodes_from(path[left:])
        else:
            graph.remove_nodes_from(path[left:-1])

    return graph


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best past between a given input path list.
    """
    w_max = max(weight_avg_list)
    index = []

    for i, j in enumerate(weight_avg_list):
        if j == w_max:
            index += [i]

    if len(index) > 1:
        l_max = max(path_length)
        shift_index = []

        for i in index:
            if path_length[i] == l_max:
                shift_index += [i]

        best_i = random.choice(shift_index)
    else:
        best_i = index[0]

    cp_path = path_list
    cp_path.pop(best_i)

    return remove_paths(graph, cp_path, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    """Return average weight of a given path.
    """
    subgraph = graph.subgraph(path).edges(data=True)
    sum = 0

    for path_unit in subgraph:
        sum += path_unit[2]["weight"]

    return sum / len(subgraph)


def solve_bubble(graph, ancestor_node, descendant_node):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    """Return all starting nodes.
    """
    node_list = []

    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            node_list += [node]

    return node_list


def get_sink_nodes(graph):
    """Return all ending nodes.
    """
    node_list = []

    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            node_list += [node]

    return node_list


def get_contigs(graph, starting_nodes, ending_nodes):
    """Give all possible contigs.
    """
    list_contigs = []

    for start in starting_nodes:
        for end in ending_nodes:
            path = nx.all_simple_paths(graph, start, end)

            for path_unit in path:
                contig = "".join(path_unit[::len(path_unit[0])])
                list_contigs += [(contig, len(contig))]

    return list_contigs


def save_contigs(contigs_list, output_file):
    """Save the contigs.
    """
    with open(output_file, "w", encoding="utf-8") as file:
        for i, sequence in enumerate(contigs_list):
            file.write(f">contig_{i} len={sequence[1]}\n")
            file.write(fill(sequence[0]))
            file.write("\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format.
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)

    # Draw the graph with networkx++
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt", encoding="utf-8") as save:
        pickle.dump(graph, save)


if __name__ == "__main__":
    # Each main variable start with a `m_`.
    m_args = get_arguments()

    m_kmer_dic = build_kmer_dict(m_args["fastq_file"], m_args["kmer_size"])

    m_diagram = build_graph(m_kmer_dic)

    m_starting_nodes = get_starting_nodes(m_diagram)
    m_ending_nodes = get_sink_nodes(m_diagram)

    m_contigs = get_contigs(m_diagram, m_starting_nodes, m_ending_nodes)

    save_contigs(m_contigs, "test.fasta")

    m_path_list = nx.all_simple_paths(m_diagram, m_starting_nodes[0],
                                      m_ending_nodes[0])

    m_length = [96, 69, 43, 50, 57, 30]
    m_mean = [96 / 5, 69 / 5, 43 / 5, 50 / 5, 57 / 5, 30 / 5]

    print(m_diagram)

    select_best_path(m_diagram, m_path_list, m_length, m_mean, True, False)

    print(m_diagram)

    # position = nx.spring_layout(diagram, seed=9001)

    # nodes = nx.draw_networkx_nodes(
    #     diagram,
    #     position
    # )
    # edges = nx.draw_networkx_edges(
    #     diagram,
    #     position
    # )

    # pc = matplotlib.collections.PatchCollection(edges, cmap=plt.cm.plasma)

    # ax = plt.gca()
    # ax.set_axis_off()
    # plt.colorbar(pc, ax=ax)
    # plt.show()
