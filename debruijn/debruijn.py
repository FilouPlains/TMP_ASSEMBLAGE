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
    with open(fastq_file, "r") as file:
        seq = []

        for line in file:
            if line[0] in ["A", "T", "C", "G"]:
                yield line.strip()


def cut_kmer(read, kmer_size):
    length_read = len(read)
    kmer = []

    for i in range(length_read):
        if i + kmer_size <= length_read:
            yield read[i:i + kmer_size]
        else:
            break


def build_kmer_dict(fastq_file, kmer_size):
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
    pass


def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass


def path_average_weight(graph, path):
    pass


def solve_bubble(graph, ancestor_node, descendant_node):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    node_list = []

    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            node_list += [node]

    return node_list


def get_sink_nodes(graph):
    node_list = []

    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            node_list += [node]

    return node_list


def get_contigs(graph, starting_nodes, ending_nodes):
    pass


def save_contigs(contigs_list, output_file):
    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
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
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


if __name__ == '__main__':
    args = get_arguments()

    kmer_dic = build_kmer_dict(args["fastq_file"], 5)

    diagram = build_graph(kmer_dic)
    
    starting_nodes = get_starting_nodes(diagram)
    ending_nodes = get_sink_nodes(diagram)
    contigs = get_contigs(diagram, starting_nodes, ending_nodes)

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
