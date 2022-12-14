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

# OTHER IMPORT
# [A]
import argparse
# [M]
import matplotlib
import matplotlib.pyplot as plt
# [N]
import networkx as nx

random.seed(9001)

__author__ = "ROUAUD Lucas"
__credits__ = __author__
__version__ = "1.0.0"
__maintainer__ = __author__
__email__ = "lucas.rouaud@gmail.com"
__status__ = "DONE"

matplotlib.use("Agg")


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = f"{path} is a directory."
        else:
            msg = f"{path} does not exist."
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
    """Return the standard deviation of a given list.
    """
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
    summation = 0

    for path_unit in subgraph:
        summation += path_unit[2]["weight"]

    return summation / len(subgraph)


def solve_bubble(graph, ancestor_node, descendant_node):
    """Solve a bubble in `networkx`.
    """
    w_list = []
    p_len = []

    p_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    for path in p_list:
        p_len += [len(path)]
        w_list += [path_average_weight(graph, path)]

    graph = select_best_path(graph, p_list, p_len, w_list)

    return graph


def simplify_bubbles(graph):
    """Simplify all bubble in a `networkx` network.
    """
    buble_to_clear = False

    for node in graph:
        node_n = node

        if graph.in_degree(node) >= 2:
            node_pred = []

            for node_unit in graph.predecessors(node):
                node_pred += [node_unit]

            if len(node_pred) <= 2:
                buble_to_clear = True
                break

        if buble_to_clear:
            break

    if not buble_to_clear:
        ancestor = nx.lowest_common_ancestor(graph, node_pred[0], node_pred[1])
        graph = solve_bubble(graph, ancestor, node_n)

    return graph


def solve_entry_tips(graph, starting_nodes):
    """Solve in entry.
    """
    n_list = []

    if len(starting_nodes) > 2:
        for n_i in starting_nodes:
            for n_j in starting_nodes:
                n_list += [(n_i, n_j)]
    else:
        n_list = [tuple(starting_nodes)]

    a_list = []

    for node in n_list:
        a_list += [nx.lowest_common_ancestor(graph.reverse(),
                                             node[0], node[1])]

    p_list = []
    w_list = []
    p_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, node[0], a_list[i]))[0]
        path_2 = list(nx.all_simple_paths(graph, node[1], a_list[i]))[0]

        p_list += [path_1, path_2]

        p_len += [len(path_1), len(path_2)]
        w_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, p_list, p_len, w_list,
                             delete_entry_node=True)

    return graph


def solve_out_tips(graph, ending_nodes):
    """Solve out entry.
    """
    n_list = []

    if len(ending_nodes) > 2:
        for n_i in ending_nodes:
            for n_j in ending_nodes:
                n_list += [(n_i, n_j)]
    else:
        n_list = [tuple(ending_nodes)]

    a_list = []

    for node in n_list:
        a_list += [nx.lowest_common_ancestor(graph, node[0], node[1])]

    p_list = []
    w_list = []
    p_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, a_list[i], node[0]))[0]
        path_2 = list(nx.all_simple_paths(graph, a_list[i], node[1]))[0]

        p_list += [path_1, path_2]

        p_len += [len(path_1), len(path_2)]
        w_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, p_list, p_len, w_list,
                             delete_sink_node=True)

    return graph


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True)
              if d["weight"] > 3]
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True)
              if d["weight"] <= 3]

    pos = nx.random_layout(graph)

    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color="b", style="dashed")

    plt.savefig(graphimg_file)


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


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt", encoding="utf-8") as save:
        pickle.dump(graph, save)


if __name__ == "__main__":
    # Each main variable start with a `m_`.
    m_args = get_arguments()

    # `k_mer` dic.
    m_kmer_dic = build_kmer_dict(m_args["fastq_file"], m_args["kmer_size"])

    # Create the diagram.
    m_diagram = build_graph(m_kmer_dic)

    # Cleaning.
    m_diagram = simplify_bubbles(m_diagram)

    m_starting_nodes = get_starting_nodes(m_diagram)
    m_diagram = solve_entry_tips(m_diagram, m_starting_nodes)

    m_ending_nodes = get_sink_nodes(m_diagram)
    m_diagram = solve_out_tips(m_diagram, m_ending_nodes)

    # Getting contigs.
    m_contigs = get_contigs(m_diagram, m_starting_nodes, m_ending_nodes)
    save_contigs(m_contigs, "test.fasta")

    if m_args["graphimg_file"] is not None:
        draw_graph(m_diagram, m_args["graphimg_file"])

    if m_args["output_file"] is not None:
        save_graph(m_diagram, m_args["output_file"])
