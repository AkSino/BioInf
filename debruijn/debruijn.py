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

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import pylab
from Bio import SeqIO

__author__ = "Aurelien B"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Aurelien B"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurelien B"
__email__ = "burieaurel@eisti.eu"
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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
	with open(fastq_file,"rt") as monfic:
		for line in monfic:
			yield next(monfic)
			next(monfic)
			next(monfic)

def fonc(fastq_file):
	for i in read_fastq(fastq_file):
		print(i)


def cut_kmer(read, kmer_size):
	for i in range(len(read)-kmer_size):
		yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
	kmer_dict = {}
	list_kmer =[]
	sequence = read_fastq(fastq_file)
	kmer = cut_kmer(next(sequence), kmer_size)
	for i in kmer:
		list_kmer.append(i)
	for kmers in list_kmer:
		kmer_dict[kmers] = list_kmer.count(kmers)
	print(kmer_dict)
	return(kmer_dict)


def build_graph(kmer_dict,draw=False):
	G = nx.DiGraph()
	for i in kmer_dict:
		G.add_edge(i[:-1], i[1:], weight=kmer_dict[i])
	edge_labels=dict([((u,v,),d['weight'])
		             for u,v,d in G.edges(data=True)])
	red_edges = [(i[:-1], i[1:]) for i in kmer_dict if kmer_dict[i]>1]
	edge_colors = ['black' if not edge in red_edges else 'red' for edge in G.edges()]
	if draw:
		pos=nx.spring_layout(G)
		node_labels = {node:node for node in G.nodes()}
		nx.draw_networkx_labels(G, pos, labels=node_labels)
		nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)
		nx.draw(G,pos, node_size=1500,edge_color=edge_colors,edge_cmap=plt.cm.Reds)
		pylab.show()
	return G


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
	g = graph
	input_nodes = [u for u, deg in g.in_degree() if deg==0]
	return(input_nodes)

def get_sink_nodes(graph):
	g = graph
	output_nodes = [u for u, deg in g.out_degree() if deg==0]
	return(output_nodes)

def get_contigs(graph, starting_nodes, ending_nodes):
	list_contigs = []
	for i in starting_nodes:
		for j in ending_nodes:
			path = nx.shortest_path(graph,source=i,target=j)
			print(path)
			contig = path[0]
			for k in range(1,len(path)):
				contig += path[k][-1]
			length_contig = len(contig)
			list_contigs.append((contig,length_contig))
	return list_contigs
	
def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
	fichier = open(output_file, "a")
	i = 1
	for contigs in contigs_list:
		fichier.write(">contig_{} len={}\n".format(i,contigs[1]))
		i+=1
		fichier.write("{}\n".format(contigs[0]))
	fichier.close()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

#==============================================================
# Main program
#==============================================================
def main():
	"""
	Main program function
	"""
	# Get arguments
	args = get_arguments()
	dict_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
	graph = build_graph(dict_kmer)
	sink = get_sink_nodes(graph)
	starting = get_starting_nodes(graph)
	print(starting)
	print(sink)
	p = get_contigs(graph,starting,sink)
	print(p)
	save_contigs(p,args.output_file)
	
if __name__ == '__main__':
    main()
