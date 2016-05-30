#! /usr/local/bin/python
# -*- coding: utf-8 -*-

import argparse
import sys
from collections import defaultdict, namedtuple
import networkx as nx
from networkx import Graph

__author__ = "Sergey Aganezov"
__email__ = "aganezov@gwu.edu"

AssPoint = namedtuple("AssPoint", field_names=["scaf1_name", "scaf2_name", "scaf1_orient", "scaf2_orient"])
graphVertex = namedtuple("graphVertex", field_names=["scaf_name", "extremity"])
chainEntry = namedtuple("chainEntry", field_names=["scaf", "orient"])


def create_graph_vertices(scaffold1, orientation1, scaffold2, orientation2):
    v1 = graphVertex(scaffold1, "h" if orientation1 == "+" else "t")
    v2 = graphVertex(scaffold2, "t" if orientation2 == "+" else "h")
    return v1, v2


def compliment_vertex(vertex):
    return graphVertex(vertex.scaf_name, "t" if vertex.extremity == "h" else "h")


def construct_graph(assembly_points):
    result_graph = Graph()
    for ass_pair in assembly_points:
        vertex1, vertex2 = create_graph_vertices(scaffold1=ass_pair.scaf1_name, orientation1=ass_pair.scaf1_orient,
                                                 scaffold2=ass_pair.scaf2_name, orientation2=ass_pair.scaf2_orient)
        vertex1_comp, vertex2_comp = compliment_vertex(vertex1), compliment_vertex(vertex2)
        result_graph.add_edge(u=vertex1, v=vertex2)
        result_graph.add_edge(u=vertex1, v=vertex1_comp)
        result_graph.add_edge(u=vertex2, v=vertex2_comp)
    return result_graph


def get_chain_from_path(edge_path):
    result = []
    current_scaffold = None
    for graph_vertex in edge_path:
        if current_scaffold != graph_vertex.scaf_name:
            current_scaffold = graph_vertex.scaf_name
            orientation = "" if graph_vertex.extremity == "t" else "-"
            result.append(chainEntry(current_scaffold, orientation))
    return result


def remove_random_adjacency_edge(connected_component):
    node = connected_component.nodes()[0]
    other_scaf_extremity = [n for n in connected_component.neighbors(node) if n.scaf_name != node.scaf_name][0]
    connected_component.remove_edge(node, other_scaf_extremity)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transforming datasets form scaffolds assembled in pairs into GRIMM format")
    parser.add_argument("input", nargs="?", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("output", nargs="?", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--error_stream", type=argparse.FileType("wt"), default=sys.stderr)
    arguments = parser.parse_args()
    error = arguments.error_stream
    output = arguments.output
    per_species = defaultdict(list)
    try:
        next(arguments.input)
        for line in arguments.input:
            if len(line.strip()) == 0:
                continue
            ass_points = line.strip().split("\t")
            try:
                species_name, scaf1, scaf2, scaf1_or, scaf2_or = ass_points[0], ass_points[1], ass_points[2], ass_points[3], ass_points[4]
            except IndexError as er:
                error.write("WARNING: couldn't find first 5 entries in \\t-separated string: '{string}'. "
                            "Check that you use tabs vs spaces. Skipping...\n"
                            "".format(string=line.strip()))
                continue
            if scaf1_or == "?" or scaf2_or == "?":
                text = "WARNING: Ambiguity in relative scaffolds orientation:" + \
                       " in '{org_name}' scaffold '{scaf_name}' has orientation '?'. Skipping a pair '{scaf1}' -- '{scaf2}' ...\n"
                error.write(text.format(org_name=species_name, scaf_name=scaf1 if scaf1_or == "?" else scaf2, scaf1=scaf1, scaf2=scaf2))
                error.flush()
                continue
            per_species[species_name].append(AssPoint(scaf1, scaf2, scaf1_or, scaf2_or))
    except IOError:
        error.write("WARNING: Encountered uncaught IOException processing input. Quiting.")
        error.close()
        output.close()
        sys.exit(1)
    finally:
        arguments.input.close()

    for organism, ass_points in per_species.items():
        output.write(">{organism}\n".format(organism=organism))
        graph = construct_graph(assembly_points=ass_points)
        for cc in nx.connected_component_subgraphs(graph):
            circular = False
            try:
                long_circular = all(map(lambda node: cc.degree(node) == 2, cc))
                if long_circular:
                    remove_random_adjacency_edge(cc)
                circular = long_circular or len(cc.nodes()) == 2
                source, sink = filter(lambda node: cc.degree(node) == 1, cc)
            except ValueError as ex:
                error.write("WARNING: Issues processing following scaffolds: ")
                err_scaf = [gv.scaf_name for gv in filter(lambda node: cc.degree(node) not in [1, 2], cc)]
                error.write("{scafs} in '{organism}';"
                            " these scaffolds are adjacent by the same extremity to more, that one other scaffold.\n"
                            "".format(scafs=err_scaf, organism=organism))
                error.flush()
                continue
            path = list(nx.all_simple_paths(G=cc, source=source, target=sink))[0]
            chain = get_chain_from_path(edge_path=path)
            chain = [chainEntry(ce.scaf, "-" if ce.orient == "" else "") for ce in chain[::-1]] if chain[0].orient == "-" else chain
            output.write(" ".join(map(lambda chain_entry: chain_entry.orient + chain_entry.scaf, chain)) +
                         (" @" if circular else " $") +
                         "\n")
            output.flush()
    output.close()
