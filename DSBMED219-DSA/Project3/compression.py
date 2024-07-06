from collections import defaultdict
from Bio import SeqIO
import networkx as nx
import sys 
import random

from collections import defaultdict, deque

sys.setrecursionlimit(100000)

def read_fasta_file(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def generate_kmers(string, k):
    kmers = []
    for key, value in string.items():
        for i in range(len(value) - k + 1):
            kmers.append((value[i:i+k], key, i))
    return kmers

def generate_index_table(string, k):
    kmers = generate_kmers(string, k)
    index_table = defaultdict(lambda: {"read_indices": defaultdict(list), "count": 0})
    for kmer, read_id, index in kmers:
        index_table[kmer]["read_indices"][read_id].append(index)
        index_table[kmer]["count"] += 1
    
    return index_table

def filter_hash_table(hash_table, coverage):
    return {kmer: data for kmer, data in hash_table.items() if data["count"] >= coverage}

def adjust_counts(hash_table):
    for kmer, data in hash_table.items():
        if 4 <= data['count'] <= 20:
            data['count'] = 1
        elif data['count'] > 20:
            data['count'] = 2
    return hash_table

def build_de_bruijn_graph(kmers):
    G = nx.MultiDiGraph()
    
    for kmer, data in kmers.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        if prefix not in G:
            G.add_node(prefix, read_ids=set())
        if suffix not in G:
            G.add_node(suffix, read_ids=set())
        
        G.nodes[prefix]['read_ids'].update(data['read_indices'].keys())
        G.nodes[suffix]['read_ids'].update(data['read_indices'].keys())
        
        edge_data = {
            'read_ids': set(),
            'kmer': kmer,
            'weight': data['count']
        }
        G.add_edge(prefix, suffix, **edge_data)
                    
    return G


def find_start_end_nodes(graph):
    start_nodes = [node for node in graph.nodes if graph.in_degree(node) == 0]
    end_nodes = [node for node in graph.nodes if graph.out_degree(node) == 0]
    return start_nodes, end_nodes

def traverse_and_reconstruct(graph):
    start_nodes, end_nodes = find_start_end_nodes(graph)
    paths = []
    visited = set()

    for start_node in start_nodes:
        path = []
        current_node = start_node
        while current_node not in visited:
            path.append(current_node)
            visited.add(current_node)
            if graph.out_degree(current_node) == 0:
                break
            successors = list(graph.successors(current_node))
            if len(successors) == 0:
                break
            current_node = successors[0]
        paths.append(path)

    return paths

def print_reconstructed_genome(paths, graph):
    lists = []
    for path in paths:
        read_ids = set()
        for i in range(len(path) - 1):
            node1, node2 = path[i], path[i + 1]
            if graph.has_edge(node1, node2):  # Check if the edge exists
                edges = graph.get_edge_data(node1, node2)
                if edges:
                    for edge_data in edges.values():
                        read_ids.update(edge_data['read_ids'])
        lists.append(path)
    return lists


def main():
    reads = read_fasta_file("/Users/blakedickerson/Desktop/School/UCLA/Spring 2024/DSBMED219/Project3/Data/project2b_reads.fasta")
    # reads = read_fasta_file("/Users/blakedickerson/Desktop/School/UCLA/Spring 2024/DSBMED219/Project3/Data/project2_sample2_reads.fasta")
    hash_table = generate_index_table(reads, k=20)
    hash_table = filter_hash_table(hash_table, coverage=3)
    hash_table = adjust_counts(hash_table)
    G = build_de_bruijn_graph(hash_table)

    paths = traverse_and_reconstruct(G)
    test = print_reconstructed_genome(paths, G)
    
    max_val = test[0]
    for i in test:
        if len(max_val) < len(i):
            max_val = i
            
    reads_ = []
    for i in max_val:
        for key, val in reads.items():
            if i in val and key not in reads_:
                reads_.append(key)

    for i in reads_:
        print(f">{i}")
    
    remaining_reads = [key for key in reads.keys() if key not in reads_]
    random.shuffle(remaining_reads)
    for key in remaining_reads:
        print(f">{key}")
        

if __name__ == '__main__':
    main()
