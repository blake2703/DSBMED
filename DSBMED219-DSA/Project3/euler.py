from collections import defaultdict, deque
from Bio import SeqIO

def read_fasta_file(file_path):
    """
    A function that parses and reads a fasta file

    Args:
        file_path (str): location of where the fasta file is locally.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences.
    """
    sequence = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequence[record.id] = str(record.seq)
    return sequence

def extract_kmers(sequence, k, read_id):
    """
    Extracts all k-mers from a given sequence with the associated read ID
    
    Args:
        sequence (str): current spectrum sequence
        k (int): size of kmer
        read_id (str): id of sequence
        
    Returns:
        list: A list of tuples
    """
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append((kmer, read_id))
    return kmers

def build_de_bruijn_graph(kmers):
    """Builds a de Bruijn graph from a list of k-mers
    
    Args:
        kmers: A list of kmers generated from the spectrum
        
    Returns:
        dict: a de bruijn graph"""
    graph = defaultdict(list)
    for kmer, read_id in kmers:
        prefix, suffix = kmer[:-1], kmer[1:]
        graph[prefix].append((suffix, read_id))
        if suffix not in graph:
            graph[suffix] = []
    return graph

def find_eulerian_path(graph):
    """Finds an Eulerian path in a given graph
    
    Args:
        graph (dict): a de briujn graph
        
    Returns:
        list: the Euler path of the de bruijn graph
        
    """
    graph = {node: deque(edges) for node, edges in graph.items()}
    
    in_degree, out_degree = defaultdict(int), defaultdict(int)
    for node in graph:
        out_degree[node] += len(graph[node])
        for neighbor, _ in graph[node]:
            in_degree[neighbor] += 1
    
    start_node = None
    for node in set(in_degree.keys()).union(out_degree.keys()):
        if out_degree[node] > in_degree[node]:
            start_node = node
            break
    if start_node is None:
        start_node = next(iter(graph))
    
    stack = [start_node]
    path = []

    while stack:
        curr = stack[-1]
        if graph[curr]:
            next_node, read_id = graph[curr].popleft()
            stack.append(next_node)
            path.append(read_id)
        else:
            stack.pop()
    
    return path

file_path = "/Users/blakedickerson/Desktop/School/UCLA/Spring 2024/DSBMED219/Project3/Data/project2a_spectrum.fasta"
sequences = read_fasta_file(file_path)

k = 20
kmers = []
for read_id, seq in sequences.items():
    kmers.extend(extract_kmers(seq, k, read_id))

graph = build_de_bruijn_graph(kmers)

reads = find_eulerian_path(graph)

reads = [read_id for read_id in reads if read_id is not None]
for read_id in reads:
    print(f">{read_id}")
