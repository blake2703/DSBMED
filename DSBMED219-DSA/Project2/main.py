""" 
Goal: predict which read goes to which sample genome
"""

from Bio import SeqIO


def read_fasta_file(file_path):
    sequence = {}
    for i in SeqIO.parse(file_path, "fasta"):
        sequence[i.id] = str(i.seq)
    return sequence

def generate_kmers(string, k):
    kmers = []
    for i in range(len(string) - k + 1):
        kmers.append(string[i:i+k])
    return kmers

def generate_index_table(sequences, k):
    index_tables = []
    for sequence in sequences:
        index_table = {}
        for genome, genome_string in sequence.items():
            kmers = generate_kmers(genome_string, k)
            for i, kmer in enumerate(kmers):
                if kmer not in index_table:
                    index_table[kmer] = [(genome, i)]
                else:
                    index_table[kmer].append((genome, i))
        index_tables.append(index_table)
    return index_tables


NUM_SAMPLES = 100
sequences = []

for i in range(NUM_SAMPLES):
    sequence = read_fasta_file(f"/Users/blakedickerson/Desktop/School/UCLA/Spring 2024/DSBMED219/Project2/project1c_data/project1c_genome_{i}.fasta")
    sequences.append(sequence)

reads = read_fasta_file(file_path="/Users/blakedickerson/Desktop/School/UCLA/Spring 2024/DSBMED219/Project2/project1c_data/project1c_reads.fasta")

hash_tables = generate_index_table(sequences=sequences, k=10)

# Preprocess hash tables to build a dictionary that maps k-mers to values
kmer_dict = {}
for table in hash_tables:
    for kmer, value in table.items():
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [value]
        else:
            kmer_dict[kmer].append(value)

def line_up_reads(kmer_dict, reads):
    for read_name, read_sequence in reads.items():
        chunk_1 = read_sequence[:10]
        chunk_2 = read_sequence[10:20]
        chunk_3 = read_sequence[20:30]
        
        def count_genome_occurrences(chunk):
            genome_counts = {}
            if chunk in kmer_dict:
                for values in kmer_dict[chunk]:
                    for genome, count in values:
                        if genome in genome_counts:
                            genome_counts[genome] += count
                        else:
                            genome_counts[genome] = count
            return genome_counts
        
        genome_counts_1 = count_genome_occurrences(chunk_1)
        genome_counts_2 = count_genome_occurrences(chunk_2)
        genome_counts_3 = count_genome_occurrences(chunk_3)
        

        dictionaries = [genome_counts_1, genome_counts_2, genome_counts_3]

        counts = {}

        for dictionary in dictionaries:
            for key, value in dictionary.items():
                if key not in counts:
                    counts[key] = 1
                else:
                    counts[key] += 1

        # Figure out what to do if the chunks are not found??
        '''
        Maybe find the longest substring and do it this way???
        '''
        if len(counts) != 0:
            print(f">{read_name}\t{max(counts, key=counts.get) }")
        
        
line_up_reads(kmer_dict, reads)