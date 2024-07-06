from Bio import SeqIO
from collections import Counter

reference_genome_path = "/Users/blakedickerson/Downloads/project1b-s/project1b-s_reference_genome.fasta"
single_read_path = "/Users/blakedickerson/Downloads/project1b-s/project1b-s_with_error_paired_reads.fasta"


def read_fasta_file(file_path):
    sequence = {}
    for i in SeqIO.parse(file_path, "fasta"):
        sequence[i.id] = str(i.seq)
    return sequence

reference_genome = read_fasta_file(file_path=reference_genome_path)
single_read_error = read_fasta_file(file_path=single_read_path)

"""
Project 1A wants us to solve for substitutions only, which means we can just go through the genome and find out at which index there
is a substitution.

substitution: a type of mutation where one nucleotide is replaced by another nucleotide -> ACGAAT: GCGAAT
"""



"""
Step 1. Decide on a split for the read -> example is 3

Step 2. Create an index for the reference genome where the key is a n-length string and the values are 
a list of starting positions the n-length string occur in the reference genome -> hash map
    - If the length of a read if 30 and our split is 3 then we will have a n-length string of 10 as the
    keys in our hashmap

Step 3. Loop over all the reads
    - 3.1 Take the first n-length substring of the current read and compare it with the keys in the hashmap
        - 3.1.2 Once the key is found line up the whole read with all starting positions at that key
        - 3.1.3 If any reads are below the threshold (perfect match or only a few character differences)
        then we align the read at that starting_index:starting_index + length of entire read
        - 3.1.4 Else skip and go to the next position
        
    - 3.2 Take the second n-length substring of current read and compare it with the keys in the hashmap
        - 3.1.2 Once the key is found line up the n-length substring with all the starting positions at
        that key
        - 3.1.3 Take the first n-length substring and line up with previous n-length characters at the 
        starting position
        - 3.1.4 Take the last n-length substring and line up with next n-length characters at the
        starting position -> starting_index - n_length:starting_index + n_length
        - 3.1.5 If any reads are below threshold then align read
        - 3.1.6 Else skip and go to next position
    
    - 3.3 Same idea as above, but take the previous n-length * 2 and compare since this is the end

Step 4. Now we know where all of our reads go we can store each read in an aligned position in comparison
with the reference genome
    - AAAAGTGTCG
       AAAG
         AGTG
      AAA

Step 5. Now that all the reads are aligned we can now look for indels (substitutions, deletions, insertions)
    - 5.1 Substitutions
        - We can say our coverage is 3 since we have a split of 3 from our original read, therefore we make a 
        prediction if 2 out of 3 reads have an error in the sampe place (subject to change)
"""



# step 1 generate a kmer hash table for every 10mer in the genome
def generate_kmers(string, k):
    kmers = []
    for i in range(len(string) - k + 1):
        kmers.append(string[i:i+k])
    return kmers

def generate_index_table(string, k):
    kmers = generate_kmers(string, k)
    index_table = {}
    for i, kmer in enumerate(kmers):
        if kmer not in index_table:
            index_table[kmer] = [i]  # Initialize with a list containing the first position
        else:
            index_table[kmer].append(i)  # Append the position to the existing list
    return index_table

hash_table = generate_index_table(string=reference_genome['genome'],
                                  k=10)

recoveredSequence = dict()

def compareReadToGenome(hash_table, n_length=20):    
    def compareChunkToGenome(chunk,read_sequence):
        # first chunk
        if chunk in hash_table:
            start_indicies = hash_table[chunk]  #get all possible start indicies 

            for start_index in start_indicies:
                count = 0
                # check if substring matches big string
                genomeSubsequence = reference_genome["genome"][start_index: start_index + len(read_sequence)]
                
                if read_sequence == genomeSubsequence:
                    # update a string at the correct alignment 
                    recoveredSequence[read_sequence] = (start_index, start_index + len(read_sequence))
                
                elif len(read_sequence) == len(genomeSubsequence):
                    for i in range(len(read_sequence)):
                        if read_sequence[i] != genomeSubsequence[i]:
                            count += 1
                    if count < 3: #arbitrary value
                        recoveredSequence[read_sequence] = (start_index, start_index + len(read_sequence))
                else:
                    pass
    
    for read_name, read_sequence in single_read_error.items():
        ch1 = read_sequence[:10] 
        compareChunkToGenome(ch1, read_sequence)


compareReadToGenome(hash_table, n_length=20)  
recoveredSequence = sorted(recoveredSequence.items(), key=lambda x: x[1][0])


bad_positions = []
for key, (start, end) in recoveredSequence:
    curr = reference_genome["genome"][start:end]
    for index, value in enumerate(curr):
        if curr[index] != key[index]:
            bad_positions.append((key, start+index, end))

second_value_counts = {}
for t in bad_positions:
    second_value = t[1]
    if second_value in second_value_counts:
        second_value_counts[second_value] += 1
    else:
        second_value_counts[second_value] = 1

# Filter the list to keep only those tuples with a second value that has more than one occurrence
filtered_tuples = [t for t in bad_positions if second_value_counts[t[1]] > 1]

seen_second_values = {}
unique_by_second_value = []
for t in filtered_tuples:
    second_value = t[1]
    if second_value not in seen_second_values:
        seen_second_values[second_value] = True
        unique_by_second_value.append(t)


for item in unique_by_second_value:
    key = item[0]
    index = item[1]
    end = item[2]
    end_pos = end - index
    
    print(f">S{index} {reference_genome['genome'][index]} {key[-end_pos]}")


# second is to align my reads
recoveredSequence = dict()

def compareReadToGenome(hash_table, n_length=20):
    insertion_positions = []
    deletion_positions = [] 
   
    # third is to break read into thirds
    for read_name, read_sequence in single_read_error.items():
        ch1 = read_sequence[:10] 
        ch2 = read_sequence[10: 10*2] 
        ch3 = read_sequence[-10:]
        
        
        # fourth is to match up the first and third chunk to check to see if there is indel in middle
        if ch1 in hash_table and ch3 in hash_table:
            # fifth grab positions of ch1 and ch3
            curr_ch1 = hash_table[ch1]
            curr_ch3 = hash_table[ch3]
            # print(len(read_sequence), curr_ch1, curr_ch3)
            if abs(curr_ch1[0] - curr_ch3[0]) == 39:
                curr_index = 0
                # insertion_positions.append((curr_ch1[0], curr_ch3[0]))
                curr_slice = reference_genome["genome"][curr_ch1[0]:curr_ch3[0]]
                
                # print(comparator, curr_slice)\
                for i in range(len(curr_slice)):
                    if read_sequence[i] != curr_slice[i]:
                        
                        # insert a dash at position in reference genome
                        curr_slice = curr_slice[:i] + '-' + curr_slice[i:]
                        curr_index = curr_ch1[0] + i
                
                dash_count = 0
                for i in range(len(curr_slice)):
                    if curr_slice[i] == '-':
                        dash_count += 1
                
                if dash_count == 1:
                    for i in range(len(curr_slice)):
                        if curr_slice[i] == '-':
                            insertion_positions.append((read_sequence[i], curr_index-1))
                    
                    
            elif abs(curr_ch1[0] - curr_ch3[0]) == 41:
                curr_index = 0
                curr_slice = reference_genome["genome"][curr_ch1[0]:curr_ch3[0]]
                
                for i in range(len(curr_slice)):
                    if curr_slice[i] != read_sequence[i]:
                        read_sequence = read_sequence[:i] + '-' + read_sequence[i:]
                        curr_index = curr_ch1[0] + i
                
                dash_count = 0
                for i in range(len(read_sequence)):
                    if read_sequence[i] == '-':
                        dash_count += 1
                
                if dash_count == 1:
                    for i in range(len(read_sequence)):
                        if read_sequence[i] == '-':
                            deletion_positions.append((curr_slice[i], curr_index-1))

    return insertion_positions, deletion_positions
        
insertions, deletions = compareReadToGenome(hash_table, n_length=20)


# Set to store unique tuples
unique_tuples_del = set()
for item in deletions:
    unique_tuples_del.add(item)
dels = list(unique_tuples_del)

unique_tuples_ins = set()
for item in insertions:
    unique_tuples_ins.add(item)
ins = list(unique_tuples_ins)

for i in dels:
    nucleotide = i[0]
    index = i[1]
    print(f">D{index} {nucleotide}")

for i in ins:
    nucleotide = i[0]
    index = i[1]
    print(f">I{index} {nucleotide}")


