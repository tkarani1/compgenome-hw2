def parse_fastq(fh):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    reads = []
    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        reads.append((name, seq, qual))
    return reads

def parse_fasta(fh):
    """ Parse reads from a FASTA filehandle.  We
        return one string of the entire dna text """
    
    fh.readline() # ignore first line
    T = ''
    while True:
        line = fh.readline()
        if len(line) == 0:
            break  # end of file
        line_str = line.rstrip()
        T += line_str
    return T

# calculate hamming
def hamming_dist(p, t, maxHammingDistance=1):
    """ Compares phrase and text and returns hamming distance. Assumes p and t are same length """
    h_dist = 0
    for i in range(len(t)): 
        if p[i] != t[i]:
            h_dist += 1
        if (h_dist > maxHammingDistance):
            h_dist = -1
            break 
    return h_dist

def make_hamtable(p, keys, maxHammingDistance): 
    p_hamtable = {}
    for key in keys:
        d = hamming_dist(p, key, maxHammingDistance)
        if (d > -1):
            p_hamtable[key] = d
    return p_hamtable

def phred33_to_q(qual):
  """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
  return ord(qual)-33

def q_to_phred33(Q):
  """ Turn Phred-scaled integer into Phred+33 ASCII-encoded quality """
  return chr(Q + 33)

import sys
import numpy as np

# get command line inputs
program_name = sys.argv[0]
arguments = sys.argv[1:]
count = len(arguments)

fasta_file = arguments[0]
fastq_file = arguments[1]
output_file = arguments[2]

# extract info from files
fasta_file = open(fasta_file, "r")
fastq_file = open(fastq_file, "r")
T = parse_fasta(fasta_file)
reads = parse_fastq(fastq_file)
fasta_file.close()
fastq_file.close()

# create k-mer hashmap for T
T_len = len(T)
k = 6
k_mer_table = {T[0:k] : [0]}

for c in range(1, T_len - k + 1):
    substr = T[c:c+k]
    if (k_mer_table.get(substr)): 
        k_mer_table.get(substr).append(c)
    else: 
        k_mer_table[substr] = [c]

# create summary 
num_reads = len(reads)
summary = np.empty((num_reads, 5), dtype='uint8')
matches = np.empty((num_reads, 5), dtype='object')
matches = [[[] for x in range(5)] for y in range(num_reads)] 

for r in range(num_reads):
        # 5 offsets
    for i in range(5): 
        query = reads[r][1][(i * 6):(i * 6 + 6)]
        if (k_mer_table.get(query)): 
            hits = k_mer_table[query]
            summary[r][i] = len(hits)
        else: 
            hits = None
            summary[r][i] = 0
    
    # create hamming dist tables for every phrase offset
    p1 = reads[r][1][0:6]
    p2 = reads[r][1][6:12]
    p3 = reads[r][1][12:18]
    p4 = reads[r][1][18:24]
    p5 = reads[r][1][24:30]

    k_mer_keys = k_mer_table.keys()

    p1_hamtable = make_hamtable(p1, k_mer_keys, 4)
    p2_hamtable = make_hamtable(p2, k_mer_keys, 4)
    p3_hamtable = make_hamtable(p3, k_mer_keys, 4)
    p4_hamtable = make_hamtable(p4, k_mer_keys, 4)
    p5_hamtable = make_hamtable(p5, k_mer_keys, 4)


    p1_keys = p1_hamtable.keys()
    for key in p1_keys: 
        key_offsets = k_mer_table[key]
        for o in key_offsets: 
            ham_sum = p1_hamtable[key]

            # ANOTHER OPTION: check if any of the keys for p2 have offsets at the next 6
            if p2_hamtable.get(T[o+6*1:o+6*2]) is None:
                break
            else: 
                ham_sum += p2_hamtable.get(T[o+6*1:o+6*2])

            if p3_hamtable.get(T[o+6*2:o+6*3]) is None: 
                break
            else: 
                ham_sum += p3_hamtable.get(T[o+6*2:o+6*3])
            
            if p4_hamtable.get(T[o+6*3:o+6*4]) is None: 
                break
            else: 
                ham_sum += p4_hamtable.get(T[o+6*3:o+6*4])
            
            if p5_hamtable.get(T[o+6*4:o+6*5]) is None: 
                break
            else: 
                ham_sum += p5_hamtable.get(T[o+6*4:o+6*5])

            if ham_sum <= 4:
                matches[r][ham_sum].append(o)

# create table 'variants' to hold information on all of the positions of T
# 0: char in T, 
# 1: number of exact matches, 
# 2: number of mismatches, 
# 3: list of match touples (read, nt, qual), 
# 4: list of mismatch touples (read, nt, qual)
# 5: total weight

variants = np.empty((T_len, 6), dtype='object')  
variants = [[[] for x in [3,4]] for y in range(T_len)]
variants = [[0 for x in [1,2,5]] for y in range(T_len)]
variants = [['' for x in [0]] for y in range(T_len)]   
greater_than_20_idx = []
greatest_weight = ()
second_greatest_weight = ()

# print(reads[0])
# print(reads[0][1][4])

for r in range(num_reads):
    for m in range(5): # diff types of mismatches
        if len(matches[r][m]) == 0: 
            continue
        
        print(matches[r][m])

        for i in matches[r][m]: 
            
            T_nt = T[i + 0]
            read_nt = reads[r][1][i + 0]
            read_qual = phred33_to_q(reads[r][2][i + 0])
            variants[i + 0][0] = T_nt

            print(T_nt)
            print(read_nt)
            print(read_qual)

            if (T_nt == read_nt):  # bases match
                print(variants[i + 0][1])
                variants[i + 0][1] += 1
                variants[i + 0][3].append((r, read_nt, read_qual))
            else:  # bases don't match
                variants[i + 0][2] += 1
                variants[i + 0][4].append((r, read_nt, read_qual))
            
            variants[i + 0][5] += read_qual
            idx_weight = variants[i + 0][5]

            if idx_weight > 20: 
                greater_than_20_idx.append(i)

                if idx_weight > second_greatest_weight[2] & idx_weight <= greatest_weight[2]:
                    greatest_weight = (second_greatest_weight[0], second_greatest_weight[1], second_greatest_weight[2])
                    second_greatest_weight = (i + 0, T_nt, idx_weight)
                
                elif idx_weight > greatest_weight[2]: 
                    second_greatest_weight = (greatest_weight[0], greatest_weight[1], greatest_weight[2])
                    greatest_weight = (i + 0, T_nt, idx_weight)









