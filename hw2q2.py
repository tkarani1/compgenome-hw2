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

def parse_fasta_1(fh):  # WHY DOESN'T THIS IGNORE NEW LINES
    """ Parse reads from a FASTA filehandle.  We
        return one string of the entire dna text """
    
    fh.readline() # ignore first line
    T =  fh.read().rstrip()
    return T

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

for c in range(1, T_len - k):
    substr = T[c:c+k]
    if (k_mer_table.get(substr)): 
        k_mer_table.get(substr).append(c)
    else: 
        k_mer_table[substr] = [c]

# create summary 
num_reads = len(reads)
summary = np.empty((num_reads, 5), dtype='object')
all_matches = np.empty((num_reads, 1), dtype = 'object')
all_fruitless = np.empty((num_reads, 1), dtype = 'object')

print(T)
print(k_mer_table)

for r in range(num_reads): 
    query = reads[r][1][0:k]
    
    # (a) The number of index hits, i.e. the total number of times its leftmost 6-mer occurs in the genome T
    if (k_mer_table.get(query)): 
        hits = k_mer_table[query]
        summary[r][0] = len(hits)
    else: 
        hits = None
        summary[r][0] = 0
    
    # determine matches and fruitless hits 
    
    matches = []
    fruitless = []
    for h in hits:  
        P_str = reads[r][1][k:]
        T_str = T[h+k:h+k+len(P_str)]
        if P_str == T_str:
            matches.append(h)
        else: 
            fruitless.append(h)

    summary[r][1] = len(matches)
    summary[r][2] = len(fruitless)
    all_matches[r][0] = matches
    all_fruitless[r][0] = fruitless
    summary[r][3] = matches
    summary[r][4] = fruitless

print(summary)
np.savetxt(output_file, summary, fmt = '%0d %0d %0d %s %s', delimiter=' ')   
# np.savetxt(output_file, summary[:][0:2], fmt = ('%0d %0d %0d'), delimiter=' ')   

