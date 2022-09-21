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
summary = np.zeros((num_reads, 5), dtype = 'uint8')

for r in range(num_reads): 
    query = reads[r][1][-6:]
    print(query)

