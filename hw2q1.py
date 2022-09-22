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

def phred33_to_q(qual):
  """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
  return ord(qual)-33

def q_to_phred33(Q):
  """ Turn Phred-scaled integer into Phred+33 ASCII-encoded quality """
  return chr(Q + 33)


from io import StringIO
import sys
import numpy as np

# get command line inputs
program_name = sys.argv[0]
arguments = sys.argv[1:]
count = len(arguments)

reads_file = arguments[0]
summary_file = arguments[1]

reads_file = open(reads_file, "r")
reads = parse_fastq(reads_file)
reads_file.close()
num_positions = len(reads[0][1]) # assume that all reads are the same length
num_reads = len(reads[:][0])

# TRY WITH EMPTY ARRAY
summary = np.zeros((num_positions, 5), dtype = 'uint8')

# calculate summary stats for each read
for r in range(num_reads): 
    name, seq, qual = reads[r]
    q_string = list(map(phred33_to_q, qual))

    for p in range(num_positions):

        if r == 0: 
           summary[p][0] = q_string[p]
           summary[p][1] = q_string[p]
        else: 
            # min quality 
            if q_string[p] < summary[p][0]: 
                summary[p][0] = q_string[p]
        
            # max quality 
            if q_string[p] > summary[p][1]: 
                summary[p][1] = q_string[p]

        # num qualities < 10
        if q_string[p] < 10: 
            summary[p][2] += 1

        # num qualities >= 30
        if q_string[p] >= 30: 
            summary[p][3] += 1

        # num chars not ACTG
        if seq[p] not in ['A', 'C', 'T', 'G']:
            summary[p][4] += 1

np.savetxt(summary_file, summary, fmt = '%0d', delimiter=' ')

# iterate through each read
# get the line of q for each read
# that gives you the string, iterate through the string and fill in the summary array


