import sys, re

'''
Read a genome fasta and find motif in all chromsome.

Usage: python motif_finder.py genome.fa [0|1] output_file_pre
[0|1] mean you chose 0 to find ACTN4ACT, chose 1 to find ACTN6ACT.

Example: 
python motif_finder.py tair10_chr_all.fas 0 motif_ACTN4ACT
python motif_finder.py tair10_chr_all.fas 1 motif_ACTN6ACT

If output_file_pre was not given, it will be set as "motif_finder".
'''

def read_fasta(file_fasta):
    '''
    read fasta file and iter generate (seqid, seq).
    Usage: for seqid, seq in read_fasta(filein):
    '''
    
    def describ2id(seq_describ):
        return seq_describ.split()[0]
        
    seq_id = ""
    seq_describ = ""
    seq = ""
    for line_num, l in enumerate(open(file_fasta)):
        l = l.rstrip()
        if l.startswith(">"):
            if seq_id:
                if not seq:
                    raise Exception("Line %s: No sequence for '%s'!" % (line_num, seq_id))
                yield( (seq_id, seq) )
            seq = ""
            seq_describ = l[1:]
            seq_id = describ2id(seq_describ)
        else:
            seq += l
    if seq_id:
        if not seq:
            raise Exception("Line %s: No sequence for '%s'!" % (line_num, seq_id))
        yield( (seq_id, seq) )
    

def find_motif(filein, pattern1, pattern2):
    
    def revcom(seq):
        d = {"A":"T","C":"G","T":"A","G":"C"}
        return ''.join([d[s] for s in seq[::-1]])
        
    data = {}
    for seqid, seq in read_fasta(filein):
        chr_data = []
        for match in pattern1.finditer(seq):
            chr_data.append([match.start()+1, "+", match.group()])
        for match in pattern2.finditer(seq):
            chr_data.append([match.start()+1, "-", revcom(match.group())])
        chr_data.sort(key = lambda x: x[0])
        data[seqid] = [seq, chr_data]
    return data

def write_motif(data, output_pre):
    
    summary_file = output_pre + ".summary.txt"
    with open(summary_file, 'w') as sumo:
        sumo.write("chr_name\tchr_length\tmotif_num\n")
        for seqid, (seq, chr_data) in data.items():
            with open('{}.{}.txt'.format(output_pre, seqid), 'w') as o:
                o.write("chr_name\tstart_pos\tstrand\tmotif_seq\n")
                for start, strand, s in chr_data:
                    start = str(start)
                    o.write('\t'.join([seqid, start, strand, s]) + "\n")
                sumo.write("\t".join([seqid, str(len(seq)), str(len(chr_data))]) + "\n")

def main():
    
    filein = sys.argv[1]
    try:
        motif = int(sys.argv[2])
    except:
        motif = 0
    try:
        output_file_pre = sys.argv[3]
    except:
        output_file_pre = "motif_finder"

    pattern1 = re.compile("ACT[ACTG][ACTG][ACTG][ACTG]ACT")
    pattern2 = re.compile("AGT[ACTG][ACTG][ACTG][ACTG]AGT")
    #for ACTN6ACT
    pattern3 = re.compile("ACT[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]ACT")
    pattern4 = re.compile("AGT[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]AGT")
    
    if motif:
        motif_data = find_motif(filein, pattern3, pattern4)
    else:
        motif_data = find_motif(filein, pattern1, pattern2)
        
    write_motif(motif_data, output_file_pre)
    
if __name__ == "__main__":
    main()   
    