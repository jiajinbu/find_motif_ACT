# find_motif_ACT

Read a genome fasta and find motif in all chromsome.

# Usage

python motif_finder.py genome.fa [0|1] output_file_pre
[0|1] mean you chose 0 to find ACTN4ACT, chose 1 to find ACTN6ACT.

# Example
 
python motif_finder.py tair10_chr_all.fas 0 motif_ACTN4ACT
python motif_finder.py tair10_chr_all.fas 1 motif_ACTN6ACT

If output_file_pre was not given, it will be set as "motif_finder".