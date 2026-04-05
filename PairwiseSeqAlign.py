#Pairwise Sequence Alignment
#Global and Local
#change global to local for local alignment

from Bio import pairwise2 #for pairwise alignment
from Bio.Align import substitution_matrices #for importing scoring matrices
alignments = pairwise2.align.localxx("GAATTC" , "GCCTTA")
#alignments = pairwise2.align.globalxx("seq1" , "seq2") first x is for match and second is for gap penalty , defailt 1 and 0 respectively
print("Default scoring and gap penalties")
for alignment in alignments: 
    print(pairwise2.format_alignment(*alignment))

alignments = pairwise2.align.localmx("GAATTC" , "GCCTTA" , match=2 , mismatch=-1)
#alignments = pairwise2.align.globalmx("seq1" , "seq2") first m is for custom score and x is for default gap penalty , set score for m , default x is 0 for gap penalty
print("Custom matching score and default gap penalty")
for alignment in alignments:
    print(pairwise2.format_alignment(*alignment))

alignments = pairwise2.align.localxs("GAATTC" , "GCCTTA" , open=-2 , extend=-1)
#alignments = pairwise2.align.globalxs("seq1" , "seq2") first x is for default pairing score and s is for custom gap penalty , set score for s , default x is 0 for gap penalty
print("Default matching score and custom gap penalty")
for alignment in alignments:
    print(pairwise2.format_alignment(*alignment))

matrix = substitution_matrices.load("BLOSUM50")
alignments = pairwise2.align.localdx("GAATTC" , "GCCTTA" , match_dict=matrix)
#alignments = pairwise2.align.global("seq1" , "seq2") first d is for custom pairing score from matrix and x is for default gap penalty , set matrix for d , default x is 0 for gap penalty
print("Matching score using substitution matrix and default gap penalty")
for alignment in alignments:
    print(pairwise2.format_alignment(*alignment))

#The first alignmnt is the best alignment coming off the loop. To show only that access thr outpot list using alignments[0] instead of for loop




#Using Bio.Align.PairwiseAligner
from Bio import Align

aligner = Align.PairwiseAligner
#Setting scores is optional, else will be in default
aligner.match_score = 1.0 # default 1
aligner.mismatch_score = -1.0 # default 0 
aligner.open_gap_score = -2.0 # default 0
aligner.extend_gap_score = -0.5 # default 0
#aligner.substitution_matrix = "BLOSUM62" #IF REQUIRED INSTEAD OF SCORES
aligner.mode = "local"  # or "local"

alignments = aligner.align("GAATTC" , "GCCTTA")
for alignment in alignments:
    print(alignment)
