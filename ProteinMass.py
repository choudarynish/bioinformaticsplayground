from Bio.SeqUtils.ProtParam import ProteinAnalysis

query = input("Enter protein sequence: ")

molwt = ProteinAnalysis(query , monoisotopic=True).molecular_weight()
print(f"{molwt}")




