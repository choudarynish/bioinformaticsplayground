from Bio.Seq import Seq
query = Seq(input("Enter RNA Sequence: "))

query_aa = Seq(query).translate(to_stop=True)
print("\n")
print(query_aa)
