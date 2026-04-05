from Bio.Seq import Seq

query = Seq(input("Enter sequence: "))

#Result variables
ac=0
gc=0
tc=0
cc=0

#Loop for calculating ATGC
for base in query:
	if base=="A":
		ac+=1
	elif base=="G":
		gc+=1
	elif base=="T":
		tc+=1
	elif base=="C":
		cc+=1
	else:
		raise TypeError

print(ac)
print(cc)
print(gc)
print(tc)

