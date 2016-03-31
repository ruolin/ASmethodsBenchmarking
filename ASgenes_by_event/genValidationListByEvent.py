import sys

all_AS_fh = open(sys.argv[1],'r')
partial_AS_fh = open(sys.argv[2], 'r')

all_AS_set = set()
partial_AS_set = set()

with all_AS_fh as fh:
	for line in fh:
		col = line.rstrip().split()
		all_AS_set.add(col[0])

with partial_AS_fh as fh:
	for line in fh:
		col = line.rstrip().split()
		if col[0] in all_AS_set:
			print col[0] + "\t" + col[1] + "\t" + "AS"
		else:
			print col[0] + "\t" + col[1] + "\t" + "NonAS"
 
