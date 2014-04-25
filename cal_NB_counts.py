## install rpy2 first
## sudo easy_install rpy

import sys, argparse, math, random
import numpy as np
from sets import Set
from scipy.stats import nbinom
from rpy2.robjects import r
import rpy2.robjects as robjects

try:
	import HTSeq
except ImportError:
	sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
	sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
	sys.exit(1)

parser = argparse.ArgumentParser(
	description = "Generated negative binomial relationship gene-level read counts for a synthetic groups based on real experimental data",
	epilog =  "Copy right: Ruolin Liu, ISU"    
	)

parser.add_argument("gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-g1", "--group1", nargs='+', type=argparse.FileType('r'), required=True, help='First group sam files separated by space.')
parser.add_argument("-g2", "--group2", nargs='+', type=argparse.FileType('r'), required=True, help='Second group sam files separated by space.')
parser.add_argument("-n", "--num-reps", type=int, default = 3, dest='nreps',
	help = "Number of replicates. Default is 3." )
parser.add_argument("-l", "--num-target-gene", type=int, default = 3, dest='ntarg',
	help = "Number of AS genes. Default is 2000. ")
parser.add_argument("-m", "--mode", choices = ['AS-genes', 'all-genes'], default = 'AS-genes',
	help = "Choose between AS-genes or all-genes: AS-genes simulates annotated AS genes only. all_genes simulates all genes in annotation.")

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

GENE = Set(["gene", "transposable_element_gene","pseudogene"])
EXON = Set(["exon", "pseudogenic_exon"])
## parse inputs
args = parser.parse_args()
name_gtf = args.gene_model
NREPS=args.nreps
MODE=args.mode
NTARG = args.ntarg

## define global variables
def countSam(_sam_file, _genes, _dic, _idx):
	'''
	This function process a single sam file and count the reads sitting on each gene
	Results are save in the global variable.
	'''
	num_reads = 0
	for alnmt in _sam_file:
		if alnmt.aligned: 
			intersection_set = None
			#print list(_genes[alnmt.iv].steps())
			for iv1, step_set in _genes[alnmt.iv].steps():
				if intersection_set is None:
					intersection_set = step_set.copy()
				else:
					intersection_set.intersection_update(step_set)
			if len(intersection_set) == 1:
				_dic[ list(intersection_set)[0] ][_idx] += 1
		num_reads += 1
		if num_reads % 100000 == 0:
			sys.stderr.write("%d reads processed.\n" % num_reads)
	return _dic

def meanVar(_files, _gff_file , _output):


	NFILE=len(_files)
	if NFILE == 1:
		sys.stderr.write("Need at least two samples for each group.\n")
		sys.exit(1)
	#####

	_dict_counts = dict() ## dictionary of gene counts
	_genes = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	idx=0
	count = 0
	if MODE == "all-genes":
		for feature in _gff_file:
			if feature.type in GENE:
				_dict_counts[ feature.name ] = [0]*NFILE
				_genes[feature.iv] += feature.name
	if MODE == "AS-genes":
		## Bug: Does not report last gene in gff if it has at least two transcript
		transcript= set()
		cur_line = None
		for feature in _gff_file:
			if feature.type in GENE:
				if len(transcript) >1:
					_dict_counts[ cur_line.name ] = [0]*NFILE
					_genes[cur_line.iv] += cur_line.name
					count +=1
				cur_line = feature
				transcript.clear()
			if feature.type in EXON:
				transcript.add(feature.attr["Parent"])
	print "number of genes", count
	_file_raw_count = open(_output+'.rawcounts','w')
	_file_nb_count = open(_output+'.nbcounts','w')
	## This loop read through the input list and call countSam for each input file  
	for f in _files:
		sam_file=HTSeq.SAM_Reader(f)
		_dict_counts=countSam(sam_file, _genes,_dict_counts, idx)
		f.close()
		idx += 1
		sys.stderr.write("library %d has generated.\n" % idx)
	## Print raw counts in file specified by <out>
	for key, value in sorted(_dict_counts.iteritems()):
		_file_raw_count.write(key+"\t"+"\t".join(map(str,value))+"\n")
	_file_raw_count.close()
	## calculate group mean and variance
	list_mean = list()
	list_var = list()
	for key, value in sorted(_dict_counts.iteritems()):
		list_mean.append(np.mean(np.array(value)))
		list_var.append(np.var(np.array(value)))
	
	## computer loess esimates	
	## The following code is using rpy2 module
	a = robjects.FloatVector(list_mean)
	b = robjects.FloatVector(list_var)
	df = robjects.DataFrame({"mean": a, "var": b})
	non0_df=df.rx(df.rx2("mean").ro > 0, True) ## subsetting if mean > 0
	loess_fit = r.loess("var ~ mean", data=non0_df, degree=2)
	'''
	#good-of-fit test:
	variance=r.predict(loess_fit, 1000)
	print variance[0]
	print (1000*1000)/(variance[0]-1000)
	'''
	var_pred = r.predict(loess_fit, a)
	# This loop overwrite global variable dict_counts for recoding new count data
	count_idx = 0
	for key, value in sorted(_dict_counts.iteritems()):
		n = math.pow(list_mean[count_idx],2)/(var_pred[count_idx]-list_mean[count_idx])
		n = int(n) # n: number of failures
		if n<=0:
			_dict_counts[key] = [0]*NREPS
		else:
			p = n/float(n+list_mean[count_idx]) # p: prob of success
			_dict_counts[key] = nbinom.rvs(n, p, size=NREPS).tolist()
		count_idx += 1
	#var_pred = r.predict(loess_fit, a)
	for key, value in sorted(_dict_counts.iteritems()):
		_file_nb_count.write(key+"\t"+"\t".join(map(str,value))+"\n")
	_file_nb_count.close()
	_file_raw_count.close()
	return _dict_counts
def main():
	## first group
	group1_f = args.group1
	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)
	
	##### Sanity Check
	samfile=HTSeq.SAM_Reader(group1_f[0])
	is_chr_sam = None
	set_chr_gff = set()
	for almnt in samfile:
		is_chr_sam = almnt.iv.chrom
		break
	for feature in gff_file:
		set_chr_gff.add(feature.iv.chrom)
	
	if  is_chr_sam not in set_chr_gff:
		sys.stderr.write("Error: Chromosome id in SAM files and GFF file does not agree!!\n")
		sys.exit(1)
	#####

	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)	
	counts1=meanVar(group1_f, gff_file,'group1')	
	file_gtf.close()
	## second group
	group2_f = args.group2
	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)
	counts2=meanVar(group2_f, gff_file,'group2')
	
	file_gtf.close()
	merged_count = {k1: v1+v2 for (k1,v1) in counts1.iteritems() for (k2,v2) in counts2.iteritems() if k1 == k2}
	non0_exp_list = list()
	for k,v in sorted(merged_count.iteritems()):
		#print k,v
		value=np.array(v)
		if len(value[value ==0]) == 0:
			non0_exp_list.append(k)
			#print k,v

	num_non0 = len(non0_exp_list)
	sys.stderr.write("randomly choose %d out of %d genes as the final target genes that are to undergo AS\n" % (NTARG,num_non0))
	l=random.sample(xrange(num_non0), NTARG)
	out = open('AS_genes_list.txt','w')
	for i in l:
		out.write(non0_exp_list[i]+"\n")
	out.close()
if __name__ == '__main__': main()
