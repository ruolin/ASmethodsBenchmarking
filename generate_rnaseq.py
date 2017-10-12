#from genNBcounts import meanVar
import sys, optparse, random,os
from itertools import izip_longest
import numpy as np
from sets import Set
import re
import subprocess
from math import sqrt, log
from itertools import repeat
Parser = optparse.OptionParser(
	usage = "python %prog [options] <group1.nbcounts> <AS_genes_list> <flux.par> <out_prefix>",
	description = "group1.nbcounts and AS_genes_list are produced by cal_NB_counts.py. flux.par is the parameter file used by flux simulator.",
	epilog =  "Copy right: Ruolin Liu, ISU"    
	)

Parser.add_option("-p", "--percent-alt", type="float", dest="palt", default = 0.2,
	help = " The percentage of signal coming from alternate splice forms. Default is 0.2" )

Parser.add_option("-c", "--mean-base-coverage", type="int", dest="coverage", default = 25,
	help = " Mean base coverage. Default is 25")

(options,args) = Parser.parse_args()
if len(args) != 4:
	Parser.print_help()
	sys.exit(1)

'''Define Global Variable'''
(FILE_SIMU,FILE_1KAS,FILE_PAR,OUT_PREFIX) = (args[0], args[1], args[2], args[3])
FLUX_DIR = os.path.dirname(os.path.realpath(FILE_PAR))
PALT = options.palt ## percentage of alternative isoform
COVERAGE = options.coverage ## mean base coverage
READLEN = 0 ## read length
PRO_NAME = None ## flux .pro file name
TRANSCRIPTOM = dict() ### average transcripts length for each gene
list_1kAS = list() ## list of known AS genes
dup_1kAS =list() ## need to be deleted
gene2count = dict() ## gene-to-count dictionary
LIST_PAR = list() ## list object stores flux .par file

gene = re.compile("gene_id \"(.*?)\";")
iso = re.compile("transcript_id \"(.*?)\";")
class JSD(object):

	def KLdivergence(self, p, q):
		""" Compute KL divergence of two vectors, K(p || q)."""
		return sum(_p * np.log2(_p / _q) for _p, _q in zip(p, q) if _p != 0)

	def sqrtJSD(self, p, q):
		""" Returns the Jensen-Shannon divergence. """
		self.JSD = 0.0
		weight = 0.5
		average = np.zeros(len(p)) #Average
		for x in range(len(p)):
			average[x] = weight * p[x] + (1 - weight) * q[x]
		self.JSD = (weight * self.KLdivergence(np.array(p), average)) + ((1 - weight) * self.KLdivergence(np.array(q), average))
		return sqrt(self.JSD)

##par flux paramter setting file:
def callFlux2PRO(_file_par):
	global PRO_NAME, READLEN
	with open(_file_par,"r") as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			else:
				LIST_PAR.append(line)
				col = line.rstrip().split()
				if col[0] == "PRO_FILE_NAME":
					PRO_NAME = FLUX_DIR+'/'+col[1]
				elif col[0] == "READ_LENGTH":
					READLEN = int(col[1])
	f.close()
	_file_par_baseName = os.path.basename(_file_par)
	args = ['bin/flux-simulator', '-xp', _file_par_baseName] 
	process=subprocess.Popen(args,shell=False, cwd=FLUX_DIR)
	process.wait()
	if process.returncode != 0:
		sys.stderr.write("Flux simulator didn't run correctly \n")
		sys.exit(2)

def split2Frac(n):
	'''Set transcripts percentage of expression randomly'''
	random.seed(80)
	dividers = sorted(random.sample(xrange(1, 100), n - 1))
	return [(a - b)/100.0 for a, b in zip(dividers + [100], [0] + dividers)]

def setAltRatio(n):
	'''Set transcripts percentage of expression according to the alternative ratio parameter PALT'''
	if n==1:
		ret = [1]
	else:
		rest = repeat((1- PALT)/float(n-1), n-1)
		ret = list(rest)+[PALT]
	return ret	


def read1kAS(_file):
	with open (_file,'r') as f1:
		for line in f1:
			col = line.rstrip().split()
			list_1kAS.append(col[0])
			dup_1kAS.append(col[0])
	f1.close()


def calcuNB_mol(_file):
	'''_file: empirical counts from real data
	   _perThred: min counts for the known 1000 AS genes 
	'''
	with open (_file,'r') as f2:
		# read first line to gene2count
		firstline=f2.readline()
		_samples = np.array(firstline.rstrip().split()[1:],dtype=np.uint32)
		gene2count[firstline.rstrip().split()[0]] = _samples
		# read rest lines to gene2count
		for line in f2:
			col = line.rstrip().split()
			values = np.array(col[1:],dtype=np.uint32)
			gene2count[col[0]] = values
			_samples = np.vstack([_samples,values])
	average = _samples.mean(axis=0)
	sum = _samples.sum(axis=0)
	return sum
def genPRO(_count, _file_PRO, _rep):
	''' Generate .pro file based on the following parameters
	1. _count:gene2count
	2. _file_PRO: flux .PRO file
	3. _rep: number of replicates simulated
	Return total transcriptom length
	'''
	#Flag = False
	global TRANSCRIPTOM 
	#if not TRANSCRIPTOM:
		#Flag =True # make sure TRANSCRIPTOM initialize only once
	total_length = 0
	_PRO_list = list()
	num_trans = 0
	total_txs_len = 0
	
        tx_file = None
	for line in LIST_PAR:
            fields = line.strip().split()
            if fields[0] == "REF_FILE_NAME": 
                tx_file = fields[1]

        print "processing transcript gtf file: ", tx_file
        tx_to_gene = {}
        prev_loci = None
        with open(tx_file, 'r') as fh:
            for line in fh:
                if line[0] == '#':
                    continue
                fields = line.strip().split("\t")
                gret = gene.search(fields[8])
                tret = iso.search(fields[8])
                if not gret or not tret:
                    sys.stderr.write("Warning: gff line '{}' does not have gene_id or transcript_id field".format(line.strip()))
                    sys.stderr.write("\n")
                else:
                    g = gret.group(1)
                    t = tret.group(1)
                    if not prev_loci:
                        prev_loci = g
                    if t not in tx_to_gene:
                        tx_to_gene[t] = g 
        #for key,value in tx_to_gene.iteritems():
            #print key, value

	isoforms = list()
	## Update every lines from toy.pro to _PRO_list
	with open(_file_PRO,'r') as f3:
		for line in f3:
			col = line.rstrip().split()
			cur_loci = tx_to_gene[col[1]]
			if prev_loci == cur_loci:
				num_trans += 1
				#if Flag:
				total_txs_len += int(col[3])
			else:
				for i in range(len(isoforms)):
					if prev_loci in gene2count:
						if prev_loci in dup_1kAS:
							isoforms[i][5] = int(setAltRatio(num_trans)[i] * gene2count[prev_loci][_rep])
						else:
							isoforms[i][5] = int(split2Frac(num_trans)[i] * gene2count[prev_loci][_rep])
					else:
						isoforms[i][5] = 0
				for each in isoforms:
					_PRO_list.append("\t".join(str(x) for x in each)+"\n")
				del isoforms[:]
				#if Flag:
				TRANSCRIPTOM[prev_loci] = int(total_txs_len/num_trans)
				total_txs_len = int(col[3])
				num_trans = 1
				prev_loci = cur_loci
			isoforms.append(col)
		# last line
                for i in range(len(isoforms)):
                        if prev_loci in gene2count:
                                if prev_loci in dup_1kAS:
                                        isoforms[i][5] = int(setAltRatio(num_trans)[i] * gene2count[prev_loci][_rep])
                                else:
                                        isoforms[i][5] = int(split2Frac(num_trans)[i] * gene2count[prev_loci][_rep])
                        else:
                                isoforms[i][5] = 0
                for each in isoforms:
                        _PRO_list.append("\t".join(str(x) for x in each)+"\n")
                TRANSCRIPTOM[prev_loci] = int(total_txs_len/num_trans)
		#_PRO_list.append("\t".join(str(x) for x in each)+"\n")
	f3.close()
	
	## print every lines in _PRO_list to file
	pro_string = OUT_PREFIX+'_'+str(_rep+1)+'.pro'
	pro_file = open(pro_string,'w')
	for each in _PRO_list: 
		pro_file.write(each)
	pro_file.close()
	
	## return total transcriptom length
	for length in TRANSCRIPTOM.itervalues():
		total_length += length
	return total_length 	

def genPAR(_rnum, _molnum, _rep):
	'''
	Generate flux .par file for the next round of running flux simulator
	'''
	global READLEN
	par_string = OUT_PREFIX+'_'+str(_rep+1)+'.par'
	par_file = open(par_string,'w')
	for line in LIST_PAR:
		col = line.rstrip().split()
		if col[0] == "PRO_FILE_NAME":
			col[1] = OUT_PREFIX+'_'+str(_rep+1)+'.pro'
		elif col[0] == "LIB_FILE_NAME":
			col[1] = OUT_PREFIX+'_'+str(_rep+1)
		elif col[0] == "NB_MOLECULES":
			col[1] = _molnum
		elif col[0] == "READ_NUMBER":
			col[1] = str(_rnum)
		par_file.write("\t".join(col)+"\n")
	par_file.close()
	return par_string

def seperatePE(_fastq_pre):
	'''Splits the flux pair-end reads to read_1 and read_2'''
	
	forward = open(_fastq_pre+"_1.fq", 'w')
	reverse = open(_fastq_pre+"_2.fq", 'w')
	_fastq_file = _fastq_pre+".fastq"
	with open(_fastq_file, 'r') as f:
		for next_4_lines in izip_longest(*[f] * 4):
			list_of_lines = [each for each in next_4_lines]
			if len(list_of_lines[3]) != len(list_of_lines[1]): continue
			if list_of_lines[0].startswith('@'):
				if list_of_lines[0].rstrip().endswith('/1'):
					list_of_lines[0] = list_of_lines[0].rstrip()[:-4]+"/1"
					for each in list_of_lines:
						forward.write(each.rstrip()+"\n")
				elif list_of_lines[0].rstrip().endswith('/2'):
					list_of_lines[0] = list_of_lines[0].rstrip()[:-4]+"/2"
					for each in list_of_lines:
						reverse.write(each.rstrip()+"\n")
			else:
				sys.stderr.write("not a valid fastq file\n")
				sys.exit(0)
	forward.close()
	reverse.close()


def main():
	callFlux2PRO(FILE_PAR)
	J = JSD()
	control = [0.8, 0.2]
	print "sqare root JSD:", J.sqrtJSD(control,setAltRatio(2))	
	read1kAS(FILE_1KAS)
	mol_number=calcuNB_mol(FILE_SIMU)
	##debug 
	#for key,value in sorted(gene2count.iteritems()):
		#if key in dup_1kAS:
			#print key,value
	for i in xrange(len(mol_number)):
		###each iteration generates a fastq file.
		
		### GENERATE .PRO FILE
		transcriptom_len = genPRO(gene2count,PRO_NAME,i)
		### GENERATE .PAR FILE
		read_number = transcriptom_len * COVERAGE / READLEN
		print mol_number[i]
		file_par = genPAR(read_number, str(mol_number[i]), i)
                print "new par: ", file_par
		### RUN FLUX to produce fastq file
		args = [FLUX_DIR+'/bin/flux-simulator', '-lsp', file_par] 
                print " ".join(args)
		process_flux=subprocess.Popen(args,shell=False)
		process_flux.wait()
		if process_flux.returncode == 0:
			sys.stderr.write("Flux simulator run successfully on iteration %d \n" % (i+1))
			fastq_prefix = OUT_PREFIX+'_'+str(i+1)
			seperatePE(fastq_prefix)
			sys.stderr.write("split output to paired end\n")
	process_delToy = subprocess.call(["rm",PRO_NAME])
	if process_delToy == 0:
		sys.stderr.write("Delete toy.pro\n")
if __name__ == "__main__":
	main()
