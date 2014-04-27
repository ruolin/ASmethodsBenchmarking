###Prerequisites:
1. Install python 2.7+ but not 3.0+.
2. Install rpy2 (python library).
3. Install scipy (python library).
4. Install numpy (python library).
5. Install Flux Simulator. http://sammeth.net/confluence/display/SIM/2+-+Download
6. Copy myPara.par to the Flux Simulator root directory.

####Step1: generate gene-level NB counts. This scripts accepts two group of samples. Recommend for using real data.
usage: cal_NB_counts.py [-h] <annotation.gff3> -g1 control_1.sam control_2.sam -g2 treat_1.sam treat_2.sam [-n 3] [-l 2000] [-m AS-genes]

Generated negative binomial relationship gene-level read counts for a synthetic groups based on real experimental data

optional arguments:

	-h, --help            show this help message and exit

	-g1 GROUP1 [GROUP1 ...], --group1 GROUP1 [GROUP1 ...]			First group sam files separated by space.

	-g2 GROUP2 [GROUP2 ...], --group2 GROUP2 [GROUP2 ...]			Second group sam files separated by space.

	-n NREPS, --num-reps NREPS			Number of replicates. Default is 3.

	-l NTARG, --num-target-gene NTARG			Number of AS genes. Default is 2000.

	-m {AS-genes,all-genes}, --mode {AS-genes,all-genes}			Choose between AS-genes or all-genes: AS-genes simulates annotated AS genes only; all_genes simulates all genes in annotation.

####Output from setp 1
group1.nbcounts and group2.nbcounts: simulated NB fragment counts.

group1.rawcounts and group2.rawcounts: fragement counts for the input data.

AS_genes.list contains the simulated AS genes.

####Step2: Simulate differentail splicing. 
Usage: python generate_rnaseq.py [options] <group1.nbcounts> <AS_genes_list> <path_to_myPara.par> <out_prefix>

Options:

	-h, --help            show this help message and exit

	-p PALT, --percent-alt=PALT			The percentage of signal coming from alternate splice forms. Default is 0.2

	-c COVERAGE, --mean-base-coverage=COVERAGE			Mean base coverage. Default is 25

Copy right: Ruolin Liu, ISU

