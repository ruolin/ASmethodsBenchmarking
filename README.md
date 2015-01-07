###Prerequisites:
1. Install python 2.7+ but not 3.0+.
2. Install rpy2 (python library).
3. Install scipy (python library).
4. Install numpy (python library).
5. Install Flux Simulator. http://sammeth.net/confluence/display/SIM/2+-+Download
6. Copy myPara.par to the Flux Simulator root directory.

####Step1: Generate gene-level fragment counts, following Negetive Binomial distribution, for two synthetic groups based on real experimental data. One group represents control group; the other is for experiment group. Scripts accept two groups of  SAM files and output simulated gene-level fragment counts. The SAM inputs function as starting point in this pipeline and it simulate differential alternative splicing base on that. We recommend for using real data.

usage: cal_NB_counts.py [-h] <annotation.gff3> **-g1** control.1.sam control.2.sam **-g2** treat.1.sam treat.2.sam [-n 3] [-l 2000] [-m AS-genes]


Arguments:

	-h, --help            show this help message and exit

	-g1 GROUP1 [GROUP1 ...], --group1 GROUP1 [GROUP1 ...]	Require		First group sam files separated by space.

	-g2 GROUP2 [GROUP2 ...], --group2 GROUP2 [GROUP2 ...]	Require		Second group sam files separated by space.

	-n NREPS, --num-reps NREPS			Number of replicates. Default is 3.

	-l NTARG, --num-target-gene NTARG		Number of AS genes. Default is 2000.

	-m {AS-genes,all-genes}, --mode {AS-genes,all-genes}			Choose between AS-genes or all-genes: AS-genes simulates annotated AS genes only; all_genes simulates all genes in annotation.

####Output from setp 1
group1.nbcounts and group2.nbcounts: simulated NB fragment counts.

group1.rawcounts and group2.rawcounts: raw fragement counts for the input data.

AS-genes.list contains the simulated AS genes.

####Step2: Simulate differentail alternative splicing for each of the two groups individually. 

Usage: python generate_rnaseq.py groupx.nbcounts AS-genes.list path-to-myPara.par out-prefix [Options]

Required Arguments:

	groupx.nbcounts		Output from Step 1, containing the simulated fragments counts following Negetive Binomial ditribution. 
	AS-genes.list		Output from Step 1, containing the differentially AS genes.
path-to-myPara.par	Absoluate path for the file myPara.par.
out-prefix	Common prefix for all output files. E.g., control_1_1.fq, control_1.bed...
Options:
	-h, --help            show this help message and exit

	-p PALT, --percent-alt=PALT			The percentage of signal coming from alternate splice forms. Default is 0.2

	-c COVERAGE, --mean-base-coverage=COVERAGE			Mean base coverage. Default is 25

Copy right: Ruolin Liu, ISU

