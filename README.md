####Prerequisites:
1. Install python 2.7+ but not 3.0+.
2. Install rpy2 (python library).
3. Install scipy (python library).
4. Install numpy (python library).
5. Install Flux Simulator. http://sammeth.net/confluence/display/SIM/2+-+Download
6. Copy myPara.par to the Flux Simulator root directory.

####Step1: generate gene-level NB counts.
usage: cal_NB_counts.py [-h] -g1 GROUP1 [GROUP1 ...] -g2 GROUP2 [GROUP2 ...]
                        [-n NREPS] [-l NTARG] [-m {AS-genes,all-genes}]
                        gene_model

Generated negative binomial relationship gene-level read counts for a
synthetic groups based on real experimental data

positional arguments:
  gene_model            Reference annotation in GFF3 format.

optional arguments:
  -h, --help            show this help message and exit
  -g1 GROUP1 [GROUP1 ...], --group1 GROUP1 [GROUP1 ...]
                        First group sam files separated by space.
  -g2 GROUP2 [GROUP2 ...], --group2 GROUP2 [GROUP2 ...]
                        Second group sam files separated by space.
  -n NREPS, --num-reps NREPS
                        Number of replicates. Default is 3.
  -l NTARG, --num-target-gene NTARG
                        Number of AS genes. Default is 2000.
  -m {AS-genes,all-genes}, --mode {AS-genes,all-genes}
                        Choose between AS-genes or all-genes: AS-genes
                        simulates annotated AS genes only. all_genes simulates
                        all genes in annotation.


