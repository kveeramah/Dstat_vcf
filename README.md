# Dstat_vcf
This python script takes a vcf file and performs D_statistic analysis for a bunch of samples direct from a vcf

Only biallelic SNP loci are considered. When encountering a heteroygote, a random allele is chosen (i.e. pseudo haploid). 

numpy and pyvcf must be installed in the version of python used.

vcf files must be bgzipped and tabix indexed. Reference genome must also be indexed
The pop_test files must have the space seperated fields in the following order: W, X, Y, Z (see qpdstat manual in Admixtools) 

Default block size for delete m_i jacknife if 5MB and qual score threshold for including a site is 30. You can change these variables in the python script

Usage is ./Dstat.vcf.py <in.vcf.gz> <pop_test_file> <reference_genome> <out_stem> <nb_threads>

Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)
