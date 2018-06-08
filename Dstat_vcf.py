#!/usr/bin/env python
# -*- coding: ASCII -*-

#####This python script takes a vcf file and performs D_statistic analysis for a bunch of samples direct from a vcf
#####Only biallelic SNP loci are considered. When encountering a heteroygote, a random allele is chosen (i.e. pseudo haploid). 
#####numpy and pyvcf must be installed in the version of python used.
#####vcf files must be bgzipped and tabix indexed. Reference genome must also be indexed
#####The pop test files must have the tab seperated fields in the following order: sister1 sister2 introgress_test outgroup 
#####Default block size for delete m_i jacknife if 5MB and qual score threshold for including a site is 30. You can change these variables in the python script
#####Usage is ./Dstat.vcf.py <in.vcf.gz> <pop_test_file> <reference_genome> <out_stem> <nb_threads>
#####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)


usable_alleles=['A','C','G','T']

###import libraries
import string
import vcf
import pysam
from sys import argv
from random import randint
import numpy as np
import numpy.ma as ma
import multiprocessing as mp

###Input arguments
VCFin=argv[1]# tabindex and bgzipped vcf file
poptest_file=argv[2] #must have the tab seperated fields in the following order chromosome, start position, end position (zero_based)
ref_file=argv[3] #reference genome
filenameout=argv[4] #stem of output file
nbthreads=int(argv[5])

bl=5000000  ###block bootstrap size
qual=30  ###variant quality filter


###Read pop_test_file
file=open(poptest_file,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

pop_test=[]
samp_try=[]
for g in range(len(data)):
    k=string.split(data[g])
    pop_test.append(k)
    for gg in range(len(k)):
        if k[gg] not in samp_try:
            samp_try.append(k[gg])

###Load up vcf file
vcf_reader = vcf.Reader(open(VCFin, 'r'))


###Load up reference genome
ref=pysam.FastaFile(ref_file)

contigs=ref.references
contig_lengths=ref.lengths

usable_chroms=[]
usable_chrom_lengths=[]
for g in range(len(contigs)):
    try:
        if contig_lengths[g]>bl:
            vcf_reader.fetch(contigs[g], 0, contig_lengths[g])
            usable_chroms.append(contigs[g])
            usable_chrom_lengths.append(contig_lengths[g])
    except:
        print 'failed to find snps on '+contigs[g]


###Get sample names
samples=vcf_reader.samples

samp_use=[]
samp_dic={}
for g in range(len(samp_try)):
    if samp_try[g] not in samples:
        print samp_try[g]+' not found in vcf file. Any tests involving this sample will not work'
    else:
        samp_use.append(samp_try[g])
        samp_dic[samp_try[g]]=samples.index(samp_try[g])

nb_samps=len(samples)
nb_samp_use=len(samp_use)

###need this to calculate block sizes
fileout_snp=open(filenameout+'.snp','w')
    
print 'counting usable biallelic loci'
count_pos=0

for g in range(len(usable_chroms)):
    for record in vcf_reader.fetch(usable_chroms[g], 0, usable_chrom_lengths[g]):
        pos=record.POS
        alleles=[record.REF]
        for gg in range(len(record.ALT)):
            alleles.append(str(record.ALT[gg]))
        if (len(alleles)==2) and (record.QUAL>qual):
            if (alleles[0] in usable_alleles) and (alleles[1] in usable_alleles):
                out=usable_chroms[g]+'\t'+str(pos)+'\t'+alleles[0]+'\t'+alleles[1]+'\n'
                fileout_snp.write(out)
                count_pos+=1
               
        if count_pos%1000==0:
            print usable_chroms[g],pos,count_pos

fileout_snp.close()
nb_snps=count_pos

###array for jacknife blocks
blck_array=np.zeros(nb_snps,dtype='int32')


####block estimation

chrom_ref=0
chrom_start=0

file=open(filenameout+'.snp','r')
x=file.readline()

count_pos=0
while x<>'':
    chrom,pos,all1,all2=string.split(x[:-1],'\t')
    if chrom<>chrom_ref:
        chrom_ref=chrom
        chrom_start=np.max(blck_array)+1
    blck_array[count_pos]=(int(pos)/bl)+chrom_start
    count_pos+=1
    x=file.readline()
    
nb_blocks=np.max(blck_array)


###arrays for individual genotypes
ind_geno={}
for g in range(len(samp_use)):
    ind_geno[samp_use[g]]=np.zeros(nb_snps,dtype='int32')
    ind_geno[samp_use[g]][:]=-9


count_pos=0
print 'Creating frequency matrices for downstream calculation'
for g in range(len(usable_chroms)):
    for record in vcf_reader.fetch(usable_chroms[g], 0, usable_chrom_lengths[g]):
        pos=record.POS
        alleles=[record.REF]
        for gg in range(len(record.ALT)):
            alleles.append(str(record.ALT[gg]))
        if (len(alleles)==2) and (record.QUAL>qual):
            if (alleles[0] in usable_alleles) and (alleles[1] in usable_alleles):
                for gg in range(len(samp_use)):
                    gt=record.samples[samp_dic[samp_use[gg]]]['GT']
                    if gt=='0/0':
                        ind_geno[samp_use[gg]][count_pos]=0
                    elif gt=='1/1':
                        ind_geno[samp_use[gg]][count_pos]=1
                    elif gt=='0/1':
                        ind_geno[samp_use[gg]][count_pos]=randint(0,1)
                count_pos+=1
               
        if count_pos%1000==0:
            print usable_chroms[g],pos,count_pos


batches=[]
for g in range(0,len(pop_test),nbthreads):
    batches.append(pop_test[g:g+nbthreads])



fileout=open(filenameout+'_f4','w')
out='samp1\tsamp2\tsamp3\tsamp4\tD_stat\tZ\tnbSNP\n'
fileout.write(out)
fileout.close()

def f4_mult(x,ind_geno,batch,nb_blocks,blck_array,output):
    
    samp1,samp2,samp3,samp4=batch[x]
    
    try:
        

        nb_snp=len(ind_geno[samp1])

        mask=np.zeros(nb_snp,dtype='int32')
        mask[np.where(ind_geno[samp1]==-9)]=1
        mask[np.where(ind_geno[samp2]==-9)]=1
        mask[np.where(ind_geno[samp3]==-9)]=1
        mask[np.where(ind_geno[samp4]==-9)]=1

        use_sites=np.where(mask==0)[0]

        m=np.zeros((nb_blocks),dtype='float32')
        for gg in range(1,nb_blocks+1):
            m[gg-1]=len(np.where(blck_array[use_sites]==gg)[0])  ##snp number per window

        n=len(use_sites)
        
        h=n/m
        no_sites=np.where(m==0)[0]

        G=len(m)-len(no_sites)

        blck_mask=np.zeros(nb_blocks,dtype='int32')
        blck_mask[np.where(m==0)[0]]=1

        ##calculate dstat
        num=(ind_geno[samp1]-ind_geno[samp2])*(ind_geno[samp3]-ind_geno[samp4])
        den=(ind_geno[samp1]+ind_geno[samp2]-2*ind_geno[samp1]*ind_geno[samp2])*(ind_geno[samp3]+ind_geno[samp4]-2*ind_geno[samp3]*ind_geno[samp4])
        
        a=np.sum(ma.masked_array(num,mask=mask),dtype='float64')   
        b=np.sum(ma.masked_array(den,mask=mask),dtype='float64') 

        Dstat=a/b
            
        D_jk=np.zeros((nb_blocks),dtype='float32')
        
        for ggg in range(1,nb_blocks+1):
            #print x,pw,ref_index,ggg
            mask_jk=np.copy(mask)
            mask_jk[np.where(blck_array==ggg)[0]]=1
            
            a1=np.sum(ma.masked_array(num,mask=mask_jk),dtype='float64')   
            b1=np.sum(ma.masked_array(den,mask=mask_jk),dtype='float64') 

            D_jk[ggg-1]=a1/b1

        ###delete-m estimator
        Dstat_Jmj=(Dstat*G) - np.sum(ma.masked_array((1-m/n)*D_jk,mask=blck_mask))

        ###delete-m variance
        Dstat_var=np.sum(ma.masked_array((1/(h-1)) * ((h*Dstat - (h-1)*D_jk - Dstat_Jmj)**2),mask=blck_mask))/G


        Z=Dstat/np.sqrt(Dstat_var)
     
        out=samp1+'\t'+samp2+'\t'+samp3+'\t'+samp4+'\t'+str(Dstat)+'\t'+str(Z)+'\t'+str(n)+'\n'

    except:
        out=samp1+'\t'+samp2+'\t'+samp3+'\t'+samp4+'\tNA\tNA\tNA\n'
        
    output.put([x,out])


for g in range(len(batches)):
    nbthreads2=len(batches[g])

    ###queue for parallelism output
    output = mp.Queue()

    # Setup a list of processes
    processes = [mp.Process(target=f4_mult, args=(x,ind_geno,batches[g],nb_blocks,blck_array,output)) for x in range(nbthreads2)]

    # Run processes
    for p in processes:
        p.start()
     
    # Exit the completed processes
    for p in processes:
        p.join()

    results = [output.get() for p in processes]

    results.sort()

    out=''
    for gg in range(len(results)):
        out=out+results[gg][1]
    
    fileout=open(filenameout+'_f4','a')
    fileout.write(out)
    fileout.close()

