#
# Run Julia on the Colonial One Cluster, using Slurm as a Job Manager
# Example Slurm commands to allocate nodes on the cluster
# sinfo # returns nodes, including those that are free
# salloc --time=4:00:00 -N 4 -p 128gb # Allocate 4 nodes on the 128gb cluster for four hours
# We got a jobid, start an interactive shell on the first node (this is purely a convention
# srun --jobid 734682 -N 1 -w node041 --pty /bin/bash                          

# Julia routines for Slurm processing
require ("slurm_utility.jl")
using SlurmUtility

#       Add cluster nodes
nodelist = slurm_nodelist_for_addprocs() # Tell julia to spin up one process for each core per machine
@time rp = convert(Array{Int,1}, addprocs(nodelist))

# Note: I don't think AppUtility needs to be @everywhere
@time require("app_utility.jl"); @everywhere using AppUtility

# define include_dir variable on all processes
@time sendto(rp, include_dir=include_dir)
@time @everywhere include(joinpath(include_dir, "StatGenDataDbootGGanova.jl"))
@time @everywhere using StatGenDataD

#this make an instance of dGenDat which distributes pieces fo the data in GenDat types on each processor from the 
#az1000snp.bed , az1000snp.bim & az1000snp.fam files
@time kdat=dGenDat(joinpath(data_dir, "smallAZdatasets/az12000snp"))
phecorefile = joinpath(data_dir, "smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt")

#this joins the phenotype data with the genotype data in the GenDat type on each process in the .fam field
@time addphe!(phecorefile,kdat)

#this updates the allele and genotype counts after the merge of snp and phenoytpe data
#this can change the .snp field in GenDat
@time updatecounts!(kdat)

#this applies a missing threshold of 5% for each snp i.e., each snp with >5% missing is no included
#this can change the .snp field in GenDat
missingthreshhold!(0.05,kdat)

#this applies a minor allele frequency threshhold to each snp, those less than 1% MAF are not included
#this can change the .snp field in GenDat
MAFthreshhold!(0.01,kdat)

#6000 snps are sent to each of the two processes
#note that after the missing and MAF thresholds the 12000 snps
#it keeps 4818 of 6000 on the first and 3143 snp on the second leaving 7961 snps to be run

#this is specific to this data
#the Series variable needs to be treated as factors & the CDR12 variable need to have 1 subtracted from each value
for i=1:length(kdat.refs)
	@spawnat kdat.refs[i].where fetch(kdat.refs[i]).fam[:Series]=PooledDataArray(fetch(kdat.refs[i]).fam[:Series])
	@spawnat kdat.refs[i].where (for k=1:1:size(fetch(kdat.refs[i]).fam,1) fetch(kdat.refs[i]).fam[:CDR12].data[k] -=1 end)
end

#defines a series of linear model formulas to be use all of the variables exist in the phenotype file except snp
#which is proxy for each snp
form_tau_ab42=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp

form_ptau_ab42=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp

form_cdr_ab42=CDR12~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp

#runs this a test of that formula (form_ptau_ab42) for each snp distributed across the processes
#in this case it tests the last term in the model (lsubAb42&snp) which is the interaction term
#it uses the data from kdat, it is a linear model, by default I am treating the snp as additve which means
#I parameterize it as a single 0,1,2 index with 1 degree of freedom, otherwise I could specify  asfactor=true
#and then the snp genotypes would be treated as factors and if 3 genotypes then it will have 2 df.
#the results are stored in a dGWAS type spread across the processes
@time ptau_ab42add12000=gwLM(form_ptau_ab42,1,kdat,responsetype=:linear)

#this take the results from each process and writes them to a file.
writeresults("ptau_ab42addrqtl12000.txt",ptau_ab42add12000)
#if I want to I could put it all in one dataframe on the main process by:
ptau_ab42add12000results=getresults(ptau_ab42add12000)
#if on the main process, I could filter or do whatever I want with it.

#this is a model using logistic regression (much slower than LM)
#@time cdr_ab42add=gwLMp(form_cdr_ab42,1,kdat,responsetype="logistic");
@time cdr_ab42add12000=gwLM(form_cdr_ab42,1,kdat,responsetype=:logistic)

@time writeresults("cdr_ab42addrqtl12000.txt",cdr_ab42add12000)
@time aa=getresults(cdr_ab42add12000)
if sum(aa[:log10pval].>7.3)>0 writetable("sig7.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>7.3,:],separator='\t') end
if sum(aa[:log10pval].>6.3)>0 writetable("sig6.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>6.3,:],separator='\t') end




#example of running LRTmv test (which also does LRTv at same time)
#these are the formulas
form_lsubtau=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
form_lsubptau=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
form_lsubab42=lsubAb42~age+gender+Series+PC1+PC2+APOE2+APOE4+snp



#uses formula, data in kdat, defines whether snp is additve or genotype, sets genotype limit per snp
@time vlsubab42add12000=gwLRTmv(form_lsubab42,kdat,AddorGen=:Add,genlimit=20)
@time writeresults("lsubab42addvqtl12000.txt",vlsubab42add12000)
aa=getresults(vlsubab42add12000)
#this can be used to pull of snps with pvalue (in -log10pval) greater than some threshold
if sum(aa[:LRTmvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>7.3,:],separator='\t') end
if sum(aa[:LRTmvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>6.3,:],separator='\t') end
if sum(aa[:LRTvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>7.3,:],separator='\t') end
if sum(aa[:LRTvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>6.3,:],separator='\t') end

#form_ptau_ab42=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
#need to pull down into dataframe and 
ParBootLM(form_ptau_ab42,EAd[idx],[10],200000000)[1]




