#start julia with more than one processor with the -p flag
#julia -p 2
#this is path to the CodeandSampleData folder which has the code in the current folder and data in the smallAZdatasets folder
dir="/Users/madmax/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData/"
dir="/home/mhimmelstein/CodeandSampleData/" # On Colonial One

#this tell julia to load all of the data in the StatGenDataDboot2.jl file on all processes
#@everywhere include("$(dir)current/StatGenDataDboot2.jl")
#for some reasion I cannot add put a variable in the the include function (in this case joining dir and the rest of the path)
#you will need to open StatGenDataDboot2.jl and change the cdir variable on line 132 to your path with /CodeandSampleData/current/
#StatGenDataDboot2.jl  defines the module, what packages to use and the functions to export and at the bottom what files to load for it.

# On Colonial One Cluster, for now
@everywhere include("/home/mhimmelstein/CodeandSampleData/current/StatGenDataDboot2.jl")
@everywhere include("/Users/madmax/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData/current/StatGenDataDboot2.jl")
#this tells it to use the model in StatGenDataDboot.jl called StatGenDataD onto all process (this takes a while)
@everywhere using StatGenDataD

#this make an instance of dGenDat which distributes pieces fo the data in GenDat types on each processor from the 
#az1000snp.bed , az1000snp.bim & az1000snp.fam files
#in this case I can use the dir variable to make a new path
@time kdat=dGenDat("$(dir)smallAZdatasets/az12000snp");
#this names the phenotype data file
phecorefile ="$(dir)smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt"

#this joins the phenotype data with the genotype data in the GenDat type on each process in the .fam field
@time addphe!(phecorefile,kdat);
#this updates the allele and genotype counts after the merge of snp and phenoytpe data
#this can change the .snp field in GenDat
@time updatecounts!(kdat);
#this applies a missing threshold of 5% for each snp i.e., each snp with >5% missing is no included
#this can change the .snp field in GenDat
missingthreshhold!(0.05,kdat);
#this applies a minor allele frequency threshhold to each snp, those less than 1% MAF are not included
#this can change the .snp field in GenDat
MAFthreshhold!(0.01,kdat);

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
@time ptau_ab42add12000=gwLM(form_ptau_ab42,1,kdat,responsetype=:linear);
#this take the results from each process and writes them to a file.
writeresults("ptau_ab42addrqtl12000.txt",ptau_ab42add12000)
#if I want to I could put it all in one dataframe on the main process by:
ptau_ab42add12000results=getresults(ptau_ab42add12000);
#if on the main process, I could filter or do whatever I want with it.

#this is a model using logistic regression (much slower than LM)
#@time cdr_ab42add=gwLMp(form_cdr_ab42,1,kdat,responsetype="logistic");
@time cdr_ab42add12000=gwLM(form_cdr_ab42,1,kdat,responsetype=:logistic);
writeresults("cdr_ab42addrqtl12000.txt",cdr_ab42add12000)
aa=getresults(cdr_ab42add12000)
if sum(aa[:log10pval].>7.3)>0 writetable("sig7.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>7.3,:],separator='\t') end
if sum(aa[:log10pval].>6.3)>0 writetable("sig6.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>6.3,:],separator='\t') end




#example of running LRTmv test (which also does LRTv at same time)
#these are the formulas
form_lsubtau=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
form_lsubtau=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
form_lsubab42=lsubAb42~age+gender+Series+PC1+PC2+APOE2+APOE4+snp



#uses formula, data in kdat, defines whether snp is additve or genotype, sets genotype limit per snp
@time vlsubab42add12000=gwLRTmv(form_lsubab42,kdat,AddorGen=:Add,genlimit=20);
writeresults("lsubab42addvqtl12000.txt",vlsubab42add12000);
aa=getresults(vlsubab42add12000);
#this can be used to pull of snps with pvalue (in -log10pval) greater than some threshold
if sum(aa[:LRTmvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>7.3,:],separator='\t') end
if sum(aa[:LRTmvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>6.3,:],separator='\t') end
if sum(aa[:LRTvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>7.3,:],separator='\t') end
if sum(aa[:LRTvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>6.3,:],separator='\t') end




