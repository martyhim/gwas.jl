# Interactive session for running on Colonial One
#
# Read configuration file to get directory locations to add to LOAD_PATH
# Other startup parameters might also be necessary. The top level variable "ini",
# defined in .juliarc is accessible to the REPL

# Based on the values in the application conf file, LOAD_PATH has several directories
# added to it in .juliarc.jl. (.juliarc.jl is run at Julia startup.)
# The conf file to use is specified in the "JULIA_APP_BOOTSTRAP",
# which is set in .bashrc. .juliarc sets a variable, ini, which is available to the REPL.
# ini contains the configuration settings set in the JULIA_APP_BOOTSTRAP conf file.


# Julia routines for Slurm processing
require ("slurm_utility.jl")
require ("mpstat_parse.jl")
using SlurmUtility

#
#       Add processors
nodelist = slurm_nodelist_for_addprocs()
@time rp = convert(Array{Int,1}, addprocs(nodelist))

require("app_utility.jl"); @everywhere using AppUtility

@time sendto(rp, include_dir=include_dir)

@time @everywhere include(joinpath(include_dir, "StatGenDataDbootGGanova.jl"))
@time @everywhere using StatGenDataD

@time kdat=dGenDat("$(data_dir)smallAZdatasets/az12000snp");
@time kdat=dGenDat(joinpath(data_dir, "smallAZdatasets/az12000snp"))
phecorefile = joinpath(data_dir, "smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt")

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

-----
@everywhere using StatGenDataD
@time kdat=dGenDat("$(data_dir)smallAZdatasets/az12000snp")

phecorefile ="$(data_dir)smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt"
#this joins the phenotype data with the genotype data in the GenDat type on each process in the .fam field
@time addphe!(phecorefile,kdat);
@time updatecounts!(kdat);
missingthreshhold!(0.05,kdat);
MAFthreshhold!(0.01,kdat);

for i=1:length(kdat.refs)
	@spawnat kdat.refs[i].where fetch(kdat.refs[i]).fam[:Series]=PooledDataArray(fetch(kdat.refs[i]).fam[:Series])
	@spawnat kdat.refs[i].where (for k=1:1:size(fetch(kdat.refs[i]).fam,1) fetch(kdat.refs[i]).fam[:CDR12].data[k] -=1 end)
end

form_tau_ab42=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
form_ptau_ab42=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
form_cdr_ab42=CDR12~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp

@time ptau_ab42add12000=gwLM(form_ptau_ab42,1,kdat,responsetype=:linear);

writeresults("$(app_dir)ptau_ab42addrqtl12000.txt",ptau_ab42add12000)