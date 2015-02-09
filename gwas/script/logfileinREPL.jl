Lily:~ lily$ cd /Users/lily/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData2/
Lily:CodeandSampleData2 lily$ cd examplescriptandresults/
Lily:examplescriptandresults lily$ julia -p 2
               _
   _       _ _(_)_     |  A fresh approach to technical computing
  (_)     | (_) (_)    |  Documentation: http://docs.julialang.org
   _ _   _| |_  __ _   |  Type "help()" for help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 0.3.5 (2015-01-08 22:33 UTC)
 _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
|__/                   |  x86_64-apple-darwin13.4.0

julia> #start julia with more than one processor with the -p flag
       #julia -p 2
       #this is path to the CodeandSampleData folder which has the code in the current folder and data in the smallAZdatasets folder
       dir="/Users/lily/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData2/"
"/Users/lily/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData2/"

julia> #this tell julia to load all of the data in the StatGenDataDbootGGanova.jl file on all processes
       #@everywhere include("$(dir)current/StatGenDataDbootGGanova.jl")
       #for some reasion I cannot add put a variable in the the include function (in this case joining dir and the rest of the path)
       #you will need to open StatGenDataDbootGGanova.jl and change the cdir variable on line 151to your path with /CodeandSampleData/current/
       #StatGenDataDbootGGanova.jl  defines the module, what packages to use and the functions to export and at the bottom what files to load for it.

       @everywhere include("/Users/lily/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData2/current/StatGenDataDbootGGanova.jl")

Warning: could not import Base.add! into NumericExtensions
Warning: could not import Base.add! into NumericExtensions
Warning: could not import Base.add! into NumericExtensions

julia> #this tells it to use the model in StatGenDataDboot.jl called StatGenDataD onto all process (this takes a while)

julia> @everywhere using StatGenDataD

julia> #this make an instance of dGenDat which distributes pieces fo the data in GenDat types on each processor from the 
       #az1000snp.bed , az1000snp.bim & az1000snp.fam files
       #in this case I can use the dir variable to make a new path
       @time kdat=dGenDat("$(dir)smallAZdatasets/az12000snp");
There are 12000 SNPs and 2886 individualselapsed time: 0.095871579 seconds (5204184 bytes allocated)

julia> #this names the phenotype data file

julia> phecorefile ="$(dir)smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt"
"/Users/lily/Desktop/desktopSep10/Houston/GWU/GWUtranfer/GWUBudget/MartyHimmelstein/CodeandSampleData2/smallAZdatasets/CSFSep06_2013_v1.1coreNAapo.txt"

julia> #this joins the phenotype data with the genotype data in the GenDat type on each process in the .fam field
       @time addphe!(phecorefile,kdat);
elapsed time: 0.02842389 seconds (1099696 bytes allocated)

julia> #this updates the allele and genotype counts after the merge of snp and phenoytpe data
       #this can change the .snp field in GenDat

julia>  From worker 2:  There are 6000 SNPs and 2886 individuals2164 out of 2886 ID's are left.
  From worker 3:  There are 6000 SNPs and 2886 individuals2164 out of 2886 ID's are left.
julia> 

julia> @time updatecounts!(kdat);
elapsed time: 0.004997049 seconds (78328 bytes allocated)

julia> #this applies a missing threshold of 5% for each snp i.e., each snp with >5% missing is no included
       #this can change the .snp field in GenDat
       missingthreshhold!(0.05,kdat);

julia> #this applies a minor allele frequency threshhold to each snp, those less than 1% MAF are not included
       #this can change the .snp field in GenDat
       MAFthreshhold!(0.01,kdat);

julia> #6000 snps are sent to each of the two processes
       #note that after the missing and MAF thresholds the 12000 snps
       #it keeps 4818 of 6000 on the first and 3143 snp on the second leaving 7961 snps to be run

       #this is specific to this data
       #the Series variable needs to be treated as factors & the CDR12 variable need to have 1 subtracted from each value
       for i=1:length(kdat.refs)
        @spawnat kdat.refs[i].where fetch(kdat.refs[i]).fam[:Series]=PooledDataArray(fetch(kdat.refs[i]).fam[:Series])
        @spawnat kdat.refs[i].where (for k=1:1:size(fetch(kdat.refs[i]).fam,1) fetch(kdat.refs[i]).fam[:CDR12].data[k] -=1 end)
       end

julia>  From worker 3:  4818 out of 6000 SNPs were kept.
  From worker 2:  3143 out of 6000 SNPs were kept.
  From worker 2:  3143 out of 3143 SNPs were kept.
  From worker 3:  4818 out of 4818 SNPs were kept.
julia> 

julia> #defines a series of linear model formulas to be use all of the variables exist in the phenotype file except snp
       #which is proxy for each snp
       form_tau_ab42=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
Formula: lsubtau ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + lsubAb42 + snp + lsubAb42 & snp

julia> form_ptau_ab42=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
Formula: lsubptau ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + lsubAb42 + snp + lsubAb42 & snp

julia> form_cdr_ab42=CDR12~age+gender+Series+PC1+PC2+APOE2+APOE4+lsubAb42+snp+lsubAb42&snp
Formula: CDR12 ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + lsubAb42 + snp + lsubAb42 & snp

julia> #runs this a test of that formula (form_ptau_ab42) for each snp distributed across the processes
       #in this case it tests the last term in the model (lsubAb42&snp) which is the interaction term
       #it uses the data from kdat, it is a linear model, by default I am treating the snp as additve which means
       #I parameterize it as a single 0,1,2 index with 1 degree of freedom, otherwise I could specify  asfactor=true
       #and then the snp genotypes would be treated as factors and if 3 genotypes then it will have 2 df.
       #the results are stored in a dGWAS type spread across the processes

julia> @time ptau_ab42add12000=gwLM(form_ptau_ab42,1,kdat,responsetype=:linear);

time started
Thu Feb  5 10:28:40 MST 2015
time finished
Thu Feb  5 10:28:40 MST 2015

elapsed time: 14.461988529 seconds (2502200 bytes allocated)

julia> #this take the results from each process and writes them to a file.
       writeresults("ptau_ab42addrqtl12000.txt",ptau_ab42add12000)

julia> #if I want to I could put it all in one dataframe on the main process by:
       ptau_ab42add12000results=getresults(ptau_ab42add12000);

julia> #if on the main process, I could filter or do whatever I want with it.

julia> #this is a model using logistic regression (much slower than LM)
       #@time cdr_ab42add=gwLMp(form_cdr_ab42,1,kdat,responsetype="logistic");
       @time cdr_ab42add12000=gwLM(form_cdr_ab42,1,kdat,responsetype=:logistic);

time started
Thu Feb  5 10:30:10 MST 2015
time finished
Thu Feb  5 10:30:10 MST 2015

elapsed time: 29.747354862 seconds (48160 bytes allocated)

julia> writeresults("cdr_ab42addrqtl12000.txt",cdr_ab42add12000)

julia> aa=getresults(cdr_ab42add12000)
7961x8 DataFrame
| Row  | chr | loc     | snpID         | minor | MAF       | chi      | df  |
|------|-----|---------|---------------|-------|-----------|----------|-----|
| 1    | 0   | 0       | "rs1004236"   | "A"   | 0.0940823 | 0.800825 | 1.0 |
| 2    | 0   | 0       | "rs10155688"  | "A"   | 0.141701  | 6.54523  | 1.0 |
| 3    | 0   | 0       | "rs1045074"   | "A"   | 0.391636  | 0.182868 | 1.0 |
| 4    | 0   | 0       | "rs10502267"  | "A"   | 0.13338   | 0.89605  | 1.0 |
| 5    | 0   | 0       | "rs1062731"   | "A"   | 0.223532  | 1.8506   | 1.0 |
| 6    | 0   | 0       | "rs1079841"   | "A"   | 0.331331  | 0.418143 | 1.0 |
| 7    | 0   | 0       | "rs10836370"  | "C"   | 0.23287   | 0.497493 | 1.0 |
| 8    | 0   | 0       | "rs11090516"  | "T"   | 0.221939  | 2.14776  | 1.0 |
⋮
| 7953 | 1   | 6490225 | "rs113421613" | "A"   | 0.0704225 | 0.825721 | 1.0 |
| 7954 | 1   | 6490487 | "rs79830917"  | "C"   | 0.0706573 | 0.84267  | 1.0 |
| 7955 | 1   | 6490566 | "rs112138943" | "A"   | 0.0706573 | 0.84267  | 1.0 |
| 7956 | 1   | 6490592 | "rs113046942" | "A"   | 0.0706573 | 0.84267  | 1.0 |
| 7957 | 1   | 6491249 | "rs61780707"  | "T"   | 0.0424683 | 3.12551  | 1.0 |
| 7958 | 1   | 6491575 | "rs111347638" | "T"   | 0.0706573 | 0.84267  | 1.0 |
| 7959 | 1   | 6492333 | "rs4908893"   | "T"   | 0.0729995 | 1.58225  | 1.0 |
| 7960 | 1   | 6493261 | "rs12079984"  | "C"   | 0.0422337 | 3.19246  | 1.0 |
| 7961 | 1   | 6493792 | "rs113453591" | "T"   | 0.0704192 | 0.683972 | 1.0 |

| Row  | log10pval |
|------|-----------|
| 1    | 0.430805  |
| 2    | 1.97813   |
| 3    | 0.174626  |
| 4    | 0.46364   |
| 5    | 0.760167  |
| 6    | 0.285783  |
| 7    | 0.318213  |
| 8    | 0.845339  |
⋮
| 7953 | 0.43948   |
| 7954 | 0.445348  |
| 7955 | 0.445348  |
| 7956 | 0.445348  |
| 7957 | 1.11308   |
| 7958 | 0.445348  |
| 7959 | 0.681024  |
| 7960 | 1.13089   |
| 7961 | 0.389103  |

julia> if sum(aa[:log10pval].>7.3)>0 writetable("sig7.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>7.3,:],separator='\t') end

julia> if sum(aa[:log10pval].>6.3)>0 writetable("sig6.3cdr_ab42addrqtl12000.txt",aa[aa[:log10pval].>6.3,:],separator='\t') end

julia> 

julia> #example of running LRTmv test (which also does LRTv at same time)
       #these are the formulas
       form_lsubtau=lsubtau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
Formula: lsubtau ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + snp

julia> form_lsubptau=lsubptau~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
Formula: lsubptau ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + snp

julia> form_lsubab42=lsubAb42~age+gender+Series+PC1+PC2+APOE2+APOE4+snp
Formula: lsubAb42 ~ age + gender + Series + PC1 + PC2 + APOE2 + APOE4 + snp

julia> #uses formula, data in kdat, defines whether snp is additve or genotype, sets genotype limit per snp
       @time vlsubab42add12000=gwLRTmv(form_lsubab42,kdat,AddorGen=:Add,genlimit=20);

elapsed time: 2401.297506928 seconds (491436 bytes allocated)

julia> writeresults("lsubab42addvqtl12000.txt",vlsubab42add12000);

julia> aa=getresults(vlsubab42add12000);

julia> #this can be used to pull of snps with pvalue (in -log10pval) greater than some threshold
       if sum(aa[:LRTmvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>7.3,:],separator='\t') end

julia> if sum(aa[:LRTmvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlMV12000.txt",aa[aa[:LRTmvPval].>6.3,:],separator='\t') end

julia> if sum(aa[:LRTvPval].>7.3)>0 writetable("sig7.3lsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>7.3,:],separator='\t') end

julia> if sum(aa[:LRTvPval].>6.3)>0 writetable("sig6.3vlsubab42addvqtlV12000.txt",aa[aa[:LRTvPval].>6.3,:],separator='\t') end

       