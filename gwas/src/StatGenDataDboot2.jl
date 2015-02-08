using DataFrames
using GLM
using Distributions
using NumericExtensions
using NumericFuns
using NLopt

module StatGenDataD

	# mwh - conf file! Change 6
	using IniFile
	ini = Inifile()
	read(ini, ENV["JULIA_APP_BOOTSTRAP"]) # conf file name
	include_dir=get(ini, "data", "include_dir") # julia source code

	using DataArrays
	using DataFrames
	import Base: getindex, size
    using DataFrames
	using GLM
	import DataFrames.DataFrameRegressionModel
	import StatsBase.sample
	import GLM.LinearModel
	import GLM.GeneralizedLinearModel
	#import GLM.LmMod
	#import GLM.GlmMod
	import GLM.Binomial
	import GLM.residuals
	#using Distributions  #this causes troubles with input to fit(LmMod,formlula,df) could not figure out
	import Distributions.ccdf
	import Distributions.Chisq
	import Distributions.FDist
	using NumericExtensions
	using NumericFuns   #this has Abs2Fun()
	import GLM.deviance
	import GLM.coef
	import DataFrames.allvars
	using NLopt


#from dGenDat2.jl
export dGenDat  #type and function for distributed GenDat, also uses GenDat type
export GenDat  #type and function
export numind   #get #of individuals in dataset
export numsnp		#get #of snps in dataset
export famname
export ival
export GenDatpart
export readbedmempart
export readbedmem
export bedfreq
export getMAF!
export updatecounts!
export dfctype!  #change type of dataarray in a df by colname of colnumber, or vector of colnames or numbers

#from dModifyGenDat2.jl
export addphe!
export pickSNPs!
export pickChroms!
export pickLocRange!
export MAFthreshhold!
export genCountthreshhold!
export pickIDs!
export HWthreshhold!
export HWtest
export HW!
export missingthreshhold!
export updateHW!

#from btGen5.jl
export ival
export btGen
export btGenFac
export btGenDat
export getmissbool
export gencount!
export getgentable
export getsnp  #for dGenDat and GenDat
export getsnpcounts
export putsnp  #retrieves given snp and places it in the fam dataframe on each processor

#from dgenglm6.jl
export gwLM
export dv
export lmtest
export effects
export anovatable
export lmssh
export logistictest
export gencut
export dGWAS   #type for distributed gwas
export getsigs  #brings all to main processor
export getresults   #brings all to main processor
export writeresults  #prints to same file without moving to main processor
##used in gwLM in dgenglm6.jl and in gwLRTmv dgwLRTmv2.jl
export allvars2   #fixes formula issue in DataFrames.allvar  where has :macrocall instead of :call

#types from  LRTtests.jl
export MeanVar
export mvinfo
export mvmod
export lrt   #type and function


#functions from  LRTtests.jl
export LRTmv
export LRTm!
export cutgen!
export logliknull
export lmssh
export loglikMV
export obj
export objsigma!

#functions from dgwLRTmv2.jl
export gwLRTmv
#have allvars2 exported above and is in dgenglm6.jl file
##used in gwLM in dgenglm6.jl and in gwLRTmv dgwLRTmv2.jl
#export allvars2   #fixes formula issue in DataFrames.allvar  where has :macrocall instead of :call

#types & functions from BootstrapLRTmv3.jl
export BootPiece
export dParBoot
export BootPiece
export ParBootLRT
export permboot!
export ParbootLRTmv
export bootlogliknull
export bootloglikMV
export ParbootLRTv

#types & functions from ParBootLM.jl
export BootPieceLM
export dParBootLM
export BootsLM
export ParBootLM

# include_dir has julia gwas source code
#fixed using symbols in DataFrames
include(joinpath(include_dir, "dGenDat3.jl"))         # has primary type for importing genetic data
include(joinpath(include_dir, "dModifyGenDat3.jl"))   # functions for importing phenotype data, filters, counts
include(joinpath(include_dir, "btGen6.jl"))           # functions for converting from binary to integers genotype data
include(joinpath(include_dir, "dgenglm9.jl"))         # functions for doing Linear model & logistic regression gwas
include(joinpath(include_dir, "LRTtests4.jl"))        # functions for the LRTmv, LRTv, & LRTm tests
include(joinpath(include_dir, "dgwLRTmv4.jl"))        # applies LRT tests to gwas data
include(joinpath(include_dir, "ParBootLRTmv4.jl"))    # functions to do Parametric bootstrap for LRTmv & LRTv
include(joinpath(include_dir, "ParBootLM2.jl"))       # functions to do Parametric bootstrap for linear models
include(joinpath(include_dir, "ParBootLogistic2.jl")) # functions to do Parametric bootstrap for logistic regression, have not tested well

end   # module StatGenDataD