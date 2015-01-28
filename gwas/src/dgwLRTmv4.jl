#@everywhere include("/Users/madmax/Desktop/desktopSep10/newdocsG5/oldG4/papers_and_stats/Talks_Papers/grants/ExonEpiGrant/EpiExonNovember/PengWeiYingCao/MeanVarMod.jl")
#using MeanVarMod
#df=DataFrame(randn(1000,5));  #make random dataset
#df["snp"]= float64([rep(0:2,333),1]); #make fake snp #rep comes from DataFrames
#@time cc=LRTmv(:(x1~x2+x3+snp),df);
#@time dd=LRTmv(:(x1~x2+x3+snp),df,:Gen);

#@everywhere include("/Users/madmax/StatGenData.jl/src/StatGenData3.jl")
#using StatGenData
#@time kk=GenDat("FTOdataAA");
##dd=ReadBed("/Users/madmax/FTOdataAA");
#phefile2="FTOpheAA.txt"
#addphe!(phefile2,kk);
#MAFthreshhold!(0.05,kk)


#@time hh=gwLRTmv(:(hdlsiu02~v1age01+sex_1+bmi01+snp),kk);
#@time ii=gwLRTmvp(:(hdlsiu02~v1age01+sex_1+bmi01+snp),kk);
##@time hh=gwLRTmv(:(hdl01~v1age01+sex_1+bmi01+snp),);


#function to work with distributed GenDat
function gwLRTmv(form::Formula,gen::dGenDat; AddorGen::Symbol=:Add,genlimit::Int=1,mafcut::FloatingPoint=0.0)
    resref=Array(RemoteRef,length(gen.refs))
    for i=1:length(gen.refs) resref[i] = @spawnat gen.refs[i].where gwLRTmv(form,fetch(gen.refs[i]),AddorGen=AddorGen,genlimit=genlimit,mafcut=mafcut) end
    #res=fetch(resref[1]);
    #for i=2:length(gen.refs) res=vcat(res,fetch(resref[i])) end
    #return res
   	if AddorGen==:Add
    	asfactor=false
    else
    	asfactor=true
    end
    return dGWAS(resref,:LRTmv,form,genlimit,asfactor,1)    
end

#function to work on single GenDat
function gwLRTmv(form::Formula,gen::GenDat; AddorGen::Symbol=:Add,genlimit::Int=1,mafcut::FloatingPoint=0.0)
	if AddorGen == :Add	#if :Gen then return snp as PooledDataArray otherwise sent out as continuous DataArray
		togen=btGenDat
	elseif AddorGen == :Gen
		togen=btGenFac
	else
		error("must be :Add or :Gen")
	end
	vars=DataFrames.allvars(form)
	#vars=allvars(form)
	vars=vars[vars.!=:snp]
	mdat=gen.fam[:,eval(vars)];  #pulls of variables from dataframe
	covmiss=complete_cases(mdat);
	mdat=mdat[covmiss,:];   #this is base dataframe to use
	mdat[:snp]=ones(size(mdat,1));  #place holder for snp
	famindskeep=gen.fam[:indx].data[covmiss]
	results=ones(Float64,numsnp(gen),6);
	bpt=convert(Ptr{Uint8},gen.gendat)
    m=size(gen.gendat,1)
    ct=[0,0,0,0];
    nsnp=numsnp(gen)
	for i=1:numsnp(gen)
		if gen.snp[i,:MAF]>=mafcut
			p=bpt + (gen.snp[:indx][i]-1)*m
		#function btGenFac{T <: Integer}(d::GenDat,snpind::Int,inds::Array{T,1},whichallele::T=1)
		#cc=btGen(jj,1,jj.fam["indx"].data,jj.snp["direction"][1])
			gencount!(ct,p,famindskeep);
			usetable, dotest =getgentable(ct,gen.snp[:direction][i],genlimit);  #gives gentable to use and dotest which is t/f to run or not
			if dotest
				#do mdat with new gen table cc[1]
                mdat[:snp]=togen(p,famindskeep,usetable)
                res= try 
						LRTmv(form,mdat,AddorGen);
				 	#catch y
				 		#if isa(y, DomainError)
				 		#[1.0,1.0,1.0,1.0,1.0,1.0]
				 		#end
				 	  catch
				 		zeros(Float64,6)
				 	end
				if typeof(res)== MeanVar
					results[i,:]=[res.LRTmv.chi,res.LRTmv.df,res.LRTmv.pval,res.LRTv.chi,res.LRTv.df,res.LRTv.pval]
				else
					results[i,:]=res
					#println("snp col $i failed")
				end
			else
				results[i,:]=zeros(Float64,6)
			end
			ct=[0,0,0,0];
		end
		if mod(i,100000)==0  println("$i done out of $nsnp") end
	end
	results=convert(DataFrame,results)
	names!(results,[:LRTmvChi,:LRTmvDF,:LRTmvPval,:LRTvChi,:LRTvDF,:LRTvPval])
	res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]],results)
	#res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]]
	#	,DataFrame(results,[:LRTmvChi,:LRTmvDF,:LRTmvPval,:LRTvChi,:LRTvDF,:LRTvPval]))
	#make sure minor allele is corrects
    for i = 1:numsnp(gen) 
        if gen.snp[:direction][i] res[:allele1][i] = gen.snp[:allele2][i] end 
    end
    #changes name from allele 1 to minor allele which genotype coding is # of minor alleles, 0,1,2
    names!(res,[:chr,:loc,:snpID,:minor,:MAF,:LRTmvChi,:LRTmvDF,:LRTmvPval,:LRTvChi,:LRTvDF,:LRTvPval])
    return res
end