 # get effects for model terms from vector for linear models Q'y = effects matrix
#effects(mod::LinearModel)=(mod.pp.qr[:Q]'*mod.rr.y)[1:size(mod.pp.X.m,2)]
effects(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64})=(mod.model.pp.qr[:Q]'*mod.model.rr.y)[1:size(mod.model.pp.X,2)]
effects(mod::LinearModel{DensePredQR{Float64}})=(mod.pp.qr[:Q]'*mod.rr.y)[1:size(mod.pp.X,2)]



#function anovatable(mod::LinearModel)
#function anovatable(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64})
#function anovatable(mod::DataFrameRegressionModel{LinearModel})
#function anovatable(mod)
function anovatable(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64})
	eff=effects(mod)
	#termnames=mod.fr.terms.terms
    termnames=mod.mf.terms.terms
	termnames=[:Intercept,termnames,:Residuals]
	#terms=[unique(mod.pp.X.assign),-1]
    terms=[unique(mod.mm.assign),-1]
	DF=zeros(Int,length(terms))
	SS=zeros(Float64,length(terms))
	MS=zeros(Float64,length(terms))
	fstat=zeros(Float64,length(terms))
	pval=zeros(Float64,length(terms))
	#DF[end]=df_residual(mod.pp)
    DF[end]=df_residual(mod.model.pp)
	#SS[end]=deviance(mod.rr)
    SS[end]=deviance(mod.model.rr)
	MS[end]=SS[end]/DF[end]

	for i=1:(length(terms)-1)
		#jj=eff[mod.pp.X.assign.==terms[i]]
        jj=eff[mod.mm.assign.==terms[i]]
		DF[i]=length(jj)
		#SS[i]=sum(Abs2(),jj)
        SS[i]=sumsq(jj)
		MS[i]=SS[i]/DF[i]
		fstat[i]=MS[i]/MS[end]
		pval[i]=ccdf(FDist(DF[i], DF[end]), fstat[i])
	end
	table=DataFrame(
		Terms = termnames[2:end],
		DF = DF[2:end],
		SS = SS[2:end],
		MS = MS[2:end],
		F = fstat[2:end],
		pval = pval[2:end]
		)
end

function logistictest(form::Formula,nterm::Int,df::AbstractDataFrame)
    #does full vs reduced model ftest
    #nterm is the number terms (not necessarily DFs) to test as the end of the model
    #for instance y~a+b+c+d+c&d  is want to test c&d int then nterm=1, if want test of
    #c+d+c&d then nterm=3   it could be more than 3 DF if any term has factors
    #requires using DataFrame and using GLM
    a=fit(GeneralizedLinearModel,form,df,Binomial())
    reduce=maximum(a.mm.assign)-(nterm-1)  #gets location of cut of for reduced model
    devFull=deviance(a)
    dfFull=df_residual(a.model.pp)
    #cut terms out of matrix model
    a.model.pp.X=a.model.pp.X[:,a.mm.assign.<reduce]
    ared=fit(GeneralizedLinearModel,a.model.pp.X,model_response(a.mf),Binomial())
    dfH=df_residual(ared.pp)-dfFull   #df of hypotheses = dfReduced-dfHypothesis
    chistat=deviance(ared)-devFull  #chi stat = deviance of reduced - deviance of full
    return [chistat,dfH, -log10(ccdf(Chisq(dfH), chistat))] #Chi square,dfH,-log10(pval)
end


#function lmssh(mod::LinearModel,nterm::Int)
function lmssh(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64},nterm::Int)
	#reduce=maximum(mod.pp.X.assign)-(nterm)
    reduce=maximum(mod.mm.assign)-(nterm)
	eff=effects(mod)  #get effect for each coefficient
	#eff=eff[mod.mm.X.assign.>reduce]  #get effects for terms not in null model
    eff=eff[mod.mm.assign.>reduce]  #get effects for terms not in null model
	return (sum(Abs2Fun(),eff),length(eff))  #returns tuple (ssh,dfh)
end

#these come from GLM package  can using GLM.dv
dv(da::DataArray) = da.data # return the data values - move this to DataFrames?
dv(pda::PooledDataArray) = int64(pda.refs)
dv{T<:Number}(vv::Vector{T}) = vv

function lmtest(form::Formula,nterm::Int,df::AbstractDataFrame)
	#does full vs reduced model ftest
	#nterm is the number terms (not necessarily DFs) to test as the end of the model
	#for instance y~a+b+c+d+c&d  is want to test c&d int then nterm=1, if want test of
	#c+d+c&d then nterm=3   it could be more than 3 DF if any term has factors
	#requires using DataFrame and using GLM
	#a=lm(form,df)
    a=fit(LinearModel,form,df)
	SSH,dfH = lmssh(a,nterm)
	#SSE=deviance(a.rr)
    SSE=deviance(a.model.rr)
	#dfE=df_residual(a)
    dfE=df_residual(a.model.pp)
	fstat=(SSH/dfH)/(SSE/dfE)
	return [fstat,dfH,dfE, -log10(ccdf(FDist(dfH, dfE), fstat))] #fstat,dfH,dfE,-log10pval
end

type dGWAS
    resref::Array{RemoteRef,1}  #remote refs to dataframe results split onto each
    testtype::Symbol  #so far either :linear, :logistic, or :LRTmv
    form::Formula    # formula
    genlimit::Int  #if >0 then any genotype with count<genlimit is made NA, if <2 geno
    asfactor::Bool  #if snp is treated as default additive (false) or as genotypic factors (true)
    nterm::Int   # # of terms from end of formula tested
    snpidx::Array{Int,2}  #c1=#snps for each partition c2=index start,c3=index finish for that partition
    function dGWAS(resref::Array{RemoteRef,1},testtype::Symbol,form::Formula,genlimit::Int,asfactor::Bool,nterm::Int)
        #c1=#snps for each partition c2=index start,c3=index finish for that partition
        snpidx=zeros(Int,length(resref),3)
        snpidx[1,1]=fetch(@spawnat resref[1].where size(fetch(resref[1]),1))
        snpidx[1,2]=1
        snpidx[1,3]=snpidx[1,1]
        for i=2:length(resref)
            snpidx[i,1]=fetch(@spawnat resref[i].where size(fetch(resref[i]),1))
            snpidx[i,2]=snpidx[i-1,3]+1
            snpidx[i,3]=snpidx[i-1,3]+snpidx[i,1]
        end
        new(resref,testtype,form,genlimit,asfactor,nterm,snpidx)
    end
end

#function to work with distributed GenDat
function gwLM(form::Formula,nterm::Int,gen::dGenDat;responsetype::Symbol=:linear,gencut::Int=0,asfactor::Bool=false,mafcut::FloatingPoint=0.0)	
    println("time started")
    run(`date`)
    nsnp=numsnp(gen)  #get actual number of snps
    resref=Array(RemoteRef,length(gen.refs))
    for i=1:length(gen.refs) resref[i] = @spawnat gen.refs[i].where gwLM(form,nterm,fetch(gen.refs[i]),responsetype=responsetype,gencut=gencut,asfactor=asfactor,mafcut=mafcut) end
    #probably not the best way to do stack dataframe results on top of each other
    #res=fetch(resref[1])
    #for i=2:length(gen.refs) res=vcat(res,fetch(resref[i])) end    
    #return res
    println("time finished")
    run(`date`)
    return dGWAS(resref,responsetype,form,gencut,asfactor,nterm)
end

#fixes weid issue with :macrocall change in formula
##used in gwLM in dgenglm6.jl and in gwLRTmv dgwLRTmv2.jl
function allvars2(ex::Expr)
    #if ex.head != :call error("Non-call expression encountered") end
    [[DataFrames.allvars(a) for a in ex.args[2:end]]...]
end

#function to work on GenDat
function gwLM(form::Formula,nterm::Int,gen::GenDat;responsetype::Symbol=:linear,gencut::Int=0,asfactor::Bool=false,mafcut::FloatingPoint=0.0)
    #formula is same as would use for LM
    #in formula put :snp where you want the snp to occur
    #nterm is how many terms you want to test for siginficance (place the terms to be tested at end of formula)
    #nterm like snp + snp&aging, nterm=2 , #nterms is not always same as df for test if one if factor
    if asfactor     #if asfactor then return snp as PooledDataArray otherwise sent out as continuous DataArray
        togen=btGenFac
    else
        togen=btGenDat
    end
    
    if responsetype==:linear
        testfunction = lmtest
        resnum=4
    elseif responsetype==:logistic
        testfunction = logistictest
        resnum=3
    else
        error("Response type was neither linear or logistic")
    end
    
    vars=DataFrames.allvars(form)
    #vars=allvars(form)
    vars=vars[vars.!=:snp]
    mdat=gen.fam[:,eval(vars)];  #pulls of variables from dgataframe
    covmiss=complete_cases(mdat);
    mdat=mdat[covmiss,:];   #this is base dataframe to use
    mdat[:snp]=ones(size(mdat,1));  #place holder for snp
    famindskeep=gen.fam[:indx].data[covmiss]
    results=ones(Float64,numsnp(gen),resnum);   #4 columns for linear, 3 for logistic
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
            usetable, dotest =getgentable(ct,gen.snp[:direction][i],gencut);  #gives gentable to use and dotest which is t/f to run or not
            if dotest
                #do mdat with new gen table cc[1]
                mdat[:snp]=togen(p,famindskeep,usetable)
                #mdat["snp"]=togen(p,faminds,direction[j],gencut);  #figure out
                #results[i,:]=testfunction(form,nterm,mdat);
                results[i,:]=   try
                                    testfunction(form,nterm,mdat);   
                                catch
                                    zeros(Float64,resnum);
                                end
            else
                results[i,:]=zeros(Float64,resnum);
            end
            ct=[0,0,0,0];
        end
        if mod(i,100000)==0  println("$i done out of $nsnp") end
    end
    if responsetype == :linear
        results=convert(DataFrame,results)
        names!(results,[:fstat,:dfH,:dfE,:log10pval])
        res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]],results)
       #res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]]
        #,DataFrame(results,[:fstat,:dfH,:dfE,:log10pval]))
    elseif responsetype == :logistic
        results=convert(DataFrame,results)
        names!(results,[:chi,:df,:log10pval])
        res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]],results)
        #res=hcat(gen.snp[[:chr,:loc,:snpID,:allele1,:MAF]]
        #,DataFrame(results,[:chi,:df,:log10pval]))
    end
    #make sure minor allele is corrects
    for i = 1:numsnp(gen) 
        if gen.snp[:direction][i] res[:allele1][i] = gen.snp[:allele2][i] end 
    end
    #changes name from allele 1 to minor allele which genotype coding is # of minor alleles, 0,1,2
    if responsetype == :linear
       names!(res,[:chr,:loc,:snpID,:minor,:MAF,:fstat,:dfH,:dfE,:log10pval])
    elseif responsetype == :logistic
        names!(res,[:chr,:loc,:snpID,:minor,:MAF,:chi,:df,:log10pval])
    end
    return res
end

function gencut(mdat::DataFrame)
    count=[0,0,0]
    gen=[0,1,2]
    for i=1:size(mdat,1)
        if !mdat[:snp].na[i] count[mdat[:snp].data[i]+1]+=1 end
    end
end

numsnp(gw::dGWAS)=sum(gw.snpidx[:,1]);  #the number of snps analysed

function getsigs(siglevel::FloatingPoint,gw::dGWAS)
    #gets a dataframe for all sig snps from all processors
    #siglevel is for -log10pval
    if gw.testtype==:LRTmv
        #sigres=@fetchrom gw.resref[1].where fetch(gw.resref[1])[:((LRTmvPval.>siglevel) | (LRTvPval.>siglevel) ),:])
        #cannot get the expression conditional things to work in remote ref so have to do LRTmvPval & LRTvPval separate and stack on each other
        sigres=@fetchfrom gw.resref[1].where fetch(gw.resref[1])[fetch(gw.resref[1])[:LRTmvPval].>siglevel,:]
        sigres=vcat(sigres,@fetchfrom gw.resref[1].where fetch(gw.resref[1])[fetch(gw.resref[1])[:LRTvPval].>siglevel,:])

    else
        #sigres=fetch(@spawnat gw.resref[1].where fetch(gw.resref[1])[:(log10pval.>siglevel),:])
        sigres=@fetchfrom gw.resref[1].where fetch(gw.resref[1])[fetch(gw.resref[1])[:log10pval].>siglevel,:]
    end
    for i=2:length(gw.resref)
        if gw.testtype==:LRTmv
            #sigres=vcat(sigres, @fetchfrom gw.resref[i].where fetch(gw.resref[i])[:((LRTmvPval.>siglevel) | (LRTvPval.>siglevel) ),:])
            sigres=vcat(sigres,@fetchfrom gw.resref[i].where fetch(gw.resref[i])[fetch(gw.resref[i])[:LRTmvPval].>siglevel,:])
            sigres=vcat(sigres,@fetchfrom gw.resref[i].where fetch(gw.resref[i])[fetch(gw.resref[i])[:LRTvPval].>siglevel,:])
        else
            #sigres=vcat(sigres,fetch(@spawnat gw.resref[i].where fetch(gw.resref[i])[:(log10pval.>siglevel),:]))
            sigres=vcat(sigres,@fetchfrom gw.resref[i].where fetch(gw.resref[i])[fetch(gw.resref[i])[:log10pval].>siglevel,:])
        end
    end
    return sigres
end

#this brings all of results to the main processor
function getresults(gw::dGWAS)
    res=fetch(gw.resref[1])
    for i=2:length(gw.resref) res=vcat(res,fetch(gw.resref[i])) end
    return res   
end

#takes pieces results from each piece of dGWAS and appends them to files
#function writeresults(outfile::String,gw::dGWAS)
#    @spawnat gw.resref[1].where writetable(outfile,fetch(gw.resref[1]),separator='\t')
#    #for i=2:length(gw.resref) @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i])) end
#    for i=2:length(gw.resref) timedwait(@spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i])),10) end
#end

#takes pieces results from each piece of dGWAS and appends them to files
#function writeresults(outfile::String,gw::dGWAS)
#    @spawnat gw.resref[1].where timedwait(writetable(outfile,fetch(gw.resref[1]),separator='\t'),10.0)
#    #for i=2:length(gw.resref) @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i])) end
#    for i=2:length(gw.resref) @spawnat gw.resref[i].where timedwait(writeresults(outfile,fetch(gw.resref[i])),10) end
#end

#takes pieces results from each piece of dGWAS and appends them to files
#function writeresults(outfile::String,gw::dGWAS)
#    @spawnat gw.resref[1].where timedwait(writetable(outfile,fetch(gw.resref[1]),separator='\t'),10.0)
#    sleep(5)
#    #for i=2:length(gw.resref) @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i])) end
#    for i=2:length(gw.resref) 
#        @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i]))
#    sleep(5)
#    end
#end




function writeresults(outfile::String,gw::dGWAS)
    @spawnat gw.resref[1].where writetable(outfile,fetch(gw.resref[1]),separator='\t')
    nn=gw.snpidx[1,1]
    sleep(10)
    for i=2:length(gw.resref)
       if countlines(outfile)==nn+1
            @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i]))
            nn=sum(gw.snpidx[1:i,:])
        else
            sleep(10)
            @spawnat gw.resref[i].where writeresults(outfile,fetch(gw.resref[i]))
            nn=sum(gw.snpidx[1:i,:])
        end
    end
end

#used in writeresults for dGWAS
function writeresults(outfile::String,df::DataFrame)
    io=open(outfile, "a")
    printtable(io,df,separator='\t',quotemark='"',header=false)
    close(io)
    return true
end



