import Distributions.Bernoulli

type BootPieceLog
	nullprobs::Array{Bernoulli,1}	
	yperm::Array{Float64,1}
	fullX::Array{Float64,2}
	nullX::Array{Float64,2}
	dfH::Int
	ostat::Float64
end

type dParBootLog
	refs::Array{RemoteRef,1}  #array of RemoteRefs for BootPiece types on each processor each holding subset of bootstraps
	nboot::Int
	pieces::Array{Int,2}  #number of bootstraps to do on each piece #c1 = #boots per ref, c2=begin position in array c3=end
	datsend::BootPieceLog
end

type BootsLog
	dist::Array{Float64,1}
	nboot::Int   #number of bootstraps tried
	nworked::Int  #number of bootstraps worked
	nfail::Int  #number of bootstraps failed
end


function ParBootLog(form::Formula,df1::AbstractDataFrame,terms::Array{Int,1},nboot::Int=10000;
					cuttoff::Array{Float64,1}=[-1.0],gencut::Int=0,
					asfactor::Bool=false,mafcut::FloatingPoint=0.0)	

	uterms=DataFrames.allvars(form)  #gets unique terms in formula
	df=df1[uterms] # makes new dataframe with only cols needed for the formula
	complete_cases!(df)
	fullmod=fit(GeneralizedLinearModel,form,df,Binomial());


	fullX=fullmod.model.pp.X

	ind=findin(fullmod.mm.assign,terms)
	keep=trues(length(fullmod.mm.assign));
	for i=1:length(ind) keep[ind[i]]=false end
	nullX=deepcopy(fullmod.model.pp.X);
	nullX=nullX[:,keep];   #cuts out last col or cols associated with snp at the end
	nullmod=fit(GeneralizedLinearModel,nullX,fullmod.model.rr.y);

	fitvals=predict(nullmod);
	nullprobs=Array(Bernoulli,length(fitvals));
	for i=1:length(fitvals) nullprobs[i]=Bernoulli(fitvals[i]) end
	yperm=zeros(Int,length(fitvals));


	chistat=deviance(nullmod)-deviance(fullmod)
	dfH= df_residual(nullmod.model.pp) - df_residual(fullmod.model.pp)
	ostat= -log10(ccdf(Chisq(dfH), chistat))

	datsend=BootPieceLM(nullprobs,yperm,fullX,nullX,dfH,ostat)
	dparts=dParBootLog(datsend,nboot);
	bootrefs=Array(RemoteRef,length(dparts.refs))
	for i=1:length(dparts.refs) bootrefs[(i)]= @spawnat dparts.refs[i].where ParBootLog(fetch(dparts.refs[i]),dparts.pieces[i,1]) end

	nwork=0
	nfail=0
	for i=1:length(bootrefs) nwork+=fetch(@spawnat bootrefs[i].where fetch(bootrefs[i]).nworked) end
	for i=1:length(bootrefs) nwork+=fetch(@spawnat bootrefs[i].where fetch(bootrefs[i]).nfail) end

	#make total distributions and retrieve them from each processor
	tdist=zeros(Float64,nboot);
	for i=1:length(bootrefs) tdist[dparts.pieces[i,2]:dparts.pieces[i,3]]=fetch(@spawnat bootrefs[i].where fetch(bootrefs[i]).dist) end

	sort!(tdist,rev=true);
	tdist=tdist[1:nwork];
	#different pvalues to get threshold values in the distribution
	#the last on is where we are comfortable calling a pvalue given the # of bootstraps
	pvalcuts=[0.05,0.01,(1/nwork)*100]

	#if any specific pvalues are imputed, it is put at the end
	if cuttoff[1]!=-1.0 pvalcuts=[pvalcuts,cuttoff] end


	vals=tdist[floor(pvalcuts.*nwork)]

	#vals=quantile(tdist,1.-pvalcuts);
	vals=hcat(pvalcuts,vals)
	pval=sum(ostat.<=tdist)/nwork;
	return(pval,vals,nboot,nwork,nfail)
end


function dParBootLog(datsend::BootPieceLog,nboot::Int)
	npuse=nworkers()
	if npuse==1 error("Only made for multiple processors") end
	refs=Array(RemoteRef,npuse)
	pieces=rep(int(floor(nboot/npuse)),npuse)
	pieces[end]=nboot-(pieces[1]*(npuse-1))
	idx=zeros(Int,npuse,2);
	idx[1,1]=1
	idx[1,2]=pieces[1]
	for i=2:npuse
		idx[i,1]=idx[i-1,2]+1
		idx[i,2]=idx[i-1,2]+pieces[i]
	end
	#c1 = #boots per ref, c2=begin position in array c3=end
	pieces=hcat(pieces,idx);
	for i in workers() refs[(i-1)]= @spawnat i datsend end
	return dParBootLog(refs,nboot,pieces,datsend);
end

function ParBootLog(dat::BootPieceLog,nboot::Int)
	ngreater=0
	dist=zeros(Float64,nboot);	
	nfail=0  #number of bootstraps that failed
	for i=1:nboot
		yperm!(dat)  #gets new yperm
		dist[i]=   try
                      ParBootLog(dat);   
                    catch
                      -1.0                  
                    end
        if dist[i]== -1.0 nfail+=1 end
	end
	nwork=nboot-nfail  # number of bootstraps that actually worked
	return BootsLog(dist,nboot,nwork,nfail)
end

yperm!(dat::BootPieceLog)= for i=1:length(dat.yperm) dat.yperm[i]=rand(dat.nullprobs[i]) end
yperm!(nullprobs::Array{Bernoulli,1},yerm::Array{Int,1}) = for i=1:length(yperm) yperm[i]=rand(nullprobs[i]) end


function ParBootLog(dat::BootPieceLog)
	full=fit(GeneralizedLinearModel,dat.fullX,dat.yperm,Binomial())
	null=fit(GeneralizedLinearModel,dat.nullX,dat.yperm,Binomial())
	chistat=deviance(null)-deviance(full)
	return -log10(ccdf(Chisq(dat.dfH), chistat))
end

