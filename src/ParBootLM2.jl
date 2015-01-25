type BootPieceLM
	nullres::Array{Float64,1}
	nullfit::Array{Float64,1}
	yperm::Array{Float64,1}
	fullX::Array{Float64,2}
	dfH::Int
	ostat::Float64
end

type dParBootLM
	refs::Array{RemoteRef,1}  #array of RemoteRefs for BootPiece types on each processor each holding subset of bootstraps
	nboot::Int
	pieces::Array{Int,2}  #number of bootstraps to do on each piece #c1 = #boots per ref, c2=begin position in array c3=end
	datsend::BootPieceLM
end

type BootsLM
	dist::Array{Float64,1}
	nboot::Int   #number of bootstraps tried
	nworked::Int  #number of bootstraps worked
	nfail::Int  #number of bootstraps failed
end


function ParBootLM(form::Formula,df1::AbstractDataFrame,terms::Array{Int,1},nboot::Int=10000;
					cuttoff::Array{Float64,1}=[-1.0],gencut::Int=0,
					asfactor::Bool=false,mafcut::FloatingPoint=0.0)	
	#uterms=DataFrames.Terms(form).eterms  #gets unique terms in formula
	#uterms=convert(Array{Symbol,1},uterms)
	uterms=DataFrames.allvars(form)  #gets unique terms in formula
	df=df1[uterms] # makes new dataframe with only cols needed for the formula
	complete_cases!(df)
	fullmod=fit(LinearModel,form,df);
	ind=findin(fullmod.mm.assign,terms)
	keep=trues(length(fullmod.mm.assign));
	for i=1:length(ind) keep[ind[i]]=false end
	cc=deepcopy(fullmod.model.pp.X);
	cc=cc[:,keep];   #cuts out last col or cols associated with snp at the end
	nullmod=fit(LinearModel,cc,fullmod.model.rr.y);
	nullfit=deepcopy(nullmod.rr.mu)  #get fitted values for null model for LRTmv (only covariates no snp or vars)
	nullres=residuals(nullmod)  #gets resids for LRTmv null
	SSE=deviance(fullmod.model.rr)
	SSR=deviance(nullmod.rr)
	SSH=SSR-SSE
	dfE=df_residual(nullmod.pp)
	dfH=size(fullmod.model.pp.X,2)-size(cc,2)
	fstat=(SSH/dfH)/(SSE/dfE)
	ostat= -log10(ccdf(FDist(dfH, dfE),fstat))
	yperm=deepcopy(fullmod.model.rr.y)
	datsend=BootPieceLM(nullres,nullfit,yperm,fullmod.model.pp.X,dfH,ostat)
	dparts=dParBootLM(datsend,nboot);
	bootrefs=Array(RemoteRef,length(dparts.refs))
	for i=1:length(dparts.refs) bootrefs[(i)]= @spawnat dparts.refs[i].where ParBootLM(fetch(dparts.refs[i]),dparts.pieces[i,1]) end

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
	pvalcuts=[0.05,0.01,(1/nwork)*100,(1/nwork)*10]

	#if any specific pvalues are imputed, it is put at the end
	if cuttoff[1]!=-1.0 pvalcuts=[pvalcuts,cuttoff] end


	vals=tdist[floor(pvalcuts.*nwork)]

	#vals=quantile(tdist,1.-pvalcuts);
	vals=hcat(pvalcuts,vals)
	pval=sum(ostat.<=tdist)/nwork;
	return(pval,vals,nboot,nwork,nfail)
end

function dParBootLM(datsend::BootPieceLM,nboot::Int)
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
	return dParBootLM(refs,nboot,pieces,datsend);
end

function ParBootLM(dat::BootPieceLM,nboot::Int)
	dist=zeros(Float64,nboot);	
	nfail=0  #number of bootstraps that failed
	for i=1:nboot
		permboot!(dat)  #gets new yperm
		dist[i]=   try
                      ParBootLM(dat);   
                    catch
                      -1.0                  
                    end
        if dist[i]== -1.0 nfail+=1 end       
	end
	nwork=nboot-nfail  # number of bootstraps that actually worked
	return BootsLM(dist,nboot,nwork,nfail)
end

#shuffles residuals, adds them to fit to make new yperm
function permboot!(dat::BootPieceLM)
	dat.nullres=sample(dat.nullres,length(dat.nullres),replace=false);
	for i=1:length(dat.yperm) dat.yperm[i]=dat.nullfit[i]+dat.nullres[i] end
end

function ParBootLM(dat::BootPieceLM)
	mod=fit(LinearModel,dat.fullX,dat.yperm)
	dfE=df_residual(mod.pp)
	MSE=deviance(mod.rr)/dfE
	MSH=sum(Abs2Fun(),effects(mod)[(end-(dat.dfH-1)):end])/dat.dfH  #this takes last dfH cols effects and get SS/dfH which = MSH with dfH
	permpval=-log10(ccdf(FDist(dat.dfH, dfE), MSH/MSE))
	return permpval
end