#get snp
#get covar & main

#define model get null and 

#start by getting original model
#get fit and resids for null model
#create new by by shuffling resids and adding back to fitted values
#for LRTmv that means mNull
#for LRTv that means mM

#MeanVar type has minfo with with mod for linear model with mean effects
#can use this to get both

type BootPiece
	dd::MeanVar
	ostat::Float64
	nullmod::LinearModel{DensePredQR{Float64}}
	nullfit::Array{Float64,1}
	nullres::Array{Float64,1}
	yperm::Array{Float64,1}  #has original y then is place holder for permuted y's
	test::Symbol
	nsnpcoef::Int  #number of columns for snp test (1 if additive up to two if genotypic)
end

type dParBoot
	refs::Array{RemoteRef,1}  #array of RemoteRefs for BootPiece types on each processor each holding subset of bootstraps
	test::Symbol  #whether LRTmv or LRTv
	nboot::Int
	pieces::Array{Int,2}  #number of bootstraps to do on each piece #c1 = #boots per ref, c2=begin position in array c3=end
	datsend::BootPiece
end

type Boots
	dist::Array{Float64,1}
	nboot::Int   #number of bootstraps tried
	nworked::Int  #number of bootstraps worked
	nfail::Int  #number of bootstraps failed
end

function ParBootLRT(form::Formula,df::AbstractDataFrame,nboot::Int=10000; test::Symbol=:LRTmv,cuttoff::Array{Float64,1}=[-1.0],AddorGen::Symbol=:Add,genlimit::Int=1,algorith::Symbol=:LN_BOBYQA)
	#do full model get LRTmv & LRTv to start
	dd=LRTmv(form,df,AddorGen,genlimit,algorith);
	cc=deepcopy(dd.minfo.mod.model.pp.X);  #copies X mat than includes snp col or cols at end
	cc=cc[:,dd.minfo.mod.mm.assign.!=maximum(dd.minfo.mod.mm.assign)];   #cuts out last col or cols associated with snp at the end
	yperm=deepcopy(dd.minfo.mod.model.rr.y);
	
	nullmod=fit(LinearModel,cc,yperm);  #fits absolute null model (just covariate no snp or vars)
	if test==:LRTmv
		ostat=dd.LRTmv.pval
		nullfit=deepcopy(nullmod.rr.mu)  #get fitted values for null model for LRTmv (only covariates no snp or vars)
		nullres=residuals(nullmod)  #gets resids for LRTmv null
	elseif test==:LRTv
		ostat=dd.LRTv.pval
		nullfit=deepcopy(dd.minfo.mod.model.rr.mu)  #get fitted values for cov+snp no var
		nullres=residuals(dd.minfo.mod)   #get residuals values for cov+snp no var
	elseif test==:LRTm
		LRTm!(dd);  #do LRTm test
		ostat=dd.LRTm.pval
		nullfit=mod.pp.X*dd.mV.betas  #fitted for mV, multiplies betas from mV to X mat from model with just covs
		nullres=yperm-nullfit  #residuals for mV (covariates + vars)
	else
		error("test can only take the form of :LRTmv, :LRTv, or :LRTm")
	end
	
	ntestcols=size(dd.minfo.mod.model.pp.X,2)-size(cc,2)
	datsend=BootPiece(dd,ostat,nullmod,nullfit,nullres,yperm,test,ntestcols);
	dparts=dParBoot(datsend,nboot,test);
	bootrefs=Array(RemoteRef,length(dparts.refs))
	for i=1:length(dparts.refs) bootrefs[(i)]= @spawnat dparts.refs[i].where ParBootLRT(fetch(dparts.refs[i]),dparts.pieces[i,1]) end

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




function dParBoot(datsend::BootPiece,nboot::Int,test::Symbol)
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
	return dParBoot(refs,test,nboot,pieces,datsend);
end

function ParBootLRT(dat::BootPiece,nboot::Int)
	#if dat.test!=:LRTmv error("not LRTmv test in BootPiece") end
	if dat.test==:LRTmv 
		Parbootfun=ParbootLRTmv
	elseif dat.test==:LRTv
	 	Parbootfun=ParbootLRTv
	else
		Parbootfun=ParbootLRTm
	end
	ngreater=0
	dist=zeros(Float64,nboot);	
	nfail=0  #number of bootstraps that failed
	for i=1:nboot
		permboot!(dat)  #gets new yperm
		dist[i]=   try
                      Parbootfun(dat);   
                    catch
                      -1.0                  
                    end
        if dist[i]== -1.0 nfail+=1 end
        if dist[i]>=dat.ostat ngreater+=1 end
	end
	nwork=nboot-nfail  # number of bootstraps that actually worked
	return Boots(dist,nboot,nwork,nfail)
end


#shuffles residuals, adds them to fit to make new yperm
function permboot!(dat::BootPiece)
	dat.nullres=sample(dat.nullres,length(dat.nullres),replace=false);
	for i=1:length(dat.yperm) dat.yperm[i]=dat.nullfit[i]+dat.nullres[i] end
end

function ParbootLRTmv(dat::BootPiece)
	#mod=lm(form,df)
	mod=fit(LinearModel,dat.dd.minfo.mod.model.pp.X,dat.yperm);  #fit with new yperm and original X, cov+snp no var
	#covM=scale(mod,true)
	covM=scale(mod,true);
	nullloglik,nparnull = bootlogliknull(mod,dat.nsnpcoef)  #get logliklihood for null model  (just covariates and no var het)
	MVloglik, nparMV=bootloglikMV(mod,dat.dd.minfo.gt,dat.dd.minfo.algorith)
	chi=nullloglik - MVloglik  #the likelihood -0.5(nlog(2pi)+sum(log(sigma))+ res'*1/sigma*res )
	df= nparMV - nparnull
	return -log10(ccdf(Chisq(df),chi))
end

function ParbootLRTv(dat::BootPiece)
	#mod=lm(form,df)
	mod=fit(LinearModel,dat.dd.minfo.mod.model.pp.X,dat.yperm)  #fit with new yperm and original X
	#covM=scale(mod,true)
	covM=scale(mod,true)
	#likM=size(mod.pp.X.m,1)*log(covM) + deviance(mod.rr)/covM
	Mloglik=size(mod.pp.X,1)*log(covM) + deviance(mod.rr)/covM
	#mM=mvmod(likM,coef(mod),[covM],true)   #get mvmod type for mean model (but no var het)
	nparM=size(mod.pp.X,2)+1  # # of parameters for model with snp but not vars = number of coefs + 1 var
	MVloglik, nparMV=bootloglikMV(mod,dat.dd.minfo.gt,dat.dd.minfo.algorith)
	chi=Mloglik - MVloglik  #the likelihood -0.5(nlog(2pi)+sum(log(sigma))+ res'*1/sigma*res )
	df= nparMV - nparM
	return -log10(ccdf(Chisq(df),chi))	
end

function ParbootLRTm(dat::BootPiece)
	#mod=lm(form,df)
	mod=fit(LinearModel,dat.dd.minfo.mod.pp.X,dat.yperm)  #fit with new yperm and original X
	#covM=scale(mod,true)
	covM=scale(mod,true)
	#likM=size(mod.pp.X.m,1)*log(covM) + deviance(mod.rr)/covM
	MVloglik, nparMV=bootloglikMV(mod,dat.dd.minfo.gt,dat.dd.minfo.algorith)
	Vloglik, nparV=bootloglikV(mod,dat)
	chi=Vloglik - MVloglik  #the likelihood -0.5(nlog(2pi)+sum(log(sigma))+ res'*1/sigma*res )
	df= nparMV - nparV
	return -log10(ccdf(Chisq(df),chi))	
end

function bootlogliknull(mod::LinearModel{DensePredQR{Float64}},nsnpcoef::Int)
	eff=effects(mod)  #get effect for each coefficient
	eff=eff[(size(mod.pp.X,2)-(nsnpcoef-1)):end]  #gets effects for cols related to snp
	ssh=sum(Abs2Fun(),eff)
	dfh=length(eff)	
	nullsse=deviance(mod.rr)+ssh
	nullvar= nullsse/(df_residual(mod.pp)+dfh)
	liknull=size(mod.pp.X,1)*log(nullvar) + nullsse/nullvar
	return (liknull, ((size(mod.pp.X,2)+1)-nsnpcoef)) #returns liklihook for null model and #parameters ncoef+var for null
end

# get effects for model terms from vector for linear models Q'y = effects matrix
#effects(mod::LinearModel{DensePredQR{Float64}})=(mod.pp.qr[:Q]'*mod.rr.y)[1:size(mod.pp.X,2)]

#function loglikMV(mod::LinearModel,gt::Vector,algorith::Symbol)
function bootloglikMV(mod::LinearModel{DensePredQR{Float64}},gt::Vector,algorith::Symbol)
	#y=mod.rr.y   #get y
	#xmat=mod.pp.X.m   #get model matrix (includes vector of 1's)
	xmat=mod.pp.X   #get model matrix (includes vector of 1's)
	varinit=scale(mod,true)  #get sigma^2
	ugen=sort(unique(gt))
	ngen=length(ugen)
	parinit=[coef(mod),rep(varinit,ngen)]
	varvec= rep(1.0,length(gt))  #this vector gets updated in the obj function
	k=length(parinit)  #number of parameters
	opt = Opt(algorith,k)
	ftol_abs!(opt, 1e-6)    # criterion on deviance changes
	#xtol_abs!(opt, 1e-6)    # criterion on all parameter value changes
	lower_bounds!(opt, [rep(-Inf,size(xmat,2)),rep(0.0,ngen)])
	#maxeval!(opt,200)
	function obj(x::Vector, g::Vector)
    	if length(g) > 0 error("gradient evaluations are not provided") end
    	#changes sigma from higher level
    	#uses xmat, ngen, ugen, and mod from higher level
   		#global count
   		#count += 1
   		#println("f_$count($x)")
   		#objsigma!(sigma)   #changes sigma from higher level so don't have t reallocate each time
   		#loglik=sum(log(sigma))+res'*diagm(1/sigma)*res
   		#uses numeric extensions
   		#    sum(log(sigma))   +      1/sigma  the r'*1/sigma*res
   		#loglik=sum(Log(),sigma)+wsum(rcp!(sigma),Multiply(),res,res)  #res = y-xmat*coef   	
   		for i=1:length(gt)    #changes sigma and uses gt, ngen, ugen from higher level
			for j=1:ngen
				if gt[i]==ugen[j]
					varvec[i]=x[end-(ngen-j)]
					#break
				end
			end
		end	
   		yhat=xmat*x[1:(end-ngen)]  
   		loglik=sum(LogFun(),varvec) + wsumfdiff(rcp!(varvec),Abs2Fun(),mod.rr.y,yhat)
   		return loglik
	end
	min_objective!(opt, obj)
	res = optimize(opt, parinit)
	return (res[1],length(parinit))  #returns logliklihood of MV and number of coefs (betas+vars)
end


function bootloglikV(mod::LinearModel{DensePredQR{Float64}},dat::BootPiece)
	#sends mM model (already caluculated and BootPiece did not want to recalc calculate null)
	#mod is absolute null model, just covariates no snp or vars
	xmat=dat.nullmod.pp.X   #get model matrix (includes vector of 1's)
	varinit=scale(mod,true)  #get sigma^2  this if from mM model but close enough
	gt=dat.dd.minfo.gt
	algorith=dat.dd.minfo.algorith
	ugen=sort(unique(gt))
	ngen=length(ugen)
	varvec= rep(1.0,length(gt))  #this vector gets updated in the obj function
	parinit=[coef(mod)[1:size(xmat,2)],rep(varinit,ngen)]  #these initial parameter come from null mod
	k=length(parinit)  #number of parameters
	opt = Opt(algorith,k)
	ftol_abs!(opt, 1e-6)    # criterion on deviance changes
	#xtol_abs!(opt, 1e-6)    # criterion on all parameter value changes
	lower_bounds!(opt, [rep(-Inf,size(xmat,2)),rep(0.0,ngen)])
	#maxeval!(opt,200)
	function obj(x::Vector, g::Vector)
    	if length(g) > 0 error("gradient evaluations are not provided") end
    	#changes sigma from higher level
    	#uses xmat, ngen, ugen, and mod from higher level
   		#global count
   		#count += 1
   		#println("f_$count($x)")
   		#objsigma!(sigma)   #changes sigma from higher level so don't have t reallocate each time
   		#loglik=sum(log(sigma))+res'*diagm(1/sigma)*res
   		#uses numeric extensions
   		#    sum(log(sigma))   +      1/sigma  the r'*1/sigma*res
   		#loglik=sum(Log(),sigma)+wsum(rcp!(sigma),Multiply(),res,res)  #res = y-xmat*coef   	
   		for i=1:length(gt)    #changes sigma and uses gt, ngen, ugen from higher level
			for j=1:ngen
				if gt[i]==ugen[j]
					varvec[i]=x[end-(ngen-j)]
					#break
				end
			end
		end	
   		yhat=xmat*x[1:(end-ngen)]  
   		loglik=sum(LogFun(),varvec) + wsumfdiff(rcp!(varvec),Abs2Fun(),dat.yperm,yhat)
   		return loglik
	end	
	min_objective!(opt, obj)
	res = optimize(opt, parinit)
	return (res[1],length(paranit))  #returns logliklihood of MV and number of coefs (betas+vars)
end