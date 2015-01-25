#for LRTmv and LRTm  may need to do
#some tests fail at the optimization step with a DomainError
#res= try 
#		LRTmv(form,mdat,AddorGen,genlimit);
#	catch y
#		if isa(y, DomainError)
			#[res.LRTmvChi,res.LRTmvDF,res.LRTmvPval,res.LRTvChi,res.LRTvDF,res.LRTvPval]
#			[1.0,1.0,1.0,1.0,1.0,1.0]  #make null numbers
			#println("test Failed")
#		end
#	end

#using GLM
#using NLopt
#using NumericExtensions

type mvinfo
	AddorGen::Symbol   #mean modelled as :Add or :Gen,  :Add is default, vars are always by genotype
	genlimit::Int
	algorith::Symbol  #algorithm for NLopt
	#mod::LinearModel    #original lm model with covariates & snp but not var het, related to mM
	mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64}  #original lm model with covariates & snp but not var het, related to mM
	gt::Vector    
end

type mvmod
	loglik::Float64     #log likelihood of model =  sum(log(sigma)) + res'*1/sigma*res  (excludes -1/2(nlog(2pi))
	betas::Array{Float64,1}
	vars::Array{Float64,1}
	fit::Bool
	ncoef::Int
	function mvmod(loglik::Float64,betas::Array{Float64,1},vars::Array{Float64,1},fit::Bool)
		new(loglik,betas,vars,fit,length(betas)+length(vars))
	end
end

type lrt
	chi::Float64
	df::Int
	pval::Float64   #reported as -log10pval
	function lrt(chi::Float64,df::Int)
		if chi != -1.0
		new(chi,df,-log10(ccdf(Chisq(df),chi)))
		else
		new(-1,1,1.1)  #sets to null with pval = 1.1
		end
	end
end

type MeanVar
	minfo::mvinfo
	mMV::mvmod  #produced in first go
	mM::mvmod   #produced in first go
	mNull::mvmod #produced in first go
	mV::mvmod	#needs to be invoked   #covariates with cov + var het
	LRTmv::lrt  #can do initially mMV & mNull
	LRTv::lrt   #can do initially mMV & mM
	LRTm::lrt   #need to do later
	function MeanVar(minfo::mvinfo, mMV::mvmod,mM::mvmod,mNull::mvmod)
		mV=mvmod(1.0,mM.betas[1:(end-(mM.ncoef-mNull.ncoef))],mMV.vars,false)  #get initial parameters from mM then set to false
		LRTmv=lrt(mMV,mNull)
		LRTv = lrt(mMV,mM)
		LRTm = lrt(-1.0,1)  #if first argument is -1 then sets to (-1,1,1.1)
		new(minfo,mMV,mM,mNull,mV,LRTmv,LRTv,LRTm)
	end
end

function lrt(full::mvmod,reduce::mvmod)
	#the real likelihood is -0.5(nlog(2pi)+sum(log(sigma))+ res'*1/sigma*res )
	# we use just sum(log(sigma))+ res'*1/sigma*res)  because -0.5(nlog(2pi) is the same for all models
	# as a result only need the difference  instead of -2(reduce - full)  we do full-reduce
	chi=reduce.loglik - full.loglik  #the likelihood -0.5(nlog(2pi)+sum(log(sigma))+ res'*1/sigma*res )
	df=full.ncoef - reduce.ncoef
	return lrt(chi,df)
end


#additive gen, genetic var
#function LRTmv(form::Expr,df::AbstractDataFrame, AddorGen::Symbol=:Add,genlimit::Int=1,algorith::Symbol=:LN_NEWUOA)
function LRTmv(form::Formula,df::AbstractDataFrame, AddorGen::Symbol=:Add,genlimit::Int=1,algorith::Symbol=:LN_BOBYQA)
	if genlimit>1  df[:snp]=cutgen!(df[:snp],genlimit)  end  #cuts out individuals with gencounts<genlimit
	
	#gt=df["snp"].data[complete_cases(df)]
	#get 0,1,2 vector for genotypes, all missing gone
	#if snp is no DataArray then it makes it one
	typeof(df[:snp])<:DataArray ? (gt=df[:snp].data[complete_cases(df)]) : (gt=DataArray(df[:snp]).data[complete_cases(df)])
	if AddorGen==:Gen df[:snp]=compact(PooledDataArray(df[:snp])) end	#if Gen then make pooledDataArray
	#mod=lm(form,df)
	mod=fit(LinearModel,form,df)
	#covM=scale(mod,true)
	covM=scale(mod.model,true)
	#likM=size(mod.pp.X.m,1)*log(covM) + deviance(mod.rr)/covM
	likM=size(mod.model.pp.X,1)*log(covM) + deviance(mod.model.rr)/covM
	mM=mvmod(likM,coef(mod),[covM],true)   #get mvmod type for mean model (but no var het)
	mNull = logliknull(mod,mM.ncoef)  #get mvmod type for null model  (just covariates and no var het)
	mMV=loglikMV(mod,gt,algorith)
	return MeanVar(mvinfo(AddorGen,genlimit,algorith,mod,gt),mMV,mM,mNull)
end



function LRTm!(mv::MeanVar)
	if mv.mV.fit == true
	else
	#xmat=mv.minfo.mod.pp.X.m[:,1:length(mv.mV.betas)]
	xmat=mv.minfo.mod.model.pp.X[:,1:length(mv.mV.betas)]
	gt=mv.minfo.gt
	ugen=sort(unique(gt))  #unique genotypes
	ngen=length(ugen)		# #unique genotypes
	varvec= rep(1.0,length(gt))  #this vector gets updated in the obj function
	parinit=[mv.mV.betas,mv.mV.vars]  #these initial parameter come from mM
	k=length(parinit)  #number of parameters
	opt = Opt(mv.minfo.algorith,k)
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
   		loglik=sum(LogFun(),varvec) + wsumfdiff(rcp!(varvec),Abs2Fun(),mv.minfo.mod.model.rr.y,yhat)
   		return loglik
	end	
	min_objective!(opt, obj)
	res = optimize(opt, parinit)
	mv.mV=mvmod(res[1],res[2][1:(end-ngen)],[res[2][(end-ngen+1):end]],true)
	mv.LRTm=lrt(mv.mMV,mv.mV)
	end
end


function cutgen!(df::AbstractDataFrame,genlimit::Int)
	use=complete_cases(df)
	counts=[0,0,0]
	for i=1:size(df,1)
		if use[i] counts[df[:snp][i]+1]=counts[df[:snp][i]+1]+1 end
	end
	keep=counts.>=genlimit
	if sum(keep)==3
		break
	else
		for j=1:size(df,1)
			if use[j] df[:snp].na[j]=keep[df[:snp][j]+1] end 
		end
	end
end

#
#null only excludes last term
#function logliknull(mod::LinearModel,ncoefM::Int)
function logliknull(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64},ncoefM::Int)
	ssh,dfh=lmssh(mod,1)
	#nullsse=deviance(mod.rr)+ssh
	nullsse=deviance(mod.model.rr)+ssh
	#nullvar= nullsse/(df_residual(mod.pp)+dfh)
	nullvar= nullsse/(df_residual(mod.model.pp)+dfh)
	#liknull=size(mod.pp.X.m,1)*log(nullvar) + nullsse/nullvar
	liknull=size(mod.model.pp.X,1)*log(nullvar) + nullsse/nullvar
	mNull=mvmod(liknull,[1.0],[nullvar],true)   #do not calculate beta coef, can get var
	mNull.ncoef=ncoefM-dfh
	return mNull::mvmod  #returns mvmod type for null model
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

# get effects for model terms from vector for linear models Q'y = effects matrix
#effects(mod::LinearModel)=(mod.pp.qr[:Q]'*mod.rr.y)[1:size(mod.pp.X.m,2)]
effects(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64})=(mod.model.pp.qr[:Q]'*mod.model.rr.y)[1:size(mod.model.pp.X,2)]


#function loglikMV(mod::LinearModel,gt::Vector,algorith::Symbol)
function loglikMV(mod::DataFrameRegressionModel{LinearModel{DensePredQR{Float64}},Float64},gt::Vector,algorith::Symbol)
	#y=mod.rr.y   #get y
	#xmat=mod.pp.X.m   #get model matrix (includes vector of 1's)
	xmat=mod.model.pp.X   #get model matrix (includes vector of 1's)
	varinit=scale(mod.model,true)  #get sigma^2
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
   		loglik=sum(LogFun(),varvec) + wsumfdiff(rcp!(varvec),Abs2Fun(),mod.model.rr.y,yhat)
   		return loglik
	end
	
	min_objective!(opt, obj)
	res = optimize(opt, parinit)
	mMV=mvmod(res[1],res[2][1:(end-ngen)],[res[2][(end-ngen+1):end]],true)
	return mMV::mvmod
end






#function obj(x::Vector, g::Vector)
#    	if length(g) > 0 error("gradient evaluations are not provided") end
#    	#global count
#    	#count += 1
#    	#println("f_$count($x)")
#    	sigma= rep(x[end-2],length(gt))
#    	for j=1:length(gt)
#    	   	if gt[j]==1
#    			sigma[j]=x[end-1]
#    		elseif gt[j]==2 
#    			sigma[j]==x[end]
#    		end
#    	end
#    	#loglik=sum(log(sigma))+res'*diagm(1/sigma)*res
#    	#uses numeric extensions
#    	#    sum(log(sigma))   +      1/sigma  the r'*1/sigma*res
#    	#loglik=sum(Log(),sigma)+wsum(rcp!(sigma),Multiply(),res,res)  #res = y-xmat*coef
#    	yhat=xmat*x[1:(end-3)]
#    	loglik=sum(Log(),sigma) + wsum_fdiff(rcp!(sigma),Abs2(),mod.rr.y,yhat)
#
#   	#println("loglik = $loglik")
#   	return loglik
#end
