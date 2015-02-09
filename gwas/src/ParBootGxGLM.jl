
#uses BootsLM from ParBootLM2.jl

type BootPieceGxGLM
	nullres::Array{Float64,1}
	nullfit::Array{Float64,1}
	yperm::Array{Float64,1}
	fullX::Array{Float64,2}
	nullX::Array{Float64,2}
	ostat::Float64
end

type dParBootGxGLM
	refs::Array{RemoteRef,1}  #array of RemoteRefs for BootPiece types on each processor each holding subset of bootstraps
	nboot::Int
	pieces::Array{Int,2}  #number of bootstraps to do on each piece #c1 = #boots per ref, c2=begin position in array c3=end
	datsend::BootPieceGxGLM
end


function ParBootGxGLM(form::Formula,gen::dGenDat,snpname1::Symbol,snpname2::Symbol; nboot::Int=10000,
					  celllimit::Int=5,cuttoff::Array{Float64,1}=[-1.0])
    df=fetch(@spawnat gen.refs[1].where fetch(gen.refs[1]).fam);
    dnames=names(df);
    #check to see if snp is alread in the dataframe, if not put it there
    if !in(snpname1,dnames) df[snpname1]=getsnp(string(snpname1),gen) end;
    if !in(snpname2,dnames) df[snpname2]=getsnp(string(snpname2),gen) end;
    ParBootGxGLM(form,df,snpname1,snpname2,nboot=nboot,celllimit=celllimit,cuttoff=cuttoff)
end		


function ParBootGxGLM(form::Formula,df::DataFrame,snpname1::Symbol,snpname2::Symbol;nboot::Int=10000,
					  celllimit::Int=5,cuttoff::Array{Float64,1}=[-1.0])
    vars=DataFrames.allvars(form);
    push!(vars,snpname1,snpname2);
    mdat=df[:,eval(vars)];  #pulls of variables from dataframe
    covmiss=complete_cases(mdat);
    mdat=mdat[covmiss,:];   #this is base dataframe to use
	y=convert(Array{Float64,1},mdat[form.lhs]);  #get phe (dependent variable)
    results=ones(Float64,4);   #4 columns for linear, 3 for logistic
    covmat=ModelMatrix(ModelFrame(form,mdat)).m;  #this is model matrix with just the intercept and covariates
    const twogentable=[[0,0,0,1,1,1,2,2,2] [0,1,2,0,1,2,0,1,2] [0,0,0,0,0,0,0,0,0] [0,0,0,0,0,0,0,0,0]]
    snp1=convert(DataArray,mdat[snpname1]);
    snp1=convert(Array{Int8,1},snp1);
    twogenvec=zeros(Int,size(snp1,1));
    #c1=indicator var for snp1 het,c2=indicator snp2 rare hom, c3=snp2het ind,c4=snp2rarehom ind
    snpind=zeros(Float64,size(snp1,1),4);
    #fill in for first snp, second snp later
    for i in 1:size(snp1,1)
        if snp1[i]==1
            snpind[i,1]=1.0  #make indicator for homozygote
        elseif snp1[i]==2
            snpind[i,2]=1.0  #make indicator for rare homozygote
        end
    end
    snp2=convert(DataArray,mdat[snpname2]);
    snp2=convert(Array{Int8,1},snp2);
    cutind,tsnpind,snpgen,dtest=twogencount!(snp1,snp2,snpind,twogenvec,twogentable,celllimit);   
   	tempy=deepcopy(y);
    tempy=tempy[cutind];
    tempcovmat=deepcopy(covmat)[cutind,:];
    fullmat=hcat(tempcovmat,snpgen);
    reducemat=hcat(tempcovmat,tsnpind);
    stats=lmtest(tempy,fullmat,reducemat);
    ostat=floor(stats[4],10)  #gets -log10pval rounds down at 8 digits
    nullmod=fit(LinearModel,reducemat,tempy);
    nullfit=deepcopy(nullmod.rr.mu);  #get fitted values for null model
	nullres=residuals(nullmod);  #gets resids for null model
	yperm=deepcopy(tempy);
	datsend=BootPieceGxGLM(nullres,nullfit,yperm,fullmat,reducemat,ostat);  #puts together data for each piece
	dparts=dParBootGxGLM(datsend,nboot);  #sents the data onto each ref
	bootrefs=Array(RemoteRef,length(dparts.refs))  #refs to hold distributions
	#distibute datsend (BootPieceGGLM) on to refs
	for i=1:length(dparts.refs) bootrefs[i]= @spawnat dparts.refs[i].where ParBootGxGLM(fetch(dparts.refs[i]),dparts.pieces[i,1]) end

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
	if cuttoff[1]!=-1.0 	
		if cuttoff[1]>=1 pvalcuts=[pvalcuts,10.^-cuttoff] #if input is as -log10 it converts back
		else
			pvalcuts=[pvalcuts,cuttoff]
		 end		
	end

	pvallogs=-log10(pvalcuts)


	vals=zeros(Float64,length(pvalcuts))

	for i in 1:length(vals)
		x=floor(pvalcuts[i].*nwork)
		if x>0 vals[i]=tdist[x] else vals[i]=0.0 end
	end
	

	#vals=quantile(tdist,1.-pvalcuts);
	vals=hcat(pvalcuts,pvallogs,vals)
	PBpval=sum(ostat.<=tdist)/nwork;
	return(PBpval,vals,stats,nboot,nwork,nfail)
end



function dParBootGxGLM(datsend::BootPieceGxGLM,nboot::Int)
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
	return dParBootGxGLM(refs,nboot,pieces,datsend);
end

function ParBootGxGLM(dat::BootPieceGxGLM,nboot::Int)
	dist=zeros(Float64,nboot);	
	nfail=0  #number of bootstraps that failed
	for i=1:nboot
		permboot!(dat)  #gets new yperm
		dist[i]=   try
                      ParBootGxGLM(dat);   
                    catch
                      -1.0                  
                    end
        if dist[i]== -1.0 nfail+=1 end       
	end
	nwork=nboot-nfail  # number of bootstraps that actually worked
	return BootsLM(dist,nboot,nwork,nfail)
end

#shuffles residuals, adds them to fit to make new yperm
function permboot!(dat::BootPieceGxGLM)
	dat.nullres=sample(dat.nullres,length(dat.nullres),replace=false);
	for i=1:length(dat.yperm) dat.yperm[i]=dat.nullfit[i]+dat.nullres[i] end
end

function ParBootGxGLM(dat::BootPieceGxGLM)	
	fm=fit(LinearModel,dat.fullX,dat.yperm)
	SSEfull=deviance(fm.rr)
	dfE=df_residual(fm)
	rm=fit(LinearModel,dat.nullX,dat.yperm)
	SSH=deviance(rm.rr)-SSEfull
	dfH=df_residual(rm)-dfE
	fstat= (SSH/dfH) / (SSEfull/dfE)
	permpval = -log10(ccdf(FDist(dfH, dfE), fstat))
	return permpval  #returns -log10pval
end
