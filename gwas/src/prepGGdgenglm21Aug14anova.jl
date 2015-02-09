#this does GxG with full & reduced models
#full = y=covariates + two-locus genotypes as factors
#reduce = y = covariates + snp1 + snp2
    #where snp1 & snp2 are factors for het & minor homozygote
#this takes into account missing data
#celllimit is the minimum count for a two-locus genotype
#all individuals with two-locus genotypes are counted as missing
#that means that that are between single locus effects and 1-4 interaction effects
#the full model explains the most that the two loci marginal & interaction effects can
#it does so without specifying the interaction effects
#the reduced model is just the covariates and the single locus effects

#function to work with distributed GenDat
function gwLMGxG(form::Formula,gen::dGenDat;responsetype::Symbol=:linear,celllimit::Int=5,
                  gencut::Int=0,asfactor::Bool=true,mafcut::FloatingPoint=0.0)
    #form::Formula only includes phe~covariates  
    println("time started")
    run(`date`)
    nsnp=numsnp(gen)  #get actual number of snps
    resref=Array(RemoteRef,length(gen.refs))
    for i=1:length(gen.refs) resref[i] = @spawnat gen.refs[i].where gwLMGxG(form,fetch(gen.refs[i]),responsetype=responsetype,gencut=gencut,asfactor=asfactor,mafcut=mafcut) end
    #probably not the best way to do stack dataframe results on top of each other
    #res=fetch(resref[1])
    #for i=2:length(gen.refs) res=vcat(res,fetch(resref[i])) end    
    #return res
    println("time finished")
    run(`date`)
    return dGWAS(resref,responsetype,form,gencut,asfactor,0)
end



#function to work on GenDat
function gwLMGxG(form::Formula,gen::GenDat;responsetype::Symbol= :linear,celllimit::Int=5,
                  gencut::Int=0,asfactor::Bool=true,mafcut::FloatingPoint=0.0)    
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
    mdat=gen.fam[:,eval(vars)];  #pulls of variables from dataframe
    covmiss=complete_cases(mdat);
    mdat=mdat[covmiss,:];   #this is base dataframe to use
    mdat[:snp]=ones(size(mdat,1));  #place holder for snp
    famindskeep=gen.fam[:indx].data[covmiss]
    results=ones(Float64,numsnp(gen),resnum);   #4 columns for linear, 3 for logistic
    bpt=convert(Ptr{Uint8},gen.gendat)
    m=size(gen.gendat,1)
    ct=[0,0,0,0];
    nsnp=numsnp(gen)
    

    y=convert(Array{Float64,1},mdat[form.lhs]);  #get phe (dependent variable)
    
    formcov=deepcopy(form)
    formcov.rhs.args=form.rhs.args[1:(end-3)]  #cut formula down to just covariates
    covmat=ModelMatrix(ModelFrame(formcov,mdat)).m  #this is model matrix with just the intercept and covariates

    snpnames=form.rhs.args[end].args[2:end]  #gets symbol names of interacting snps
    snpname1=snpnames[snpnames.!=:snp][1]  #gets name of interacting snp not called :snp
    #c1=snp1 genotype,c2=snp2 genotype, c3=2 loc genotype count,c4=(0=keep,1=cut)
    const twogentable=[[0,0,0,1,1,1,2,2,2] [0,1,2,0,1,2,0,1,2] [0,0,0,0,0,0,0,0,0] [0,0,0,0,0,0,0,0,0]]
    snp1=convert(DataArray,mdat[snpname1])
    snp1=convert(Array{Int8,1},snp1)
    twogenvec=zeros(Int,size(snp1,1))

    #c1=indicator var for snp1 het,c2=indicator snp2 rare hom, c3=snp2het ind,c4=snp2rarehom ind
    snpind=zeros(Float64,size(snp1,1),4)
    
    #fill in for first snp, second snp later
    for i in 1:size(snp1,1)
        if snp1[i]==1
            snpind[i,1]=1.0  #make indicator for homozygote
        elseif snp1[i]==2
            snpind[i,2]=1.0  #make indicator for rare homozygote
        end
    end

    for i=1:numsnp(gen)
        if gen.snp[i,:MAF]>=mafcut
            p=bpt + (gen.snp[:indx][i]-1)*m
            gencount!(ct,p,famindskeep);
            usetable, dotest =getgentable(ct,gen.snp[:direction][i],gencut);  #gives gentable to use and dotest which is t/f to run or not
            if dotest
                snp2=btGen(p,famindskeep,usetable);
                tsnpind=deepcopy(snpind)

               	cutind,tsnpind,snpgen,dtest=twogencount!(snp1,snp2,tsnpind,twogenvec,twogentable,celllimit);
               	if dtest
               		tempy=deepcopy(y);
                	tempy=tempy[cutind];
                	tempcovmat=deepcopy(covmat)[cutind,:];
                	fullmat=hcat(tempcovmat,snpgen);
                	reducemat=hcat(tempcovmat,tsnpind);

         			results[i,:]=   try
                                    	testfunction(tempy,fullmat,reducemat);   
                                	catch
                                    	#zeros(Float64,resnum);
                                        rep(-1.0,resnum); #not run because test failed
                                	end
                else
                	#results[i,:]=zeros(Float64,resnum);
                    results[i,:]=rep(-2.0,resnum);  #not run because of two loc counts (dtest)
                end

			else             
				#results[i,:]=zeros(Float64,resnum);
                results[i,:]=rep(-3.0,resnum);  #not run because of MAF (dotest)
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

#does not use dataframe, input y vector, and X matrices for full model and reduced models
function lmtest(y::Array{Float64,1},fullmat::Array{Float64,2},reducemat::Array{Float64,2})
	fm=fit(LinearModel,fullmat,y)
	SSEfull=deviance(fm.rr)
	dfE=df_residual(fm)
	rm=fit(LinearModel,reducemat,y)
	SSH=deviance(rm.rr)-SSEfull
	dfH=df_residual(rm)-dfE
	fstat= (SSH/dfH) / (SSEfull/dfE)
	return [fstat,dfH,dfE, -log10(ccdf(FDist(dfH, dfE), fstat))] #fstat,dfH,dfE,-log10pval
end

#does not use dataframe, input y vector, and X matrices for full model and reduced models
function logistictest(y::Array{Float64,1},fullmat::Array{Float64,2},reducemat::Array{Float64,2})
	fm=fit(GeneralizedLinearModel,fullmat,y,Binomial())  #full model
	devFull=deviance(fm)
    dfFull=df_residual(fm)
    rm=fit(GeneralizedLinearModel,reducemat,y,Binomial())  #reduced model
    dfH=df_residual(rm)-dfFull	  #df of hypotheses = dfReduced-dfHypothesis
    chistat=deviance(fm)-devFull  #chi stat = deviance of reduced - deviance of full
    return [chistat,dfH, -log10(ccdf(Chisq(dfH), chistat))] #Chi square,dfH,-log10(pval)
end


function twogencount!(snp1::Array{Int8,1},snp2::Array{Int8,1},snpind::Array{Float64,2},twogenvec::Array{Int,1},twogentable::Array{Int,2},celllimit::Int=5)
	#take two snp's coded at 0,1,2, second can have missing data (not 0,1,2)
	#gets 2 locus counts, use celllimit cut certain genotypes and individuals with those genotypes
	#mames indicator variables for het and minor hom for snp2 (snp1 is already done as input in tsnpind)
	#makes indicator for n-1 two locus genotypes with count >= celllimit
	#outputs bool array for individuals to keep (true) or cut (false), new float64 matrix indicator variable for
	#if less than 4 useable two locus genotypes then has test = false
	#it also modifies tsnpind by cutting out missing, creating genotype indicators for snp2, cuts out any indicators with count <celllimit
	cutind=trues(size(snp2,1))
	for i in 1:size(snp2,1)
		if snp2[i]==1
			snpind[i,3]=1.0  #make indicator for homozygote
        elseif snp2[i]==2
            snpind[i,4]=1.0  #make indicator for rare homozygote
        elseif snp2[i]== -1
        	cutind[i]=false  #tags individuals as missing(false)
        else 
        	snpind[i,3]=0.0  #makes sure these are zero'd out
        	snpind[i,4]=0.0
        end
   end

   for i = 1:size(twogentable,1)   #zeros out counts and cuts
		twogentable[i,3]=0
		twogentable[i,4]=0
	end

	for i in 1:size(twogenvec,1)   #makes vector with 1-9 for each twolocus genotype, 0 is missing
        twogenvec[i]=0  #zeros it out first
        for j in 1:9
            if (snp1[i]==twogentable[j,1]) && (snp2[i]==twogentable[j,2])
                twogentable[j,3]+=1
                twogenvec[i]=j
                break
            end
        end
    end

    tocut=[0]  #gets which two locus genotypes to cut
    tokeep=[0]  #gets which two locus genotypes to keep
    for j in 1:9 
        if twogentable[j,3]<celllimit 
            twogentable[j,4]=1
            push!(tocut,j)
        else
        	push!(tokeep,j)
        end
    end

	#labels all individuals with 2 loc genotypes that are cut as false 
    for i in 1:size(twogenvec,1)
        for j in 1:size(tocut,1)
            if twogenvec[i]==tocut[j] 
                cutind[i]=false
                break
            end
        end
    end

    ntwogen=sum(twogentable[:,4].==0)  #of two loc genotypes to keep
    if ntwogen>=4
    	tokeep=tokeep[3:end]  #gets rid of 0 a beginning an the first two-loc genotype to keep, leave n-1 2loc to keep
    	twogenveccut=twogenvec[cutind]  #cuts out missing
    	snpgen=zeros(Float64,size(twogenveccut,1),size(tokeep,1))
    	for i in 1:size(twogenveccut,1)
    		for j in 1:size(tokeep,1)
    			if twogenveccut[i]==tokeep[j] 
    				snpgen[i,j]=1.0
    				break
    			end
    		end
    	end
    	
    	#mcounts=vec(sum(snpind,1))
    	tsnpind=snpind[cutind,:]
        mcounts=vec(sum(tsnpind,1))
    	tsnpind=tsnpind[:,mcounts.>=celllimit]  #cuts out marginal vectors with < celllimit
    	dtest=true
    else
    	dtest=false
    	snpgen=zeros(Float64,2,2) #sends out dummy matrix
    	tsnpind=zeros(Float64,2,2) #sends out dummy matrix
    end
	return(cutind,tsnpind,snpgen,dtest)
end

function GxG(form::Formula,gen::dGenDat,snpname1::Symbol,snpname2::Symbol;
              responsetype::Symbol = :linear,celllimit::Int=5)
    df=fetch(@spawnat gen.refs[1].where fetch(gen.refs[1]).fam)
    dnames=names(df)
    #check to see if snp is alread in the dataframe, if not put it there
    if !in(snpname1,dnames) df[snpname1]=getsnp(string(snpname1),gen) end
    if !in(snpname2,dnames) df[snpname2]=getsnp(string(snpname2),gen) end
    GxG(form,df,snpname1,snpname2,responsetype= :linear,celllimit=celllimit)
end

function GxG(form::Formula,df::DataFrame,snpname1::Symbol,snpname2::Symbol;
             responsetype::Symbol = :linear,celllimit::Int=5)
    if responsetype==:linear
        testfunction = lmtest
        resnum=4
    elseif responsetype==:logistic
        testfunction = logistictest
        resnum=3
    else
        error("Response type was neither linear or logistic")
    end
    vars=DataFrames.allvars(form);
    push!(vars,snpname1,snpname2);
    mdat=df[:,eval(vars)];  #pulls of variables from dataframe
    covmiss=complete_cases(mdat);
    mdat=mdat[covmiss,:];   #this is base dataframe to use

    y=convert(Array{Float64,1},mdat[form.lhs]);  #get phe (dependent variable)
    results=ones(Float64,resnum);   #4 columns for linear, 3 for logistic
    covmat=ModelMatrix(ModelFrame(form,mdat)).m  #this is model matrix with just the intercept and covariates
    const twogentable=[[0,0,0,1,1,1,2,2,2] [0,1,2,0,1,2,0,1,2] [0,0,0,0,0,0,0,0,0] [0,0,0,0,0,0,0,0,0]]
    snp1=convert(DataArray,mdat[snpname1])
    snp1=convert(Array{Int8,1},snp1)
    twogenvec=zeros(Int,size(snp1,1))
    #c1=indicator var for snp1 het,c2=indicator snp2 rare hom, c3=snp2het ind,c4=snp2rarehom ind
    snpind=zeros(Float64,size(snp1,1),4)
    #fill in for first snp, second snp later
    for i in 1:size(snp1,1)
        if snp1[i]==1
            snpind[i,1]=1.0  #make indicator for homozygote
        elseif snp1[i]==2
            snpind[i,2]=1.0  #make indicator for rare homozygote
        end
    end
    snp2=convert(DataArray,mdat[snpname2])
    snp2=convert(Array{Int8,1},snp2)

    cutind,tsnpind,snpgen,dtest=twogencount!(snp1,snp2,snpind,twogenvec,twogentable,celllimit);
    if dtest
        tempy=deepcopy(y);
        tempy=tempy[cutind];
        tempcovmat=deepcopy(covmat)[cutind,:];
        fullmat=hcat(tempcovmat,snpgen);
        reducemat=hcat(tempcovmat,tsnpind);
        results=    try
                        testfunction(tempy,fullmat,reducemat);   
                    catch
                        #zeros(Float64,resnum);
                        rep(-1.0,resnum); #not run because test failed
                    end
    else
        #results[i,:]=zeros(Float64,resnum);
        results=rep(-2.0,resnum);  #not run because of two loc counts (dtest)
    end

    return results
end