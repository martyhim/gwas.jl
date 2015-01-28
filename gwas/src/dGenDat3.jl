type dGenDat
    snpidx::Array{Int,2}  #c1=processor #, c2=first snp, c3=last snp
    refs::Array{RemoteRef,1}  #array of RemoteRefs for GenDat types on each processor each holding subset of snps
    originaldims::Array{Int,1}  #the original number of individuals and snps, used to delimin snp length in bedfile
    genfile::ASCIIString  #original base name of plink file without the (.bed,.bim, or .fam)
    phefile::ASCIIString
end

function numsnp(g::dGenDat)
    n=0
    for i in g.refs n+= fetch(@spawnat i.where numsnp(fetch(i))) end
    return n
end

numind(g::dGenDat)=fetch(@spawnat g.snpidx[1,1] numind(fetch(g.refs[1])))  #all refs should have same fam dataframe
famname(g::dGenDat)=fetch(@spawnat g.snpidx[1,1] famname(fetch(g.refs[1])))  #all refs should have smae fam dataframe

# ival returns the 0,1,2,3 snp value at offset i from the byte pointer
#ival(p::Ptr{Uint8},i::Integer) = 0<=i ? (unsafe_load(p,i>>2+1)>>(2i&0x06))&0x03 : throw(DomainError)

ival(p::Ptr{Uint8},i::Integer) = i>0 ? (unsafe_load(p,(i-1)>>2+1)>>(2(i-1)&0x06))&0x03 : throw(DomainError)

type GenDat
  snp::DataFrame  #snps with alleles, direction,MAF,genotype counts, missing
  fam::DataFrame  #famfile from plink, adds index at end, phenotype data get put here if incorporated
  gendat::Matrix{Uint8} #points to original plink bedfile matrix
  originaldims::Array{Int,1}  #the original number of individuals and snps, used to delimin snp length in bedfile
  genfile::ASCIIString  #original base name of plink file without the (.bed,.bim, or .fam)
  phefile::ASCIIString  #original file to import phenotype data with
end

numind(gd::GenDat)=size(gd.fam,1)
numsnp(gd::GenDat)=size(gd.snp,1)
famname(gd::GenDat)=names(gd.fam)

function dGenDat(file::ASCIIString)
    bimfile="$(file).bim";
    famfile="$(file).fam";
    snpfile="$(file).bed";
    nsnp = countlines("$(file).bim"); 
    nsubj = countlines("$(file).fam");

    np=nprocs()
    if np==1 error("Only 1 processor, use the gwLM function") end
    npuse=nworkers()
    parts=int(floor(nsnp/npuse))
    idx=zeros(Int,npuse,3)    #rows are c1=beginning, c2=end of indexes to send out to each processor
    idx[:,1]=[2:np]
    idx[:,2]=[1:parts:nsnp][1:npuse]
    idx[1:(npuse-1),3]=([1:parts:nsnp].-1)[2:npuse]
    idx[npuse,3]=nsnp

    #distribute snps in different GenDat's across procs
    refs=Array(RemoteRef,npuse)
    for i=1:npuse refs[i]= @spawnat idx[i,1] GenDatpart(file,idx[i,2],idx[i,3]) end

    print("There are $nsnp SNPs and $nsubj individuals")
    return dGenDat(idx,refs,[nsubj,nsnp],file,"empty")
end



function GenDatpart(file::ASCIIString,firstsnp::Int,lastsnp::Int) # file = base file name without extension
  	bimfile="$(file).bim";
  	famfile="$(file).fam";
    snpfile="$(file).bed";
    nrowskeep=lastsnp-firstsnp+1
    #snpdat=readdlm(bimfile,'\t');
    snames=[:chr,:snpID,:gendist,:loc,:allele1,:allele2]
    snp=readtable(bimfile,separator='\t',header=false,names=snames,skipstart=firstsnp-1,nrows=nrowskeep)
    if size(snp,2) != 6
        error("file $bimfile should have 6 tab-delimited columns")
    end
    ##dfctype!(Int8,["chr","gendist"],snp)  #changes these columns to Int8
       #fam=readdlm(famfile,' ')
    dfctype!(Int8,[:chr],snp)  #changes these columns to Int8
    fnames=[:fID,:ID,:pID,:mID,:sex,:phe]
    fam=readtable(famfile,separator=' ',header=false,names=fnames)
    if size(fam,2) != 6
        error("file $famfile should have 6 space-delimited columns")
    end
    dfctype!(Int8,[:sex,:phe],fam)  #changes these columns to Int8
    nsnp = size(snp,1)
    nsubj = size(fam,1)
    gendat=readbedmempart(file,nsubj,nsnp,firstsnp,lastsnp)

    #c1=00 hom11,c2=10 miss, c3=01 het,c4=11 hom22
    counts=bedfreq(gendat,nsubj)'
    maf=zeros(Float64,nsnp)
    dir=falses(nsnp)   #if true the allele 1 is major allele 2 is minor, opposite if false
    minorhom=counts[:,1]
    for i=1:nsnp
      getMAF!(i,counts,maf,dir,minorhom)
    	#maf[i]=float64((0.5*counts[i,3]+counts[i,4])/(counts[i,2]+counts[i,3]+counts[i,4]))  #initially maf for allele 2 freq
      #if maf[i]>0.5 (maf[i]=1-maf[i]); dir[i]=false; minorhom=counts[i,1] end
    	#if maf[i]>0.5 switchcounts!(i,maf,dir,counts,snp) end
    end
    #snp names
    #"chr","snpID","gendist","loc","major","minor","direction","MAF","MajHom","Het","MinorHom","missing"
	snp[:direction]=dir
	snp[:MAF]=maf
    snp[:minorhom]=minorhom
	snp[:hom11]= counts[:,1]
	snp[:het]= counts[:,3]
	snp[:hom22]= counts[:,4]
	snp[:missing]= counts[:,2]
    snp[:indx]=1:(size(snp,1))  #add internal rowindex corresponds with columns of gendat matrix
	fam[:indx]=1:(size(fam,1))  #add internal rowindex corresponds with ival index for gendat matrix
	print("There are $(size(snp,1)) SNPs and $(size(fam,1)) individuals")
	GenDat(snp,fam,gendat,[nsubj,nsnp],file,"empty")
end

function GenDat(file::ASCIIString) # file = base file name without extension
    bimfile="$(file).bim";
    famfile="$(file).fam";
    #snpdat=readdlm(bimfile,'\t');
    snames=[:chr,:snpID,:gendist,:loc,:allele1,:allele2]
    snp=readtable(bimfile,separator='\t',header=false,names=snames)
    if size(snp,2) != 6
        error("file $bimfile should have 6 tab-delimited columns")
    end
    #dfctype!(Int8,["chr","gendist"],snp)  #changes these columns to Int8
       #fam=readdlm(famfile,' ')
    dfctype!(Int8,[:chr],snp)  #changes these columns to Int8
    fnames=[:fID,:ID,:pID,:mID,:sex,:phe]
    fam=readtable(famfile,separator=' ',header=false,colnames=fnames)
    if size(fam,2) != 6
        error("file $famfile should have 6 space-delimited columns")
    end
    dfctype!(Int8,[:sex,:phe],fam)  #changes these columns to Int8
    nsnp = size(snp,1)
    nsubj = size(fam,1)
    gendat=readbedmem(file,nsubj,nsnp)

    #c1=00 hom11,c2=10 miss, c3=01 het,c4=11 hom22
    counts=bedfreq(gendat,nsubj)'
    maf=zeros(Float64,nsnp)
    dir=falses(nsnp)   #if true the allele 1 is major allele 2 is minor, opposite if false
    minorhom=counts[:,1]
    for i=1:nsnp
      getMAF!(i,counts,maf,dir,minorhom)
        #maf[i]=float64((0.5*counts[i,3]+counts[i,4])/(counts[i,2]+counts[i,3]+counts[i,4]))  #initially maf for allele 2 freq
      #if maf[i]>0.5 (maf[i]=1-maf[i]); dir[i]=false; minorhom=counts[i,1] end
        #if maf[i]>0.5 switchcounts!(i,maf,dir,counts,snp) end
    end
    #snp names
    #"chr","snpID","gendist","loc","major","minor","direction","MAF","MajHom","Het","MinorHom","missing"
    snp[:direction]=dir
    snp[:MAF]=maf
    snp[:minorhom]=minorhom
    snp[:hom11]= counts[:,1]
    snp[:het]= counts[:,3]
    snp[:hom22]= counts[:,4]
    snp[:missing]= counts[:,2]
    snp[:indx]=1:(size(snp,1))  #add internal rowindex corresponds with columns of gendat matrix
    fam[:indx]=1:(size(fam,1))  #add internal rowindex corresponds with ival index for gendat matrix
    print("There are $(size(snp,1)) SNPs and $(size(fam,1)) individuals")
    GenDat(snp,fam,gendat,[nsubj,nsnp],file,"empty")
end


function readbedmempart(file::ASCIIString,nsubj::Int,nsnp::Int,firstsnp::Int,lastsnp::Int)
    #this reads into memory (not memory mapped) the set of snps it is asked to
    # file = basename of .bed, .bim and .fam files
    # count the lines in the .bim and .fam files to get the dimensions
    #nsnp = countlines("$(file).bim"); nsubj = countlines("$(file).fam")
    bednm = "$(file).bed"
    nrowskeep=lastsnp-firstsnp+1
    s = open(bednm)
    read(s,Uint8) == 0x6c && read(s,Uint8) == 0x1b || error("wrong magic number in file $bednm")
    read(s,Uint8) == 1 || error(".bed file, $bednm, is not in correct orientation")
    m = div(nsubj+3,4) # of bytes per col., nsubj rounded up to next multiple of 4
    skip(s,m*(firstsnp-1))
    #readbytes!(stream, b::Vector{Uint8}, nb=length(b))
    #bb=Array(Uint8,readbytes!(s,bb,m*nrowskeep))
    bb=read(s,Uint8,m,nrowskeep)  #gets Uint8 matrix with just snps for that part
    close(s)
    return bb
end

function readbedmem(file::ASCIIString,nsubj::Int,nsnp::Int)     # file = basename of .bed, .bim and .fam files
    # count the lines in the .bim and .fam files to get the dimensions
    #nsnp = countlines("$(file).bim"); nsubj = countlines("$(file).fam")
    bednm = "$(file).bed"
    s = open(bednm)
    read(s,Uint8) == 0x6c && read(s,Uint8) == 0x1b || error("wrong magic number in file $bednm")
    read(s,Uint8) == 1 || error(".bed file, $bednm, is not in correct orientation")
    m = div(nsubj+3,4) # of bytes per col., nsubj rounded up to next multiple of 4
    bb=read(s,Uint8,m,nsnp)  #gets Uint8 matrix with just snps for that part
    close(s)
    bb
    #GenData2(bb,nsubj,bedfreq(bb,nsubj)')
end


function getMAF!(i::Int,counts::Array{Int,2},maf::Array{Float64,1},dir::BitArray{1},minorhom::Array{Int,1})
  maf[i]=float64((0.5*counts[i,3]+counts[i,1])/(counts[i,1]+counts[i,3]+counts[i,4]))
  if maf[i]>0.5 (maf[i]=1-maf[i]); dir[i]=true; minorhom[i]= counts[i,4] end 
end

function bedfreq(b::Matrix{Uint8},ns::Integer)
    m,n = size(b)
    div(ns+3,4) == m || error("ns = $ns should be in [$(4m-3),$(4m)] for size(b) = $(size(b))")
    cc = zeros(Int, 4, n)
    bpt = convert(Ptr{Uint8},b)
    for j in 0:(n-1)
        pj = bpt + j*m
        ## the ival function is inlined explicitly
        for i in 1:ns cc[1+ival(pj,i),1+j] += 1 end
        #for i in 0:(ns-1) cc[1+j,1+(unsafe_load(pj,i>>2+1)>>(2i&0x06))&0x03] += 1 end
    end
    cc
end

#change type of dataarray in a df by colname of colnumber
dfctype!(dtype::DataType,colname::Symbol,df::AbstractDataFrame)= (df[colname]=convert(DataArray{dtype},df[colname]))
dfctype!(dtype::DataType,colnum::Int,df::AbstractDataFrame)= (df[:,colnum]=convert(DataArray{dtype},df[:,colnum]))

#change datatype for a vector of colnames or numbers in a dataframe
dfctype!(dtype::DataType,colnames::Array{Symbol,1},df::AbstractDataFrame)= (for i in colnames dfctype!(dtype,i,df) end)
dfctype!(dtype::DataType,colnums::Array{Int,1},df::AbstractDataFrame)= (for i in colnums dfctype!(dtype,i,df) end)


updatecounts!(g::dGenDat) = for i in g.refs @spawnat i.where updatecounts!(fetch(i)) end


function updatecounts!(g::GenDat)
  #after any snp or individual changes this recalculated maf,gencounts,direction
  #this can be done after any modifying changes in ModifyGenDat.jl
  nsnp=numsnp(g)
  nind=numind(g)
  m=size(g.gendat,1)
  bpt=convert(Ptr{Uint8},g.gendat)
  counts=zeros(Int,4,nsnp)
  minorhom=zeros(Int,nsnp)
  updatecounts!(nsnp,nind,g.snp[:indx].data,g.fam[:indx].data,counts,bpt,m)
  counts=counts'
  #for i=1:nsnp
  #  pj = bpt + (g.snp["index"][i]-1)*m
  #  for j=1:nind counts[1+ival(pj,g.fam["index"][j]),i] += 1 end
  #  getMAF!(i,g,counts,minorhom)
    #g.snp["MAF"][i]= float64((0.5*counts[i,3]+counts[i,4])/(counts[i,2]+counts[i,3]+counts[i,4]))
    #if g.snp["MAF"][i]>0.5 switchcounts!(i,count,g.snp) end
  #end
  minorhom=counts[:,1]
  dir=falses(nsnp)
  maf=zeros(Float64,nsnp)
  for i=1:nsnp
    getMAF!(i,counts,maf,dir,minorhom)
  end
  g.snp[:MAF]=maf
  g.snp[:minorhom]=minorhom
  g.snp[:direction]=dir
  g.snp[:hom11]= counts[:,1]
  g.snp[:het]= counts[:,3]
  g.snp[:hom22]= counts[:,4]
  g.snp[:missing]= counts[:,2]
end

#used in updatecounts!
function updatecounts!(nsnp::Int,nind::Int,snpind::Array{Int,1},famind::Array{Int,1},counts::Array{Int,2},bpt::Ptr{Uint8},m::Int)
  for i=1:nsnp
    pj = bpt + (snpind[i]-1)*m
    for j=1:nind counts[1+ival(pj,famind[j]),i] += 1 end
    
    #g.snp["MAF"][i]= float64((0.5*counts[i,3]+counts[i,1])/(counts[i,1]+counts[i,3]+counts[i,4]))
    #if g.snp["MAF"][i]>0.5 switchcounts!(i,count,g.snp) end
  end
end
