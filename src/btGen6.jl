#@time bb=btGen(jj,1);
#@time cc=btGen(jj,1,jj.fam["index"].data,jj.snp["direction"][1]);

#r1=00 hom11, #r2=10 = miss,#r3=01=het, r4=11=hom22

#make -1,0,1,2 for missing, majorhom,het,minorhom
const allelefor=int8([0,-1,1,2]);  #makes 00=hom11=maj=0, 10=miss=-1,01=het=1,11=hom22=minor=2
const allelerev=int8([2,-1,1,0]);  #makes 00=hom11=minor=2, 10=miss=-1,01=het=1,11=hom22=major=0

#col1 = pattern if count at position 1 is <genlimit, col2=regular, col3=if position 2 < genlimit
#col4= pattern if count at position 4 is <genlimit
#if >1 genotype <genlimit then don't do test
const alleleforcut=int8([[-1,-1,1,2] [0,-1,1,2] [0,-1,-1,2] [0,-1,1,-1]]);
const allelerevcut=int8([[-1,-1,1,0] [2,-1,1,0] [2,-1,-1,0] [2,-1,1,-1]]);

#makes 0,1,2,3 for missing, majorhom,het,minorhom
#const allelefor=int8([0,-1,1,2]);  #makes 00=hom11=maj=1, 10=miss=0,01=het=2,11=hom22=minor=3
#const allelerev=int8([2,-1,1,0]);  #makes 00=hom11=minor=3, 10=miss=0,01=het=2,11=hom22=major=1


# ival returns the 0,1,2,3 snp value at offset i from the byte pointer
#ival(p::Ptr{Uint8},i::Integer) = 0<=i ? (unsafe_load(p,i>>2+1)>>(2i&0x06))&0x03 : throw(DomainError)

#ival(p::Ptr{Uint8},i::Integer) = i>0 ? (unsafe_load(p,(i-1)>>2+1)>>(2(i-1)&0x06))&0x03 : throw(DomainError)

#bpt = convert(Ptr{Uint8},b)
#    for j in 1:n
#        ptj = bpt + (snpnum-1)*nind
#        ## the ival function is inlined explicitly
#        for i in 0:(ns-1) counts[1+(unsafe_load(ptj,i>>2+1)>>(2i&0x06))&0x03,j] += 1 end
#    end


#bpt= convert(Ptr{Uint8},b)
#ptj = bpt + (snpnum-1)*nind


# these are for a whole snp without any specification
#function btGen(d::GenDat,snpind::Int,whichallele::Bool=true)
#  bpt=convert(Ptr{Uint8},d.gendat)
#  nind=numind(d)
#  p= bpt + (snpind-1)*nind   #points to beginning of that snp
#  whichallele ? (return btGen(p,nind,allelefor)) : (return btGen(p,nind,allelerev))
#end

#function btGen(p::Ptr{Uint8},nind::Int,gentable::Array{Int8,1}=allelefor)
#  snp=Array(Int8,nind)
#  #this skips the ival function and just runs through consecutive bits
#  for i in 0:(nind-1) snp[i+1]=gentable[1+(unsafe_load(p,i>>2+1)>>(2i&0x06))&0x03] end
#  #for i=1:nind snp[i]=gentable[ival(p,i)] end
#  return snp::Array{Int8,1}
#end


#function btGenFac{T <: Integer}(d::GenDat,snpind::Int,whichallele::T=1)
#  x=btGen(p,nind,whichallele)
#  return PooledDataArray(x,getmissbool(x))
#end

#function btGenDat{T <: Integer}(d::GenDat,snpind::Int,whichallele::T=1)
#  x=btGen(p,nind,whichallele)
#  return DataArray(x,getmissbool(x))
#end

# ival returns the 0,1,2,3 snp value at offset i from the byte pointer
#ival(p::Ptr{Uint8},i::Integer) = 0<=i ? (unsafe_load(p,i>>2+1)>>(2i&0x06))&0x03 : throw(DomainError)

ival(p::Ptr{Uint8},i::Integer) = i>0 ? (unsafe_load(p,(i-1)>>2+1)>>(2(i-1)&0x06))&0x03 : throw(DomainError)


function getmissbool{T <: Integer}(x::Array{T,1})
  mis=falses(length(x))
  for i=1:length(x); if x[i]==-1 mis[i]=true; end; end
  return mis::BitArray{1}
end


##these use vector of integer indexes for individuals
function btGen{T <: Integer}(d::GenDat,snpind::Int,inds::Array{T,1},whichallele::Bool=true)
  bpt=convert(Ptr{Uint8},d.gendat)
  m=size(d.gendat,1)
  p= bpt + (snpind-1)*m  #points to beginning of that snp  
  whichallele ? (return btGen(p,inds,allelefor)) : (return btGen(p,inds,allelerev))
end

function btGen{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},gentable::Array{Int8,1}=allelefor)
  snp=Array(Int8,size(inds,1))
  for i=1:size(inds,1) snp[i]=gentable[ival(p,inds[i])+1] end
  return snp::Array{Int8,1}
end

function btGenFac{T <: Integer}(d::GenDat,snpind::Int,inds::Array{T,1},whichallele::Bool=true)
  x=btGen(d,snpind,inds,whichallele)
  return PooledDataArray(x,getmissbool(x))
end

function btGenDat{T <: Integer}(d::GenDat,snpind::Int,inds::Array{T,1},whichallele::Bool=true)
  x=btGen(d,snpind,inds,whichallele)
  return DataArray(x,getmissbool(x))
end

#to be used for gwLMp
#this one assumes you have pointer to the correct SNP already i.e. convert(Ptr{Uint8},d.gendat) + (snpind-1)-m
function btGenFac{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},whichallele::Bool=true)
  whichallele ? (x=btGen(p,inds,allelefor)) : (x=btGen(p,inds,allelerev))
  return PooledDataArray(x,getmissbool(x))
end

#to be used for gwLMp
#this one assumes you have pointer to the correct SNP already i.e. convert(Ptr{Uint8},d.gendat) + (snpind-1)-m
function btGenDat{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},whichallele::Bool=true)
  whichallele ? (x=btGen(p,inds,allelefor)) : (x=btGen(p,inds,allelerev))
  return DataArray(x,getmissbool(x))
end

#allows just specification GenDat and snpind
#btGen(d::GenDat,snpind::Int)=btGen(d::GenDat,snpind::Int,d.fam["indx"].data,d.snp["direction"][snpind])
#btGenFac(d::GenDat,snpind::Int)=btGenFac(d::GenDat,snpind::Int,d.fam["indx"].data,d.snp["direction"][snpind])
#btGenDat(d::GenDat,snpind::Int)=btGenDat(d::GenDat,snpind::Int,d.fam["indx"].data,d.snp["direction"][snpind])

btGen(d::GenDat,snprow::Int)=btGen(d,d.snp[:indx][snprow],d.fam[:indx].data,d.snp[:direction][snprow])
btGenFac(d::GenDat,snprow::Int)=btGenFac(d,d.snp[:indx][snprow],d.fam[:indx].data,d.snp[:direction][snprow])
btGenDat(d::GenDat,snprow::Int)=btGenDat(d,d.snp[:indx][snprow],d.fam[:indx].data,d.snp[:direction][snprow])


#this uses a bitarray, trues for individuals to keep
function btGen(d::GenDat,snpind::Int,inds::BitArray{1},whichallele::Bool=true)
  bpt=convert(Ptr{Uint8},d.gendat)
  m=size(d.gendat,1)
  p= bpt + (snpind-1)*nind  #points to beginning of that snp
  whichallele ? (return btGen(p,inds,allelefor)) : (return btGen(p,inds,allelerev)) 
end


function btGen(p::Ptr{Uint8},inds::BitArray{1},gentable::Array{Int8,1}=allelefor)
  snp=Array(Int8,sum(inds))
  i=0
  for j=1:size(inds,1) if inds[j] i+=1; snp[i]=gentable[ival(p,inds[j])+1] end end
  return snp::Array{Int8,1}
end

function btGenFac(d::GenDat,snpind::Int,inds::BitArray{1},whichallele::Bool=true)
  x=btGen(p,inds,whichallele)
  return PooledDataArray(x,getmissbool(x))
end

function btGenDat(d::GenDat,snpind::Int,inds::BitArray{1},whichallele::Bool=true)
  x=btGen(p,inds,whichallele)
  return DataArray(x,getmissbool(x))
end


#use these for ChunkLM
function btGen{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},gentable::Array{Int8,1})
  snp=Array(Int8,size(inds,1))
  for i=1:size(inds,1) snp[i]=gentable[ival(p,inds[i])+1] end
  return snp::Array{Int8,1}
end

function btGenFac{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},gentable::Array{Int8,1})
  x=btGen(p,inds,gentable)
  return PooledDataArray(x,getmissbool(x))
end

function btGenDat{T <: Integer}(p::Ptr{Uint8},inds::Array{T,1},gentable::Array{Int8,1})
  x=btGen(p,inds,gentable)
  return DataArray(x,getmissbool(x))
end


function gencount!(ct::Array{Int,1},p::Ptr{Uint8},faminds::Array{Int,1})
    for i in faminds ct[1+ival(p,i)] += 1 end
end

function getgentable(ct::Array{Int,1},whichallele::Bool,gencut::Int)
    #take genotype counts, returns gentable to use and t/f if run test if > 1 gen >=gencut
    uselimit=sum(ct)*.8   #if total number left is <80% then don't run test
    usesum=0
    ncut=0
    cti=0
    for i in [1,3,4]
        if ct[i]<gencut 
            cti=i
            ncut+=1
        else
            usesum+=ct[i]
        end
    end
    if usesum<uselimit return [0,0,0,0],false end

    if whichallele
        if ncut==0 
            return alleleforcut[:,2],true
        elseif ncut>1 
            return [0,0,0,0],false
        else 
            return alleleforcut[:,cti],true
        end
    else
        if ncut==0 
            return allelerevcut[:,2],true
        elseif ncut>1 
            return [0,0,0,0],false
        else 
            return allelerevcut[:,cti],true
        end
    end
end

function getsnp(snpname::ASCIIString,gen::GenDat;gencut::Int=0,asfactor::Bool=false)
  if asfactor     #if asfactor then return snp as PooledDataArray otherwise sent out as continuous DataArray
        togen=btGenFac
    else
        togen=btGenDat
    end

  tp=findin(gen.snp["snpID"],[snpname])
  if length(tp)==1
    tp=tp[1]
  elseif length(tp)==0
    print("SNP $(snpname) was not found.\n")
  else
    tp=tp[1]
    print("SNP $(snpname) is duplicated in the snpinfo list, the first was chosen.\n")
  end

  tp=findin(gen.snp[:snpID],[snpname])[1]
  sid=gen.snp[:indx][tp]
  bpt=convert(Ptr{Uint8},gen.gendat)
    m=size(gen.gendat,1)
    p=bpt + (gen.snp[:indx][tp]-1)*m
    ct=[0,0,0,0];
    gencount!(ct,p,gen.fam[:indx].data);
    usetable, dotest =getgentable(ct,gen.snp[:direction][tp],gencut);
    return togen(p,gen.fam[:indx].data,usetable)
end

function getsnpcounts(snpname::ASCIIString,gen::GenDat)
  tp=findin(gen.snp[:snpID],[snpname])
  if length(tp)==1
    tp=tp[1]
  elseif length(tp)==0
    print("SNP $(snpname) was not found.\n")
  else
    tp=tp[1]
    print("SNP $(snpname) is duplicated in the snpinfo list, the first was chosen.\n")
  end

  tp=findin(gen.snp[:snpID],[snpname])[1]
  sid=gen.snp[:indx][tp]
  bpt=convert(Ptr{Uint8},gen.gendat)
    m=size(gen.gendat,1)
    p=bpt + (gen.snp[:indx][tp]-1)*m
    ct=[0,0,0,0];
    gencount!(ct,p,gen.fam[:indx].data);
    #returns vector of counts for (major hom, het, minor hom, missing)
    gen.snp[:direction][tp] ? (return ct=ct[[1,3,4,2]]) : (return ct=ct[[4,3,1,2]])
end

#to get snp from dGenDat
function getsnp(snpname::ASCIIString,gen::dGenDat;gencut::Int=0,asfactor::Bool=false)
  np=length(gen.refs)
  for i=1:np
    tp = fetch(@spawnat gen.refs[i].where findin(fetch(gen.refs[i]).snp[:snpID],[snpname]))
    if length(tp)!=0 return fetch(@spawnat gen.refs[i].where getsnp(tp[1],fetch(gen.refs[i]),gencut=gencut,asfactor=asfactor)) end
  end
  #error("SNP not found")    
end

#used with dGendat snpind means row number on gen.snp that snp resides
function getsnp(snpind::Int,gen::GenDat;gencut::Int=0,asfactor::Bool=false)
  if asfactor     #if asfactor then return snp as PooledDataArray otherwise sent out as continuous DataArray
        togen=btGenFac
    else
        togen=btGenDat
    end
  sid=gen.snp[:indx][snpind]
  bpt=convert(Ptr{Uint8},gen.gendat)
  m=size(gen.gendat,1)
  p=bpt + (gen.snp[:indx][snpind]-1)*m
  ct=[0,0,0,0];
    gencount!(ct,p,gen.fam[:indx].data);
    usetable, dotest =getgentable(ct,gen.snp[:direction][snpind],gencut);
    return togen(p,gen.fam[:indx].data,usetable)
end

#for dGenDat
function getsnpcounts(snpname::ASCIIString,gen::dGenDat)
  np=length(gen.refs)
  for i=1:np
    tp = fetch(@spawnat gen.refs[i].where findin(fetch(gen.refs[i]).snp[:snpID],[snpname]))
    if length(tp)!=0 return fetch(@spawnat gen.refs[i] getsnpcounts(tp,fetch(gen.refs[i]))) end
  end
end

#used with dGendat snpind means row number on gen.snp that snp resides
function getsnpcounts(snpind::Int,gen::GenDat)
  sid=gen.snp[:indx][snpind]
    bpt=convert(Ptr{Uint8},gen.gendat)
    m=size(gen.gendat,1)
    p=bpt + (gen.snp[:indx][snpind]-1)*m
    ct=[0,0,0,0];
    gencount!(ct,p,gen.fam[:indx].data);
    #returns vector of counts for (major hom, het, minor hom, missing)
    gen.snp[:direction][snpind] ? (return ct=ct[[1,3,4,2]]) : (return ct=ct[[4,3,1,2]])
end

#retrieves given snp and places it in the fam dataframe for each processor which can be used
#as a "phenotype" for GxG
#function putsnp(snpname::ASCIIString,gen::dGenDat;gencut::Int=0,asfactor::Bool=false)
#  snp=getsnp(snpname,gen,gencut=gencut,asfactor=asfactor)
#  for i=1:length(gen.refs) @spawnat gen.refs[i].where fetch(gen.refs[i]).fam[symbol(snpname)]=snp end
#end

#retrieves given snp and places it in the fam dataframe for each processor which can be used
#as a "phenotype" for GxG
#if snp name is invalid name for a dataframe ie. anything with a "-" then you can rename it
function putsnp(snpname::ASCIIString,gen::dGenDat;usename=snpname,gencut::Int=0,asfactor::Bool=false)
  snp=getsnp(snpname,gen,gencut=gencut,asfactor=asfactor)
  for i=1:length(gen.refs) @spawnat gen.refs[i].where fetch(gen.refs[i]).fam[symbol(usename)]=snp end
end

function putphe(phename::Symbol,phe::AbstractDataArray{1},gen::dGenDat)
  for i=1:length(gen.refs) @spawnat gen.refs[i].where fetch(gen.refs[i]).fam[phename]=phe end
end




