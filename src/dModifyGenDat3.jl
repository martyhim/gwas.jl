

#not the most efficient way, it does it on all process
#instead of doing it on one then sending it to the others but easy to programp
function addphe!(file::String,g::dGenDat)
  for i in g.refs @spawnat i.where addphe!(file,fetch(i)) end
  g.phefile=file
end

function addphe!(file::String, gen::GenDat)
  #any phefile column with name (fID,pID,mID,sex,or phe) will be renamed with a _1 next to it, ID in the phefile
  #will be matched to ID in the fam file
  #imports phenotypes (with header)  ## SHOULD have one column labelled ID corresponding to ID in faminfo
  #merges into fam, only keeps only individiuals with phe
  #we do not modify gen.gendat because it points to full bed uint8 matrix
  #we just use gen.fam["indexfam"] as index to pull out gen.gendat with ival function
  nind=numind(gen)
  phedat=readtable(file,separator='\t')  #check with development in DataFrames, may change
  phedat[:pheindx]=1:(size(phedat,1))
  keepfam=falses(numind(gen))
  keepphe=falses(size(phedat,1))
  gen.fam=join(gen.fam,phedat, on = :ID,kind =:inner)
  gen.fam=gen.fam[sortperm(gen.fam[:indx]),:]
  keepphe[phedat[:pheindx]]=true   #T/F if keep phe in merged dataset
  keepfam[gen.fam[:indx]]=true  #T/F if keep individ in merged dataset
  gen.phefile=file
  print("$(numind(gen)) out of $nind ID's are left.\n")
end

#*****still need to write dGenDat version for pick snps


#uses findin function from base
function pickSNPs!(snplist::Array{ASCIIString,1},gen::GenDat)
  nsnp=length(snplist)
  keep=zeros(Int,nsnp)
  for i = 1:nsnp
    tp=findin(gen.snp[:snpID],[snplist[i]])
    if length(tp)==1
      keep[i]=tp[1]
    elseif length(tp)==0
      print("SNP $(snplist[i]) was not found.\n")
    else
      keep[i]=tp[1]
      print("SNP $(snplist[i]) is duplicated in the snpinfo list, the first was chosen.\n")
    end
  end
  if sum(keep.==0)==nsnp
    error("None of the SNPs were found")
  else
    keep=keep[keep.!=0]
    gen.snp=gen.snp[keep,:]
    print("$(length(keep)) of $nsnp SNPs were found.\n")
  end
end

#*****still need to write dGenDat version for pick chroms

function pickChroms!{T <: Integer}(chromlist::Array{T,1},gen::GenDat)
  nsnp=numsnp(gen)
  keep=findin(gen.snp[:chr],chromlist)
  if length(keep)==0
    print("No SNPs were found in the chromosomes from the list.\n")
    print("No Changes were made.\n")
  else
  gen.snp=gen.snp[keep,:]
  print("$(numsnp(gen)) out of $nsnp SNPs were found.\n")
  end
end

#*****still need to write dGenDat version for pickLocRange!
function pickLocRange!{T <: Integer}(chrom::T,lower,upper,gen::GenDat)
  nsnp=numsnp(gen)
  keep=findin(gen.snp[:chr],chrom)
  gen.snp=gen.snp[keep,:]
  gen.snp=gen.snp[((gen.snp[:loc].<=upper) & (gen.snp[:loc].>=lower)),:]
  print("$(numsnp(gen)) out of $nsnp SNPs were kept.\n")
end


MAFthreshhold!{T <: FloatingPoint}(minMAF::T,g::dGenDat) = for i in g.refs @spawnat i.where MAFthreshhold!(minMAF,fetch(i)) end

function MAFthreshhold!{T <: FloatingPoint}(minMAF::T,gen::GenDat)
  #this is the MAF from the original data set before the phenotype data was added
  nsnp=numsnp(gen)
  gen.snp=gen.snp[gen.snp[:MAF].>=minMAF,:]
  print("$(numsnp(gen)) out of $nsnp SNPs were kept.\n")
end


genCountthreshhold!{T <: Integer}(minGenCount::T,g::dGenDat) = for i in g.refs @spawnat i.where genCountthreshhold!(minGenCount,fetch(i)) end

function genCountthreshhold!{T <: Integer}(minGenCount::T,gen::GenDat)
  #actually just puts count threshhold on the minor allele homozygote
  nsnp=numsnp(gen)
  gen.snp=gen.snp[gen.snp[:MinorHom].>=minGenCount,:]
  print("$(numsnp(gen)) out of $nsnp SNPs were kept.\n")
end


#not the most efficient way, it does it on all process
#instead of doing it on one then sending it to the others but easy to program
pickIDs!(IDlist::Array{ASCIIString,1},g::dGenDat) = for i in g.refs @spawnat i.where pickIDs!(IDlist,fetch(i)) end


function pickIDs!(IDlist::Array{ASCIIString,1},gen::GenDat)
  keepfam=falses(numind(gen))
  nind=length(IDlist)
  keep=zeros(Int,nsnp)
  for i = 1:nind
    tp=findin(gen.fam[:ID],[IDlist[i]])
    if length(tp)==1
      keep[i]=tp[1]
    elseif length(tp)==0
      print("ID $(IDlist[i]) was not found.\n")
    else
      keep[i]=tp[1]
      print("ID $(IDlist[i]) is duplicated in the geninfo list, the first was chosen.\n")
    end
  end
  if sum(keep.==0)==nind
    error("None of the IDs were found")
  else
    keep=keep[keep.!=0]
    keepfam[keep]=true
    gen.fam=gen.fam[keepfam,:]
    print("$(numind(gen)) out of $nind ID's were found.\n")
  end
end

#get -log10 of pval for HW
function HW!(gen::GenDat)
  HWlog10pval=zeros(FloatingPoint,numsnp(gen))
  for i=1:numsnp(gen)
      HWlog10pval[i]=HWtest(i,gen.snp)
  end
  gen.snp[:HWlog10pval]=HWlog10pval
end


function HWtest(i::Int, snp::DataFrame)
  q=snp[:MAF][i]
  o=[snp[:hom11][i],snp[:het][i],snp[:hom22][i]]
  n=sum(o)
  e=[(1-q)*(1-q)*n, 2*(1-q)*p*n, q*q ]
  chi=0.0
  for j=1:3 (chi += ((o[j]-e[j])^2)/e[j]) end
  return -log10(ccdf(Chisq(1),chi))
end

HWthreshhold!{T <: FloatingPoint}(maxLog10HW::T,g::dGenDat) = for i in g.refs @spawnat i.where HWthreshhold!(maxLog10HW,fetch(i)) end

function HWthreshhold!{T <: FloatingPoint}(maxLog10HW::T,gen::GenDat)
  #need to first check if HWlog10pval has been calculated yet
  #if not then it gets calculation first
  if sum(names(gen.snp).==:HWlog10pval) == 0 HW!(gen) end
  nsnp=numsnp(gen)
  gen.snp=gen.snp[gen.snp[:HWlog10pval].<=maxLog10HW,:]
  print("$(numsnp(gen)) out of $nsnp SNPs were kept.\n")
end

updateHW!(gen::dGenDat) = for i in g.refs fetch(@spawnat i.where updateHW!(fetch(i))) end

#if change ID's like updatephe! or pickIDs! then might want to recalculate HW
updateHW!(gen::GenDat) = for i=1:numsnp(gen) gen.snp[:HWlog10pval][i]=HWtest(i,gen.snp) end

missingthreshhold!{T <: FloatingPoint}(minmissingpercent::T,g::dGenDat) = for i in g.refs @spawnat i.where missingthreshhold!(minmissingpercent,fetch(i)) end

function missingthreshhold!{T <: FloatingPoint}(minmissingpercent::T,gen::GenDat)
  nsnp=numsnp(gen)
  nind=numind(gen)
  percentmissing=zeros(FloatingPoint,numsnp(gen))
  for i=1:numsnp(gen)
      percentmissing[i]= gen.snp[:missing][i]/nind
  end
  gen.snp[:percentmissing]=percentmissing
  gen.snp=gen.snp[gen.snp[:percentmissing].<=minmissingpercent,:]
  print("$(numsnp(gen)) out of $nsnp SNPs were kept.\n")
end



