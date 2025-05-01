### This script demonstrates the logic progression in building the simulation, first simulating stem cell niche dynamics, followed by crypt cell dynamics, followed by putting of the two parts together

### A: Simulating the Stem Cell Niche ###
rm(list=ls())
# Settable variables
sim=500    # number of simulated generations 
mutrate=1  # number of mutation(s) per genome per generation
numSC=4    # number of starting stem cell(s) in niche

# Initialize a data frame to store stem cell data. One row per generation. Columns for each stem cell X 2 to allow recording of stem cell division. Each data frame cell value records the cell lineage history and accumulated mutations.  
tempdf=data.frame(SC_gen=paste("gen",(0:sim),sep=""))
tempdf[tempdf$SC_gen=="gen0",paste("StemCell",1:(2*numSC),sep="")]=c(paste("sc",1:numSC,"//",sep=""),rep(NA,numSC))
tempdf[tempdf$SC_gen=="gen0","keep"]=paste(1:numSC,collapse=",")
head(tempdf)

# start running through generations to populate tempdf
mutcount=0
gencount=0
for (i in 1:sim) {
  #i=1 i=2 2=3
  prevgen=paste("gen",gencount,sep="")
  gencount=gencount+1
  currgen=paste("gen",gencount,sep="")
  prevkeep=as.numeric(unlist(strsplit(tempdf[tempdf$SC_gen==prevgen,"keep"],",")))
  
  #cell division with mutation acquired   
  filled=0
  for (j in prevkeep) {
    #j=1 j=2 j=3
    
    # parent cell
    parent=tempdf[tempdf$SC_gen==prevgen,paste("StemCell",j,sep="")]
    parentid=sub("//\\S*","",parent)
    parentmut=sub("^[^/]+//","",parent)
    
    # daughter cell 1
    d1id=paste(parentid,".1",sep="") 
    d1mut=parentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d1mut=paste(d1mut,":m",mutcount,sep="")
    }
    d1mut=sub("^:","",d1mut)
    d1=paste(d1id,d1mut,sep="//")
    d1
    filled=filled+1
    tempdf[tempdf$SC_gen==currgen,paste("StemCell",filled,sep="")]=d1
    
    # daughter cell 2
    d2id=paste(parentid,".2",sep="")
    d2mut=parentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d2mut=paste(d2mut,":m",mutcount,sep="")
    }
    d2mut=sub("^:","",d2mut)
    d2=paste(d2id,d2mut,sep="//")
    d2
    filled=filled+1
    tempdf[tempdf$SC_gen==currgen,paste("StemCell",filled,sep="")]=d2
  }
  
  # randomly select the stem cells to keep in niche for the next generation
  tempdf[tempdf$SC_gen==currgen,"keep"]=paste(sort(sample(1:(numSC*2), numSC, replace = FALSE, prob = NULL)),collapse=",")
  if (i%%100==0) {
    cat(i," ",sep="")
  }
}
head(tempdf)
tail(tempdf)

# summarize allele frequencies
sumdf=data.frame(row.names=tempdf$SC_gen,SC_gen=tempdf$SC_gen)
for (i in 2:length(row.names(tempdf))) {
  #i=1  i=2  i=3 i=100
  tempdf[i,paste("StemCell",seq(1,(numSC*2),2),sep="")]
  muts=unlist(strsplit(sub("^\\S+//","",tempdf[i,paste("StemCell",seq(1,(numSC*2),2),sep="")]),":"))
  
  sumdf[i,"num_mut"]=length(muts)
  sumdf[i,"num_uniq_mut"]=length(unique(muts))
  sumdf[i,"unfixed_genappear"]=floor(min(as.numeric(sub("^m","",names(table(muts)[table(muts)<numSC])))/(numSC*2)))
  
  for (j in (sort(unique(table(muts)),decreasing=T))) {
    #j=1 j=4
    sumdf[i,paste("num_AF",(j/(numSC*2)),sep="")]=sum(table(muts)==j)
  }
  
  if (i%%100==0) {
    cat(i," ",sep="")
  }
}
head(sumdf,10)
tail(sumdf)

# some data exploration and visualizations
(1:length(row.names(sumdf)))-sumdf$unfixed_genappear
hist((1:length(row.names(sumdf)))-sumdf$unfixed_genappear,breaks=50,main="unfixed_genappear")
summary((1:length(row.names(sumdf)))-sumdf$unfixed_genappear)
table(sumdf$num_AF0.5)
table(sumdf$num_AF0.375)
table(sumdf$num_AF0.25)
table(sumdf$num_AF0.125)
plot(sumdf$num_AF0.5)
plot(sumdf$num_AF0.375)
plot(sumdf$num_AF0.25)
plot(sumdf$num_AF0.125)

### B: Simulating the Crypt ###
rm(list=ls())
# Settable variables
sim=100             # number of simulated generations
mutrate=1           # number of mutation(s) per genome per generation
numSC=1             # number of starting cell(s) in crypt
numcryptcells=2048  # crypt size i.e. number of cells it can hold

# Initialize a data frame to store crypt cell data. One row per generation. 2048 columns represent the number of cells the crypt can hold. Each data frame cell value records the cell lineage history and accumulated mutations.
tempdf=data.frame(Crypt_gen=paste("gen",(0:sim),sep=""))
tempdf[tempdf$Crypt_gen=="gen0",paste("CryptCell",1:2048,sep="")]=c(paste("cc",1,"//",sep=""),rep(NA,numcryptcells-numSC))
head(tempdf)
tempdf[1:5,1:5]

# start running through generations to populate tempdf
mutcount=0
gencount=0
for (i in 1:sim) {
  #i=1 i=2 2=3
  prevgen=paste("gen",gencount,sep="")
  gencount=gencount+1
  currgen=paste("gen",gencount,sep="")

  #cell division with mutation acquired   
  filled=0
  for (j in 1:numcryptcells) {
    #j=1 j=2 j=3
    
    # parent cell
    parent=tempdf[tempdf$Crypt_gen==prevgen,paste("CryptCell",j,sep="")]
    if (is.na(parent)) {
      break   
    }
    if (filled > numcryptcells) {
      break   
    }
    parentid=sub("//\\S*","",parent)
    parentmut=sub("^[^/]+//","",parent)
    
    # daughter cell 1
    d1id=paste(parentid,".1",sep="")
    d1mut=parentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d1mut=paste(d1mut,":m",mutcount,sep="")
    }
    d1mut=sub("^:","",d1mut)
    d1=paste(d1id,d1mut,sep="//")
    d1
    filled=filled+1
    if (filled <= numcryptcells) {
      tempdf[tempdf$Crypt_gen==currgen,paste("CryptCell",filled,sep="")]=d1
    }
    
    # daughter cell 2
    d2id=paste(parentid,".2",sep="")
    d2mut=parentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d2mut=paste(d2mut,":m",mutcount,sep="")
    }
    d2mut=sub("^:","",d2mut)
    d2=paste(d2id,d2mut,sep="//")
    d2
    filled=filled+1
    if (filled <= numcryptcells) {
      tempdf[tempdf$Crypt_gen==currgen,paste("CryptCell",filled,sep="")]=d2
    }
  }
  if (i%%10==0) {
    cat(i," ",sep="")
  }
}
tempdf[1:6,1:8]
tempdf[1:6,2048:2049]
tempdf[51:55,1001:1005]
apply(apply(tempdf,1,is.na),2,sum)

# summarize allele frequencies
sumdf=data.frame(row.names=tempdf$Crypt_gen,Crypt_gen=tempdf$Crypt_gen)
for (i in 2:length(row.names(tempdf))) {
  #i=1  i=2  i=3 i=10 i=11 i=12 i=13 i=20 i=40 i=60 i=100
  
  muts=unlist(strsplit(sub("^\\S+//","",tempdf[i,paste("CryptCell",1:numcryptcells,sep="")]),":"))
  muts=muts[!is.na(muts)]
  
  sumdf[i,"num_mut"]=length(muts)
  sumdf[i,"num_uniq_mut"]=length(unique(muts))
  for (j in (sort(unique(table(muts)),decreasing=T))) {
    #j=1 j=4
    sumdf[i,paste("num_AF",(j/(numcryptcells*2)),sep="")]=sum(table(muts)==j)
  }
  
  if (i%%10==0) {
    cat(i," ",sep="")
  }
}
head(sumdf)
tail(sumdf)


### C: Simulating the Stem Cell Niche and Crypt together         ###
### Assumptions:                                                 ###
###  - Generation time of niche cells double of crypt cells      ###
###  - After each generation of stem cell division, half is      ###
###    randomly chosen to be retained in niche and half will be  ###
###    pushed into crypt (i.e. symmetric).                       ###
rm(list=ls())
# Settable variables
sim=500             # number of simulated generations
mutrate=1           # number of mutation(s) per genome per generation
numstemcells=4      # number of starting stem cell(s) in niche
numcryptcells=2048  # crypt size i.e. number of cells it can hold

# Initialize a single data frame to store stem cell niche and crypt data. 
# One row per "generation". But with assumption that stem cell generation time is double that of crypt cell generation time, stem cell niche columns will be simulated/populated only at every alternate row, while crypt columns will be simulated/populated in every row.     
tempdf=data.frame()
tempdf["sim0","SC_gen"]="scgen0"
tempdf["sim0",paste("StemCell",1:(numstemcells*2),sep="")]=c(paste("sc",1:numstemcells,"//",sep=""),rep(NA,numstemcells))
tempdf["sim0","keep"]=paste(1:numstemcells,collapse=",")
tempdf["sim0","push"]=NA
tempdf["sim0","CC_gen"]="ccgen0"
tempdf["sim0",paste("CryptCell",1:numcryptcells,sep="")]=rep(NA,numcryptcells)
head(tempdf)
tempdf[,1:20]

# start running through generations to populate tempdf
mutcount=0
simgencount=0
scgencount=0
ccgencount=0
for (i in 1:sim) {
  #i=1 i=2 i=3
  simgencount=i
  currrowname=paste("sim",i,sep="")
  tempdf[currrowname,]=NA

  ### SC Niche : run every alternate row
  if (simgencount%%2==0) {
    scprevgen=paste("scgen",scgencount,sep="")
    scgencount=scgencount+1
    sccurrgen=paste("scgen",scgencount,sep="")
    tempdf[currrowname,"SC_gen"]=sccurrgen
    scprevkeep=as.numeric(unlist(strsplit(tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==scprevgen),"keep"],",")))
    
    scfilled=0
    for (j in scprevkeep) {
      #j=1 j=2 j=3
      
      scparent=tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==scprevgen),paste("StemCell",j,sep="")]
      scparentid=sub("//\\S*","",scparent)
      scparentmut=sub("^[^/]+//","",scparent)
      
      d1id=paste(scparentid,".1",sep="")
      d1mut=scparentmut
      for (k in 1:mutrate) {
        mutcount=mutcount+1;
        d1mut=paste(d1mut,":m",mutcount,sep="")
      }
      d1mut=sub("^:","",d1mut)
      d1=paste(d1id,d1mut,sep="//")
      d1
      scfilled=scfilled+1
      tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),paste("StemCell",scfilled,sep="")]=d1
      
      d2id=paste(scparentid,".2",sep="")
      d2mut=scparentmut
      for (k in 1:mutrate) {
        mutcount=mutcount+1;
        d2mut=paste(d2mut,":m",mutcount,sep="")
      }
      d2mut=sub("^:","",d2mut)
      d2=paste(d2id,d2mut,sep="//")
      d2
      scfilled=scfilled+1
      tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),paste("StemCell",scfilled,sep="")]=d2
    }
    select=sample(1:(numstemcells*2), (numstemcells*2), replace = FALSE, prob = NULL)
    tokeep=sort(select[1:numstemcells])
    topush=(select[(numstemcells+1):(numstemcells*2)]) # note: not sorted to maintain randomness in order of cells being pushed into crypt
    tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),"keep"]=paste(tokeep,collapse=",")
    tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),"push"]=paste(topush,collapse=",")
  }
  
  ### Crypt : run every row
  ccprevgen=paste("ccgen",ccgencount,sep="")
  ccgencount=ccgencount+1
  cccurrgen=paste("ccgen",ccgencount,sep="")
  tempdf[currrowname,"CC_gen"]=cccurrgen
  
  ccfilled=0 # to track the progress of the 2048 crypt positions being filled
  
  # start populating row with the stem cells to be pushed into crypt... 
  scpushes=tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),"push"]
  if (!is.na(scpushes)) {
    for (scpush in unlist(strsplit(scpushes,","))) {
      ccfilled=ccfilled+1
      tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=sub("//",".cc//",tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("StemCell",scpush,sep="")])
    }
  } 
  
  # ...followed by crypt cells from previous generation from left to right. this maintains positional crypt cell division in place
  for (j in 1:numcryptcells) {
    #j=1 j=2 j=3
    ccparent=tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==ccprevgen),paste("CryptCell",j,sep="")]
    if (is.na(ccparent)) {
      break
    }
    if (ccfilled > numcryptcells) {
      break
    }
    
    ccparentid=sub("//\\S*","",ccparent)
    ccparentmut=sub("^[^/]+//","",ccparent)
    
    d1id=paste(ccparentid,".1",sep="")
    d1mut=ccparentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d1mut=paste(d1mut,":m",mutcount,sep="")
    }
    d1mut=sub("^:","",d1mut)
    d1=paste(d1id,d1mut,sep="//")
    d1
    ccfilled=ccfilled+1
    if (ccfilled <= numcryptcells) {
      tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=d1
    }
    
    d2id=paste(ccparentid,".2",sep="")
    d2mut=ccparentmut
    for (k in 1:mutrate) {
      mutcount=mutcount+1;
      d2mut=paste(d2mut,":m",mutcount,sep="")
    }
    d2mut=sub("^:","",d2mut)
    d2=paste(d2id,d2mut,sep="//")
    d2
    ccfilled=ccfilled+1
    if (ccfilled <= numcryptcells) {
      tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=d2
    }
  }
  if (i%%100==0) {
    mutcount
    
  }
  if (i%%10==0) {
    cat(i," ",sep="")
  }
  
}
tempdf[1:5,1:16]
tempdf[11:15,2036:2060]
tempdf[101:105,1:16]
#write.table(tempdf,"tempdf.txt",quote=F,sep="\t",row.names=T,col.names=T)
#tempdf<-read.table("tempdf.txt",header=T,sep="\t",stringsAsFactors=F)

# summarize allele frequencies
sumdf=data.frame(row.names=tempdf$CC_gen,Crypt_gen=tempdf$CC_gen)
for (i in 2:length(row.names(tempdf))) {
  #i=2  i=3 i=10 i=11 i=12 i=13 i=14 i=20 i=40 i=60 i=100
  
  muts=unlist(strsplit(sub("^\\S+//","",tempdf[i,paste("CryptCell",1:numcryptcells,sep="")]),":"))
  muts=muts[!is.na(muts)]
  
  sumdf[i,"num_mut"]=length(muts)
  sumdf[i,"num_uniq_mut"]=length(unique(muts))
  for (j in (sort(unique(table(muts)),decreasing=T))) {
    #j=1 j=4
    sumdf[i,paste("num_AF",(j/(numcryptcells*2)),sep="")]=sum(table(muts)==j)
  }
  
  if (i%%10==0) {
    cat(i," ",sep="")
  }
}
for (j in 4:length(names(sumdf))) {
  #j=4  j=30
  sumdf[is.na(sumdf[,j]),j]=0
}
head(sumdf)
tail(sumdf)
#write.table(sumdf,"sumdf.txt",quote=F,sep="\t",row.names=T,col.names=T)
#sumdf<-read.table("sumdf.txt",header=T,sep="\t",stringsAsFactors=F)

# some data exploration and visualizations of the allele frequency summary table
sumdf<-read.table("sumdf.txt",header=T,sep="\t",stringsAsFactors=F)
sumdf[1:5,1:5]
names(sumdf)
sumdf[,4:length(names(sumdf))]
apply(sumdf[,4:length(names(sumdf))],1,sum,na.rm=T)
apply(sumdf[,4:length(names(sumdf))],1,sum,na.rm=T)==sumdf$num_uniq_mut
sum(apply(sumdf[,4:length(names(sumdf))],1,sum,na.rm=T)==sumdf$num_uniq_mut,na.rm=T)

#create data tab to plot histograms
sumdf[1:5,1:5]
library(tidyr)
library(reshape2)
plottab=melt(sumdf[,c(1,4:length(names(sumdf)))],"Crypt_gen")
plottab=plottab[plottab$value>0,]
as_tibble(plottab)
uncount(as_tibble(plottab),value)
plottab=as.data.frame(uncount(as_tibble(plottab),value))
plottab$variable=as.numeric(sub("num_AF","",plottab$variable))
names(plottab)[2]="AF"
plottab$gen2=as.numeric(sub("ccgen","",plottab$Crypt_gen))
plottab$Crypt_gen=sprintf("ccgen%03d", plottab$gen2)
head(plottab)
write.table(plottab,"plottab.txt",quote=F,sep="\t",row.names=T,col.names=T)
plottab<-read.table("plottab.txt",header=T,sep="\t",stringsAsFactors=F)
head(plottab)
summary(plottab$AF)
table(plottab$AF)

plottabsmall=plottab[(plottab$AF >= 0.01)&(plottab$gen2==95),]
hist(plottabsmall$AF,freq=T,breaks=seq(0,0.6,0.01),xlab="Allele Frequency",main="AF>=0.01")
mtext(paste("(n=",length(row.names(plottabsmall)),")",sep=""),3,line=-1)

library(ggplot2)
plottab<-read.table("plottab.txt",header=T,sep="\t",stringsAsFactors=F)
p1 <- ggplot(plottab[(plottab$gen2 %in% seq(210,400,10))&(plottab$AF >= 0.01),], aes(x=AF,)) 
p1 <- p1 + geom_vline(xintercept=c(0.5,0.25,0.125,0.0625,0.03125,0),lty=2,col="gray90")
p1 <- p1 + geom_histogram(binwidth=0.01)
p1 <- p1 + theme_classic()
p1 <- p1 + xlim(0,0.6) + ylim(0,100)
p1 <- p1 + facet_wrap(~Crypt_gen,  ncol=4)
p1


### D: Running the simulation at different variable settings to gather plotting data for Shiny App ###
### To reduce simulation memory and storage requirements (due to large tempdf), a different strategy for generating tempdf and sumdf is implemented. Instead of generating the full tempdf followed by sumdf, sumdf is progressively populated after every row of tempdf is generated, allowing earlier rows of tempdf to be deleted progressively as well. Therefore the entirety of tempdf is not stored.
library(tidyr)
library(reshape2)
rm(list=ls())

# List of settable variables
# sim           : number of simulated generations
# numstemcells  : number of starting stem cell(s) in niche
# numcryptcells : crypt size i.e. number of cells it can hold
# mutrate       : number of mutation(s) per genome per generation

sim=500
numcryptcells=2048
for (numstemcells in c(4,8,12,16)) { 
  #numstemcells=4  numstemcells=8
  for (mutrate in 1:3) { 
    #mutrate=1
    cat("\nmutrate=",mutrate,",numSC=",numstemcells," - ",sep="")

    # Initialize a single data frame to store stem cell niche and crypt data. 
    tempdf=data.frame()
    tempdf["sim0","SC_gen"]="scgen0"
    tempdf["sim0",paste("StemCell",1:(numstemcells*2),sep="")]=c(paste("sc",1:numstemcells,"//",sep=""),rep(NA,numstemcells))
    tempdf["sim0","keep"]=paste(1:numstemcells,collapse=",")
    tempdf["sim0","push"]=NA
    tempdf["sim0","CC_gen"]="ccgen0"
    tempdf["sim0",paste("CryptCell",1:numcryptcells,sep="")]=rep(NA,numcryptcells)
    #head(tempdf)
    #tempdf[,1:25]
    
    sumdf=data.frame()
    sumdf["ccgen0","Crypt_gen"]="ccgen0"
    
    # start running through generations to populate tempdf
    mutcount=0
    simgencount=0
    scgencount=0
    ccgencount=0
    for (i in 1:sim) {
      #i=1 i=2 i=3
      simgencount=i
      currrowname=paste("sim",i,sep="")
      tempdf[currrowname,]=NA
      
      ### SC Niche 
      if (simgencount%%2==0) {
        scprevgen=paste("scgen",scgencount,sep="")
        scgencount=scgencount+1
        sccurrgen=paste("scgen",scgencount,sep="")
        tempdf[currrowname,"SC_gen"]=sccurrgen
        scprevkeep=as.numeric(unlist(strsplit(tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==scprevgen),"keep"],",")))
        
        scfilled=0
        for (j in scprevkeep) {
          #j=1 j=2 j=3
          
          scparent=tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==scprevgen),paste("StemCell",j,sep="")]
          scparentid=sub("//\\S*","",scparent)
          scparentmut=sub("^[^/]+//","",scparent)
          
          d1id=paste(scparentid,".1",sep="")
          d1mut=scparentmut
          for (k in 1:mutrate) {
            mutcount=mutcount+1;
            d1mut=paste(d1mut,":m",mutcount,sep="")
          }
          d1mut=sub("^:","",d1mut)
          d1=paste(d1id,d1mut,sep="//")
          #d1
          scfilled=scfilled+1
          tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),paste("StemCell",scfilled,sep="")]=d1
          
          d2id=paste(scparentid,".2",sep="")
          d2mut=scparentmut
          for (k in 1:mutrate) {
            mutcount=mutcount+1;
            d2mut=paste(d2mut,":m",mutcount,sep="")
          }
          d2mut=sub("^:","",d2mut)
          d2=paste(d2id,d2mut,sep="//")
          #d2
          scfilled=scfilled+1
          tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),paste("StemCell",scfilled,sep="")]=d2
        }
        select=sample(1:(numstemcells*2), (numstemcells*2), replace = FALSE, prob = NULL)
        tokeep=sort(select[1:numstemcells])
        topush=(select[(numstemcells+1):(numstemcells*2)])
        tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),"keep"]=paste(tokeep,collapse=",")
        tempdf[(!is.na(tempdf$SC_gen))&(tempdf$SC_gen==sccurrgen),"push"]=paste(topush,collapse=",")
      }
      
      ### Crypt 
      ccprevgen=paste("ccgen",ccgencount,sep="")
      ccgencount=ccgencount+1
      cccurrgen=paste("ccgen",ccgencount,sep="")
      tempdf[currrowname,"CC_gen"]=cccurrgen
      
      ccfilled=0
      scpushes=tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),"push"]
      if (!is.na(scpushes)) {
        for (scpush in unlist(strsplit(scpushes,","))) {
          ccfilled=ccfilled+1
          tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=sub("//",".cc//",tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("StemCell",scpush,sep="")])
        }
      } 
      
      for (j in 1:numcryptcells) {
        #j=1 j=2 j=3
        ccparent=tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==ccprevgen),paste("CryptCell",j,sep="")]
        if (is.na(ccparent)) {
          break
        }
        if (ccfilled > numcryptcells) {
          break
        }
        
        ccparentid=sub("//\\S*","",ccparent)
        ccparentmut=sub("^[^/]+//","",ccparent)
        
        d1id=paste(ccparentid,".1",sep="")
        d1mut=ccparentmut
        for (k in 1:mutrate) {
          mutcount=mutcount+1;
          d1mut=paste(d1mut,":m",mutcount,sep="")
        }
        d1mut=sub("^:","",d1mut)
        d1=paste(d1id,d1mut,sep="//")
        #d1
        ccfilled=ccfilled+1
        if (ccfilled <= numcryptcells) {
          tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=d1
        }
        
        d2id=paste(ccparentid,".2",sep="")
        d2mut=ccparentmut
        for (k in 1:mutrate) {
          mutcount=mutcount+1;
          d2mut=paste(d2mut,":m",mutcount,sep="")
        }
        d2mut=sub("^:","",d2mut)
        d2=paste(d2id,d2mut,sep="//")
        #d2
        ccfilled=ccfilled+1
        if (ccfilled <= numcryptcells) {
          tempdf[(!is.na(tempdf$CC_gen))&(tempdf$CC_gen==cccurrgen),paste("CryptCell",ccfilled,sep="")]=d2
        }
      }
      if (i%%100==0) {
        mutcount
      }
      if (i%%10==0) {
        cat(i," ",sep="")
        tempdf=tempdf[row.names(tempdf) %in% paste("sim",(i-5):(i+5),sep=""),]
      }
      
      sumdf[cccurrgen,"Crypt_gen"]=cccurrgen
      if (i>1) {
        muts=unlist(strsplit(sub("^\\S+//","",tempdf[tempdf$CC_gen==cccurrgen,paste("CryptCell",1:numcryptcells,sep="")]),":"))
        muts=muts[!is.na(muts)]
        muts=muts[muts != "NA"]
        
        sumdf[cccurrgen,"num_mut"]=length(muts)
        sumdf[cccurrgen,"num_uniq_mut"]=length(unique(muts))
        for (j in (sort(unique(table(muts)),decreasing=T))) {
          #j=1 j=4
          sumdf[cccurrgen,paste("num_AF",(j/(numcryptcells*2)),sep="")]=sum(table(muts)==j)
        }
      }
    }
    #tempdf[1:6,1:16]
    #tempdf[101:105,1:16]

    for (j in 4:length(names(sumdf))) {
      #j=4  j=30
      sumdf[is.na(sumdf[,j]),j]=0
    }
    #head(sumdf)
    #tail(sumdf)
    sumoutfile=paste("sumdf_mutrate",mutrate,"_startSC",numstemcells,".txt",sep="")
    write.table(sumdf,sumoutfile,quote=F,sep="\t",row.names=T,col.names=T)
    #sumdf<-read.table(sumoutfile,header=T,sep="\t",stringsAsFactors=F)
    
    #create data tab to plot histograms
    #sumdf[1:5,1:5]
    #library(tidyr)
    #library(reshape2)
    plottab=melt(sumdf[,c(1,4:length(names(sumdf)))],"Crypt_gen")
    plottab=plottab[plottab$value>0,]
    #as_tibble(plottab)
    #uncount(as_tibble(plottab),value)
    plottab=as.data.frame(uncount(as_tibble(plottab),value))
    plottab$variable=as.numeric(sub("num_AF","",plottab$variable))
    names(plottab)[2]="AF"
    plottab$gen2=as.numeric(sub("ccgen","",plottab$Crypt_gen))
    plottab$Crypt_gen=sprintf("ccgen%03d", plottab$gen2)
    #head(plottab)
    plottab=plottab[plottab$AF>0.01,]
    plotoutfile=paste("plottab_mutrate",mutrate,"_startSC",numstemcells,".txt",sep="")
    write.table(plottab,plotoutfile,quote=F,sep="\t",row.names=T,col.names=T)
    #plottab<-read.table(plotoutfile,header=T,sep="\t",stringsAsFactors=F)
    #head(plottab)
    #summary(plottab$AF)
    #table(plottab$AF)
  }
}

# combine all the plottab files into a single plottaball file
plotfiles=dir(".",full.names=T)
plotfiles=plotfiles[grep("plottab_mutrate",plotfiles)]
plotfiles
plottaball=data.frame()
for (i in 1:length(plotfiles)) {
  #i=1
  plotfile=plotfiles[i]
  mutrate=as.numeric(sub("_\\S+$","",sub("^\\S+mutrate","",plotfile)))
  startsc=as.numeric(sub(".txt","",sub("^\\S+startSC","",plotfile)))
  plottab<-read.table(plotfile,header=T,sep="\t",stringsAsFactors=F)
  plottab$mutrate=mutrate
  plottab$numstemcells=startsc
  plottab
  if (i==1) {
    plottaball=plottab
  }else {
    plottaball=rbind(plottaball,plottab)
  }
  cat(i," ",sep="")
}
head(plottaball)
head(plottaball$AF)
head(plottaball$gen2)
head(plottaball$mutrate)
head(plottaball$numstemcells)
write.table(plottaball,"ShinyApp/plottaball.txt",quote=F,sep="\t",row.names=T,col.names=T)




