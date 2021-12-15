##import libraries
library(ape)
library(phangorn)
library(optparse)

##Define Arguments

option_list = list(
  #make_option(c("-t", "--tree"), action="store", type="character", help="Table (in comma-separated values) from parseRepeats.py [default %default]", metavar="string");
  make_option(c("-d", "--directory"), action="store", type="character", help="Table (in comma-separated values) from parseRepeats.py [default %default]", metavar="string"),
  make_option(c("-o", "--outfile"), action="store", type="character", help="Table (in comma-separated values) from parseRepeats.py [default %default]", metavar="string"));
opt = parse_args(OptionParser(option_list=option_list))

outfile <- opt$outfile
folder <- opt$directory


folder_name = strsplit(folder,"_")
date_name = strsplit(folder,"_R_")
date_name1 = strsplit(date_name[[1]],"D_")
dst <- date_name1[[1]][2]
ded <- date_name[[1]][2]
p = folder_name[[1]][1]

### Function Definitions

# This function checks if the reference outgroup subtype is monophyletic
# If the reference outgroup subtype is monophyletic it will root the tree on those labels
# Finally, the terminals (tips) containing the outgroup will be removed (dropped)
# The resulting tree can be assessed and assigned a topological class
monoRootdrop <- function(treefile){
  trs <- read.tree(treefile)
  if(is.monophyletic(trs, grep("Ref",trs$tip.label))){
  tr <- root(trs, outgroup = grep("Ref",trs$tip.label))
  tr1 <- drop.tip(tr, tip = grep("Ref",tr$tip.label))
  return(tr1)
  }
  else{
    return(trs)
  }
}
#-----------------------------------------------
# Calculates the 
# 1) # of mutations in the tree
# 2) root label 
# 3) topological class
# 4) and number of monophyletic clades of class A and B
  num.intro = function(tr, lab1, lab2, mu=0.80e-2/365, sig=0.14e-2/365, seq.len=346, gen.len=2){
    if(!is.rooted(tr)){return(NULL)}
    if(!all(tr$tip.label%in%c(lab1, lab2))){print("Tip labels not covered"); return(NULL)}
    root = Ntip(tr)+1
    ind1 = which(tr$tip.label==lab1)
    ind2 = which(tr$tip.label==lab2)
    lab.vect = c(tr$tip.label, rep(NA, Nnode(tr)))
    anc.vect = c(tr$tip.label, rep(NA, Nnode(tr)))
    for(x in (Ntip(tr)+Nnode(tr)):root){
      kids.id = tr$edge[tr$edge[,1]==x,2]
      kids.lab = lab.vect[kids.id]    
      kids.anc = anc.vect[kids.id]
      #--------------------Propagation for number monophyletic clades
      # labels are the same, propagate up
      if(kids.lab[1]==kids.lab[2]){
        lab.vect[x] = kids.lab[1]
        lab.vect[kids.id] = '?'
      }else{
        # labels are A? in either order 
        if(setequal(kids.lab, c(lab1, "?"))){lab.vect[x] = '?'}
        # labels are B? in either order 
        if(setequal(kids.lab, c(lab2, "?"))){lab.vect[x] = '?'}
        # labels are AB in either order 
        if(setequal(kids.lab, c(lab1, lab2))){lab.vect[x] = "?"}
      }
      #--------------------Propagation for root label
      # lables are the same, propagate up
      if(kids.anc[1]==kids.anc[2]){
        anc.vect[x] = kids.anc[1]
      }else{
        # labels are A? in either order 
        if(setequal(kids.anc, c(lab1, "?"))){anc.vect[x] = lab1}
        # labels are B? in either order 
        if(setequal(kids.anc, c(lab2, "?"))){anc.vect[x] = lab2}
        # labels are AB in either order 
        if(setequal(kids.anc, c(lab1, lab2))){anc.vect[x] = "?"}
      }
    }
    # All trees have one or more AB coal, if more than one AB coal->PP topo, if root AB-> MM topo else PM topo 
    if(sum(anc.vect=='?')>1){
      topo='PP'
      nd = node.depth.edgelength(tr)
      nd = nd/max(nd)
      tIA = nd[which(lab.vect==lab1)]
      tIB = nd[which(lab.vect==lab2)]
      tIA.mean = mean(tIA)
      tIB.mean = mean(tIB)
      tIA.sd = sd(tIA)
      tIB.sd = sd(tIB)
    }else{
      nd = node.depth.edgelength(tr)
      nd = nd/max(nd)
      tIA = nd[which(lab.vect==lab1)]
      tIB = nd[which(lab.vect==lab2)]
      tIA.mean = mean(tIA)
      tIB.mean = mean(tIB)
      tIA.sd = 0
      tIB.sd = 0
      if(anc.vect[root]!='?')topo='PM' else topo = "MM" 
    }
    theta = sig^2/mu
    k = mu^2/sig^2 
    evo.rate = rgamma(1, k, scale=theta)
    t.dist=sum(tr$edge.length*gen.len*evo.rate*seq.len)
    
    d = cophenetic(tr)
    ind.A = which(rownames(d)==lab1)
    ind.B = which(rownames(d)==lab2)
    tmp.AB = d[ind.A,ind.B]
    tmp.A = d[ind.A,ind.A]
    tmp.B = d[ind.B,ind.B]
    pwd.AB = mean(tmp.AB[upper.tri(tmp.AB)]*gen.len*seq.len*evo.rate)
    pwd.A = mean(tmp.A[upper.tri(tmp.A)]*gen.len*seq.len*evo.rate)
    pwd.B = mean(tmp.B[upper.tri(tmp.B)]*gen.len*seq.len*evo.rate)
    
    
    if(any(tr$edge.length<0))print("Negative branch lengths in sims")
    
    return(list(anc=anc.vect[root], IA=unname(table(lab.vect)[lab1]), IB=unname(table(lab.vect)[lab2]), topo=topo, 
                t.dist=t.dist, pwd.AB=pwd.AB, pwd.A=pwd.A, pwd.B=pwd.B, tIA.mean=tIA.mean, tIA.sd=tIA.sd, tIB.mean=tIB.mean, tIB.sd=tIB.sd))
  }

### MAIN SCRIPT

##Read Data
treeList <- list.files(path=folder, pattern="*.tre$", full.names=TRUE, recursive=TRUE)
phyList <- list.files(path=folder, pattern="*.phy$", full.names=TRUE, recursive=TRUE)

"""
x Start_Coord
x End_Coord
Outgroup_monophyly
Ancestral_state
IA
IB
t.dist
pwd.AB
pd.A
pwd.B
tIA.mean
tIA.sd
tIB.mean
tIB.sd
Topology
x Number Donor Tips
x Number Recipient Tips
Total number of tips
x Donor Date
x Recipient Date
mean blen
"""

## stats1
pair <-c()
dateStart <-c()       # "Donor Date"
dateEnd <-c()         # "Recipient Date"
treeStart <-c()       # "Start_Coord"
treeEnd <-c()         # "End_Coord"
ConsIN <-c()
HeIN <-c()
RetIN <-c()
numberRecipient <-c() # "Number Donor Tips"
numberDonor <-c()     # "Number Recipient Tips"

## stats2
outgroup <-c()
ancestral <-c()
ia <-c()
ib <-c()
tdist <-c()
pwdAB <-c()
pdA <-c()
pwdB <-c()
tIAmean <-c()
tIAstd <-c()
tIBmean <-c()
tIBstd <-c()
topology <-c()
meanBlen <-c()

for(i in 1:length(treeList)){
  pair <- c(pair,p)
  dateStart <- c(dateStart,dst)
  dateEnd <- c(dateEnd,ded)
  tname = treeList[i]
  phyname = phyList[i]
  tname_1 =strsplit(tname,"/")
  tname_split = strsplit(tname_1[[1]][3], "_")
  treeStart <- c(treeStart,tname_split[[1]][1])
  treeEnd <- c(treeEnd,tname_split[[1]][2])

  treedat <-read.tree(tname)
  phydat <-read.phyDat(phyname,format = "interleaved",type = "DNA")

  numberRecipient<-c(numberRecipient,length(grep("Recipient",treedat$tip.label)))
  numberDonor<-c(numberDonor,length(grep("Donor",treedat$tip.label)))
  
  ## consistency index
  ConsIN <- c(ConsIN, CI(treedat,phydat))
  ## homoplasy index
  HeIN <- c(HeIN, 1-CI(treedat,phydat))
  ## retention index
  RetIN <- c(RetIN, RI(treedat,phydat))

  tf=monoRootdrop(treefile = fileList[i])
  numberRecipients<-length(grep("Recipient",tf$tip.label))
  numberDonors<-length(grep("Donor",tf$tip.label))
  totaltips<-numberDonors+numberRecipients
  ablen<- mean(tf$edge.length)
  
  terminal = sapply(strsplit(tf$tip.label,"_"),"[",2)
  tf$tip.label= ifelse(terminal=="Recipient", "Recipient", "Donor")
  
  #original
  A=num.intro(tf, "Recipient", "Donor")
  
  #swap recipient with donor
  #A=num.intro(tf, "Donor", "Recipient")
  
  if(!is.null(A[4])){
  cat(c(start_coord,end_coord,"yes",A[[1]],A[[2]],A[[3]],A[[5]],A[[6]],A[[7]],A[[8]],A[[9]],A[[10]],A[[11]],A[[12]],A[[4]],numberDonors,numberRecipients,totaltips,dt[1],rt[1],ablen,"\n"),file=output,sep="\t",append = TRUE)
  }
  #instead of printing to console here, print to file
  else{
    cat(c(start_coord,end_coord,"no","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","NMO",numberDonors,numberRecipients,totaltips,dt[1],rt[1],ablen,"\n"),file=output,sep="\t",append = TRUE)
  }

}

measures.data <- data.frame(pair,dateStart,dateEnd,treeStart,treeEnd,ConsIN,HeIN,RetIN,numberDonor,numberRecipient)
write.csv(measures.data,outfile, row.names = FALSE)