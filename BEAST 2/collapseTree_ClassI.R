library(phytools)
library(phangorn)


collapse_tree <- function(tree, matchString) {

    nodes<-1:tree$Nnode+Ntip(tree) ## all nodes
    subtrees<-list()

    for(i in 1:tree$Nnode) subtrees[[i]]<-extract.clade(tree,nodes[i])
    names(subtrees)<-nodes ## all subtrees

    ## all nodes with only taxa names matching matchString as descendant
    all.foos<-nodes[sapply(subtrees,function(x) all(grepl(matchString, x$tip.label)))]
    
    if (length(all.foos) > 0) {
       foo.mrcas<-all.foos[sapply(all.foos,function(x,tree,y) !Ancestors(tree,x,"parent")%in%y, tree=tree,y=all.foos)]

       w.foos<-grep(matchString , tree$tip.label)

       collapsed<-tree
       ## iterate over all MRCAs of foo clades
       for(i in 1:length(foo.mrcas)){
           M<-matchNodes(tree,collapsed)
           nn<-M[which(M[,1]==foo.mrcas[i]),2]
           dd<-Descendants(collapsed,nn)[[1]]
           h<-sapply(dd,nodeheight,tree=collapsed)
           collapsed$tip.label[dd[1]]<-paste(length(dd),matchString)
           ind<-which(collapsed$edge[,2]==dd[1])
           collapsed$edge.length[ind]<-collapsed$edge.length[ind]+mean(h)-h[1]
           if(length(dd)>1)
             collapsed<-drop.tip(collapsed,collapsed$tip.label[dd[2:length(dd)]])
       }
    } else {
       collapsed<-tree 
    }

    return (collapsed)
}

trees <- read.nexus(file="~/ClassI_Arch_BD.trees")

trees <- trees[999:1001]

ct <- list()

cnames <- c("tyr_", "met_", "trp_", "lys_", "leu_", "val_", "ile_", "glu_", "gln_", "arg_", "cys_")

for (i in 1:length(trees)) {
	tree <- trees[[i]]
	
	# funky things!
	tree$tip.label <- gsub("try_", "tyr_", tree$tip.label)
	tree$tip.label <- gsub("lysi_", "lys_", tree$tip.label)
	tree <- drop.tip(tree, grep("gln_", tree$tip.label))
		
	coltree <- tree	
	for (j in 1:length(cnames)) {
	  print(paste0("Collapsing ",cnames[j]))
	  coltree <- collapse_tree(coltree, cnames[j])
	}
	ct[[i]] <- coltree
}

plot(ct[[1]])