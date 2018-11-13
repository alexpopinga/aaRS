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

trees <- read.nexus(file="~/NormalBEAST2_sim1_consensus.tree")

trees <- trees[0:1]

ct <- list()

cnames <- c("arg", "cys", "glu", "gln", "ile", "leu", "met", "tyr", "trp", "val", "lys", "thr", "ser", "pro", "phe", "his", "gly", "asp", "asn", "ala")

for (i in 1:length(trees)) {
	tree <- trees[[i]]
	
	# funky things!
	tree$tip.label <- gsub("AARS_R_", "arg", tree$tip.label)
	tree$tip.label <- gsub("AARS_C_", "cys", tree$tip.label)
	tree$tip.label <- gsub("AARS_E_", "glu", tree$tip.label)
	tree$tip.label <- gsub("AARS_Q_", "gln", tree$tip.label)
	tree$tip.label <- gsub("AARS_I_", "ile", tree$tip.label)
	tree$tip.label <- gsub("AARS_L_", "leu", tree$tip.label)
	tree$tip.label <- gsub("AARS_M_", "met", tree$tip.label)
	tree$tip.label <- gsub("AARS_Y_", "tyr", tree$tip.label)
	tree$tip.label <- gsub("AARS_W_", "trp", tree$tip.label)
	tree$tip.label <- gsub("AARS_V_", "val", tree$tip.label)
	tree$tip.label <- gsub("AARS_K_", "lys", tree$tip.label)
	tree$tip.label <- gsub("AARS_T_", "thr", tree$tip.label)
	tree$tip.label <- gsub("AARS_S_", "ser", tree$tip.label)
	tree$tip.label <- gsub("AARS_P_", "pro", tree$tip.label)
	tree$tip.label <- gsub("AARS_F_", "phe", tree$tip.label)
	tree$tip.label <- gsub("AARS_H_", "his", tree$tip.label)
	tree$tip.label <- gsub("AARS_G_", "gly", tree$tip.label)
	tree$tip.label <- gsub("AARS_D_", "asp", tree$tip.label)
	tree$tip.label <- gsub("AARS_N_", "asn", tree$tip.label)
	tree$tip.label <- gsub("AARS_A_", "ala", tree$tip.label)
		
	coltree <- tree	
	for (j in 1:length(cnames)) {
	  print(paste0("Collapsing ",cnames[j]))
	  coltree <- collapse_tree(coltree, cnames[j])
	}
	ct[[i]] <- coltree
}

plot(ct[[1]])