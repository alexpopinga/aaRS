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

cnames <- c("AARS_Y_", "AARS_M_", "AARS_W_", "AARS_K_", "AARS_L_", "AARS_V_", "AARS_I_", "AARS_E_", "AARS_Q_", "AARS_R_", 
"AARS_C_", "AARS_A_", "AARS_D_", "AARS_N_", "AARS_P_", "AARS_F_", "AARS_S_", "AARS_T_", "AARS_H_", "AARS_G_","Y", "M", 
"W", "K", "L", "V", "I", "E", "Q", "R", "C", "A", "D", "N", "P", "F", "S", "T", "H", "G")

for (i in 1:length(trees)) {
	tree <- trees[[i]]
	
	# funky things!
	tree$tip.label <- gsub("AARS_R_", "R", tree$tip.label)
	tree$tip.label <- gsub("AARS_Y_", "Y", tree$tip.label)
	tree$tip.label <- gsub("AARS_K_", "K", tree$tip.label)
	tree$tip.label <- gsub("AARS_M_", "M", tree$tip.label)
	tree$tip.label <- gsub("AARS_W_", "W", tree$tip.label)
	tree$tip.label <- gsub("AARS_L_", "L", tree$tip.label)
	tree$tip.label <- gsub("AARS_V_", "V", tree$tip.label)
	tree$tip.label <- gsub("AARS_I_", "I", tree$tip.label)
	tree$tip.label <- gsub("AARS_E_", "E", tree$tip.label)
	tree$tip.label <- gsub("AARS_Q_", "Q", tree$tip.label)
	tree$tip.label <- gsub("AARS_C_", "C", tree$tip.label)
	tree$tip.label <- gsub("AARS_A_", "A", tree$tip.label)
	tree$tip.label <- gsub("AARS_D_", "D", tree$tip.label)
	tree$tip.label <- gsub("AARS_N_", "N", tree$tip.label)
	tree$tip.label <- gsub("AARS_P_", "P", tree$tip.label)
	tree$tip.label <- gsub("AARS_F_", "F", tree$tip.label)
	tree$tip.label <- gsub("AARS_S_", "S", tree$tip.label)
	tree$tip.label <- gsub("AARS_T_", "T", tree$tip.label)
	tree$tip.label <- gsub("AARS_H_", "H", tree$tip.label)
	tree$tip.label <- gsub("AARS_G_", "G", tree$tip.label)
	
		
	coltree <- tree	
	for (j in 1:length(cnames)) {
	  print(paste0("Collapsing ",cnames[j]))
	  coltree <- collapse_tree(coltree, cnames[j])
	}
	ct[[i]] <- coltree
	write.tree(ct[[i]],file="~/NormalBEAST2_sim1_collapsed.tree")
}

#plot(ct[[1]])