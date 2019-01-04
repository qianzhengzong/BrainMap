library("mvMORPH")
mytree <- read.nexus("/home/qzz/Documents/brain_map_analysis/tree/consensusTree_10kTrees_Primates_Version3_consensus_tree.nex")
new.exp.norm <-exp.norm
mean.exp <- cbind(apply(new.exp.norm[,1:4],1,mean),apply(new.exp.norm[,5:7],1,mean),
                  apply(new.exp.norm[,8:10],1,mean),apply(new.exp.norm[,11:13],1,mean))
mean.exp <- cbind(mean.exp[,1],mean.exp[,4],mean.exp[,3],mean.exp[,2])
colnames(mean.exp) <- mytree$tip.label
data <- t(mean.exp)

sta <- as.vector(mytree$tip.label); names(sta) <- mytree$tip.label

tree <- make.simmap(mytree,sta,model = "ER",nsim = 1)

data.select <- data[,1:2]
data.select.result <- mvOU(tree,data.select,model = "OU1")
