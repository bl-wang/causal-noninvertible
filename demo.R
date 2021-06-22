rm(list = ls())


setwd("~/Desktop/RCode")
source(file = "pairdatagenerator.R")
source(file = "DAGdatagen.R")
source(file = "DAGupdate.R")


n = 500; nnode = 4; nedge = 3; nedge.nl = 3
node_labels = paste("X", 1:nnode, sep = "")

g = matrix(0, nnode, nnode)
colnames(g) = rownames(g) = node_labels
g[1, 3] = g[2, 3] = g[3, 4] = 1



res = NULL

seeds=12345
res$dagList = DAGdatagen(n = n, nnode = nnode, nedge = nedge, nedge.nl = nedge.nl,
                         dagEdgeList = edgeList(g), se = 0.5,
                         labels = node_labels)

# pairs(res$dagList$DAGdata, cex = 0.1)


res$resList = DAGupdate(dat = res$dagList$DAGdata, labels = node_labels, p_value = 0.01, 
            stat = "pearson", iter=100, ntau = min(n/20, 100), 
            UsePC = T, UseSparsebn = T, UseEmptyG = T, conIndTest = T)


### true dag
plot(res$dagList$dag.bn)
### true cpdag
plot(res$dagList$cpdag.bn)




### pc cpdag
plot(res$resList$res1$res_pc$G_est$cpdag.learn)
### pc + NNCL restricted cpdag
plot(res$resList$res1$res_pc$G_update$cpdagAddBN)
### NNCL+pc restricted cpdag
plot(res$resList$res2$res_pc$G_update$cpdagAddBN)


### sparsebn dag
plot(res$resList$res1$res_sbn$G_est$G.learn) 
### sparsebn cpdag
plot(res$resList$res1$res_sbn$G_est$cpdag.learn) 
### sparsebn + NNCL restricted cpdag
plot(res$resList$res1$res_sbn$G_update$cpdagAddBN) 
### NNCL + sparsebn restricted cpdag
plot(res$resList$res2$res_sbn$G_update$cpdagAddBN)

### NNCL dag
plot(res$resList$res1$res_emptyG$cpdagAddBN)




