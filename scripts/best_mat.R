# Find best performing matrices for experimen

best_mat4exp = function(roc, exp, mat) {
N=dim(roc)[1]; L=dim(roc)[2]
Roc=rep(0,N); Mat=rep(0,N); Gene_exp=rep(0,N); Gene_mat=rep(0,N);
for(i in 1:N) {n=order(roc[i,])[L];
Roc[i]=roc[i,n]; Mat[i]=mat[n,2]
Gene_exp[i]=exp[i,3]; Gene_mat[i]=mat[n,3]}
best_mat4exp=cbind(exp[,2], Mat, Gene_exp, Gene_mat, Roc)
return(best_mat4exp)
}

best_mat4gene = function(roc,exp,mat,genes2exp) {
genes=unique(genes2exp[,1])

rank_mat=roc; N=dim(roc)[1]; L=dim(roc)[2];
troc=t(roc); for(i in 1:N) {rank_mat[i,]=L-rank(troc[,i])+1}

gm_mean = function(a){exp(mean(log(a)))}
N=length(genes); rank=matrix(0,nrow=N, ncol=L)
bestmat=rep("",N); bestgmr=rep(0,N)
rownames(rank)=genes; colnames(rank)=colnames(roc)
for(i in 1:N) {
xmat=rank_mat[genes2exp[which(genes2exp[,1] == genes[i]),2],]
gm=sort(apply(xmat,2,gm_mean))
bestmat[i]=names(gm[1]); bestgmr[i]=gm[1]}
best_mat4gene=cbind(genes, bestmat, bestgmr, mat[bestmat,3])
return(best_mat4gene)
}

best_analysis = function(name,roc,exp,mat,genes2exp) {
best=best_mat4exp(roc,exp, mat)
write.table(best,paste0("../results/best4exp_",name,".txt"),
sep="\t", quote=F,col.names=F,row.names=F)
best=best_mat4gene(roc,exp, mat, genes2exp)
write.table(best,paste0("../results/best4gene_",name,".txt"),
sep="\t", quote=F,col.names=F,row.names=F)
}

read_hocomoco = function() {
options(stringsAsFactors = FALSE)
mat=          read.table("../hocomoco.txt")
rownames(mat)=mat[,1]; mat <<- mat}
read_jaspar = function() {
options(stringsAsFactors = FALSE)
mat=          read.table("../jaspar.txt")
rownames(mat)=mat[,1]; mat <<- mat}
read_cisbp = function() {
options(stringsAsFactors = FALSE)
mat=          read.table("../cisbp.txt")
rownames(mat)=mat[,1]; mat <<- mat}

read_mat_all = function() {
options(stringsAsFactors = FALSE)
mat=          read.table("../hocomoco.txt")
mat=rbind(mat,read.table("../jaspar.txt"))
mat=rbind(mat,read.table("../cisbp.txt"))
rownames(mat)=mat[,1]; mat <<- mat}

read_remap = function() {
options(stringsAsFactors = FALSE)
exp=read.table("../remap.txt"); rownames(exp)=exp[,1];
genes2exp=read.table("../genes2remap.txt")
exp <<- exp; genes2exp <<- genes2exp}

read_jolma_yang = function() {
options(stringsAsFactors = FALSE)
exp=read.table("../jolma_yang.txt"); rownames(exp)=exp[,1];
genes2exp=read.table("../genes2jolma_yang.txt")
exp <<- exp; genes2exp <<- genes2exp}

read_uniprobe2 = function() {
options(stringsAsFactors = FALSE)
exp=read.table("../uniprobe2.txt"); rownames(exp)=exp[,1];
genes2exp=read.table("../genes2uniprobe2.txt")
exp <<- exp; genes2exp <<- genes2exp}

read_exp_all = function() {
options(stringsAsFactors = FALSE)
exp=          read.table("../remap.txt")
exp=rbind(exp,read.table("../jolma_yang.txt")) 
exp=rbind(exp,read.table("../uniprobe2.txt")) 
rownames(exp)=exp[,1]
genes2exp=                read.table("../genes2remap.txt")
genes2exp=rbind(genes2exp,read.table("../genes2jolma_yang.txt")) 
genes2exp=rbind(genes2exp,read.table("../genes2uniprobe2.txt"))
exp <<- exp; genes2exp <<- genes2exp}

read_remap_all = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../remap_hocomoco_roc.txt")
roc=cbind(roc,read.table("../remap_jaspar_roc.txt"))
roc=cbind(roc,read.table("../remap_cisbp_roc.txt"))
roc <<- roc
}
read_jolma_yang10_all = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../jolma_yang_hocomoco_roc10.txt")
roc=cbind(roc,read.table("../jolma_yang_jaspar_roc10.txt"))
roc=cbind(roc,read.table("../jolma_yang_cisbp_roc10.txt"))
roc <<- roc
}

read_jolma_yang50_all = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../jolma_yang_hocomoco_roc50.txt")
roc=cbind(roc,read.table("../jolma_yang_jaspar_roc50.txt"))
roc=cbind(roc,read.table("../jolma_yang_cisbp_roc50.txt"))
roc <<- roc
}

read_uniprobe2_all = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../uniprobe2_hocomoco_cor.txt")
roc=cbind(roc,read.table("../uniprobe2_jaspar_cor.txt"))
roc=cbind(roc,read.table("../uniprobe2_cisbp_cor.txt"))
roc <<- roc
}

read_all_hocomoco = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../remap_hocomoco_roc.txt")
roc=rbind(roc,read.table("../jolma_yang_hocomoco_roc10.txt"))
roc=rbind(roc,read.table("../uniprobe2_hocomoco_cor.txt"))
roc <<- roc
}

read_all_jaspar = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../remap_jaspar_roc.txt")
roc=rbind(roc,read.table("../jolma_yang_jaspar_roc10.txt"))
roc=rbind(roc,read.table("../uniprobe2_jaspar_cor.txt"))
roc <<- roc
}

read_all_cisbp = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../remap_cisbp_roc.txt")
roc=rbind(roc,read.table("../jolma_yang_cisbp_roc10.txt"))
roc=rbind(roc,read.table("../uniprobe2_cisbp_cor.txt"))
roc <<- roc
}

read_all_all = function() {
options(stringsAsFactors = FALSE)
roc1=          read.table("../remap_hocomoco_roc.txt")
roc1=cbind(roc1,read.table("../remap_jaspar_roc.txt"))
roc1=cbind(roc1,read.table("../remap_cisbp_roc.txt"))
roc2=          read.table("../jolma_yang_hocomoco_roc10.txt")
roc2=cbind(roc2,read.table("../jolma_yang_jaspar_roc10.txt"))
roc2=cbind(roc2,read.table("../jolma_yang_cisbp_roc10.txt"))
roc3=          read.table("../uniprobe2_hocomoco_cor.txt")
roc3=cbind(roc3,read.table("../uniprobe2_jaspar_cor.txt"))
roc3=cbind(roc3,read.table("../uniprobe2_cisbp_cor.txt"))
roc=rbind(roc1,roc2,roc3)
roc <<- roc
}

read_all_hocomoco_filtered = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../filtered/remap_hocomoco_roc.txt")
roc=rbind(roc,read.table("../filtered/jolma_yang_hocomoco_roc10.txt"))
roc=rbind(roc,read.table("../filtered/uniprobe2_hocomoco_cor.txt"))
roc <<- roc
}

read_all_jaspar_filtered = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../filtered/remap_jaspar_roc.txt")
roc=rbind(roc,read.table("../filtered/jolma_yang_jaspar_roc10.txt"))
roc=rbind(roc,read.table("../filtered/uniprobe2_jaspar_cor.txt"))
roc <<- roc
}

read_all_cisbp_filtered = function() {
options(stringsAsFactors = FALSE)
roc=          read.table("../filtered/remap_cisbp_roc.txt")
roc=rbind(roc,read.table("../filtered/jolma_yang_cisbp_roc10.txt"))
roc=rbind(roc,read.table("../filtered/uniprobe2_cisbp_cor.txt"))
roc <<- roc
}

# analyze remap

exp=0; options(stringsAsFactors = FALSE)

read_remap(); read_hocomoco(); name="remap_hocomoco";
roc=read.table("../remap_hocomoco_roc.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_remap(); read_jaspar(); name="remap_jaspar";
roc=read.table("../remap_jaspar_roc.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_remap(); read_cisbp(); name="remap_cisbp";
roc=read.table("../remap_cisbp_roc.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_mat_all()
read_remap(); read_remap_all(); name="remap_all";
best_analysis(name,roc,exp,mat,genes2exp)

# analyze jolma_yang10

exp=0; options(stringsAsFactors = FALSE)

read_jolma_yang(); read_hocomoco(); name="jolma_yang10_hocomoco";
roc=read.table("../jolma_yang_hocomoco_roc10.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_jolma_yang(); read_jaspar(); name="jolma_yang10_jaspar";
roc=read.table("../jolma_yang_jaspar_roc10.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_jolma_yang(); read_cisbp(); name="jolma_yang10_cisbp";
roc=read.table("../jolma_yang_cisbp_roc10.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_mat_all()
read_jolma_yang(); read_jolma_yang10_all(); name="jolma_yang10_all";
best_analysis(name,roc,exp,mat,genes2exp)

# analyze jolma_yang50

exp=0; options(sstringsAsFactors = FALSE)

read_jolma_yang(); read_hocomoco(); name="jolma_yang50_hocomoco";
roc=read.table("../jolma_yang_hocomoco_roc50.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_jolma_yang(); read_jaspar(); name="jolma_yang50_jaspar";
roc=read.table("../jolma_yang_jaspar_roc50.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_jolma_yang(); read_cisbp(); name="jolma_yang50_cisbp";
roc=read.table("../jolma_yang_cisbp_roc50.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_mat_all()
read_jolma_yang(); read_jolma_yang50_all(); name="jolma_yang50_all";
best_analysis(name,roc,exp,mat,genes2exp)

# analyze uniprobe2  

exp=0; options(stringsAsFactors = FALSE)

read_uniprobe2(); read_hocomoco(); name="uniprobe2_hocomoco";
roc=read.table("../uniprobe2_hocomoco_cor.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_uniprobe2(); read_jaspar(); name="uniprobe2_jaspar";
roc=read.table("../uniprobe2_jaspar_cor.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_uniprobe2(); read_cisbp(); name="uniprobe2_cisbp";
roc=read.table("../uniprobe2_cisbp_cor.txt")
best_analysis(name,roc,exp,mat,genes2exp)

read_mat_all()
read_uniprobe2(); read_uniprobe2_all(); name="uniprobe2_all";
best_analysis(name,roc,exp,mat,genes2exp)

# analyze all

read_exp_all()
read_hocomoco(); read_all_hocomoco(); name="all_hocomoco";
best_analysis(name,roc,exp,mat,genes2exp)

read_exp_all()
read_jaspar(); read_all_jaspar(); name="all_jaspar";
best_analysis(name,roc,exp,mat,genes2exp)

read_exp_all()
read_cisbp(); read_all_cisbp(); name="all_cisbp";
best_analysis(name,roc,exp,mat,genes2exp)

read_exp_all()
read_mat_all(); read_all_all(); name="all_all";
best_analysis(name,roc,exp,mat,genes2exp)

# all filtered

read_exp_all()
read_hocomoco(); read_all_hocomoco();
genes2exp=read.table("../genes2all_hocomoco_filtered.txt")
best=best_mat4gene(roc,exp, mat, genes2exp)
write.table(best,"../results/best4gene_all_hocomoco_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

read_exp_all()
read_jaspar(); read_all_jaspar(); 
genes2exp=read.table("../genes2all_jaspar_filtered.txt")
best=best_mat4gene(roc,exp, mat, genes2exp)
write.table(best,"../results/best4gene_all_jaspar_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

read_exp_all()
read_cisbp(); read_all_cisbp(); 
genes2exp=read.table("../genes2all_cisbp_filtered.txt")
best=best_mat4gene(roc,exp, mat, genes2exp)
write.table(best,"../results/best4gene_all_cisbp_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F) 

read_exp_all()
read_mat_all(); read_all_all(); 
genes2exp=                read.table("../genes2all_hocomoco_filtered.txt")
genes2exp=rbind(genes2exp,read.table("../genes2all_jaspar_filtered.txt"))
genes2exp=rbind(genes2exp,read.table("../genes2all_cisbp_filtered.txt"))
best=best_mat4gene(roc,exp, mat, genes2exp)
write.table(best,"../results/best4gene_all_all_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

exp=read.table("../all_hocomoco_filtered.txt"); rownames(exp)=exp[,1]
mat=read.table("../hocomoco.txt")
read_all_hocomoco_filtered()
best=best_mat4exp(roc,exp, mat)
write.table(best,"../results/best4exp_all_hocomoco_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

exp=read.table("../all_jaspar_filtered.txt"); rownames(exp)=exp[,1]
mat=read.table("../jaspar.txt")
read_all_jaspar_filtered()
best=best_mat4exp(roc,exp, mat)
write.table(best,"../results/best4exp_all_jaspar_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

exp=read.table("../all_cisbp_filtered.txt"); rownames(exp)=exp[,1]
mat=read.table("../cisbp.txt")
read_all_cisbp_filtered()
best=best_mat4exp(roc,exp, mat)
write.table(best,"../results/best4exp_all_cisbp_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)

exp=read.table("../all_filtered.txt"); rownames(exp)=exp[,1]
mat=read_mat_all(); 
roc=          read.table("../filtered/all_hocomoco.txt")
roc=cbind(roc,read.table("../filtered/all_jaspar.txt"))
roc=cbind(roc,read.table("../filtered/all_cisbp.txt"))
best=best_mat4exp(roc,exp, mat)
write.table(best,"../results/best4exp_all_all_filtered.txt",
sep="\t", quote=F,col.names=F,row.names=F)


