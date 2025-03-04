# codes for 2-step deconvolution
#library(devtools)
#install_github("Danko-Lab/BayesPrism/BayesPrism")
library(Seurat)
library(GSVA)
library(BayesPrism)

# step1: bayesprism
bk.dat <- readRDS('../count_combat_9.rds')
bk.dat <- t(bk.dat)
sc <- readRDS('../sc/sc_wang_240712.rds')
sc.dat <- t(sc@assays$RNA@counts)
sc.dat <- as.matrix(sc.dat)
meta <- sc@meta.data

sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                  input.type="count.matrix",
                                  species="mm", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","chrX","chrY") ,
                                  exp.cells=5)

sc.dat.filtered.sig<-sc.dat.filtered[colnames(sc),]
myPrism <- new.prism(
  reference=sc.dat.filtered.sig, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = meta1$cell.type.label, 
  cell.state.labels = meta1$cell.state.label,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res <- run.prism(prism = myPrism, n.cores=20)
# get main celltype ratio
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
# get main celltype exp
Z <- get.exp(bp=bp.res,state.or.type="type")

# step2:gsva
# get the markers from scrnaseq
sc.markers=read.csv('second.deconv.markers.for.t.mye.csv') 
sc.set <- split(sc.markers$gene,sc.markers$cluster)
# calc ssgsea in bulk samples
gsva_bulk_ct <- gsva(expr = as.matrix(bk.data) , gset.idx.list = sc.set,parallel.sz=10,method='ssgsea')
# calc ssgsea in sc samples
gsva_sc_ct <- gsva(expr = as.matrix(sc@assays$RNA@data) , gset.idx.list = sc.set,parallel.sz=10,method='ssgsea')
gsva_sc <- aggregate(gsva_sc_ct,list(sc$subtype),mean)
rownames(gsva_sc)=gsva_sc$Group.1
# divide by mean activity of sc
gsva_bulk <- gsva_bulk_ct/diag(as.matrix(gsva_sc[rownames(gsva_bulk_ct),rownames(gsva_bulk_ct)]))
gsva_bulk <- as.data.frame(t(gsva_bulk))
gsva_bulk$myeloid <- rowSums(gsva_bulk[,c("monocyte","macrophage","neutrophol")])
gsva_bulk$tcell <- rowSums(gsva_bulk[,c("Cd4+ T cell","Cd8+ T cell")])


theta$Cd4T <- theta$`T cell`*gsva_bulk$`Cd4+ T cell`/gsva_bulk$tcell
theta$Cd8T <- theta$`T cell`*gsva_bulk$`Cd8+ T cell`/gsva_bulk$tcell
theta$Neu <- theta$Myeloid*gsva_bulk$neutrophol/gsva_bulk$myeloid
theta$Mono <- theta$Myeloid*gsva_bulk$monocyte/gsva_bulk$myeloid
theta$Macro<- theta$Myeloid*gsva_bulk$macrophage/gsva_bulk$myeloid


