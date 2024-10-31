#load pockage
library(Seurat)

suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(DoubletFinder))

#Integrate Data
ctrl<-readRDS("ctrl_doublet_mt10%.rds")
scop <- readRDS("scop_doublet_mt10%.rds")
tam <- readRDS("tam_doublet_mt10%.rds")
ctrl_singlet <- subset(ctrl, subset = doublet == "Singlet")
scop_singlet <- subset(scop, subset = doublet == "Singlet")
tam_singlet <- subset(tam, subset = doublet == "Singlet")
ctrl_singlet$Sample <- "Ctrl"
scop_singlet$Sample <- "SCOP"
tam_singlet$Sample <- "TAM"
object.list <- c(ctrl_singlet,scop_singlet,tam_singlet)
for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- FindNeighbors(combined, dims = 1:13)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunTSNE(combined,dims = 1:13)

#select cd45+cell
FeaturePlot(combined,features = c("Ptprc"),min.cutoff = 'q9', order = T, reduction = "tsne")
Idents(combined)<-combined$seurat_clusters
VlnPlot(combined,features = "Ptprc")
cd45_com<-subset(combined,idents=c("7","11","13"))
#Integrate again
Idents(cd45_com)<-cd45_com$Sample
cd45_con<-subset(cd45_com,idents="Ctrl")
cd45_scop<-subset(cd45_com,idents="SCOP")
cd45_tam<-subset(cd45_com,idents="TAM")
object.list <- c(cd45_con,cd45_scop,cd45_tam)

for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], nfeatures = 2000, verbose = FALSE)
}

anchors <- FindIntegrationAnchors(object.list = object.list,anchor.features = 2000)

cd45_com <- IntegrateData(anchorset = anchors)
dim(cd45_com)
cd45_com <- ScaleData(cd45_com)
cd45_com <- RunPCA(cd45_com, features = VariableFeatures(cd45_com))
ElbowPlot(cd45_com, ndims = 20)
DimHeatmap(cd45_com, dims = 1:20, cells = 500, balanced = TRUE)
cd45_com <- FindNeighbors(cd45_com, dims = 1:15)
cd45_com <- FindClusters(cd45_com, resolution = 0.5)
cd45_com <- RunTSNE(cd45_com,dims = 1:15)
#delete doublet 
cd45_com_del<-subset(cd45_com,idents=c("0","3","7"),invert=TRUE)
cd45_com_del <- FindNeighbors(cd45_com_del, dims = 1:15)
cd45_com_del <- FindClusters(cd45_com_del, resolution = 0.3)
cd45_com_del <- RunTSNE(cd45_com_del,dims = 1:15)
cd45_com<-cd45_com_del
DefaultAssay(cd45_com)<-"RNA"
feature<-c("Themis","Cd3e","Il7r",
           "Klra8","Klrc2","Gzma",
           "Ebf1","Ighm","Ms4a1",
           "Jchain","Igkc","Igha",
           "Wdfy4","Cbfa2t3","Fn1","Plxnc1","Samsn1","Relb",
           "Ms4a7","Tgfbr1","Hpgds","Cd163")
DotPlot(cd45_com,features = feature)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,face = "italic"))+
  scale_color_gradientn(colours = colorRampPalette(colors = c("#FEC173","#A8DBA8",'#DF044E'))(100))
Idents(cd45_com)<-cd45_com$seurat_clusters
cd45_com<-RenameIdents(cd45_com,"0"="Macrophages")
cd45_com<-RenameIdents(cd45_com,"1"="T cells")
cd45_com<-RenameIdents(cd45_com,"2"="Monocyte/DCs")
cd45_com<-RenameIdents(cd45_com,"3"="Macrophages")
cd45_com<-RenameIdents(cd45_com,"4"="Plasma")
cd45_com<-RenameIdents(cd45_com,"5"="B cells")
cd45_com<-RenameIdents(cd45_com,"6"="NK cells")
cd45_com<-RenameIdents(cd45_com,"7"="T cells")
cd45_com<-RenameIdents(cd45_com,"8"="T cells")
cd45_com<-RenameIdents(cd45_com,"9"="Monocyte/DCs")
Idents(cd45_com)<-factor(Idents(cd45_com),levels = c("T cells","NK cells","B cells","Plasma",
                                                     "Monocyte/DCs","Macrophages"))

p1<-DimPlot(cd45_com,label = T,label.size = 6,cols = c("#D53E4F","#FDAE61","#9E0142","#66C2A5","#F46D43","#3288BD"))
pdf("cd45_com_umap.pdf", width = 8, height = 6)
p1
dev.off()