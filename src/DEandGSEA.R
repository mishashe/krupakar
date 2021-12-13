# export OPENBLAS_NUM_THREADS=1
# export GOTO_NUM_THREADS=1
# export OMP_NUM_THREADS=1


# For github
# cd /home/m.sheinman/Development/krupakar
# git add ./src/*.*
# git commit --all -m "added flexGSEA"
# git push -u origin --all


setwd("/home/m.sheinman/Development/krupakar")
library(readxl)
library(readr)
library(edgeR)
library(flexgsea)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
require(doParallel)
registerDoParallel(40)
library(foreach)
library(ComplexHeatmap)
library(globaltest)
############################################## import clinical data ###########################################################################################
SamplesInfo <- as.data.frame(read_tsv("./data/sample_annotation.tsv", col_names = TRUE,show_col_types = FALSE))
SamplesInfo <- SamplesInfo[SamplesInfo$material_type=="RNA",]
SamplesInfo <- SamplesInfo[SamplesInfo$event_type=="PRI",]
SamplesInfo <- SamplesInfo[SamplesInfo$experimental_technique=="RNA Sequencing" | SamplesInfo$sample_used_for=="RNA Sequencing",]

SamplesInfo$experimental_technique_column <- 0
SamplesInfo <- unique(SamplesInfo)
SamplesInfo <- SamplesInfo[!is.na(SamplesInfo$sample_id),]
rownames(SamplesInfo) <- SamplesInfo$sample_id

############################################## import PCs ###########################################################################################
PCs <- as.data.frame(read_excel("./data/Scores_PCA_meanspectra_NKI.xlsx"))
colnames(PCs) <- paste0("PC_",colnames(PCs))
colnames(PCs)[1] <- "patient"
PCs$patient <- paste0("PRE_NKI",sapply(1:length(PCs$patient),function(i){strsplit(PCs$patient[i],"_")[[1]][3]}))

############################################## import expression levels ###########################################################################################
xNKI1 <- as.matrix(read.table("./data/gene-expression-nki.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("./data/excluded_samples_nki.tsv",show_col_types = FALSE))
xNKI1 <- xNKI1[,!(colnames(xNKI1) %in% Excluded$precision_sample)]; 
xNKI1 <- xNKI1[,(colnames(xNKI1) %in% rownames(SamplesInfo))]; 
SamplesInfoNKI1 <- SamplesInfo[colnames(xNKI1),]
SamplesInfoNKI1$Set <- "NKI1"
xNKI1 <- xNKI1[which(rowMeans(xNKI1)>1),]

xNKI2 <- as.matrix(read.table("./data/gene-expression-nki2.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("./data/excluded_samples_nki2.tsv",show_col_types = FALSE))
xNKI2 <- xNKI2[,!(colnames(xNKI2) %in% Excluded$precision_sample)]
xNKI2 <- xNKI2[,(colnames(xNKI2) %in% rownames(SamplesInfo))]; 
SamplesInfoNKI2 <- SamplesInfo[colnames(xNKI2),]
SamplesInfoNKI2$Set <- "NKI2"
xNKI2 <- xNKI2[which(rowMeans(xNKI2)>1),]

############################################## combine datasets ############################################## 
Genes <- intersect(rownames(xNKI1),rownames(xNKI2))
x <- cbind(xNKI1[Genes,],xNKI2[Genes,])
SamplesInfo <- rbind(SamplesInfoNKI1,SamplesInfoNKI2)

############################################## keep samples with PCs ############################################## 
x <- x[,SamplesInfo$site_accession %in% PCs$patient]
SamplesInfo <- SamplesInfo[SamplesInfo$site_accession %in% PCs$patient,]
PCs <- PCs[PCs$patient %in% SamplesInfo$site_accession,]
rownames(PCs) <- PCs$patient
PCs <- PCs[SamplesInfo$site_accession,]

############################################## change genes names ############################################## 
GenesEnsembl <- sapply(rownames(x),function(g){strsplit(g,split="[.]")[[1]][1]})
GenesNames <- ensembldb::select(EnsDb.Hsapiens.v86, key=GenesEnsembl,columns=c("SYMBOL"),keytype="GENEID")
biotypes <- ensembldb::select(EnsDb.Hsapiens.v86, key=GenesEnsembl,columns=c("TXBIOTYPE"),keytype="GENEID")

rownames(GenesNames) <- GenesNames[,1]

# biotypes <- biotypes[biotypes$TXBIOTYPE=="protein_coding",]
# rownames(biotypes) <- biotypes[,1]

rownames(x) <- sapply(1:nrow(x),function(i){strsplit(rownames(x)[i],"\\.")[[1]][1]})

# x <- x[rownames(x) %in% rownames(biotypes),]
GenesNames <- GenesNames[!duplicated(GenesNames[,2]),]
x <- x[rownames(x) %in% rownames(GenesNames),]
GenesNames[duplicated(GenesNames[,2]),2] <- paste0(GenesNames[duplicated(GenesNames[,2]),2],"_2")

rownames(x) <- GenesNames[rownames(x),2]



############################################## limma ############################################## 
pdf(paste0("./plots/pValHist.pdf"))
for (iPC in 1:25)
{
  print(iPC)
  PC <- PCs[,paste0("PC_",iPC)]
  Batch <- SamplesInfo$Set
  mm <- model.matrix(~ PC + Batch)
  pdf(paste0("./plots/voom.pdf"))
  xx <- voom(x, mm, plot = T)
  dev.off()
  dat <- cbind(t(xx$E),data.frame(Set=SamplesInfo$Set))
  
  fit <- lmFit(xx, mm)
  tmp <- contrasts.fit(fit, coef = 2) # test "HR" coefficient
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  top.table$ID <- rownames(top.table)
  top.table$cor <- cor(t(xx$E[top.table$ID,]),PC)
  
  write.table(top.table,
              file = paste0("./plots/limma_PC",iPC,".csv"),
              append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
  gt_i <- gt(PCs[,paste0("PC_",iPC)], t(xx$E), null = ~ Set, data=dat)
  hist(top.table$P.Value,plot=TRUE,main=paste0("GT p-val = ",summary(gt_i)[1]))
}
dev.off()

##############################################  flexGSEA ############################################## 
flexgsea_limmaContrast <- function (x,y, abs=T) {
  fit <- lmFit(x, y)
  tmp <- contrasts.fit(fit, coef = 2) # test "HR" coefficient
  tmp <- eBayes(tmp)
  t_stat <- tmp$t
  if (abs) {
    abs(t_stat)
  } else{
    t_stat
  }
}
for (iPC in 1:25)
{
  print(iPC)
  PC <- PCs[,paste0("PC_",iPC)]
  Batch <- SamplesInfo$Set
  mm <- model.matrix(~ PC + Batch)
  pdf(paste0("./plots/voom.pdf"))
  xx <- voom(x, mm, plot = T)
  dev.off()
  
  # flexgsea
  GMT <- "h.all.v7.4.symbols.gmt"
  SetFile <- flexgsea::read_gmt(file=paste0("./data/",GMT))
  GeneSets <- flexgsea(x=xx,y=mm, 
                       gene.score.fn = flexgsea_limmaContrast,
                       es.fn = flexgsea_weighted_ks, 
                       sig.fun = flexgsea_calc_sig,
                       gene.names = NULL, 
                       nperm = 40*100, gs.size.min = 2,
                       gene.sets=SetFile,
                       gs.size.max = 5000, verbose = TRUE, block.size = 1,
                       parallel = TRUE, abs = FALSE, return_values = character())
  GeneSets <- GeneSets$table$PC
  GeneSets$p.adjust <- p.adjust(GeneSets$p, method = "fdr")
  write.table(GeneSets[order(GeneSets[,4],decreasing=FALSE),],
              file = paste0("./plots/flexGSEA_PC",iPC,".csv"),
              append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
}


GeneSets <- as.data.frame(read_tsv(paste0("./plots/flexGSEA_PC1.csv"), col_names = TRUE,show_col_types = FALSE))
SetsMatrix <- matrix(0,nrow=nrow(GeneSets),ncol=25)
rownames(SetsMatrix) <- GeneSets$GeneSet
colnames(SetsMatrix) <- paste0("PC",1:ncol(SetsMatrix))
for (iPC in 1:ncol(SetsMatrix))
{
  GeneSets <- as.data.frame(read_tsv(paste0("./plots/flexGSEA_PC",iPC,".csv"), col_names = TRUE,show_col_types = FALSE))
  rownames(GeneSets) <- GeneSets$GeneSet
  SetsMatrix[,iPC] <- -log10(GeneSets[rownames(SetsMatrix),]$fdr)
}

library(circlize)
col_fun = colorRamp2(c(min(SetsMatrix), max(SetsMatrix)), c("white", "red"))
pdf(paste0("./plots/",GMT,".pdf"),height=8)
p <- Heatmap(SetsMatrix, col = col_fun,
             column_names_gp = grid::gpar(fontsize = 6),
             row_names_gp = grid::gpar(fontsize = 6),
             show_row_names = T,
             show_column_names = T,
             cluster_columns = F, 
             cluster_rows = F)
print(p)
dev.off()  







