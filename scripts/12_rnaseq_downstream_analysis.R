# RNA-sequencing differential expression analysis


library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(enrichplot)
library(patchwork)


#--------------------------------------------------------------------------------
# Part 1. Exploratory data analysis
#--------------------------------------------------------------------------------
# step 1. create the DESeqDataSet object
# modifiy the column name of count.txt
count_data <- read.table("D:/Bioinformatics/rna_seq/counts_modified.txt",header = TRUE, row.names = 1)
original_colnames <- colnames(count_data)
print(original_colnames)
# use the sample name as the new column name
new_colnames <- c("SRR7821918","SRR7821919","SRR7821920",
                  "SRR7821921","SRR7821922","SRR7821923",
                  "SRR7821924","SRR7821925","SRR7821927",
                  "SRR7821937","SRR7821938","SRR7821939",
                  "SRR7821940","SRR7821941","SRR7821942")
colnames(count_data) <- new_colnames
print(colnames(count_data))
#create the colData for DESeqDataSetFromMatrix()
sample_info <- data.frame(
  condition = factor(c(rep("Case",9),rep("Control",6)),levels = c("Control","Case")),
  genotype = factor(c(rep("WT",5),rep("DKO",4),rep("WT",3),rep("DKO",3)),levels = c("WT","DKO"))
)

rownames(sample_info) <- new_colnames
sample_info$condition <- factor(sample_info$condition,levels = c("Control","Case"))
sample_info$condition
sample_info$genotype <- factor(sample_info$genotype, levels = c("WT","DKO"))
sample_info$genotype




# transfer the count_data to matrix format, it can use as countData for DESeqDataSetFromMatrix()
count_matrix <- as.matrix(count_data)
count_matrix # Gene_ID do not need for matrix format, and for the downstream analysis
table(duplicated(rownames(count_matrix)))
# step 2. read counts_modified.txt

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition * genotype
)

dds
# Step 3. execute DESeq() to get the analysis result, use the same variable name dds
dds <- DESeq(dds)
# Step 4. use rlog to transform the dataset (rlog is more robust than vsd)
rlog <- rlog(dds, blind = TRUE) # blind = TRUE means that do not consider the condition 
rlog
# Step 5. check the PCA
# draw PCA plot
width_inch <- 8
height_inch <- 6
dpi <- 600
png(filename = "PCA_plot.png",width = width_inch * dpi,
    height = height_inch * dpi,
    units = "px",res = dpi,
    pointsize = 12)
pca_data <- plotPCA(rlog,intgroup = c("condition","genotype"),returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(x = PC1, y = PC2,
                     color = condition,
                     shape = genotype)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  show.legend = FALSE
            ) +
  
  xlab(paste0("PC1:", percent_var[1], "% variance"))+
  ylab(paste0("PC2:", percent_var[2],"% variance"))+
  ggtitle("PCA Plot with sample labels")+
  theme_bw()+
  theme(legend.position = "right")
dev.off()
#--------------------------------------------------------------------------------------
# Part 2. get the DESeq2 results of target comparison
#---------------------------------------------------------------------------------------
resultsNames(dds)
# DKO vs WT in case and control
res <- results(dds,name = "conditionCase.genotypeDKO")
#--------------------------------
# volcano_plot function
#--------------------------------
volcano_plot<- function(result,title,p_value=0.01,threshold=1,l2f_cutoff = 3, p_cutoff = 3){
  # do some statistic(sum of the up and down genes)
  sum_up <- sum(result$padj < p_value & result$log2FoldChange > threshold, na.rm = TRUE)
  sum_down <- sum(result$padj < p_value & result$log2FoldChange < -threshold, na.rm = TRUE)
  # draw the volcano plot,up genes are red, down genes are blue
  plot(result$log2FoldChange, -log10(result$padj),
       pch = 16, cex = 0.5,
       col = ifelse(result$padj < p_value & result$log2FoldChange > threshold, "red",
                    ifelse(result$padj < p_value & result$log2FoldChange < -threshold, "blue", "gray")),
       xlab = "log2 Fold Change",
       ylab = "-log10(padj)"
  )
  # add title and subtitle
  title(main = title,line = 2,cex.main = 1.2)
  mtext(paste0("sum of up regulate genes:" ,sum_up,"; sum of down regulate genes:", sum_down ),
        side = 3, line = 0.4, cex = 0.6, col = "black")
  # filter significant fold change genes 
  idx <- which(abs(result$log2FoldChange) > l2f_cutoff & -log10(result$padj) > p_cutoff)
  
  
  # check if have significant fold change genes(if have, transfer the gene_id)
  if(length(idx)>0){
    symbol_idx <- mapIds(org.Mm.eg.db,
                         keys = rownames(result)[idx],
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    points(result$log2FoldChange[idx], -log10(result$padj[idx]), pch = 16, col = "darkgreen",cex = 0.4)
    text(result$log2FoldChange[idx], -log10(result$padj[idx]), 
         labels = symbol_idx, pos = 3, cex = 0.4,offset = 0.2)
  }else{}
  
  # add threshold line
  abline(v = c(-threshold, threshold), h = -log10(p_value), lty = 2, col = "gray")
  # add legend with counts on respective sides
  legend("topright",legend = paste("Up:",sum_up),col = "red",pch = 16, bty = "n",cex = 0.8)
  legend("topleft", legend = paste("Down:", sum_down), col = "blue", pch = 16, bty = "n", cex = 0.8)
  
} 

#generate the  volcano plot
width_inch <- 8
height_inch <- 6
dpi <- 600
png(filename = "volcano_plot.png",width = width_inch * dpi,
    height = height_inch * dpi,
    units = "px",res = dpi,
    pointsize = 12)
volcano_plot(result = res, title = "DKO vs WT effects compared across conditions")

dev.off()


#----------------------------------------------------------------------------------------------------------
# Part3. GO analysis function 
#----------------------------------------------------------------------------------------------------------
# Step 1. Function of GO analysis and execute function to get enrichment result
go_enrich_analysis <- function(result,lfc_threshold = 1){
  # prepare the gene list
  de_genes <- subset(result,padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > lfc_threshold)
  de_ensembl <- rownames(de_genes)
  all_genes <- subset(result,!is.na(padj))
  all_ensembl <- rownames(all_genes)
  # run GO enrich analysis
  GO_result <- enrichGO(gene = de_ensembl,
                        universe = all_ensembl,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENSEMBL",
                        ont = "BP")
  return(GO_result)
}
# execute the function and get the result of GO analysis
go_result <- go_enrich_analysis(res)




# Step 2. Enrich result visualization
# 1st plot: Dotplot
width_inch <- 8
height_inch <- 6
dpi <- 600
png(filename = "dot_plot.png",width = width_inch * dpi,
    height = height_inch * dpi,
    units = "px",res = dpi,
    pointsize = 12)
# draw dotplot 
dotplot(go_result, x = "GeneRatio",color = "p.adjust",
        showCategory = 10,title = "DKO vs WT effects compared across conditions")
dev.off()

# 2nd plot: cnet plot function and draw cnet plot
cnet_plot <- function(result,go_result,threshold = 1){
  de_genes <- subset(result,padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > threshold)
  de_genes_ID <- rownames(de_genes)
  genelist <- setNames(de_genes$log2FoldChange,de_genes_ID)
  symbol <- mapIds(org.Mm.eg.db,
                   keys = names(genelist),
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")
  genelist_symbol <- setNames(genelist,symbol)
  geneList_symbol <- na.omit(genelist_symbol)
  go_result@result <- go_result@result[order(go_result$p.adjust), ]
  go_result_ID <- setReadable(go_result,"org.Mm.eg.db","ENSEMBL")
  
  p_cnet <- cnetplot(go_result_ID,foldChange = geneList_symbol,showCategory = 10,
                     size_item = 1, color_edge = "gray",font.size = 0.5)
  return(p_cnet)
}
pdf("cnet_plot.pdf",width = 12,height = 9,
    pointsize = 14,paper = "special",
    onefile = FALSE)
cnet_plot(result = res,go_result = go_result)

dev.off()


# check the enrichment pathway of significant genes 
go_result_ID <- setReadable(go_result,'org.Mm.eg.db',"ENSEMBL")

pathway_wars1 <- subset(go_result_ID@result, grepl("Wars1",geneID))
print(pathway_wars1)
pathway_batf2 <- subset(go_result_ID@result, grepl("Batf2",geneID))
print(pathway_batf2)
pathway_Irgm1 <- subset(go_result_ID@result, grepl("Irgm1",geneID))
print(pathway_Irgm1)
pathway_Igtp <- subset(go_result_ID@result, grepl("Igtp",geneID))
print(pathway_Igtp)
pathway_Chil4 <- subset(go_result_ID@result, grepl("Chil4",geneID))
print(pathway_Chil4)


#--------------------------------------------------------------------------------------------
# Part 4. Visualization of the target genes (get the genes for specific comparison(volcano plot)) 
#--------------------------------------------------------------------------------------------
#------------------------------------------------------
# step 1. Transform the data format(all data )
#-------------------------------------------------------
# get the counts after normalized
normalized_counts <- counts(dds,normalized= TRUE)


#transfer Ensembl IDto geneID 
gene_ID <- mapIds(org.Mm.eg.db,
                  keys = rownames(normalized_counts),
                  column = "SYMBOL",
                  keytype = "ENSEMBL")
gene_ID
# transfer the Ensemble to synbol
rownames(normalized_counts) <- gene_ID
normalized_counts
# transfer the data format to long format
# transposition the matrix
expr_df <- as.data.frame(t(normalized_counts))
# add sample name to a new column
expr_df$sample <- rownames(expr_df)
# add sample information
expr_df$condition <- sample_info$condition
expr_df$genotype <- sample_info$genotype
# transfer to long format
expr_long <- melt(expr_df,
                  id.vars = c("sample","condition","genotype"),
                  variable.name = "gene",
                  value.name = "expression")





#-------------------------------------
# Step 2. draw point plot
#-------------------------------------
# prepare the data
selected_data <- function(condition_para,genotype_para,gene_para){
  selected_data <- subset(expr_long,
                          condition == condition_para &
                            genotype == genotype_para &
                            gene %in% gene_para)
  return(selected_data)
}

case_wt <- selected_data(condition_para = "Case", genotype_para = "WT", 
                         gene_para = c("Wars1","Batf2","Igtp","Irgm1"))
case_DKO <- selected_data(condition_para = "Case", genotype_para = "DKO", 
                          gene_para = c("Wars1","Batf2","Igtp","Irgm1"))
control_wt <- selected_data(condition_para = "Control", genotype_para = "WT", 
                            gene_para = c("Wars1","Batf2","Igtp","Irgm1"))
control_DKO <- selected_data(condition_para = "Control", genotype_para = "DKO", 
                             gene_para = c("Wars1","Batf2","Igtp","Irgm1"))


#---------------------------------------------------------------------------------------------------    
# draw point plot
points_plot <- function(gene_name){

gene_case_wt <- subset(case_wt, gene == gene_name)
gene_case_DKO <- subset(case_DKO, gene == gene_name)
gene_control_wt <- subset(control_wt, gene == gene_name)
gene_control_DKO <- subset(control_DKO, gene == gene_name)

plot(0,0,type = "n",
     xlim = c(0.5,4.5),ylim = range(c(gene_case_wt$expression,gene_case_DKO$expression,
                                      gene_control_wt$expression,gene_control_DKO$expression)),
     xlab = "Experimental Groups",ylab = paste(gene_name, "Expression"),
     main = "DKO vs WT in Control and Case Conditions",
     xaxt = "n", cex.lab = 1.2, cex.main = 1.3)
# add x axis
axis(1, at = 1:4, labels = c("Case\nWT","Case\nDKO",
                             "Control\nWT","Control\nDKO"), padj = 0.5)

# add points
points(jitter(rep(1,5),amount = 0.2),gene_case_wt$expression,pch = 16, col = "blue", cex = 1.2)
points(jitter(rep(2,4),amount = 0.2),gene_case_DKO$expression, pch = 16, col = "red", cex = 1.2)
points(jitter(rep(3,3),amount = 0.2),gene_control_wt$expression, pch = 16, col = "blue", cex = 1.2)
points(jitter(rep(4,3),amount = 0.2),gene_control_DKO$expression, pch = 16, col = "red", cex = 1.2)

# add legend
legend("topright", legend = c("WT","DKO"),
       col = c("blue","red"), pch = 16,
       title = "Genotype", title.font = 2,cex = 0.7,
       pt.cex = 1, inset = 0.02, bty = "n")
}
#-----------------------------------------------------------------------------------------------------    

pdf("point_plot_gene.pdf",width = 12,height = 9,
    onefile = FALSE)
layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))


points_plot("Irgm1")

points_plot("Igtp")

dev.off()

