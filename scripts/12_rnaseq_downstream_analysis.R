# RNA-sequencing differential expression analysis

library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(ggpubr)
library(rstatix)
library(enrichplot)
library(patchwork)


#--------------------------------------------------------------------------------
# Part 1. Exploratory data analysis
#--------------------------------------------------------------------------------
# step 1. create the DESeqDataSet object
# modifiy the column name of count.txt
count_data <- read.table("counts_modified.txt",header = TRUE, row.names = 1)
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

sample_info$group <- factor(paste(sample_info$condition, sample_info$genotype, sep = "_"),
                            levels = c("Control_WT", "Control_DKO", "Case_WT", "Case_DKO"))
sample_info$group


# transfer the count_data to matrix format, it can use as countData for DESeqDataSetFromMatrix()
count_matrix <- as.matrix(count_data)
count_matrix # Gene_ID do not need for matrix format, and for the downstream analysis

# step 2. read counts_modified.txt

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ group
)

dds
# Step 3. execute DESeq() to get the analysis result, use the same variable name dds
dds <- DESeq(dds)
# Step 4. use rlog to transform the dataset (rlog is more robust than vsd)
rlog <- rlog(dds, blind = TRUE) # blind = TRUE means that do not consider the condition 
rlog
# Step 5. check the PCA
# draw PCA plot
pca_data <- plotPCA(rlog,intgroup = c("condition","genotype"),returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$hjust <- ifelse(pca_data$PC1 > median(pca_data$PC1), -0.2, 1.2)
pca_data$vjust <- ifelse(pca_data$PC2 > median(pca_data$PC2), -0.2, 1.2)
ggplot(pca_data, aes(x = PC1, y = PC2,
                     color = condition,
                     shape = genotype)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = name, hjust = hjust,
                      vjust = vjust),size = 3,
                      check_overlap = FALSE) +
             
  xlab(paste0("PC1:", percent_var[1], "% variance"))+
  ylab(paste0("PC2:", percent_var[2],"% variance"))+
  ggtitle("PCA Plot with sample labels")+
  theme_bw()+
  theme(legend.position = "right")
#--------------------------------------------------------------------------------------
# Part 2. get the DESeq2 results of each comparison
#---------------------------------------------------------------------------------------
resultsNames(dds)
# Case vs Control in WT
res1 <- results(dds,contrast = c("group","Case_WT","Control_WT"))
# DKO vs WT in Control
res2 <- results(dds, contrast = c("group","Control_DKO","Control_WT"))
# Case vs Control in DKO
res3 <- results(dds, contrast = c("group","Case_DKO","Control_DKO"))
# DKO vs WT in Case
res4 <- results(dds,contrast = c("group","Case_DKO","Case_WT"))
#--------------------------------
# volcano_plot function
#--------------------------------
volcano_plot<- function(result,title,p_value=0.05,threshold=2,l2f_cutoff = 5, p_cutoff = 150){
  # do some statistic(sum of the up and down genes)
  sum_up <- sum(result$padj < p_value & result$log2FoldChange > threshold, na.rm = TRUE)
  sum_down <- sum(result$padj < p_value & result$log2FoldChange < -threshold, na.rm = TRUE)
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
  # add specific genes and gene name
  idx <- which(abs(result$log2FoldChange) > l2f_cutoff & -log10(result$padj) > p_cutoff)
  # check whether have the specific gene
  if(length(idx)>0){
    symbol_idx <- mapIds(org.Mm.eg.db,
                         keys = rownames(result)[idx],
                         column = "SYMBOL",
                         keytype = "ENSEMBL")
    points(result$log2FoldChange[idx], -log10(result$padj[idx]), pch = 17, col = "darkgreen",cex = 0.4)
    text(result$log2FoldChange[idx], -log10(result$padj[idx]), 
         labels = symbol_idx, pos = 3, cex = 0.4,offset = 0.2)
  }else{}
  # geneID transfer
  
  # add threshold line
  abline(v = c(-threshold, threshold), h = -log10(p_value), lty = 2, col = "gray")
 
} 

#generate the 4 volcano plot
layout(matrix(c(1,2,3,4),nrow = 2, byrow = TRUE))
p1 <- volcano_plot(result = res1, title = "Case vs Control in WT")
p2 <- volcano_plot(result = res2, title = "DKO vs WT in Control")
p3 <- volcano_plot(result = res3, title = "Case vs Control in DKO")
p4 <- volcano_plot(result = res4, title = "DKO vs WT in Case")


#----------------------------------------------------------------------------------------------------------
#  GO analysis function 
#----------------------------------------------------------------------------------------------------------
go_enrich_analysis <- function(result,dds,lfc_threshold = 1){
  # prepare the gene list
  de_genes <- subset(result,padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > lfc_threshold)
  de_ensembl <- rownames(de_genes)
  all_ensembl <- rownames(dds)
  # run GO enrich analysis
  GO_result <- enrichGO(gene = de_ensembl,
                        universe = all_ensembl,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENSEMBL",
                        ont = "BP",
                        readable = TRUE)
  return(GO_result)
}



go_result_4 <- go_enrich_analysis(res4, dds)


# enrich result visualization
# barplot function
bar_plot <- function(go_result){
  enrich_p1 <- mutate(go_result, qscore = -log(p.adjust, base=10)) |>
  barplot(x = "qscore", 
          title = "GO Enrichment",
          showCategory = 10)
  enrich_p2 <- mutate(go_result, qscore = -log(p.adjust, base=10)) |>
  barplot(x = "Count", 
          title = "GO Enrichment",
          showCategory = 10)
  enrich_p1 + enrich_p2
}
bar_plot(go_result_4)

# cnet plot function
cnet_plot <- function(result,go_result,threshold = 1){
  de_genes <- subset(result,padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > threshold)
  de_genes_ID <- rownames(de_genes)
  de_genes_symbol <- mapIds(org.Mm.eg.db,
                            keys = de_genes_ID,
                            column = "SYMBOL",
                            keytype = "ENSEMBL")  
  genelist <- de_genes$log2FoldChange
  names(genelist) <- de_genes_symbol
  go_result_ID <- setReadable(go_result, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
  p_cnet <- cnetplot(go_result_ID,foldChange = genelist,showCategory = 5)
}

cnet_4 <- cnet_plot(result = res4,go_result = go_result_4)


  
  
  
#------------------------------------------------------
# Transform the data format(all data )
#-------------------------------------------------------
# get the counts after rlog
rld <- rlog(dds,blind = FALSE)
rld_counts <- assay(rld)
rld_counts
#transfer Ensembl IDto geneID 
gene_ID <- mapIds(org.Mm.eg.db,
                         keys = rownames(rld_counts),
                         column = "SYMBOL",
                         keytype = "ENSEMBL")
gene_ID
# transfer the Ensemble to synbol
rownames(rld_counts) <- gene_ID
rld_counts
# transfer the data format to long format
# transposition the matrix
expr_df <- as.data.frame(t(rld_counts))
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


 
  
 
#--------------------------------------------------------------------------------------------
# Visualization of the target genes # get the genes for specific comparison
#--------------------------------------------------------------------------------------------
# draw box plot

# prepare the data
selected_data <- function(condition_para,genotype_para,gene_para){
  selected_data <- subset(expr_long,
                          condition == condition_para &
                            genotype %in% genotype_para &
                            gene %in% gene_para)
  return(selected_data)
}

select_data <- selected_data(condition_para = "Case", genotype_para = c("DKO","WT"), gene_para = c("Gm7449","Becn2","Uxs1"))


# draw the plot

plot_gene_stats <- function(gene_name,select_data) {
  # abstract data
  data_gene <- subset(select_data, gene == gene_name)
  y_min <- min(data_gene$expression)
  y_max <- max(data_gene$expression)
  # plot
  boxplot(expression ~ genotype, data = data_gene, 
          main = paste(gene_name, "Expression"),
          col = c("skyblue", "lightgreen"),
          ylim = c(y_min, y_max * 1.2))
  
  # p value
  p_val <- t.test(expression ~ genotype, data = data_gene)$p.value
  
  # significance
  sig <- if(p_val < 0.001) "***" else
    if(p_val < 0.01) "**" else
      if(p_val < 0.05) "*" else "ns"
  
  
  y_max <- max(data_gene$expression)
  lines(c(1, 2), rep(y_max * 1.05, 2), col = "gray")
  
  # add p value and significance
  text(1.5, y_max * 1.04, paste("p =", round(p_val, 4)), col = "red",cex = 0.5)
  text(1.5, y_max * 1.06, sig, col = "red", cex = 0.5)
}
plot_gene_stats("Uxs1",select_data)
#---------------------------------------------------------------------------------------------------    
# draw point plot


    
#-----------------------------------------------------------------------------------------------------    

              
    
    
  

  








