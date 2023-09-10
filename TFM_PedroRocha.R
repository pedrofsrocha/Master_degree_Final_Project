####Load libraries####
install.packages("OlinkAnalyze")
library(OlinkAnalyze)
library(tidyverse)
library(readxl)
library(openxlsx)
library(ggrepel)
library(ggpubr)
library(stringr)
library(atable)
library(reshape2)
library(pheatmap)
library(data.table)
library(kableExtra)
library(knitr)
library(rmarkdown)
install.packages("compareGroups")
library(compareGroups)
library(limma)
library(VennDiagram)
library(ggvenn)
library(DT)
library(RColorBrewer)
library(ggplot2)
library(Hmisc)
library(atable)
library(gridExtra)
library(maftools)
library('NMF')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)

####Load clinical data####
load("clinical_data.Rdata")
head(clinical_data)
####Load Proteomics data####
load("olink_qc_data_clin.Rdata")
head(data_qc)

####Proteomics - Subset of samples from Timepoint 1.####
TP1_data <- data_qc %>% 
  filter(Timepoint == "TP1")

####Proteomics - Subset of samples from Timepoint 2####
TP2_data <- data_qc %>% 
  filter(Timepoint == "TP2")
#####Proteomics - USING LIMMA  TIMEPOINT 1 - OS_ICI_36months #####
head(TP1_data)
colnames(TP1_data)
#Create a matrix for Timepoint 1.
matrix_tp1_1 <- data.table(subset(TP1_data, select = c(1,6,13,64)))
matrix_tp1_2 <- dcast.data.table(matrix_tp1_1, formula = Assay~ SampleID + OS_ICI_36months, value.var = c('NPX'))
colnames(matrix_tp1_2) <- sub("-", ".", colnames(matrix_tp1_2))
pheno <- colnames(matrix_tp1_2)[2:length(colnames(matrix_tp1_2))] #here I get the names to create a pheno table to input in the matrix.
pheno <- as.data.frame(pheno)
#final pheno table to input in the matrix design in limma.
pheno <- separate(pheno, pheno, into = c("SampleID", "OS_ICI_36months"),"_")
#This is the final matrix:
matrix_tp1_2 <- as.data.frame(matrix_tp1_2)
rownames(matrix_tp1_2) <- matrix_tp1_2$Assay
matrix_tp1_3 <- matrix_tp1_2[,2:37] #This is the final  dataframe to be used.
head(matrix_tp1_3)


#Apply limma:
cond <- pheno$OS_ICI_36months
design <- model.matrix(~0 + cond)

design <- model.matrix(~0+cond)
rownames(design) <- pheno$OS_ICI_36months
colnames(design) <- gsub("cond", "", colnames(design))


# model fit
fit <- lmFit(matrix_tp1_3, design) 

# contrasts
contrast.matrix <- makeContrasts(con1=LTR-STR,
                                 levels = design) 

# contrasts fit and Bayesian adjustment
fit2 <- contrasts.fit(fit, contrast.matrix)
fite <- eBayes(fit2)

#number of proteins DE using adjuste p-values.
summary(decideTests(fite, method = "separate")) # zero
# in this case we cannot adjust for multiple comparisons.
summary(decideTests(fite, adjust.method = "none", method = "separate"))

top.table <- topTable(fite, number = Inf, adjust = "fdr")
head(top.table)
hist(top.table$P.Value, breaks = 100, main = "Results p-values")

#Significant genes table
signif.genes_limma <- top.table[(top.table$P.Value < 0.05), ]
signif.genes_limma <- top.table[(top.table$P.Value < 0.05  & abs(top.table$logFC)>0.5), ]
signif.genes_limma <- signif.genes_limma[1:39,]
signif.genes_limma$UP_DOWN <- ifelse(signif.genes_limma$logFC < 0, "DOWN", "UP")
head(signif.genes_limma)

#Save the table as an excel file.
signif.genes_limma$Protein <- rownames(signif.genes_limma)
write.xlsx(signif.genes_limma, file = "Comparison_NR_LTR_tp1_limma_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#HEATMAP USING LIMMA  TIMEPOINT 1 - OS_ICI_36months.
#list of genes DE between LTR and NR.
head(signif.genes_limma)
proteins_hm_tp1 <-  signif.genes_limma$Protein
#Matrix for heatmap.
hm_limma_tp1 <- matrix_tp1_3[rownames(matrix_tp1_3) %in% proteins_hm_tp1,]

#Annotations for heatmap.
#Annotations - should be a dataframe with SampleID and Response variable:
annotations <- data_qc %>% 
  filter(Timepoint == "TP1")
annotations1 <- annotations[,c(1,28,31,32,34,35,38,39,40,62,64,53,66, 44)]
annotations1 <- annotations1 %>% 
  distinct(SampleID, .keep_all = TRUE)
annotations1$SampleID <- gsub("-", ".", annotations1$SampleID)
rownames(annotations1) <- annotations1$SampleID
head(annotations1)
annotations2 <- annotations1[,c(11,13,10,12,14,8,9,5,3,2)]

colnames(hm_limma_tp1) <- rownames(annotations2)
#Using correlation as a cluster method.
pheatmap(hm_limma_tp1,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         border_color = "white",
         annotation_col = annotations2
)

#Using Euclidean as a cluster method.
pheatmap(hm_limma_tp1,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         border_color = "beige",
         annotation_col = annotations2
)

#VOLCANO PLOT.
volcano_data <- top.table
volcano_data$Protein <- rownames(volcano_data)
#Add a new column for the proteins with pvalue<0.05 and FC 0.5.
volcano_data$Different_prot <- "No"
# if estimate > 0.5 and pvalue < 0.05, set as "UP" 
volcano_data$Different_prot[volcano_data$logFC >0.5 & volcano_data$P.Value < 0.05] <- "Up"
# if estimate < -0.5 and pvalue < 0.05, set as "DOWN"
volcano_data$Different_prot[volcano_data$logFC < -0.5 & volcano_data$P.Value < 0.05] <- "Down"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
volcano_data$labels <- NA
volcano_data$labels[volcano_data$Different_prot != "No"] <- volcano_data$Protein[volcano_data$Different_prot != "No"]

#Volcano plot.
ggplot(volcano_data, aes(x=logFC, y=-log10(P.Value),col=Different_prot, label = labels)) +
  geom_point() +
  theme_classic2() +
  geom_vline(xintercept=c(-0.5, 0.5), col="gray40",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="gray40",linetype="dashed") +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  scale_color_manual(values=c("cornflowerblue", "gray90", "salmon")) +
  xlim(-1.7,1.7) +
  labs(title = "Serum baseline comparison LTR vs STR")

#Save as an excel file the data from the volcanoplot TP1.
write.xlsx(volcano_data, file = "Comparison_NR_LTR_tp1_limma_OS_ICI_36months_volcano_data.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#Boxplots - LTR Vs STR in TP1.
#Do a for loop to print all the plots at once.
proteins_hm_tp1

for (protein in proteins_hm_tp1) {
  # Filter the data for the current gene and select only the necessary columns
  protein_subset <- data_qc %>%
    filter(OS_ICI_36months=="LTR" | OS_ICI_36months == "STR") %>% 
    filter(Timepoint == "TP1") %>% 
    filter(Assay %in% protein)
  
  # Create the boxplot using ggplot
  p <- ggplot(data = protein_subset, aes(x = OS_ICI_36months, y = NPX)) +
    geom_boxplot(aes(fill=OS_ICI_36months)) +
    scale_fill_manual(values=c("#30b9f0", "#d6a02b"))+
    #stat_compare_means(paired=FALSE,label.y = -3.5,size=4) +
    theme_classic()+
    labs(title = paste("Boxplot for", protein),
         x = "",
         y = "NPX")
  
  # Save the plot as a PDF file
  pdf_filename <- paste("Timepoint 1 - Boxplot_for_", protein, ".pdf", sep = "")
  ggsave(filename = pdf_filename, plot = p)
}

#MSigDb analysis

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(hallmark_genes)

#####Proteomics - USING LIMMA  TIMEPOINT 2 - OS_ICI_36months #####
head(TP2_data)
colnames(TP2_data)
#Create a matrix for Timepoint 1.
# matrix_tp1 <- TP1_data %>% 
#   filter(Exceptional_responders_cufoff_TTF_24months =="NR" | Exceptional_responders_cufoff_TTF_24months == "LTR")
matrix_tp2_1 <- data.table(subset(TP2_data, select = c(1,6,13,64)))
matrix_tp2_2 <- dcast.data.table(matrix_tp2_1, formula = Assay~ SampleID + OS_ICI_36months, value.var = c('NPX'))
colnames(matrix_tp2_2) <- sub("-", ".", colnames(matrix_tp2_2))
pheno_tp2 <- colnames(matrix_tp2_2)[2:length(colnames(matrix_tp2_2))] #here I get the names to create a pheno table to input in the matrix.
pheno_tp2 <- as.data.frame(pheno_tp2)
#final pheno table to input in the matrix design in limma.
pheno_tp2 <- separate(pheno_tp2, pheno_tp2, into = c("SampleID", "OS_ICI_36months"),"_")
#This is the final matrix:
matrix_tp2_2 <- as.data.frame(matrix_tp2_2)
rownames(matrix_tp2_2) <- matrix_tp2_2$Assay
matrix_tp2_3 <- matrix_tp2_2[,2:38] #This is the final  dataframe to be used.
head(matrix_tp2_3)


#Apply limma:
cond_tp2 <- pheno_tp2$OS_ICI_36months
design_tp2 <- model.matrix(~0 + cond_tp2)

design_tp2 <- model.matrix(~0+cond_tp2)
rownames(design_tp2) <- pheno_tp2$OS_ICI_36months
colnames(design_tp2) <- gsub("cond_tp2", "", colnames(design_tp2))


# model fit
fit_tp2 <- lmFit(matrix_tp2_3, design_tp2) 

# contrasts
contrast.matrix_tp2 <- makeContrasts(con1=LTR-STR,
                                     levels = design_tp2) 

# contrasts fit and Bayesian adjustment
fit2_tp2 <- contrasts.fit(fit_tp2, contrast.matrix_tp2)
fite_tp2 <- eBayes(fit2_tp2)

#number of proteins DE using adjuste p-values.
summary(decideTests(fite_tp2, method = "separate")) # zero
# in this case we cannot adjust for multiple comparisons.
summary(decideTests(fite_tp2, adjust.method = "none", method = "separate"))

top.table_tp2 <- topTable(fite_tp2, number = Inf, adjust = "fdr")
head(top.table_tp2)
hist(top.table_tp2$P.Value, breaks = 100, main = "Results p-values")

#Significant genes table
signif.genes_limma_tp2 <- top.table_tp2[(top.table_tp2$P.Value < 0.05 & abs(top.table_tp2$logFC)>0.5), ]
signif.genes_limma_tp2 <- signif.genes_limma_tp2[1:8,]
signif.genes_limma_tp2$UP_DOWN <- ifelse(signif.genes_limma_tp2$logFC < 0, "DOWN", "UP")
head(signif.genes_limma_tp2)
dim(signif.genes_limma_tp2)
range(signif.genes_limma_tp2$logFC)

#Save the table as an excel file.
signif.genes_limma_tp2$Protein <- rownames(signif.genes_limma_tp2)
write.xlsx(signif.genes_limma_tp2, file = "Comparison_NR_LTR_tp2_limma_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)


#HEATMAP USING LIMMA  TIMEPOINT 2 - OS_ICI_36months.
#list of genes DE between LTR and NR.
head(signif.genes_limma_tp2)
proteins_hm_tp2 <-  signif.genes_limma_tp2$Protein
#Matrix for heatmap.
hm_limma_tp2 <- matrix_tp2_3[rownames(matrix_tp2_3) %in% proteins_hm_tp2,]
colnames(hm_limma_tp2) <- gsub("_LTR|_STR","", colnames(hm_limma_tp2))
#Annotations for heatmap.
#Annotations - should be a dataframe with SampleID and Response variable:
annotations_tp2 <- data_qc %>% 
  filter(Timepoint == "TP2")
annotations1_tp2 <- annotations_tp2[,c(1,28,31,32,34,35,37,38,39,40,62,64)]
annotations1_tp2 <- annotations1_tp2 %>% 
  distinct(SampleID, .keep_all = TRUE)
annotations1_tp2$SampleID <- gsub("-", ".", annotations1_tp2$SampleID)
rownames(annotations1_tp2) <- annotations1_tp2$SampleID
head(annotations1_tp2)
annotations2_tp2 <- annotations1_tp2[,c(5,11,12)]

#Using correlation as a cluster method.
pheatmap(hm_limma_tp2,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         border_color = "beige",
         annotation_col = annotations2_tp2
)

#Using Euclidean as a cluster method.
pheatmap(hm_limma_tp2,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         border_color = "beige",
         annotation_col = annotations2_tp2
)

#VOLCANO PLOT.
volcano_data_tp2 <- top.table_tp2
volcano_data_tp2$Protein <- rownames(volcano_data_tp2)
#Add a new column for the proteins with pvalue<0.05 and FC 0.5.
volcano_data_tp2$Different_prot <- "No"
# if estimate > 0.5 and pvalue < 0.05, set as "UP" 
volcano_data_tp2$Different_prot[volcano_data_tp2$logFC >0.5 & volcano_data_tp2$P.Value < 0.05] <- "Up"
# if estimate < -0.5 and pvalue < 0.05, set as "DOWN"
volcano_data_tp2$Different_prot[volcano_data_tp2$logFC < -0.5 & volcano_data_tp2$P.Value < 0.05] <- "Down"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
volcano_data_tp2$labels <- NA
volcano_data_tp2$labels[volcano_data_tp2$Different_prot != "No"] <- volcano_data_tp2$Protein[volcano_data_tp2$Different_prot != "No"]

#Volcano plot.
ggplot(volcano_data_tp2, aes(x=logFC, y=-log10(P.Value),col=Different_prot, label = labels)) +
  geom_point() +
  theme_classic2() +
  geom_vline(xintercept=c(-0.5, 0.5), col="gray40",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="gray40",linetype="dashed") +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  scale_color_manual(values=c("cornflowerblue", "gray90", "salmon")) +
  xlim(-1.7,1.7) +
  labs(title = "Timepoint 2 comparison LTR vs STR")


#Boxplots - LTR Vs STR in TP2.
#Do a for loop to print all the plots at once.
proteins_hm_tp2

for (protein in proteins_hm_tp2) {
  # Filter the data for the current gene and select only the necessary columns
  protein_subset <- data_qc %>%
    filter(OS_ICI_36months=="LTR" | OS_ICI_36months == "STR") %>% 
    filter(Timepoint == "TP2") %>% 
    filter(Assay %in% protein)
  
  # Create the boxplot using ggplot
  p <- ggplot(data = protein_subset, aes(x = OS_ICI_36months, y = NPX)) +
    geom_boxplot(aes(fill=OS_ICI_36months)) +
    scale_fill_manual(values=c("#30b9f0", "#d6a02b"))+
    #stat_compare_means(paired=FALSE,label.y = -3.5,size=4) +
    theme_classic()+
    labs(title = paste("Boxplot for", protein),
         x = "",
         y = "NPX")
  
  # Save the plot as a PDF file
  pdf_filename <- paste("Timepoint 2 - Boxplot_for_", protein, ".pdf", sep = "")
  ggsave(filename = pdf_filename, plot = p)
}


#####Proteomics - USING LIMMA  TIMEPOINT 2 - TIMEPOINT 1 - OS_ICI_36months     #####
TP2_TP1_data <- data_qc %>% 
  filter(Timepoint!="TP3")

head(TP2_TP1_data)
dim(TP2_TP1_data)
colnames(TP2_TP1_data)
#Transform the data frame in order to subtract TP2-TP1:
TP_dif <- subset(TP2_TP1_data, select = c(1,3,4,6,13,17,19,61,62,64,65,66))
TP_dif_dt <- data.table(TP_dif) #convert to a data.table.
#Unmelt the table and create two new columns for TP1 and TP2 for each ID and for each Assay:
TP_dif2 <- dcast.data.table(TP_dif_dt, formula = ID + Assay ~ Timepoint, value.var = c('NPX'))#Create a table that allows the computation of TP2-TP1
TP_dif2[,'TP2_TP1':=TP2-TP1] #Create a new column that computes TP2-TP1.

#Here I am repllicating the previous code but adding the OS_ICI_36months variable.
TP_dif3 <- dcast.data.table(TP_dif_dt, formula = ID + Assay + OS_ICI_36months + OlinkID ~ Timepoint, value.var = c('NPX'))#Create a table that allows the computation of TP2-TP1
TP_dif3[,'TP2_TP1':=TP2-TP1] #Create a new column that computes TP2-TP1.

head(TP_dif3)# This is the final table where we have the values for the substraction of TP2-TP1:
dim(TP_dif3)
TP_dif4 <- dcast.data.table(TP_dif3, formula = Assay ~ ID, value.var = c('TP2_TP1'))
TP_dif4  <-  as.data.frame(TP_dif4)
rownames(TP_dif4) <- TP_dif4$Assay
TP_dif4  <-TP_dif4[,2:38]
dim(TP_dif4)
head(TP_dif4)#This is the matrix to be used to calculate DEG.

#To create a pheno dataframe we will use the previous dataframe TP_dif3 and remove the duplicates by ID.
# Create a new data frame without duplicated cases based on a specific column
pheno_tp2_tp1 <- TP_dif3 %>% distinct(ID, .keep_all = TRUE)
pheno_tp2_tp1 <- pheno_tp2_tp1[,c(1,3)]


#Apply limma:
cond_tp2_tp1 <- pheno_tp2_tp1$OS_ICI_36months
design_tp2_tp1 <- model.matrix(~0 + cond_tp2_tp1)

rownames(design_tp2_tp1) <- pheno_tp2_tp1$OS_ICI_36months
colnames(design_tp2_tp1) <- gsub("cond_tp2_tp1", "", colnames(design_tp2_tp1))


# model fit
fit_tp2_tp1 <- lmFit(TP_dif4, design_tp2_tp1) 

# contrasts
contrast.matrix_tp2_tp1 <- makeContrasts(con1=LTR-STR,
                                         levels = design_tp2_tp1) 

# contrasts fit and Bayesian adjustment
fit_tp2_tp1_1 <- contrasts.fit(fit_tp2_tp1, contrast.matrix_tp2_tp1)
fite_tp2_tp1 <- eBayes(fit_tp2_tp1_1)

#number of proteins DE
summary(decideTests(fite_tp2_tp1, method = "separate")) # zero
# in this case we cannot adjust for multiple comparisons.
summary(decideTests(fite_tp2_tp1, adjust.method = "none", method = "separate"))

top.table_tp2_tp1 <- topTable(fite_tp2_tp1, number = Inf, adjust = "fdr")
head(top.table_tp2_tp1)
hist(top.table_tp2_tp1$P.Value, breaks = 100, main = "Results p-values")

#Significant genes table
signif.genes_limma_tp2_tp1 <- top.table_tp2_tp1[(top.table_tp2_tp1$P.Value < 0.05 & abs(top.table_tp2_tp1$logFC)>0.5), ]
signif.genes_limma_tp2_tp1 <- signif.genes_limma_tp2_tp1[1:12,]
head(signif.genes_limma_tp2_tp1)

#Save the table as an excel file.
signif.genes_limma_tp2_tp1$Protein <- rownames(signif.genes_limma_tp2_tp1)
write.xlsx(signif.genes_limma_tp2_tp1, file = "Comparison_NR_LTR_tp2_tp1_limma_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)


#HEATMAP USING LIMMA  TIMEPOINT 2 - TIMEPOINT 1 - OS_ICI_36months.
#list of genes DE between LTR and NR.
head(signif.genes_limma_tp2_tp1)
proteins_hm_tp2_tp1 <-  signif.genes_limma_tp2_tp1$Protein
#Matrix for heatmap.
hm_limma_tp2_tp1 <- TP_dif4[rownames(TP_dif4) %in% proteins_hm_tp2_tp1,]
#need to remove sample 56.
hm_limma_tp2_tp1 <- hm_limma_tp2_tp1[,-25]

#Annotations for heatmap.
#Annotations - should be a dataframe with SampleID and Response variable:
annotations_tp2_tp1 <- data_qc %>% 
  filter(Timepoint == "TP2")
annotations_tp2_tp1 <- annotations_tp2_tp1[,c(17,28,31,32,34,35,37,38,39,40,62,64,66)]
annotations_tp2_tp1 <- annotations_tp2_tp1 %>% 
  distinct(ID, .keep_all = TRUE)
rownames(annotations_tp2_tp1) <- annotations_tp2_tp1$ID
head(annotations_tp2_tp1)
annotations_tp2_tp1 <- annotations_tp2_tp1[,c(7,11,13,12)]

#Using correlation as a cluster method.
pheatmap(hm_limma_tp2_tp1,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         border_color = "black",
         annotation_col = annotations_tp2_tp1
)

#Using Euclidean as a cluster method.
pheatmap(hm_limma_tp2_tp1,scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "euclidean",
         clustering_method = "ward.D",
         border_color = "beige",
         annotation_col = annotations_tp2_tp1
)

#VOLCANO PLOT.
volcano_data_tp2_tp1 <- top.table_tp2_tp1
volcano_data_tp2_tp1$Protein <- rownames(volcano_data_tp2_tp1)
#Add a new column for the proteins with pvalue<0.05 and FC 0.5.
volcano_data_tp2_tp1$Different_prot <- "No"
# if estimate > 0.5 and pvalue < 0.05, set as "UP" 
volcano_data_tp2_tp1$Different_prot[volcano_data_tp2_tp1$logFC >0.5 & volcano_data_tp2_tp1$P.Value < 0.05] <- "Up"
# if estimate < -0.5 and pvalue < 0.05, set as "DOWN"
volcano_data_tp2_tp1$Different_prot[volcano_data_tp2_tp1$logFC < -0.5 & volcano_data_tp2_tp1$P.Value < 0.05] <- "Down"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
volcano_data_tp2_tp1$labels <- NA
volcano_data_tp2_tp1$labels[volcano_data_tp2_tp1$Different_prot != "No"] <- volcano_data_tp2_tp1$Protein[volcano_data_tp2_tp1$Different_prot != "No"]

#Volcano plot.
ggplot(volcano_data_tp2_tp1, aes(x=logFC, y=-log10(P.Value),col=Different_prot, label = labels)) +
  geom_point() +
  theme_classic2() +
  geom_vline(xintercept=c(-0.5, 0.5), col="gray80") +
  geom_hline(yintercept=-log10(0.05), col="gray80") +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  scale_color_manual(values=c("cornflowerblue", "gray80", "salmon")) +
  xlim(-1.7,1.7) +
  labs(title = "TP2-TP1 - LTR vs STR - OS_ICI 36months")



#Boxplots - LTR Vs STR in TP2 - TP1.
#Do a for loop to print all the plots at once.
proteins_hm_tp2_tp1

for (protein in proteins_hm_tp2_tp1) {
  # Filter the data for the current gene and select only the necessary columns
  protein_subset <- TP_dif3 %>%
    filter(OS_ICI_36months=="LTR" | OS_ICI_36months == "STR") %>% 
    filter(Assay %in% protein)
  
  # Create the boxplot using ggplot
  p <- ggplot(data = protein_subset, aes(x = OS_ICI_36months, y = TP2_TP1)) +
    geom_boxplot(aes(fill=OS_ICI_36months)) +
    scale_fill_manual(values=c("lightblue1", "salmon"))+
    stat_compare_means(paired=FALSE,label.y = -3.5,size=4) +
    theme_classic()+
    labs(title = paste("Boxplot for", protein),
         x = "Group",
         y = "Protein Expression")
  
  # Save the plot as a PDF file
  pdf_filename <- paste("Boxplot_for_", protein, ".pdf", sep = "")
  ggsave(filename = pdf_filename, plot = p)
}

####Proteomics - Limma TIMESERIES ####

cond5 <- paste0(cond,"_TP1")
cond6 <- paste0(cond_tp2, "_TP2")
cond_tutorial <- c(cond5, cond6)


design_tutorial <- model.matrix(~0+cond_tutorial)

all.equal(rownames(matrix_tp1_3), rownames(matrix_tp2_3))

matrix_tutorial <- cbind(matrix_tp1_3, matrix_tp2_3)

ID <- gsub("\\..*", "", colnames(matrix_tutorial))

cor <- duplicateCorrelation(matrix_tutorial, design_tutorial, block=ID)
cor$consensus.correlation

fit_tutorial <- lmFit(matrix_tutorial, design_tutorial, block=ID,
                      correlation=cor$consensus.correlation)

fit_tutorial$coefficients

contrasts_tutorial <- makeContrasts(
  LTR_T2vsT1= cond_tutorialLTR_TP2-cond_tutorialLTR_TP1, 
  STR_T2vsT1=cond_tutorialSTR_TP2-cond_tutorialSTR_TP1, 
  LTRvsSTR=(cond_tutorialLTR_TP2-cond_tutorialLTR_TP1)-(cond_tutorialSTR_TP2-cond_tutorialSTR_TP1), 
  levels=colnames(design_tutorial)) 

contrasts_tutorial

fit_tutorial2 <- contrasts.fit(fit_tutorial, contrasts_tutorial)

fit_tutorial2 <- eBayes(fit_tutorial2)

#number of proteins DE
summary(decideTests(fit_tutorial2, method = "separate")) # zero
# in this case we cannot adjust for multiple comparisons.
summary(decideTests(fit_tutorial2, adjust.method = "none", method = "separate", lfc = 0.585))
vennDiagram(decideTests(fit_tutorial2, adjust.method = "none", method = "separate", lfc = 0.585))

top.table_tutorial <- topTable(fit_tutorial2, number = Inf, adjust = "fdr")
top.table_tutorial1 <- topTable(fit_tutorial2, number = Inf, adjust = "fdr",coef = 1)
top.table_tutorial2 <- topTable(fit_tutorial2, number = Inf, adjust = "fdr",coef = 2)
top.table_tutorial3 <- topTable(fit_tutorial2, number = Inf, adjust = "fdr",coef = 3)
head(top.table_tutorial1)


#Significant genes table - TO CREATE A TABLE AND AN EXCEL FILE FOR ALL THE COEFFICIENTS.
signif.genes_tutorial <- top.table_tutorial[(top.table_tutorial$P.Value < 0.05), ]
signif.genes_tutorial <- signif.genes_tutorial[1:18,]
head(signif.genes_tutorial)

#Save the table as an excel file.
signif.genes_tutorial$Protein <- rownames(signif.genes_tutorial)
write.xlsx(signif.genes_tutorial, file = "Comparison_LTR_timeseries_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)


#Significant genes table - TO CREATE A TABLE AND AN EXCEL FILE FOR LTR TP2 Vs TP1.
signif.genes_tutorial1 <- top.table_tutorial1[(top.table_tutorial1$P.Value < 0.05), ]
signif.genes_tutorial1 <- signif.genes_tutorial1[1:11,]
head(signif.genes_tutorial1)

#Save the table as an excel file.
signif.genes_tutorial1$Protein <- rownames(signif.genes_tutorial1)
write.xlsx(signif.genes_tutorial1, file = "Comparison_LTR_timeseries_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)


#Significant genes table - TO CREATE A TABLE AND AN EXCEL FILE FOR STR TP2 Vs TP1.
signif.genes_tutorial2 <- top.table_tutorial2[(top.table_tutorial2$P.Value < 0.05), ]
signif.genes_tutorial2 <- signif.genes_tutorial2[1:46,]
head(signif.genes_tutorial2)

#Save the table as an excel file.
signif.genes_tutorial2$Protein <- rownames(signif.genes_tutorial2)
write.xlsx(signif.genes_tutorial2, file = "Comparison_STR_timeseries_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)

#Significant genes table - TO CREATE A TABLE AND AN EXCEL FILE FOR LTR Vs STR.
signif.genes_tutorial3 <- top.table_tutorial3[(top.table_tutorial3$P.Value < 0.05), ]
signif.genes_tutorial3 <- signif.genes_tutorial3[1:12,]
head(signif.genes_tutorial3)

#Save the table as an excel file.
signif.genes_tutorial3$Protein <- rownames(signif.genes_tutorial3)
write.xlsx(signif.genes_tutorial3, file = "Comparison_LTR_STR_timeseries_OS_ICI_36months.xlsx", sheetName = "Sheet1", rowNames = FALSE)

####Genomic data####
PGDx_data <- read_excel("PGDx_data.xlsx")
PGDx_data <- as.data.frame(PGDx_data)
head(PGDx_data)
####Genomic data - Create  a MAFTOOLLS FILE TO PROCEED WITH THE ANALYSIS####
#load excell file.
mutation_data <- read_excel("PGDx_data.xlsx", 
                            sheet = "Sequence Alterations")
mutation_data <- as.data.frame(mutation_data)
head(mutation_data)
#Separate Nucleotide_Position_(Genomic_hg19) into Chromossome, Start and End position.
mutation_data <- separate(mutation_data, 'Nucleotide_Position_(Genomic_hg19)', into = c("Chromosome", "Start_Position", "End_Position"), sep = "[:-]")

#Changing "SBS" with "SNP" in the variant_type column.
mutation_data <- mutation_data %>% 
  mutate(Variant_Type = ifelse(Variant_Type == "SBS", "SNP", Variant_Type))

#Change Variant_Classification as provided in Default uses Variant Classifications with High/Moderate variant consequences. 
#https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html: "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
#"Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"
mutation_data <- mutation_data %>% 
  mutate(Variant_Classification = case_when(
    Variant_Classification == "Frameshift" ~ "Frame_Shift_Del",
    Variant_Classification == "In-frame Deletion" ~ "In_Frame_Del",
    Variant_Classification == "In-frame Insertion" ~ "In_Frame_Ins",
    Variant_Classification == "Missense" ~ "Missense_Mutation",
    Variant_Classification == "Nonsense" ~ "Nonsense_Mutation",
    Variant_Classification == "Splice site acceptor" ~ "Splice_Site",
    Variant_Classification == "Splice Site Acceptor" ~ "Splice_Site",
    TRUE ~ Variant_Classification  # Keep other values unchanged
  ))

#Generate MAFTOOLS file.
maf <- maftools::read.maf(mutation_data_LTR, clinicalData = clinical_data)

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

#Summary plot of the mutation data available.
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

#Oncoplot including the the 40 most frequent mutations:
oncoplot(maf = maf, top = 40, removeNonMutated = FALSE, clinicalFeatures = c("Histology"),colors = colors_coBarplot, pwLineCol = "white",
         draw_titv = TRUE, titv_col = NULL, bgCol = "gray90",
         borderCol = "white",annoBorderCol = "white",
         sepwd_genes =1,
         annotationColor = oncoplotcolors, gene_mar = 6, fontSize = 0.7)

oncoplotcolors = c("#9ED2BA", "#EBCB29", "#BF4525")
names(oncoplotcolors) = c("Adenocarcinoma", "Squamous cell carcinoma", "NOS")
oncoplotcolors = list(Histology = oncoplotcolors)


annotationColor_oncoplot <- c("Adenocarcinoma" = "#9ED2BA",
                              "Squamous cell Carcinoma" = "#EBCB29",
                              "NOS"= "#BF4525",
                              "LTR" = "#FF9C62",
                              "STR" = "#00DBD3")

#Oncoplot including the the 40 most frequent mutations:
clinical_data$OS_ICI_36months <- factor(clinical_data$OS_ICI_36months)
oncoplot(maf = maf, top = 40, removeNonMutated = FALSE, clinicalFeatures = c("OS_ICI_36months"), pwLineCol = "white",
         draw_titv = TRUE, titv_col = NULL, bgCol = "gray90",
         borderCol = "white",annoBorderCol = "white",
         sepwd_genes =1, gene_mar = 6, fontSize = 0.7,
         annotationColor = oncoplotcolors_rx)

#By response:
oncoplotcolors_rx = c("#FF9C62","#00DBD3")
names(oncoplotcolors_rx) = c("LTR", "STR")
oncoplotcolors_rx = list(OS_ICI_36months = oncoplotcolors_rx)
oncoplotcolors = c("#9ED2BA", "#EBCB29", "#BF4525")
names(oncoplotcolors_rx) = c("Adenocarcinoma", "Squamous cell carcinoma", "NOS")
oncoplotcolors = list(Histology = oncoplotcolors)

#titv function classifies SNPs into transitions (Ti) and transversions (Tv).
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)


#Compare with other TCGA datasets - TYPICAL PLOT of TMB.
nsclc.mutload = tcgaCompare(maf = maf, cohortName = 'NSCLC-HMAR', logscale = TRUE, capture_size = 2.1, bg_col = c("#EDF8B1", "#2C7FB8"))
#Important to double if the capture size is the same compared to the "panel size".
nsclc.mutload = tcgaCompare(maf = maf, cohortName = 'NSCLC-HMAR', logscale = TRUE, capture_size = 2.1, bg_col = c("white", "gray5"))

#Exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))

#Survival analysis.
maf@clinical.data$last_FU_status <- factor(maf@clinical.data$last_FU_status)#last_FU_status should be a logical variable.

#Survival analysis based on grouping of gene mutation status
#Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset = survGroup(maf = maf, top = 20, geneSetSize = 2, time = "OS_months", Status = "last_FU_status", verbose = FALSE)

mafSurvGroup(maf = maf, geneSet = c("TP53", "KMT2D"), time = "OS_months", Status = "last_FU_status")
mafSurvGroup(maf = maf, geneSet = c("KEAP1","KRAS"), time = "OS_months", Status = "last_FU_status")

####RNA seq analysis - LIMMA ONLY BASELINE SAMPLES BASED ON QC FROM MS####
head(rnaseq_metadata_ms)
head(rna_qc_ms)
#Limma
cond_rna_ms <- as.factor(rnaseq_metadata_ms$response_MS)

# design matrix
design_rna_ms <- model.matrix(~0+cond_rna_ms)
rownames(design_rna_ms) <- rnaseq_metadata_ms$n_bx
colnames(design_rna_ms) <- gsub("cond_rna_ms", "", colnames(design_rna_ms))

voom.res_rna_ms <- voom(rna_qc_ms, design_rna_ms, plot = T) 

fit_rna_ms <- lmFit(voom.res_rna_ms, design_rna_ms) 

# contrasts
contrast.matrix_rna_ms <- makeContrasts(con1=LR-NR,
                                        levels = design_rna_ms) 

# contrasts fit and Bayesian adjustment
fit2_rna_ms <- contrasts.fit(fit_rna_ms, contrast.matrix_rna_ms)
fite_rna_ms <- eBayes(fit2_rna_ms)

# summary 
summary(decideTests(fite_rna_ms, method = "separate"))

summary(decideTests(fite_rna_ms, adjust.method = "none", method = "separate")) 

top.table_rna_ms <- topTable(fite_rna_ms, number = Inf, adjust = "fdr")
head(top.table_rna_ms)
hist(top.table_rna_ms$P.Value, breaks = 100, main = "Results p-values")

#Significant genes table
signif.genes_rna_ms <- top.table_rna_ms[(top.table_rna_ms$P.Value < 0.05), ]
head(signif.genes_rna_ms)

signif.genes_rna_ms_adj <- top.table_rna_ms[(top.table_rna_ms$adj.P.Val < 0.05), ]

##HEATMAP.
#list of genes DE between LTR and NR.
head(signif.genes_rna_ms_adj)
geneshm_rna_ms_adj<-  rownames(signif.genes_rna_ms_adj)

#Matrix for heatmap.
hm_rna_ms_adj <- countsTMM[rownames(countsTMM) %in% geneshm_rna_ms_adj,]

#Annotations for heatmap.
#Annotations - should be a dataframe with SampleID and Response variable:
head(rnaseq_metadata_ms)
colnames(rnaseq_metadata_ms)
#Order the 
annotations_hm_rna_ms <- rnaseq_metadata_ms[,c(3,4,7)]
rownames(annotations_hm_rna_ms) <- rnaseq_metadata_ms$n_bx
head(annotations_hm_rna_ms)

#Heatmap using the genes with adj-p-value <0.05.
pheatmap(hm_rna_ms_adj, scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         border_color = "black",
         annotation_col = annotations_hm_rna_ms
)

#Heatmap using the genes with unadj-p-value <0.05.
pheatmap(hm_rna_ms, scale="row",
         cluster_cols = TRUE, clustering_distance_cols = "correlation",
         clustering_method = "ward.D",
         border_color = "black",
         annotation_col = annotations_hm_rna_ms
)

#VOLCANO PLOT - for adjusted p-value <0.05
#Volcano plot only with genes with adjusted p-value.
volcano_data_rna_ms <- top.table_rna_ms
volcano_data_rna_ms$Gene <- rownames(volcano_data_rna_ms)
#Add a new column for the proteins with pvalue<0.05 and FC 0.5.
volcano_data_rna_ms$Different_gene <- "No"
# if estimate > 0.5 and pvalue < 0.05, set as "UP" 
volcano_data_rna_ms$Different_gene[volcano_data_rna_ms$logFC >2 & volcano_data_rna_ms$adj.P.Val < 0.05] <- "Up"
# if estimate < -0.5 and pvalue < 0.05, set as "DOWN"
volcano_data_rna_ms$Different_gene[volcano_data_rna_ms$logFC < -2 & volcano_data_rna_ms$adj.P.Val< 0.05] <- "Down"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
volcano_data_rna_ms$labels <- NA
volcano_data_rna_ms$labels[volcano_data_rna_ms$Different_gene != "No"] <- volcano_data_rna_ms$Gene[volcano_data_rna_ms$Different_gene != "No"]


ggplot(volcano_data_rna_ms, aes(x=logFC, y=-log10(adj.P.Val),col=Different_gene, label = labels)) +
  geom_point() +
  theme_classic2() +
  geom_vline(xintercept=c(-2, 2), col="gray40",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="gray40",linetype="dashed") +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  scale_color_manual(values=c("gray", "salmon")) +
  labs(title = "Tissue RNAseq baseline comparison LTR Vs STR")

#VOLCANO PLOT - for unadjusted p-value <0.05
volcano_data_rna_ms <- top.table_rna_ms
volcano_data_rna_ms$Gene <- rownames(volcano_data_rna_ms)
#Add a new column for the proteins with pvalue<0.05 and FC 0.5.
volcano_data_rna_ms$Different_gene <- "No"
# if estimate > 0.5 and pvalue < 0.05, set as "UP" 
volcano_data_rna_ms$Different_gene[volcano_data_rna_ms$logFC >2 & volcano_data_rna_ms$P.Value < 0.05] <- "Up"
# if estimate < -0.5 and pvalue < 0.05, set as "DOWN"
volcano_data_rna_ms$Different_gene[volcano_data_rna_ms$logFC < -2 & volcano_data_rna_ms$P.Value < 0.05] <- "Down"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
volcano_data_rna_ms$labels <- NA
volcano_data_rna_ms$labels[volcano_data_rna_ms$Different_gene != "No"] <- volcano_data_rna_ms$Gene[volcano_data_rna_ms$Different_gene != "No"]

#Volcano plot.
ggplot(volcano_data_rna_ms, aes(x=logFC, y=-log10(P.Value),col=Different_gene, label = labels)) +
  geom_point() +
  theme_classic2() +
  geom_vline(xintercept=c(-2, 2),  col="gray40",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05),  col="gray40",linetype="dashed") +
  geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  scale_color_manual(values=c("cornflowerblue", "gray80", "salmon")) +
  labs(title = "Tissue RNAseq baseline comparison LTR Vs STR")