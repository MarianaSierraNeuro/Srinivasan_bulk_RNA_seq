###following reading_count_data_optional

#Title: "SCA1 Excercise vs Control Mice"
#Authors: Sharan Srinivasan, Mariana Sierra
#date: 05/20/2024

#Install Packages
install.packages(c('dyplr','ggplots','ggplot2','ggrepel','RColorBrewer','airway'))
BiocManager::install(c('limma','DESeq2','AnnotationDbi','org.Mm.eg.db','ReportingTools','Go.db','GOstats','pathview','gage','gageData','select'))

#install libraries 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(airway)
library(tidyverse)


###0. Treatment Information/Metadata ----

#we need to make the genotype and treatment groups as two new columns in our dataset so we can have them as a metadata available to use for DESeq to actually differentiate our groups
#right now I have them set randomly so this will need to be changed - I'm setting them as vectors as well to work with them downstream our you're not going to be able to do much with the dataframe
genotype = c('wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','sca1')
treatment = c('sedentary','sedentary','sedentary','sedentary','sedentary','sedentary','early start','early start','early start','early start','early start','early start','late start','late start','late start','late start','late start')

#here I am appending the datasets I made to a new dataframe with the sample, geneotype, and treament
colData = as.data.frame(cbind(colnames(counts.data),genotype,treatment))

###1. DESeq2 Quick analysis ----

dds = DESeqDataSetFromMatrix(countData = counts.data,
                             colData = colData,
                             design = ~genotype+treatment+genotype*treatment) #when we do the * we are saying we are interested in seeing the relationship between the genotype and exercise that would then be amplified by one or the other (we can always just comapre one group or the other but since they might have different expression profiles due to the genotype this is the safest bet to start)
dds = DESeq(dds)
nrow(dds) #we expect the same number of rows as we have genes available - in our case all good :)

#optional: but we can remove rows with very low gene expression, I am opting for it since we have some genes that are not expressed in any of the samples so it doesn't make sense to keep them
dds = dds[rowSums(counts(dds))>5,]
nrow(dds) #we've now cut down significantly :)

#Size factor: The DESeq size factor in differential expression analysis tells you about the normalization factor applied to the raw counts of each sample to account for differences in sequencing depth or library size between samples.
sizeFactors(dds)
### 1.1 PCA Plot ----
#PCA plot for individual treamtent
rld = rlog(dds)
plotPCA(rld, intgroup="genotype")
plotPCA(rld, intgroup="treatment")
#PCA plot for both genotype and treament on the same plot
plotPCA(rld, intgroup= c("genotype","treatment"))

### 1.2 Heat Maps ----

#will look at our data at tell us which sample belongs to which group
detectGroups <- function (x) { #in this code x is equvalent to the column names 
  tem <- gsub("","",x) #removes all numbers from the end 
  #tem = gsub("_Rep/_rep/REP","",tem)
  tem <- gsub("_$","",tem); #remove _ from the end
  tem <- gsub("_Rep$","",tem); #remove _Rep from the end
  tem <- gsub("_rep$","",tem); #remove _rep from the end
  tem <- gsub("_REP$","",tem); #remove _REP from the end
  return(tem)
}

detectGroups(colnames(counts.data))

#caluclates distance between all of our groups
dist2 <- function(x,...){
  as.dist(1-cor(t(x),method = "pearson"))
}

#how similar a genes is between treamtents 
hclust2 <- function(x,method='average',...){
  hclust(x,method=method,...)
}

n = 50 #number of genes top 50 genes

x=assay(rld) 
if(n>dim(x)[1]) n = dim(x)[1] # max	as data

x = x[order(apply(x,1,sd),decreasing=TRUE),]  # sort genes by standard deviation

x = x[1:n,]   # only keep the n genes

# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
x[x< cutoff] <- cutoff

groups = detectGroups(colnames(x) )
groups.colors = rainbow(length(unique(groups) ) )


lmat = rbind(c(5,4),c(0,1),c(3,2))
lwid = c(1.5,4)
lhei = c(1,.2,4)


heatmap.2(x, distfun = dist2,hclustfun=hclust2,
          col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
          ,key=T, symkey=F
          ,ColSideColors=groups.colors[ as.factor(groups)]
          ,margins=c(8,12)
          ,cexRow=1
          ,srtCol=45
          ,cexCol=1.  # size of font for sample names
          ,lmat = lmat, lwid = lwid, lhei = lhei
)  

###1.3 Heatmap to compare samples to each other ----

#make variance stabalized transformation 
vsd = vst(dds,blind=FALSE)
#to generate distance matrix (without the function above - I'm trying someting different)
sampleDists = dist(t(assay(vsd)))
sampleDistMartix = as.matrix(sampleDists)
colnames(sampleDistMartix)

#setting color scheme
colors = colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

#sample heatmap
pheatmap(sampleDistMartix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col=colors)
         
#generating the heatmap


###1.4 Volcano Plot ----

# Set reference levels for genotype and treatment
dds$genotype <- relevel(dds$genotype, ref = "wt")
dds$treatment <- relevel(dds$treatment, ref = "sedentary")

plotMA(dds)




  
