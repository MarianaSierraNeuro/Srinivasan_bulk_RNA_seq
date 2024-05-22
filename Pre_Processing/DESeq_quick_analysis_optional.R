###following reading_count_data_optional
#____
#Title: "SCA1 Excercise vs Control Mice"
#Authors: Sharan Srinivasan, Mariana Sierra
#date: "05/20/2024"
#output = html_document
#____

#'''{r setup, include = FALSE}
#knit::opts_chunk$set(echo=TRUE)
#false = FALSE
#true = TRUE
#'''

#Install Packages
#'''{r, eval=false}
install.packages(c('dyplr','ggplots','ggplot2','ggrepel'))
BiocManager::install(c('limma','DESeq2','AnnotationDbi','org.Mm.eg.db','ReportingTools','Go.db','GOstats','pathview','gage','gageData','select'))
#'''
#'''{r}
library(DESeq2)
library(dplyr)
#'''

###0. Treatment Information/Metadata ----

genotype = c('wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','wt','sca1','sca1')

treatment = c('sedentary','sedentary','sedentary','sedentary','sedentary','sedentary','early start','early start','early start','early start','early start','early start')


