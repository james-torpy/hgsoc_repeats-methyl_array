### diff_exp_of_meth_regions_from_array.R ###

# This script takes methylation beta-values from 450k EPIC array data and
# calculates differential expression of regions between samples and 
# controls, outputting Gviz methylation plots of DE regions containing
# repeats

library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
#library(missMethyl)
#library(matrixStats)
#library(minfiData)
#library(Gviz)
library(DMRcate)
library(stringr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "/projects/hgsoc_repeats/methyl_array/")
result_dir <- paste0(project_dir, "/results/")
Robject_dir <- paste0(project_dir, "/Robjects/")

system(paste0("mkdir -p " result_dir))
system(paste0("mkdir -p " Robject_dir))


##############################################################################
### 0. Load in data and convert to m-values ###
##############################################################################

# load in beta-values:
b_values <- read.csv(file = paste0(result_dir, 
  "/bowtell_trimmed_methylation_array_b-values.csv"), header = T)

# make first column rownames and remove:
rownames(b_values) <- b_values[,1]
b_values <- subset(b_values, select = -X)

b_IDs <- gsub("\\ICGC.*_", "", colnames(b_values))

# generate vector of values with duplicated AOCS IDs (keep earliest):
rm(list = c("rem", "rec", "dupes", "keep"))
rem <- c()
for ( i in 1:length(b_IDs) ) {
  if ( length(grep("ICGC", colnames(b_values)[i])) > 0 ) {
    print(i)
    if (i==1) {
      rec <- c(b_IDs[i])

    } else {
      id <- gsub("\\_[0-9][0-9][0-9]$", "", b_IDs[i])

      if ( length(grep(id, rec)) > 0 ) {
        rec <- append(rec, b_IDs[i])
        dupes <- rec[grep(id, rec)]
        keep <- dupes[grep(min(gsub("^.*_", "", dupes)), dupes)]
        rem <- append(rem, dupes[!(dupes %in% keep)])
      } else {
        rec <- append(rec, b_IDs[i])
      }
    }
  }  
}

# remove duplicated AOCS ID samples:
b_temp <- b_values[, !(b_IDs %in% rem)]

# convert all NAs to the value '2' (all other values are in range 0-1):
b_temp[is.na(b_temp)] <- 2

# convert all beta-values to m_values - '2' will automatically be converted
# to NA:
m_values <- log2(b_temp/(1-b_temp))

saveRDS(m_values, file = paste0(Robject_dir, 
        "/bowtell_trimmed_methylation_array_m-values.rds"))


##############################################################################
### 1. Annotate data and generate bedgraph file ###
##############################################################################

m_mat <- as.matrix(m_values)

# create factor vector of categories for samples:
cats <- rep(NA, ncol(m_values))
cats[grep("ICGC", colnames(m_values))] <- "HGSOC"
cats[grep("ICGC", colnames(m_values), invert = T)] <- "FT"
cats <- as.factor(cat)
levels(cats) = c("FT", "HGSOC")

# create design for differential expression:
design <- model.matrix(~0+cats)
colnames(design) <- c("FT", "HGSOC")
cont_matrix <- makeContrasts(HGSOCvsFT=HGSOC-FT, levels=design)

m_annot <- cpg.annotate(
  object = m_mat, 
  datatype = "array", 
  what = "M",
  analysis.type = "differential", 
  design = design,
  contrasts = TRUE,
  cont.matrix = cont_matrix,
  coef = "HGSOCvsFT", 
  arraytype = "450K",
  fdr = 0.1
)







