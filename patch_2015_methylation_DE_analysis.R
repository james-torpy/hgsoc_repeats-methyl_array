library(tibble)
library(dplyr)
library(qvalue)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

home_dir <- "/Users/jamestorpy/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "/projects/hgsoc_repeats/methyl_array/")
result_dir <- paste0(project_dir, "/results/")
Robject_dir <- paste0(project_dir, "/Robjects/")
ref_dir <- paste0(project_dir, "/refs/"
plot_dir <- paste0(result_dir, "plots/")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", Robject_dir))


##############################################################################
### 0. Load in data, filter and create m-values ###
##############################################################################

# load in beta-values:
b_values <- read.csv(file = paste0(result_dir, 
  "/bowtell_trimmed_methylation_array_b-values.csv"), header = T)

# make first column rownames and remove:
rownames(b_values) <- b_values[,1]
b_values <- subset(b_values, select = -X)

# create function to split values into control and treatment data frames and
# filter for QC:
split_and_filter <- function(x, min_rows) {
  # split into HGSOC and FT values:
  FT <- x[,grep("ICGC", colnames(x), invert = T)]
  HGSOC <- x[,grep("ICGC", colnames(x))]
  
  # remove any rows with less than 4 values per group:
  print(paste0("No. rows in FT before filtering is: ", nrow(FT)))
  FT <- FT %>%
    rownames_to_column('probe_id') %>%
    dplyr::filter((rowSums(!is.na(FT)) >= min_rows)) %>%
    column_to_rownames('probe_id')
  print(paste0("No. rows in FT after filtering is: ", nrow(FT)))
  
  print(paste0("No. rows in HGSOC before filtering is: ", nrow(HGSOC)))
  HGSOC <- HGSOC %>%
    rownames_to_column('probe_id') %>%
    dplyr::filter((rowSums(!is.na(HGSOC)) >= min_rows)) %>%
    column_to_rownames('probe_id')
  print(paste0("No. rows in HGSOC after filtering is: ", nrow(HGSOC)))
  
  return(list(FT, HGSOC))
}

if ( !file.exists(paste0(Robject_dir, 
                         "/bowtell_trimmed_methylation_array_b-values.rds")) ) {
  # split and filter beta values:
  b_list <- split_and_filter(b_values, 4)
  saveRDS(b_list, file = paste0(Robject_dir, 
                                "/bowtell_trimmed_methylation_array_b-values.rds"))
  
} else {
  b_list <- readRDS(paste0(Robject_dir, 
                           "/bowtell_trimmed_methylation_array_b-values.rds"))
}

# assign control and treatment beta value data frames:
FT_b <- b_list[[1]]
HGSOC_b <- b_list[[2]]


if ( !file.exists(paste0(Robject_dir, 
                         "/bowtell_trimmed_methylation_array_m-values.rds")) ) {
  # replace NA beta-values with '2':
  b_temp <- b_values
  b_temp[is.na(b_temp)] <- 2
  # convert all beta-values to m_values - '2' will automatically be converted
  # to NA:
  m_values <- log2(b_temp/(1-b_temp))
  m_list <- split_and_filter(m_values, 4)
  saveRDS(m_list, file = paste0(Robject_dir, 
                                "/bowtell_trimmed_methylation_array_m-values.rds"))
  
} else {
  m_list <- readRDS(paste0(Robject_dir, 
                           "/bowtell_trimmed_methylation_array_m-values.rds"))
}

# assign control and treatment m-value data frames:
FT_m <- m_list[[1]]
HGSOC_m <- m_list[[2]]


##############################################################################
### 1. Apply statistical tests to beta and m-values for HGSOC vs FT DE 
# analysis ###
##############################################################################

apply_test <- function(ctl, treat, test_type) {
	if ( test_type == "t-test" ) {
		i=1
		pvals <- apply(ctl, 1, function(x) {
			print(i)
		
			tres <- tryCatch(
				t.test(as.numeric(x), as.numeric(treat[i,]), conf.level = 0.9),
				error=function(err) NA
			)
		
			if ( !is.na(tres) ) {
				res <- tres$p.value
			} else {
				res <- NA
			}
		
			i <<- i+1
			return(res)
		})

	} else if ( test_type == "wilcox") {
		i=1
		pvals <- apply(ctl, 1, function(x) {
			print(i)
		
			wres <- tryCatch(
				wilcox.test(as.numeric(x), as.numeric(treat[i,]), conf.level = 0.9),
				error=function(err) NA
			)
		
			if ( !is.na(wres) ) {
				res <- wres$p.value
			} else {
				res <- NA
			}
		
			i <<- i+1
			return(res)
		})
	}
	return(pvals)
}

if ( !file.exists(paste0(Robject_dir, "/bval_ps.rds")) ) {
	# apply unpaired Wilcox test to beta values:
	bval_ps <- apply_test(FT_b, HGSOC_b, "wilcox")
	names(bval_ps) <- rownames(FT_b)
	saveRDS(bval_ps, file = paste0(Robject_dir, "/bval_ps.rds"))
} else {
  bval_ps <- readRDS(file = paste0(Robject_dir, "/bval_ps.rds"))
}


if ( !file.exists(paste0(Robject_dir, "/mval_ps.rds")) ) {
	# apply unpaired Wilcox test to beta values:
	# apply unpaired t-test to beta values:
	mval_ps <- apply_test(FT_m, HGSOC_m, "t-test")
	names(mval_ps) <- rownames(FT_m)
	saveRDS(mval_ps, file = paste0(Robject_dir, "/mval_ps.rds"))
} else {
  mval_ps <- readRDS(file = paste0(Robject_dir, "/mval_ps.rds"))
}


##############################################################################
### 2. Make median beta and m-value bedGraphs files for IGV ###
##############################################################################

# load EPIC 450K annotation:
annot <- read.csv(paste0(ref_dif, "/methyl450k_annot.csv", fill=T, header=F)


med_FT_b <- apply(FT_b, 1, median)
med_HGSOC_b <- apply(HGSOC_b, 1, median)



##############################################################################
### 3. Calculate q-values, median difference between HGSOC and FT 
# methylation ###
##############################################################################

# check p-value distributions:
hist(bval_ps, nclass=20)
pdf(paste0(plot_dir, "/bval_ps_hist.pdf"))
hist(bval_ps, nclass=20)
dev.off()


hist(mval_ps, nclass=20)
pdf(paste0(plot_dir, "/mval_ps_hist.pdf"))
hist(mval_ps, nclass=20)
dev.off()

# calculate qvalues:
bval_qs <- qvalue(bval_ps)$qvalues
rownames(bval_qs) <- rownames(bval_ps)

mval_qs <- qvalue(mval_ps)$qvalues
rownames(mval_qs) <- rownames(mval_ps)

# limit results to those with q-values < 0.1
bval_sig_q <- bval_qs[bval_qs < 0.1]

mval_sig_q <- mval_qs[mval_qs < 0.1]

print(paste0("64,313 probes were reported as differentially",
	"methylated between HGSOC and FT with q < 0.1, compared", 
	"to ", length(bval_sig_q), " probes reported as DE using", 
	"beta-values, and ", length(mval_sig_q), " probes reported",
	"using m-values"))

# create function to add medians for control and treatment samples:
add_medians <- function(qvals, ctl, treat) {
	ctl <- ctl[rownames(ctl) %in% rownames(qvals),]
	ctl_med <- median(ctl)

	treat <- treat[rownames(treat) %in% rownames(qvals),]
	treat_med <- median(treat)

	return(data.frame(qvals, ctl_med, treat_med))
}

# calculate and add medians to bval_qs as data frame:
b_stats <- add_medians(bval_sig_q, FT_b, HGSOC_b)
colnames(m_stats) <- c("qvalue", "FT_median", "HGSOC_median")
rownames(b_stats) <- rownames(bval_sig_q)

# calculate differences between medians:
b_stats$median_diffs <- b_stats$HGSOC_median - b_stats$FT_median
b_stats$direction <- NA
b_stats$direction[b_stats$median_diffs < 0] <- "less_methylated"
b_stats$direction[b_stats$median_diffs > 0] <- "more_methylated"


# calculate and add medians to bval_qs as data frame:
m_stats <- add_medians(mval_sig_q, FT_m, HGSOC_m)
colnames(m_stats) <- c("qvalue", "FT_median", "HGSOC_median")
rownames(m_stats) <- rownames(mval_sig_q)

# calculate differences between medians:
m_stats$median_diffs <- m_stats$HGSOC_median - m_stats$FT_median
m_stats$direction <- NA
m_stats$direction[m_stats$median_diffs < 0] <- "less_methylated"
m_stats$direction[m_stats$median_diffs > 0] <- "more_methylated"













#############################################
i=1
meth_ps <- apply(FT_m, 1, function(x) {
	print(i)

	tres <- tryCatch(
		t.test(x, HGSOC_m[i,], conf.level = 0.9),
		error=function(err) NA
	)

	if ( !is.na(tres) ) {
		res <- tres$p.value
	} else {
		res <- NA
	}

	i <<- i+1
	return(res)
})

[41790:41795,]

names(meth_ps) <- rownames(FT_m)


