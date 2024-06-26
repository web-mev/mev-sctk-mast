suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("optparse")))

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
  make_option(
    c('-f','--input_file'),
    help='Path to the count matrix input.'
  ),
  make_option(
    c('-o','--output_file_prefix'),
    help='The prefix for the output file'
  ),
  make_option(
    c('-a','--experimental_samples'),
    help='The experimental samples, provided as a comma-delimited string.'
  ),
  make_option(
    c('-b','--base_samples'),
    help='The base samples, provided as a comma-delimited string. Optional.'
  ),
  make_option(
    c('--experimental_group_name'),
    help='The name of the experimental sample group'
  ),
  make_option(
    c('--base_group_name'),
    help='The name of the base sample group. Optional, and only relevant if --base_samples given.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}

# No matter what, we need at least the name of the experimental group and which samples
# are in that group. 
EXP_CONDITION_NAME <- opt$experimental_group_name
if (is.null(EXP_CONDITION_NAME)){
    message('Improper call. Need to specify the name of the experimental group via the --experimental_group_name arg.')
    quit(status=1)
}
EXPERIMENTAL_CONDITION_SAMPLES <- opt$experimental_samples
if (is.null(EXPERIMENTAL_CONDITION_SAMPLES)){
    message('Improper call. Need to specify the names of the samples in the group of interest.')
    quit(status=1)
}
exp_samples <- make.names(strsplit(EXPERIMENTAL_CONDITION_SAMPLES, ',')[[1]])

if (is.null(opt$base_samples)) {

    # a flag for later use
    ONE_VS_ALL = TRUE

    # Creates a 1-vs-all contrast
    contrast_str = paste0(EXP_CONDITION_NAME, '_vs_all')
    BASE_CONDITION_NAME <- 'other'
    base_samples <- c() 

} else { 

    # a flag for later use
    ONE_VS_ALL = FALSE

    # create a string to identify the contrast:
    BASE_CONDITION_NAME <- opt$base_group_name
    if (is.null(BASE_CONDITION_NAME)){
        message('Improper call. Need to specify the name of the base group via the --base_group_name arg.')
        quit(status=1)
    }

    contrast_str = paste0(EXP_CONDITION_NAME, '_vs_', BASE_CONDITION_NAME)

    # to enter this block, we had to have opt$base_samples as non-null, so don't check again
    BASE_CONDITION_SAMPLES <- opt$base_samples
    
    # the sample names are given as a comma-delimited string. Split them
    base_samples <- make.names(strsplit(BASE_CONDITION_SAMPLES, ',')[[1]])
    exp_samples <- make.names(strsplit(EXPERIMENTAL_CONDITION_SAMPLES, ',')[[1]])
    intersection_list = intersect(base_samples, exp_samples)

    # Checks if there is any overlap of samples in the two conditions
    # Returns an error code if TRUE
    if (length(intersection_list) > 0){
        sample_list = paste0(intersection_list, collapse=',')
        message(paste(
        'The following samples were in both contrast groups. Fix this and try again: ',
        sample_list
        ))
        quit(status=1)
    }

    # a list of all the samples we need. Note that these could be completely bogus and
    # not correspond to anything in the matrix.
    all_samples <- c(base_samples, exp_samples)
}


# Import counts as a data.frame
cnts <- read.table(
    file = opt$input_file,
    sep = "\t",
    row.names = 1,
    header=T
)

# ensure the sample names are kosher and subset the matrix if necessary.
if (ONE_VS_ALL) {
    # only need to check the experimental samples in this case.
    check_samples = exp_samples
} else {
    check_samples = all_samples
}
 
# check the set difference. If there are samples in `check_samples` 
# that are NOT columns of the matrix, then error out immediately.
sd <- setdiff(check_samples, colnames(cnts))
if ( length(sd) > 0 ){
    sample_list = paste0(sd, collapse=',')
    message(paste(
        'The following samples were not in your matrix: ',
        sample_list
    ))
    quit(status=1)
}

# If we are doing a A vs. B contrast, subset to cut down the matrix
if (!ONE_VS_ALL){
    cnts = cnts[,all_samples]
}

# change to a sparse matrix representation, necessary for SCE
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# Preprocess the counts with log-normalization
sce <- scaterlogNormCounts(sce, assayName = "logcounts")

# Create the labels list
labels <- c()
# A more elegant solution should be implemented.
# This is "slow".
for (sample in colnames(sce)) {
    if (sample %in% exp_samples) {
        labels <- c(labels, EXP_CONDITION_NAME)
    } else {
        labels <- c(labels, BASE_CONDITION_NAME)
        if (ONE_VS_ALL) {
            # if a one-vs-all biomarker analysis, no 'base_samples'
            # were defined. However, we will need this list explicitly
            # later on. 
            base_samples <- c(base_samples, sample)
        }
    }
}

if(is.null(base_samples)){
  if(length(exp_samples) == dim(sce)[2]){
    message('All the samples were contained in one contrast group and there is nothing to compare against. Check your inputs and try again.')
  } else {
    message('There were no samples to compare to. Please check your inputs.')
  }
  quit(status=1)
}

# Add the labels to the SCE object
colLabels(sce) <- as.factor(labels)

# Run the DGE
sce <- runDEAnalysis(
    inSCE = sce, 
    method = "MAST", 
    useAssay = "logcounts",
    class = "label", 
    classGroup1 = EXP_CONDITION_NAME, 
    classGroup2 = BASE_CONDITION_NAME,
    groupName1 = "exp", 
    groupName2 = "base", 
    analysisName = "exp_VS_base",
    log2fcThreshold = NULL, # no filters
    fdrThreshold = NULL
)

# Output results to dataframe
df.results <- as.data.frame(metadata(sce)$diffExp$exp_VS_base$result)
rownames(df.results) <- df.results$Gene
df.results <- subset(df.results, select=-Gene)

# Summarize the counts
counts <- assay(sce, 'counts')
cs = Matrix::rowSums(t(counts)) # gets the column sums
rescaled = 1e6*sweep(counts,2,cs, '/')

expCounts <- as.matrix(rescaled[, exp_samples])
baseCounts <- as.matrix(rescaled[, base_samples])

rqExp <- as.data.frame(rowQuantiles(expCounts))
rqBase <- as.data.frame(rowQuantiles(baseCounts))
colnames(rqExp) <- c('minExp', 'q1Exp', 'medianExp', 'q3Exp', 'maxExp')
colnames(rqBase) <- c('minBase', 'q1Base', 'medianBase', 'q3Base', 'maxBase')
rqExp$iqrExp <- rqExp$q3 - rqExp$q1
rqBase$iqrBase <- rqBase$q3 - rqBase$q1

# merge the summary info
rq <- merge(rqExp, rqBase, by='row.names')
rownames(rq) <- rq$Row.names
rq <- subset(rq, select=-Row.names)

# merge the two dataframes:
merged_data <- merge(df.results, rq, by='row.names')
rownames(merged_data) <- merged_data$Row.names
merged_data <- subset(merged_data, select=-Row.names)

# rename the Log2_FC and Pvalue columns so it works with other WebMeV tools like gsea
# which expect to see columns with particular names
names(merged_data)[names(merged_data) == 'Log2_FC'] <- 'log2FoldChange'
names(merged_data)[names(merged_data) == 'Pvalue'] <- 'pvalue'
names(merged_data)[names(merged_data) == 'FDR'] <- 'padj'

merged_data <- merged_data[order(merged_data$pvalue),]

# Write results to file
output_filename <- paste(opt$output_file_prefix, contrast_str, 'tsv', sep='.')
write.table(
    merged_data, 
    output_filename, 
    sep='\t', 
    quote=F, 
    row.names = TRUE
)
