suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
BASE_CONDITION_SAMPLES <- args[2]
EXPERIMENTAL_CONDITION_SAMPLES <- args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
OUTPUT_MAST_FILE_BASE <- 'mast_results'

# I split this from the other conditional later so that we can get the checks 
# out of the way before the time intensive import of data, etc.
if (CONDITION_B == "FAKE") {# some conditional for the 1-vs-all = TRUE
    CONDITION_A <- "background" # overwrite / define the CONDITION_A
    contrast_str = paste0(CONDITION_B, '_vs_', CONDITION_A)
} else { # 
    # create a string to identify the contrast:
    contrast_str = paste0(CONDITION_B, '_vs_', CONDITION_A)


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
    all_samples <- c(base_samples, exp_samples)
}


# Import counts as a data.frame
cnts <- read.table(
    file = RAW_COUNT_MATRIX,
    sep = "\t",
    row.names = 1
)
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# Preprocess the counts with normalization
sce <- seuratSCTransform(
    inSCE = sce, 
    normAssayName = "SCTCounts", 
    useAssay = "counts"
)

# Create the labels list
if (CONDITION_B == "FAKE") { # some conditional for the 1-vs-all = TRUE
    labels <- c()
    # A more elegant solution should be implemented.
    # This is "slow".
    for (sample in colnames(sce)) {
        if (sample %in% exp_samples) {
            labels <- c(labels, CONDITION_B)
        } else {
            labels <- c(labels, CONDITION_A)
        }
    }
} else { # else, do a-vs-b differential expression
    labels <- c()
    for (sample in colnames(sce)) {
        if (sample %in% exp_samples) {
            labels <- c(labels, CONDITION_B)
        } else if (sample %in% base_samples) {
            labels <- c(labels, CONDITION_A)
        } else {
            labels <- c(labels, "bobby_drop_tables") # We need some protected NULL type here
        }
    }
}

# Add the labels to the SCE object
colLabels(sce) <- as.factor(labels)

# Run the DGE
sce <- runDEAnalysis(
    inSCE = sce, 
    method = "MAST", 
    useAssay = "SCTCounts",
    class = "label", 
    classGroup1 = CONDITION_B, 
    classGroup2 = CONDITION_A,
    groupName1 = "exp", 
    groupName2 = "base", 
    analysisName = "exp_VS_base",
    log2fcThreshold = NULL, # no filters
    fdrThreshold = NULL
)

# Output results to dataframe
df.results <- metadata(sce)$diffExp$exp_VS_base$result

# Write results to file
output_filename <- paste(OUTPUT_MAST_FILE_BASE, contrast_str, 'tsv', sep='.')
write.table(
    df.results, 
    output_filename, 
    sep='\t', 
    quote=F, 
    row.names = FALSE
)