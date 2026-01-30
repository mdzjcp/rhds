## -------------------------------------------- ##
## predict-proteins.r
## -------------------------------------------- ##
library(meffonym)

my.write.table <- function(x, filename) {
  cat("saving", basename(filename), "...\n")
  write.table(x, file = filename, row.names = T, col.names = T, sep = "\t")
}

args <- commandArgs(trailingOnly = T)
datadir <- args[1]
resultsdir <- args[2]

methylation.file <- file.path(datadir, "methylation-clean-score-sites.csv.gz")

## read dnam file
data <- as.data.frame(data.table::fread(methylation.file))
rownames(data) <- data[,1]
data <- as.matrix(data[,-1])

## check number of rows missing per sample
# miss <- apply(data, 2, function(i) table(is.na(i)), simplify=F)
# miss.df <- as.data.frame(do.call(rbind, miss))
# summary(miss.df$"TRUE")

## Before the all na row drop above ~90k observations
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  89519   89566   89645   89790   89817   95558

## get gadd et all episcores models
models <- subset(
  meffonym.models(full = T),
  grepl("^episcores", filename)
)

# get list of proteins to estimate
proteins <- models$name

# apply protein abundance coefs to dna methylation
pred.proteins <- sapply(
  proteins,
  function(model) {
    cat(date(), model, " ")
    ret <- meffonym.score(data, model)
    cat(
      " used ", length(ret$sites), "/",
      length(ret$vars), "sites\n"
    )
    ret$score
  }
)
rownames(pred.proteins) <- colnames(data)
pred.proteins <- scale(pred.proteins)
colnames(pred.proteins) <- make.names(colnames(pred.proteins))

## export results
my.write.table(
  pred.proteins,
  file.path(resultsdir, "predicted-proteins.txt")
)

## -------------------------------------------- ##
## combine.r
## -------------------------------------------- ##
args <- commandArgs(trailingOnly = T)
datadir <- args[1]
resultsdir <- args[2]

my.write.table <- function(x, filename) {
  cat("saving", basename(filename), "...\n")
  write.table(x, file = filename, row.names = T, col.names = T, sep = "\t")
}

## The format of sample identifiers/barcodes is described here:
## https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
##
## Here is a summary:
## e.g. TCGA-3C-AAAU-01A-11D-A41Q-05
##   project TCGA
##   tissue source site 3C
##   participant AAAU
##   sample 01 (01-09 tumor, 10-19 normal, 20-29 controls)
##   vial A
##   portion 11
##   analyte D (as in DNA)
##   plate A41Q
##   analysis center 05
##
## The following function extracts the participant identifier
## from a sample id/barcode.

extract.participant <- function(id) {
  sub("TCGA-[^-]+-([^-]+)-.*", "\\1", id)
}

extract.tissue <- function(id) {
  sub("TCGA-[^-]+-[^-]+-([0-9]+)[^-]+-.*", "\\1", id)
}

pred.protein.filename <- file.path(resultsdir, "predicted-proteins.txt")
clinical.filename <- file.path(resultsdir, "clinical-clean.txt")

pred.proteins <- read.table(pred.protein.filename,
  header = T, sep = "\t", stringsAsFactors = F
)

## extract participant tissue information
tissues <- data.frame(
  participant = extract.participant(rownames(pred.proteins)),
  tissue = extract.tissue(rownames(pred.proteins)),
  participant.tissue = paste(extract.participant(rownames(pred.proteins)),
    extract.tissue(rownames(pred.proteins)),
    sep = "-"
  )
)
tissues <- subset(tissues, tissue != "06" & tissue != "V582")

## update pred.proteins to use participant.tissue rownames
samples <- rownames(pred.proteins)
rownames(pred.proteins) <- paste(extract.participant(samples),
  extract.tissue(samples),
  sep = "-"
)

## get cleaned clinical data
clinical <- read.table(clinical.filename,
  header = T, sep = "\t", stringsAsFactors = F
)

## combine with participant tissue info from predicted protein dataset
clinical <- merge(clinical, tissues, by.x = "participant")
clinical$tumor.or.normal <- ifelse(as.numeric(clinical$tissue) < 9, "tumor", "normal")
clinical$tumor <- sign(clinical$tumor.or.normal == "tumor")


table(rownames(pred.proteins) %in% clinical$participant.tissue)

## combine the clinical info with the methylation predicted protein abundances
out <- cbind(
  clinical,
  pred.proteins[match(clinical$participant.tissue, rownames(pred.proteins)), ]
)

## export results
my.write.table(
  out,
  file.path(resultsdir, "combined-clin-pred-proteins.txt")
)
