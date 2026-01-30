## -------------------------------------------- ##
## extract-data.r
## -------------------------------------------- ##

args <- commandArgs(trailingOnly = TRUE)

datadir <- args[1]
resultsdir <- args[2]

## function for extracting tcga tar.gz's to named output
extract.file <- function(tar.file, extract.file, new.file, resultsdir) {
  # get file path to extracted file
  x.file <-
    grep(extract.file,
      untar(tar.file, list = T),
      value = T
    )
    
  # extract the tar file
  cat("Extracting", tar.file, "to", new.file, "\n")
  untar(tar.file, exdir=resultsdir, extras="--no-same-owner")
  x.file = file.path(resultsdir,x.file)

  # move the data to named output
  file.copy(x.file, new.file)

  # remove untared directory
  unlink(dirname(x.file), recursive = TRUE)
}

#######################
## extract the clinical data
clinical.file <- file.path(resultsdir, "clinical.txt")
if (!file.exists(clinical.file)) {
  extract.file(
    tar.file =
      file.path(
        datadir,
        grep(".*_HNSC\\..*_Clinical\\.Level_1\\..*\\.tar\\.gz$",
          list.files(datadir),
          value = T
        )
      ),
    extract.file = "HNSC.clin.merged.txt",
    new.file = clinical.file,
    resultsdir = resultsdir
  )
}


########################
## extract the protein data
protein.file <- file.path(resultsdir, "protein.txt")
if (!file.exists(protein.file)) {
  extract.file(
    tar.file =
      file.path(
        datadir,
        grep("*_protein_normalization__data.Level_3.*.tar.gz$",
          list.files(datadir),
          value = T
        )
      ),
    extract.file = "data.txt",
    new.file = protein.file,
    resultsdir = resultsdir
  )
}
## clean protein output:
## 	- remove 2nd row
lines <- readLines(protein.file)[-2]
writeLines(lines, file.path(resultsdir, "protein-clean.txt"))


########################
## methylation data is pre-extracted into the 'methylation-clean-score-sites.csv' file

## -------------------------------------------- ##
## clean-clinical.r 
## -------------------------------------------- ##

args <- commandArgs(trailingOnly = T)
datadir <- args[1]
resultsdir <- args[2]

clinical.filename <- file.path(resultsdir, "clinical.txt")
pan.cancer.filename <- file.path(
  resultsdir,
  "TCGA-CDR-SupplementalTableS1.txt"
)
output.filename <- file.path(resultsdir, "clinical-clean.txt")

cat(
  "extract-clinical.r",
  "\n ", clinical.filename,
  "\n ", pan.cancer.filename,
  "\n ", output.filename, "\n"
)

raw <- readLines(clinical.filename)
raw <- strsplit(raw, "\t")
raw <- sapply(raw, function(sample) sample)
colnames(raw) <- raw[1, ]
raw <- raw[-1, ]
raw <- as.data.frame(raw, stringsAsFactors = F)

clinical <- data.frame(
  participant = sub("[^-]+-[^-]+-", "", raw$patient.bcr_patient_barcode),
  stringsAsFactors = F
)
clinical$participant <- toupper(clinical$participant)

clinical$female <- raw$patient.gender == "female"
clinical$histology <- raw$patient.tumor_samples.tumor_sample.tumor_histologies.tumor_histology.histological_type
clinical$age.at.diagnosis <- as.numeric(raw$patient.age_at_initial_pathologic_diagnosis)
clinical$estrogen.receptor.status <- raw$patient.breast_carcinoma_estrogen_receptor_status
clinical$progesterone.receptor.status <- raw$patient.breast_carcinoma_progesterone_receptor_status
clinical$her2.status <- raw$patient.lab_proc_her2_neu_immunohistochemistry_receptor_status
clinical$ethnicity <- raw$patient.ethnicity
clinical$race <- raw$patient.race_list.race
clinical$positive.lymphnodes <- as.numeric(raw$patient.number_of_lymphnodes_positive_by_he)
clinical$stage <- raw$patient.stage_event.pathologic_stage
clinical$tnm.m.category <- raw$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
clinical$tnm.n.category <- raw$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
clinical$tnm.t.category <- raw$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t
clinical$lymphocyte.infiltration <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_lymphocyte_infiltration)
clinical$monocyte.infiltration <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_monocyte_infiltration)
clinical$neutrophil.infiltration <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_neutrophil_infiltration)
clinical$necrosis.percent <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_necrosis)
clinical$normal.cells.percent <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_normal_cells)
clinical$stromal.cells.percent <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_stromal_cells)
clinical$tumor.cells.percent <- as.numeric(raw$patient.samples.sample.portions.portion.slides.slide.percent_tumor_cells)

clinical$stage[clinical$stage == "stage x"] <- NA

clinical$tnm.m.category <- factor(
  as.character(clinical$tnm.m.category),
  levels = c("m0", "m1")
)

clinical$tnm.t.category[clinical$tnm.t.category == "tx"] <- NA
clinical$tnm.t.category[grepl("t1", clinical$tnm.t.category)] <- "t1"
clinical$tnm.t.category[grepl("t2", clinical$tnm.t.category)] <- "t2"
clinical$tnm.t.category[grepl("t3", clinical$tnm.t.category)] <- "t3"
clinical$tnm.t.category[grepl("t4", clinical$tnm.t.category)] <- "t4"

clinical$tnm.n.category[clinical$tnm.n.category == "nx"] <- NA
clinical$tnm.n.category[grepl("n0", clinical$tnm.n.category)] <- "n0"
clinical$tnm.n.category[grepl("n1", clinical$tnm.n.category)] <- "n1"
clinical$tnm.n.category[grepl("n2", clinical$tnm.n.category)] <- "n2"
clinical$tnm.n.category[grepl("n3", clinical$tnm.n.category)] <- "n3"


clinical.pan <- read.table(pan.cancer.filename, header = T, sep = "\t", stringsAsFactors = F)
clinical.pan <- clinical.pan[which(clinical.pan$type == "HNSC"), ]
clinical.pan$participant <- sub("[^-]+-[^-]+-", "", clinical.pan$bcr_patient_barcode)
clinical.pan <- clinical.pan[match(clinical$participant, clinical.pan$participant), ]
clinical$pfi <- clinical.pan$PFI
clinical$pfi.time <- clinical.pan$PFI.time
clinical$dfi <- clinical.pan$DFI
clinical$dfi.time <- clinical.pan$DFI.time

write.table(clinical, file = output.filename, row.names = F, col.names = T, sep = "\t")

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
## -------------------------------------------- ##
## analysis.r
## -------------------------------------------- ##
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = T)
datadir <- args[1]
resultsdir <- args[2]

combined.filename <- file.path(resultsdir, "combined-clin-pred-proteins.txt")
data <- read.table(combined.filename,
  header = T, sep = "\t", stringsAsFactors = F
)

protein.names <-
  subset(
    meffonym::meffonym.models(full = T),
    grepl("^episcores", filename)
  )$name
protein.names <- make.names(protein.names)

table(protein.names %in% colnames(data))

## Run glms for association between proteins levels and tissue type (tumor vs. normal)

## define glm formulae with pred.proteins as predictors of 'tumor.or.normal'
## tissue i.e. tumor.or.normal ~ pred.protein
formulae <- sapply(protein.names, function(i) {
  reformulate(i, response = "tumor")
}, simplify = F)

# run glms
fit <- sapply(formulae, function(i) {
  glm(i, data = data, family = binomial())
}, simplify = F)

fit.summary <- sapply(fit, function(i) {
  out <- summary(i)$coefficients
  out[, "Estimate"] <- out[, "Estimate"]
  out
}, simplify = F)

fit.coefs <- sapply(fit.summary, function(i) {
  i[2, c("Estimate", "Pr(>|z|)")]
}, simplify = F)
fit.coefs <- {
  x <- do.call(rbind, fit.coefs)
  data.frame(
    pred.protein = rownames(x),
    coef = x[, "Estimate"],
    p.value = x[, "Pr(>|z|)"]
  )
}

bonferroni <- -log10(0.05 / length(fit))

### Visualize results

fit.coefs |>
  ggplot(aes(x = pred.protein, y = -log10(p.value))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text_repel(
    data = fit.coefs[which((-log10(fit.coefs$p.value) > bonferroni)), ],
    aes(x = pred.protein, y = -log10(p.value), label = pred.protein)
  ) +
  geom_hline(
    yintercept = bonferroni,
    linetype = "dashed"
  )

fit.coefs |>
  ggplot(aes(x = coef, y = -log10(p.value))) +
  geom_point() +
  geom_text_repel(
    data = fit.coefs[which((-log10(fit.coefs$p.value) > bonferroni)), ],
    aes(x = coef, y = -log10(p.value), label = pred.protein)
  ) +
  geom_hline(
    yintercept = bonferroni,
    linetype = "dashed"
  )


## Run glms for association between proteins levels and progression free interval (PFI)

# This analysis should be restricted to measurements taken from tumor samples

tumor.data <- subset(data, tumor == 1)

## define glm formulae with pred.proteins as predictors of pfi
## tissue i.e. pfi ~ pred.protein
formulae <- sapply(protein.names, function(i) {
  reformulate(i, response = "pfi")
}, simplify = F)

# run glms
fit <- sapply(formulae, function(i) {
  glm(i, data = tumor.data, family = binomial())
}, simplify = F)

fit.summary <- sapply(fit, function(i) {
  out <- summary(i)$coefficients
  out[, "Estimate"] <- out[, "Estimate"]
  out
}, simplify = F)

fit.coefs <- sapply(fit.summary, function(i) {
  i[2, c("Estimate", "Pr(>|z|)")]
}, simplify = F)
fit.coefs <- {
  x <- do.call(rbind, fit.coefs)
  data.frame(
    pred.protein = rownames(x),
    coef = x[, "Estimate"],
    p.value = x[, "Pr(>|z|)"]
  )
}

### Visualize results

fit.coefs |>
  ggplot(aes(x = pred.protein, y = -log10(p.value))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text_repel(
    data = fit.coefs[which((-log10(fit.coefs$p.value) > bonferroni)), ],
    aes(x = pred.protein, y = -log10(p.value), label = pred.protein)
  ) +
  geom_hline(
    yintercept = bonferroni,
    linetype = "dashed"
  )

fit.coefs |>
  ggplot(aes(x = coef, y = -log10(p.value))) +
  geom_point() +
  geom_text_repel(
    data = fit.coefs[which((-log10(fit.coefs$p.value) > bonferroni)), ],
    aes(x = coef, y = -log10(p.value), label = pred.protein)
  ) +
  geom_hline(
    yintercept = bonferroni,
    linetype = "dashed"
  )
  
