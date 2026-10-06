# Prepare local caches for the published pancreas examples (single-cell Jamboree
# revision below; credit Peter Carbonetto and the study authors, see CREDITS.md).
# Explicit, network-using step run outside package builds and R CMD check. Upstream
# data are downloaded to a local directory and are not redistributed in the package.
# Usage: FLASHIER_UTILS_PANCREAS_DIR=<dir> Rscript .github/scripts/prepare-pancreas.R
library(Matrix)
pancreas_dir <- Sys.getenv("FLASHIER_UTILS_PANCREAS_DIR")
stopifnot(nzchar(pancreas_dir))
dir.create(pancreas_dir,recursive=TRUE,showWarnings=FALSE)
revision <- "6cf9aa720cc16dc543e93b36d57bc1a85555860d"
base_url <- paste0("https://raw.githubusercontent.com/stephenslab/",
                   "single-cell-jamboree/",revision,"/")
files <- c("data/pancreas.RData", "output/pancreas_celseq2_factors.RData")
for (file in files) {
  destination <- file.path(pancreas_dir,basename(file))
  if (!file.exists(destination))
    download.file(paste0(base_url,file),destination,mode="wb")
}
original <- new.env(parent=emptyenv())
coordinates <- new.env(parent=emptyenv())
load(file.path(pancreas_dir,"pancreas.RData"),envir=original)
load(file.path(pancreas_dir,"pancreas_celseq2_factors.RData"),envir=coordinates)
metadata <- original$sample_info
stopifnot(identical(rownames(original$counts),as.character(metadata$id)))
reference <- coordinates$fl_nmf_ldf
signed <- coordinates$fl_snmf_ldf
genes <- rownames(reference$F)
stopifnot(identical(rownames(signed$F),genes))
cell_types <- c("acinar","ductal","alpha","beta","delta","gamma")
select_cells <- function(technology,n) {
  unlist(lapply(cell_types,function(type)
    head(sort(as.character(metadata$id[metadata$tech==technology &
                                       metadata$celltype==type])),n)),use.names=FALSE)
}
reference_ids <- select_cells("celseq2",8)
query_ids <- select_cells("smartseq2",4)
# Match the published gene filter before calculating library sizes.
reference_totals <- Matrix::rowSums(original$counts[
  metadata$id[metadata$tech=="celseq2"],genes,drop=FALSE])
reference_library_mean <- mean(reference_totals)
shifted_log <- function(ids) {
  counts <- original$counts[ids,genes,drop=FALSE]
  totals <- Matrix::rowSums(counts)
  stopifnot(all(totals>0))
  Y <- Matrix::Diagonal(x=reference_library_mean/totals) %*% counts
  Y <- log1p(Y)
  dimnames(Y) <- list(ids,genes)
  Y
}
cache <- list(reference=reference,signed=signed,
  reference_expression=shifted_log(reference_ids),
  query_expression=shifted_log(query_ids),
  metadata=metadata[match(c(reference_ids,query_ids),metadata$id),],
  revision=revision,reference_library_mean=reference_library_mean,
  source_files=files,
  source_sha256=vapply(basename(files),function(file)
    digest::digest(file=file.path(pancreas_dir,file),algo="sha256"),character(1)))
saveRDS(cache,file.path(pancreas_dir,"pancreas-worked-example.rds"))

# Second cache: balanced cells from four technologies, raw counts on the published
# gene set, for cross-technology (cohort), holdout and embedding examples.
technologies <- c("celseq2","smartseq2","fluidigmc1","inDrop3")
n_per_group <- 40L
cohort_ids <- unlist(lapply(technologies,function(technology)
  unlist(lapply(cell_types,function(type)
    head(sort(as.character(metadata$id[metadata$tech==technology &
                                       metadata$celltype==type])),n_per_group)),
  use.names=FALSE)),use.names=FALSE)
cells <- list(counts=original$counts[cohort_ids,genes,drop=FALSE],
  metadata=metadata[match(cohort_ids,metadata$id),],
  reference_library_mean=reference_library_mean,
  technologies=technologies,cell_types=cell_types,n_per_group=n_per_group,
  revision=revision,source_files=files,source_sha256=cache$source_sha256)
saveRDS(cells,file.path(pancreas_dir,"pancreas-cells.rds"))
