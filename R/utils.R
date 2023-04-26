
# This function saves a dgcMatrix in a 10X-like format
# It outputs 3 files int the output_dir:
# - matrix.mtx: the coo format matrix, read in R using readMM, read in python using scipy.io.mmread
# - features.tsv
# - barcodes.tsv
# The output is 10X-like, note it's not strictly 10X format because the features.tsv file only has one column of the gene name.
# If want to read back in, use "Read10X(..., gene.column = 1)" from Seurat, the default is "gene.column=2" because the features.tsv from 10X has multiple columns and gene symbol is usually the 2nd column
write_count_matrix_sparse <- function(exp_mtx, output_dir){
    dir.create(output_dir, showWarnings = F)

    genes <- row.names(exp_mtx)
    cells <- colnames(exp_mtx)

    filename <- file.path(output_dir, paste0("matrix.mtx"))
    Matrix::writeMM(obj = exp_mtx, file = filename)
    R.utils::gzip(filename, overwrite = T, remove = T)

    filename <- file.path(output_dir, "features.tsv")
    utils::write.table(data.frame(gene = genes), filename, sep = "\t", row.names = F, col.names = F, quote = F)
    R.utils::gzip(filename, overwrite = T, remove = T)

    filename <- file.path(output_dir, "barcodes.tsv")
    utils::write.table(data.frame(cell = cells), filename, sep = "\t", row.names = F, col.names = F, quote = F)
    R.utils::gzip(filename, overwrite = T, remove = T)
}



# Filter genes by number of cells detected and total UMI count
#' Title
#'
#' @param exp_mtx
#' @param proportion_cells_detected
#' @param min_cell_number
#' @param total_umi_cutoff
#' @param gene_names_remove_pattern
#'
#' @return
#' @export
#'
#' @examples
filter_exp_mtx_genes <- function(exp_mtx,
                                 proportion_cells_detected = 5e-3,
                                 min_cell_number = 3,
                                 total_umi_cutoff = "default",
                                 gene_names_remove_pattern = "^MT-"){

    # Filter for number of cells detected
    number_cells_detected <- max(ncol(exp_mtx) * proportion_cells_detected, min_cell_number)

    # Filter for total number of UMI
    if (is.null(total_umi_cutoff)){
        total_umi_cutoff <- 0
    }else if (is.numeric(total_umi_cutoff)){
        total_umi_cutoff <- total_umi_cutoff
    }else if (total_umi_cutoff == "default"){
        total_umi_cutoff <- ncol(exp_mtx) * proportion_cells_detected * 2
    }

    # Filter
    genes.select <- row.names(exp_mtx)[(Matrix::rowSums(exp_mtx > 0) >= number_cells_detected) &
                                           (Matrix::rowSums(exp_mtx) > total_umi_cutoff)]

    # Filter out genes by gene names
    genes.select <- genes.select[!grepl(gene_names_remove_pattern, genes.select)]

    return(exp_mtx[genes.select, ])
}


# columns are neighbor indices for each cell
#' Title
#'
#' @param data_multi
#' @param k_smooth
#' @param cells
#'
#' @return
#' @export
#'
#' @examples
neighbor_matrix <- function(data_multi, k_smooth, cells = colnames(data_multi)){
    temp <- data_multi@neighbors$weighted.nn@nn.idx[, 1:k_smooth, drop = F]
    temp_cells <- data_multi@neighbors$weighted.nn@cell.names
    neighbor_mtx <- Matrix::sparseMatrix(i = rep(seq(1, nrow(temp)), each = ncol(temp)),
                                 j = as.vector(t(temp)),
                                 dims = c(nrow(temp), nrow(temp)),
                                 dimnames = list(temp_cells, temp_cells))
    neighbor_mtx <- t(as.matrix(neighbor_mtx)[cells, cells, drop = F])
    return(neighbor_mtx)
}

# Extend Granges object
# Adapted from Signac::Extend
# Support extension length being a vector of the same length of the GRanges object
# Suppress warning for strand
#' Title
#'
#' @param x
#' @param upstream
#' @param downstream
#' @param from.midpoint
#'
#' @return
#' @export
#'
#' @examples
Extend_Granges <- function (x,
                            upstream = 0,
                            downstream = 0,
                            from.midpoint = FALSE){
    # if (any(strand(x = x) == "*")) {
    #     warning("'*' ranges were treated as '+'")
    # }
    on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
    if (from.midpoint) {
        midpoints <- GenomicRanges::start(x = x) + (GenomicRanges::width(x = x)/2)
        new_start <- midpoints - ifelse(test = on_plus, yes = upstream, no = downstream)
        new_end <- midpoints + ifelse(test = on_plus, yes = downstream, no = upstream)
    }
    else {
        new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
        new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
    }
    IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
    x <- GenomicRanges::trim(x = x)
    return(x)
}



#' Title
#'
#' @param x
#' @param y
#' @param sep
#' @param sep_x
#' @param sep_y
#' @param keep_name
#' @param invert
#'
#' @return
#' @export
#'
#' @examples
subsetByOverlaps_GRstr <- function(x, y, sep = c("-", "-"),
                                   sep_x = NULL, sep_y = NULL,
                                   keep_name = F,
                                   invert = F){

    if (length(x) == 0){ return(NULL) }
    if (length(y) == 0){ return(NULL) }
    if (length(sep_x) == 0){sep_x = sep}
    if (length(sep_y) == 0){sep_y = sep}
    if (keep_name){
        temp_name <- names(x)
        x <- Signac::StringToGRanges(x, sep = sep_x)
        x$name <- temp_name
    }else{
        x <- Signac::StringToGRanges(x, sep = sep_x)
    }
    y <- Signac::StringToGRanges(y, sep = sep_y)
    x <- GenomicRanges::subsetByOverlaps(x, y, invert = invert)

    if(length(x) > 0){
        if (keep_name){
            temp_name <- x$name
            x <- Signac::GRangesToString(x, sep = sep_x)
            names(x) <- temp_name
        }else{
            x <- Signac::GRangesToString(x, sep = sep_x)
        }

    }else{
        x <- character(0)
    }
    return(x)
}
