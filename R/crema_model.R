

#' Get access matrix for each motif
#'
#' For each TF, get the ATAC accessibility matrix of each motif loccations of the TF
#'
#' @param regions
#' @param motifs
#' @param fragment_object The fragment object of a `Seurat`/`Signac` object.
#' @param tfs
#' @param cells
#' @param site_extension
#' @param genome
#' @param p.cutoff
#' @param bg
#' @param region_extension
#'
#' @return
#' @export
#'
#' @examples
get_motif_access_mtx <- function(regions,
                                 motifs,
                                 fragment_object,
                                 tfs = NULL,
                                 cells = NULL,
                                 site_extension = 200,
                                 genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                                 p.cutoff = 5e-5,
                                 bg = "even",
                                 region_extension = 0){

    # Filter motifs by the selected TFs
    motifs <- filter_motifs_by_tf_names(motifs, tfs = tfs)

    # Motif site positions -----
    motif_ix_gr <- scan_TF_positions(regions = regions,
                                     motifs = motifs,
                                     genome = genome,
                                     p.cutoff = p.cutoff,
                                     bg = bg,
                                     region_extension = region_extension)

    # Set cells to calculate access matrix for
    if (is.null(cells)){
        cells <- unique(unlist(lapply(fragment_object, function(x){Signac::GetFragmentData(object = x, slot = "cells")})))
    }

    # Extend the motif range -----
    motif_ix_gr_extended <- Extend_Granges(motif_ix_gr,
                                           upstream = site_extension,
                                           downstream = site_extension,
                                           from.midpoint = T)

    # Access matrix -----
    motif_sites_access_mtx <- Signac::FeatureMatrix(fragments = fragment_object,
                                                    features = motif_ix_gr_extended,
                                                    sep = c("-","-"),
                                                    cells = cells,
                                                    process_n = 20000,
                                                    verbose = F)

    # Convert the rownames of the access matrix to un-extended ranges
    motif_ix_gr_str <- Signac::GRangesToString(motif_ix_gr, sep = c("-", "-"))
    row.names(motif_sites_access_mtx) <- motif_ix_gr_str

    # Split the access matrix by the TF names
    motif_ix_gr_str_list <- split(x = motif_ix_gr_str, f = motif_ix_gr$TF_processed, drop = F)
    motif_sites_access_mtx_list <- lapply(motif_ix_gr_str_list, function(x){motif_sites_access_mtx[x, , drop = F]})

    return(motif_sites_access_mtx_list)
}



#' ATAC weighted RNA matrix for the TFs
#'
#' @param target
#' @param TFs
#' @param regions_str
#' @param motifs
#' @param exp_mtx
#' @param fragment_object
#' @param cells
#' @param genome
#' @param motif_match_p_cutoff
#' @param site_extension
#' @param smooth
#' @param smooth_cutoff
#' @param neighbor_mtx
#'
#' @return
#' @export
#'
#' @examples
ATAC_weighted_RNA_matrix <- function(target,
                                     TFs,
                                     regions_str,
                                     motifs,
                                     exp_mtx,
                                     fragment_object,
                                     cells,
                                     genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                                     motif_match_p_cutoff = 5e-5,
                                     site_extension = 200,
                                     smooth = F, smooth_cutoff = 1, neighbor_mtx = NULL){

    # Get candidate tfs -----
    temp_tfs <- TFs
    temp_tfs <- temp_tfs[!is.na(temp_tfs)]
    temp_tfs <- intersect(temp_tfs, row.names(exp_mtx))
    temp_tfs <- setdiff(temp_tfs, target)
    if(length(temp_tfs) == 0){return(matrix(nrow = 0, ncol = 0))}

    # Select regions -----
    regions_str <- regions_str[!is.na(regions_str)]
    if (length(regions_str) == 0){ return(matrix(nrow = 0, ncol = 0)) }
    regions_gr <- Signac::StringToGRanges(regions_str, sep=c("-","-"))


    # Obtain match position of each TF's motifs -----
    # tf_access_mtx_list <- TF_positions_region(peaks = regions_gr,
    #                                           motifs = motifs,
    #                                           tfs = temp_tfs,
    #                                           fragment_object = fragment_object,
    #                                           cells = cells,
    #                                           genome = genome,
    #                                           site_extension = site_extension)
    tf_access_mtx_list <- get_motif_access_mtx(regions = regions_gr,
                                               motifs = motifs,
                                               fragment_object = fragment_object,
                                               tfs = temp_tfs,
                                               cells = cells,
                                               genome = genome,
                                               p.cutoff = motif_match_p_cutoff,
                                               site_extension = site_extension)
    tf_access_mtx_list <- tf_access_mtx_list[!is.na(tf_access_mtx_list)]
    if (length(tf_access_mtx_list) == 0){ return(matrix(nrow = 0, ncol = 0)) }

    # # only keep sites that's open in at least n cells
    # n_cells <- 10
    # tf_access_mtx_list <- lapply(tf_access_mtx_list,
    #                              function(x){x[Matrix::rowSums(x>0) > n_cells, , drop=F]})
    # tf_access_mtx_list <- tf_access_mtx_list[sapply(tf_access_mtx_list, nrow) > 0]
    # if (length(tf_access_mtx_list) == 0){ return(NA) }
    #
    # # Weigh TF expression by the site accessibility
    # tf_site_access_mtx <- do.call(rbind, tf_access_mtx_list)
    # tf_site_access_mtx <- (tf_site_access_mtx > 0) + 0
    # tf_mtx <- exp_mtx[rep(names(tf_access_mtx_list),
    #                       times = sapply(tf_access_mtx_list, nrow)), , drop = F]
    # tf_mtx <- tf_mtx * tf_site_access_mtx
    # row.names(tf_mtx) <- paste0(row.names(tf_mtx), "_", row.names(tf_site_access_mtx))

    # Construct binarized access matrix for each site -----
    tf_site_access_mtx <- do.call(rbind, tf_access_mtx_list)
    tf_site_access_mtx <- (tf_site_access_mtx > 0) + 0
    tf_site_access_mtx_row_tfs <- rep(names(tf_access_mtx_list), times = sapply(tf_access_mtx_list, nrow))

    # Smoothing (optional) -----
    if (smooth){
        tf_site_access_mtx <- (tf_site_access_mtx + tf_site_access_mtx %*% neighbor_mtx) >= smooth_cutoff
        tf_site_access_mtx <- tf_site_access_mtx + 0
    }

    # Filter sites: only keep sites that's open in at least n cells -----
    n_cells_open <- 10
    temp_inds <- Matrix::rowSums(tf_site_access_mtx) > n_cells_open
    if (sum(temp_inds) == 0) { return(matrix(nrow = 0, ncol = 0)) }
    tf_site_access_mtx <- tf_site_access_mtx[temp_inds, , drop = F]
    tf_site_access_mtx_row_tfs <- tf_site_access_mtx_row_tfs[temp_inds]

    # Filter sites: only keep sites that's open in at least n cells where the target gene is positive -----
    n_cells_pos <- 5
    cells_pos <- colnames(exp_mtx)[exp_mtx[target, ] > 0]
    temp_inds <- Matrix::rowSums(tf_site_access_mtx[, cells_pos, drop = F]) > n_cells_pos
    if (sum(temp_inds) == 0) { return(matrix(nrow = 0, ncol = 0)) }
    tf_site_access_mtx <- tf_site_access_mtx[temp_inds, , drop = F]
    tf_site_access_mtx_row_tfs <- tf_site_access_mtx_row_tfs[temp_inds]

    # Weigh TF expression by the site accessibility -----
    tf_mtx <- exp_mtx[tf_site_access_mtx_row_tfs, cells, drop = F]
    tf_mtx <- tf_mtx * tf_site_access_mtx
    row.names(tf_mtx) <- paste0(row.names(tf_mtx), "_", row.names(tf_site_access_mtx))

    # Remove duplicate row names -----
    # TFs having multiple motifs with the same location mapping
    tf_mtx <- tf_mtx[!duplicated(row.names(tf_mtx)), , drop = F]

    # Filter for informative TF-site combinations -----
    # tf_mtx <- tf_mtx[Matrix::rowSums(tf_mtx) > 0 & matrixStats::rowVars(tf_mtx) > 0, , drop = F]
    tf_mtx <- tf_mtx[Matrix::rowSums(tf_mtx) > 0, , drop = F]
    if (nrow(tf_mtx) == 0) { return(matrix(nrow = 0, ncol = 0)) }
    tf_mtx <- as.matrix(tf_mtx)

    return(tf_mtx)
}


#' Title
#'
#' @param x
#' @param y
#' @param cell_grouping
#' @param cell_grouping_clusters
#' @param regression_method
#'
#' @return
#' @export
#'
#' @examples
regression_model <- function(x, y,
                             cell_grouping = F, cell_grouping_clusters = NULL,
                             regression_method = c("ols", "poisson")){
    if (cell_grouping){
        cells_keep_ind <- Reduce("|", list(!duplicated(cell_grouping_clusters),
                                           !duplicated(x),
                                           !duplicated(y)))
        y <- y[cells_keep_ind, , drop = F]
        x <- x[cells_keep_ind, , drop = F]
    }

    tryCatch({
        if (regression_method == "ols"){
            temp_lm <- summary(stats::lm(y ~ x))
        }else if (regression_method == "poisson"){
            temp_lm <- summary(stats::glm(y ~ x, family = stats::poisson(link = "log")))
        }
        temp_lm <- as.data.frame(stats::coef(temp_lm))
        return(temp_lm[2, , drop = F])
    },
    error = function(e){return(NULL)})
}



#' CREMA model for selected gene
#'
#' Run the crema model on one gene.
#'
#' @param target the target gene of the regulatory circuit.
#' @param TFs TFs to include in the modeling of the regulatory circuit.
#' @param regions_str the genomic region to search for motifs in a "chr-start-end" string.
#' @param motifs the motifs object provided from the CREMA package.
#' @param exp_mtx the RNA count matrix
#' @param fragment_object the fragments object as described in the `Signac` package. Ideally, use the processed fragments object of the cutsites instead of the ATAC fragments.
#' @param cells cells to include in the modeling. If `NULL` (default), use all the cells in the fragments object.
#' @param genome the BSgenome object of the the corresponding species.
#' @param motif_match_p_cutoff the p value cutoff for matching motifs, i.e. the `p.cutoff` argument for the `motifmatchr::matchMotifs` function.
#' @param site_extension the extension window surrounding the motif site for calculating chromatin accessibility (ATAC Tn5 insert sites)
#' @param regression_method the regression method for the linear model. Currently supports "ols" (ordinary least squares) and "poisson" (Poisson regression). Note that Poisson regression is much slower.
#' @param smooth (Not fully tested yet) whether to smooth the accessibility matrix of the TF binding sites.
#' @param smooth_cutoff (Not fully tested yet) if `smooth` set to `TRUE`, the cutoff to binarize the accessibility matrix of the TF binding sites.
#' @param neighbor_mtx (Not fully tested yet) if `smooth` set to `TRUE`, the neighbor matrix of cells for smoothing the accessibility matrix of the TF binding sites.
#' @param cell_grouping (Not fully tested yet) whether to group cells together into "meta cells" for modeling. May improve speed without sacrificing on accuracy.
#' @param cell_grouping_clusters (Not fully tested yet) is `cell_grouping` set to `TRUE`, a vector (the same length of the number of cells) of grouping memberships of the cells.
#' @param return_val the format of the return value. See below for details.
#'
#' @return
#' * list: a list containing "regression_method", "site_extension", a vector of p values and a vector of t values of all the TF-site regressors.
#' * vector: a vector of p values and t values of all the TF-site regressors.
#' * df: a dataframe containing the p values and t values of all the TF-site regressors.
#'
#' @export
#'
#' @examples
#'
ATAC_weighted_tf_model_highres <- function(target,
                                           TFs,
                                           regions_str,
                                           motifs,
                                           exp_mtx,
                                           fragment_object,
                                           cells = NULL,
                                           genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                                           motif_match_p_cutoff = 5e-5,
                                           site_extension = 200,
                                           return_val = c("list", "vector", "df"),
                                           regression_method = c("ols", "poisson"),
                                           smooth = F, smooth_cutoff = 1, neighbor_mtx = NULL,
                                           cell_grouping = F, cell_grouping_clusters = NULL){

    return_val <- match.arg(return_val)
    regression_method <- match.arg(regression_method)

    # Get ATAC weighted TF matrix -----
    tf_mtx <- ATAC_weighted_RNA_matrix(target = target,
                                       TFs = TFs,
                                       regions_str = regions_str,
                                       motifs = motifs,
                                       exp_mtx = exp_mtx,
                                       fragment_object = fragment_object,
                                       cells = cells,
                                       genome = genome,
                                       motif_match_p_cutoff = motif_match_p_cutoff,
                                       site_extension = site_extension,
                                       smooth = smooth, smooth_cutoff = smooth_cutoff, neighbor_mtx = neighbor_mtx)
    if (nrow(tf_mtx) == 0){ return(NA) }

    # Model -----
    # if (regression_method == "ols"){
    #     tryCatch({
    #         temp_lm <- sapply(row.names(tf_mtx),
    #                           function(x){summary(lm(t(exp_mtx[target, , drop = F]) ~
    #                                                      t(tf_mtx[x, , drop = F])))},
    #                           USE.NAMES = T, simplify = F)
    #         temp_lm <- lapply(temp_lm, function(x){as.data.frame(coef(x))})
    #         temp_p <- sapply(temp_lm, function(x){x[2, "Pr(>|t|)"]})
    #         temp_t <- sapply(temp_lm, function(x){x[2, "t value"]})
    #
    #         if (return_val == "vector"){
    #             # return a named vector for each TF_site combination's p value
    #             return(temp_p)
    #         }else if (return_val == "list"){
    #             # return a list named by TFs, with each element being a named vector of this TF's p value at each site
    #             # message("Returning a list not Implemented yet")
    #             return(list(regression_method = regression_method,
    #                         site_extension = site_extension,
    #                         p = temp_p,
    #                         t = temp_t))
    #         }
    #     },
    #     error = function(e){return(NA)})
    #
    # }else if (regression_method == "poisson"){
    #     tryCatch({
    #         temp_lm <- sapply(row.names(tf_mtx),
    #                           function(x){summary(glm(t(exp_mtx[target, , drop = F]) ~
    #                                                       t(tf_mtx[x, , drop = F]),
    #                                                   family = poisson(link = "log")))},
    #                           USE.NAMES = T, simplify = F)
    #         temp_lm <- lapply(temp_lm, function(x){as.data.frame(coef(x))})
    #         temp_p <- sapply(temp_lm, function(x){x[2, "Pr(>|z|)"]})
    #         temp_t <- sapply(temp_lm, function(x){x[2, "z value"]})
    #
    #         if (return_val == "vector"){
    #             return(temp_p)
    #         }else if (return_val == "list"){
    #             return(list(regression_method = regression_method,
    #                         site_extension = site_extension,
    #                         p = temp_p,
    #                         t = temp_t))
    #         }
    #     },
    #     error = function(e){return(NA)})
    # }

    tryCatch({
        temp_lm <- sapply(row.names(tf_mtx),
                          function(x){regression_model(x = t(tf_mtx[x, , drop = F]), y = t(exp_mtx[target, , drop = F]),
                                                       cell_grouping = cell_grouping, cell_grouping_clusters = cell_grouping_clusters,
                                                       regression_method = regression_method)},
                          USE.NAMES = T, simplify = F)
        temp_regressor_names <- names(temp_lm)
        temp_lm <- do.call("rbind", temp_lm)

        if (return_val == "df"){
            return(temp_lm)
        }

        if (regression_method == "ols"){
            temp_p <- temp_lm$`Pr(>|t|)`
            temp_t <- temp_lm$`t value`
        }else if (regression_method == "poisson"){
            temp_p <- temp_lm$`Pr(>|z|)`
            temp_t <- temp_lm$`z value`
        }
        names(temp_p) <- temp_regressor_names
        names(temp_t) <- temp_regressor_names

        if (return_val == "vector"){
            return(temp_p)
        }else if (return_val == "list"){
            return(list(regression_method = regression_method,
                        site_extension = site_extension,
                        p = temp_p,
                        t = temp_t))
        }
    },
    error = function(e){return(NA)})
}




