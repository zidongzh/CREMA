# TF scan ==========


#' filter motifs by the selected TFs
#'
#' motifs names are stored as "TF_motifID"
#'
#' @param motifs
#' @param tfs
#'
#' @return
#' @export
#'
#' @examples
filter_motifs_by_tf_names <- function(motifs, tfs = NULL){
    if (!is.null(tfs)){
        motif_tfs <- sapply(names(motifs), function(x){strsplit(x, split="_")[[1]][1]})
        motifs <- motifs[motif_tfs %in% tfs]
    }
    if (length(motifs) == 0){
        stop("Selected TFs not in the motif database")
    }
    return(motifs)
}



#' Find the matching positions of the motifs
#'
#' A wrapper for the `motifmatchr::matchMotifs` function to find motif match (potential binding sites) locations in a given region
#'
#' @param regions a GRanges object of regions to search for motifs. Will be used as the `subject` argument for the `matchMotifs` function
#' @param motifs a list of position weight matrices. Will be use as the `pwm` argument for the `matchMotifs` function
#' @param genome a BSgenome object. Will be used as the `genome` argument for the `matchMotifs` function
#' @param p.cutoff the `p.cutoff` argument for the `matchMotifs` function
#' @param bg the `bg` argument for the `matchMotifs` function
#' @param region_extension whether to extend the region to search for motifs.
#' @param reduce (Not implemented yet) whether to combine / reduce binding sites for the same TF if they are close to each other
#' @param reduce_gap the gap allowed for reducing the binding sites.
#'
#' @return a Granges object of the positions of all the mapped motif (motif with at least one match). The GRanges object has metadata column for the TF names
#' @export
#'
#' @examples
scan_TF_positions <- function(regions,
                              motifs,
                              genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                              p.cutoff = 5e-5,
                              bg = "even",
                              region_extension = 0,
                              reduce = F,
                              reduce_gap = 0){

    # TF binding sites by motifmatchr
    # Return a list of GRanges object for each motif
    motif_ix <- motifmatchr::matchMotifs(pwms = motifs,
                                         subject = Extend_Granges(regions,
                                                                  upstream = region_extension,
                                                                  downstream = region_extension),
                                         out = "positions",
                                         genome = genome,
                                         p.cutoff = p.cutoff,
                                         bg = bg)

    # Keep motifs with at least one match
    motif_ix_each_length <- sapply(motif_ix, length)

    # Simplify the list of GRanges to one combined GRanges
    motif_ix_gr <- unlist(methods::as(motif_ix, "GRangesList"))

    if (length(motif_ix_gr) == 0){ return(NA) }

    # Add motif ID to the netadata of the combined GRanges object
    motif_ix_gr$ID <- rep(sapply(motifs, TFBSTools::ID), times = motif_ix_each_length)
    motif_ix_gr$Name <- rep(sapply(motifs, TFBSTools::name), times = motif_ix_each_length)
    motif_ix_gr$TF_processed <- rep(sapply(strsplit(names(motifs), split = "_"), function(x){x[1]}), times = motif_ix_each_length)
    motif_ix_gr$ID_processed <- rep(sapply(strsplit(names(motifs), split = "_"), function(x){paste0(x[-1], collapse = "_")}), times = motif_ix_each_length)

    # # Combine / Reduce binding sites for the same TF if they are close to each other
    # # Not Implemented Yet
    # if (reduce){}

    return(motif_ix_gr)
}







