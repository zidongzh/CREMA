


#' Select candidate region for each gene
#'
#' This function selects the cnadidate region for each gene to run CREMA model on. Deefault is +/- 100kb around the TSS
#'
#' @param genes list of genes to run CREMA on
#' @param gene_body_gr a GRanges object of the genomic positions of each gene. The gene names are in the "gene_name" metadata
#' @param window_up search window in the upstream of the TSS
#' @param window_down search window in the downstream of the TSS
#' @param manual_truncate default FALSE. if the supplied gene_body_gr doesn't have seqinfo, need to manually truncate the regions so they don't surpass the boundaries of chromosomes
#' @param genome is manual_truncate = TRUE, supply the genome to get the length of each chromosome
#' @param gene_region_inclusion see details
#' @param self_region_inclusion default TRUE. Where to include the coding region of the gene itself.
#'
#' @return a list of GRanges objects named by the gene names. Each GRanges object correspond to the selected candidate region of the gene.
#' @export
#'
#' @examples
select_proximal_regions <- function(genes,
                                    gene_body_gr,
                                    window_up = 100000,
                                    window_down = window_up,
                                    manual_truncate = F,
                                    genome = NULL,
                                    gene_region_inclusion = c("all", "all_noncoding", "closest_noncoding"),
                                    self_region_inclusion = T){

    gene_region_inclusion <- match.arg(gene_region_inclusion)

    # Get the gene body regions of selected genes
    gene_regions_select_gr <- gene_body_gr[gene_body_gr$gene_name %in% genes]
    gene_regions_select_gr$TSS <- ifelse(test = (GenomicRanges::strand(gene_regions_select_gr) == "+" | GenomicRanges::strand(gene_regions_select_gr) == "*"),
                                         yes = GenomicRanges::start(gene_regions_select_gr), no = GenomicRanges::end(gene_regions_select_gr))

    gene_TSS_select_gr <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(gene_regions_select_gr),
                                  ranges = IRanges::IRanges(start = gene_regions_select_gr$TSS,
                                                   end = gene_regions_select_gr$TSS),
                                  strand = GenomicRanges::strand(gene_regions_select_gr),
                                  symbol = gene_regions_select_gr$symbol,
                                  TSS = gene_regions_select_gr$TSS,
                                  gene_strand = as.character(GenomicRanges::strand(gene_regions_select_gr)),
                                  seqinfo = GenomeInfoDb::seqinfo(gene_body_gr))

    # convert to list of ranges
    gene_regions_select_list <- split(gene_regions_select_gr, gene_regions_select_gr$symbol)
    gene_TSS_select_list <- split(gene_TSS_select_gr, gene_TSS_select_gr$symbol)


    # Obtain the candidate regulatory regions -----
    if (gene_region_inclusion %in% c("all", "all_noncoding")){

        # simply extend to the specified window
        reg_regions_select_gr <- Signac::Extend(gene_TSS_select_gr, upstream = window_up, downstream = window_down, from.midpoint = F)
        GenomicRanges::strand(reg_regions_select_gr) <- "*"

        # truncate the region by the chromosome length
        # Not necessary if gene_body_gr contains appropriate seqinfo metadata
        if (manual_truncate){
            # genome maximum length
            genome_lengths <- GenomeInfoDb::seqlengths(genome)
            # truncate the regulatory region
            GenomicRanges::start(reg_regions_select_gr) <- pmax(GenomicRanges::start(reg_regions_select_gr), 0)
            GenomicRanges::end(reg_regions_select_gr) <- pmin(GenomicRanges::end(reg_regions_select_gr),
                                                              genome_lengths[as.character(GenomicRanges::seqnames(reg_regions_select_gr))])
        }

        if (gene_region_inclusion == "all"){

            # all ---
            # simply extend to the specified window
            reg_regions_select_list <- split(reg_regions_select_gr, reg_regions_select_gr$symbol)

        }else if (gene_region_inclusion == "all_noncoding"){

            # all noncoding ---
            # extend to the specified window, but removing other gene regions from the specified window
            # i.e. only keep un-transcribed regions in this window

            # # Use setdiff to subtract gene regions
            # # (~ 10min for ~2500 genes)
            # reg_regions_select_list <- lapply(reg_regions_select_list, function(x){setdiff(x, gene_body_gr, ignore.strand = T)})

            # Use psetdiff to subtract gene regions
            # (faster than setdiff, ~30s for ~2500 genes)
            reg_regions_select_list <- GenomicRanges::psetdiff(reg_regions_select_gr,
                                                               rep(GenomicRanges::GRangesList(gene_body_gr), length(reg_regions_select_gr)),
                                                               ignore.strand = T)
            names(reg_regions_select_list) <- reg_regions_select_gr$symbol

            #TODO
            # # New function `subtract()` in GenomicRange_1.48, should be a lot faster
            # reg_regions_select_list <- subtract(reg_regions_select_list, gene_body_gr)


            # Add back the queried gene region themselves
            # Decide whether to include the self-region in the last step
            reg_regions_select_list <- lapply(reg_regions_select_list, function(x){GenomicRanges::reduce(c(x, gene_regions_select_list[[x$symbol[1]]]))})
        }

    }else if (gene_region_inclusion == "closest_noncoding"){

        # closest noncoding ---
        # extend until encountering the region of another gene, or reaching the specified window, whichever comes first
        # i.e. only selected the closest continuous un-transcribed region in the window

    }

    # Remove 0-length GRanges
    reg_regions_select_list <- reg_regions_select_list[sapply(reg_regions_select_list, length) > 0]

    # Whether to include the gene body of the queried gene itself into the regulatory regions of this gene
    if (!(self_region_inclusion)){
        reg_regions_select_list <- lapply(names(reg_regions_select_list),
                                          function(x){GenomicRanges::setdiff(reg_regions_select_list[[x]], gene_regions_select_list[[x]], ignore.strand = T)})
    }

    return(reg_regions_select_list)
}









