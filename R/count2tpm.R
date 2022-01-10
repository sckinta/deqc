#' normalize count to TPM
#'
#' @param count a matrix or data.frame with unique id as row names and unique sample name as colnames required
#' @param geneLength a named vector with unique id (required)
#' @param pseudoCount the integer pseudoCount added to count matrix/data.frame. default = 0
#' @param log a boolean value to determine whether the output tpm is Logarithm. default = FALSE. If is TRUE, use pseudoTPM.
#' @param pseudoTPM a integer pseudoCount added to TPM matrix/data.frame after normalization for Logarithm. default = 1
#'
#' @return a data.frame of normalized count as TPM
#' @importFrom stats setNames
#' @export
#'
#'
count2tpm <- function(count, geneLength, pseudoCount=0, log=F, pseudoTPM=1){
        # count is a matrix or data.frame with unique gene id as rownames
        # unique sample name as colnames. geneLength is a gene id named vector

        stopifnot(any(class(count) %in% c("data.frame","matrix")))
        stopifnot(all(rownames(count) %in% names(geneLength)))

        if (!"data.frame" %in% class(count)){
                count <- count %>%
                        as.data.frame()
        }

        samples <- colnames(count)

        count <- count + pseudoCount

        rpk <- count %>%
                dplyr::mutate_at(samples, ~.x/(geneLength[rownames(count)] / 1000))

        scale_factor <- apply(rpk, 2, sum) / 10^6

        tpm <- samples %>%
                setNames(samples) %>%
                purrr::map_dfc(~rpk[,.x]/scale_factor[.x]) %>%
                as.data.frame()

        attr(tpm, "row.names") <- attr(rpk, "row.names")

        if(log) {
                tpm <- tpm %>%
                        dplyr::mutate_at(samples, ~log2(.x+pseudoTPM))
        }

        tpm
}
