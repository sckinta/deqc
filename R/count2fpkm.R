#' normalize count to FPKM
#'
#' @param count a matrix or data.frame with unique id as row names and unique sample name as colnames required
#' @param geneLength a named vector with unique id (required)
#' @param pseudoCount the integer pseudoCount added to count matrix/data.frame. default = 0
#' @param log a boolean value to determine whether the output fpkm is Logarithm. default = FALSE. If is TRUE, use pseudoFPKM.
#' @param pseudoFPKM a integer pseudoCount added to FPKM matrix/data.frame after normalization for Logarithm. default = 1
#'
#' @return a data.frame of normalized count as FPKM
#' @importFrom stats setNames
#' @export
#'
#'
count2fpkm <- function(count, geneLength, pseudoCount=0, log=F, pseudoFPKM=1){
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

        scale_factor <- apply(count, 2, sum)/10^6

        rpm <- samples %>%
                setNames(samples) %>%
                purrr::map_dfc(~count[,.x]/scale_factor[.x]) %>%
                as.data.frame()

        attr(rpm, "row.names") <- attr(count, "row.names")


        rpkm <- rpm %>%
                dplyr::mutate_at(samples, ~.x/(geneLength[rownames(rpm)] / 1000))

        if(log){
                rpkm <- rpkm %>%
                        dplyr::mutate_at(samples, ~log2(.x+pseudoRPKM))
        }
        rpkm
}
