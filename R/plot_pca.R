#' plot PCA with given normalized count matrix and sample annotation
#'
#' @param tpm named count matrix/data.frame with column name representing sample names and row name representing gene name/id (required)
#' @param col_anno the column (sample) annotation data.frame with row name match the tpm matrix column name,
#' column name represents the factor to label for PCA plot (required)
#' @param color_col a column name from col_anno that be color aesthetics for PCA (required)
#' @param shape_col a column name from col_anno that be shape aesthetics for PCA. Default = NULL
#' @param text_col a column name from col_anno that be text label aesthetics for PCA. Default = NULL
#' @param point_size sample point size. Default = 5
#' @param plot_file save file to pdf. Default = ""
#'
#' @return a list of 4 plots, PCA matrix and PC variance summary. list(pc=pc, pca_summary=pca_summary, plots=plots)
#' @importFrom stats prcomp
#' @importFrom grDevices dev.off pdf
#' @import ggplot2 gridExtra
#' @export

plot_pca <- function(tpm, col_anno, color_col, shape_col=NULL, text_col=NULL, point_size=5, plot_file=""){

        m <- tpm[apply(tpm, 1, sd) > 0, ]

        pca_m <- prcomp(t(m),center = TRUE, scale. = TRUE)

        pc <- pca_m$x %>%
                as.data.frame() %>%
                tibble::rownames_to_column("sample") %>%
                dplyr::as_tibble()

        pc <- pc %>%
                dplyr::left_join(
                        col_anno %>%
                                tibble::rownames_to_column("sample")
                )
        pca_summary <- tibble::tibble(
                pc = paste0("PC",1:length(pca_m$sdev)),
                eigs = pca_m$sdev^2
        ) %>%
                dplyr::mutate(
                        sd=sqrt(eigs),
                        prop=eigs/sum(eigs)
                ) %>%
                dplyr::mutate(cumprop=cumsum(prop))

        pc3_labs <- purrr::map_chr(
                paste0("PC",1:3),
                ~pca_summary %>%
                        dplyr::filter(pc==.x) %>%
                        dplyr::pull(prop) %>%
                        scales::percent(accuracy=0.01)
        )
        pc3_labs <- purrr::map_chr(1:3, ~glue::glue("PC{.x} ({pc3_labs[.x]})"))
        names(pc3_labs) <- paste0("PC",1:3)

        .pca_ggplot <- function(x,y){
                p1 <- ggplot(pc, aes_string(x, y)) +
                        geom_point(aes_string(col=color_col, shape=shape_col), size=point_size) +
                        labs(x=pc3_labs[x], y=pc3_labs[y]) +
                        theme_bw()

                if(!is.null(text_col)){
                        p1 <- p1+
                                ggrepel::geom_text_repel(aes_string(label=text_col))
                }
                p1
        }

        plots <- tidyr::expand_grid(x=1:3,y=1:3) %>%
                dplyr::filter(x < y) %>%
                dplyr::mutate_all(~glue::glue("PC{.x}")) %>%
                dplyr::mutate(p = purrr::map2(x, y, ~.pca_ggplot(.x, .y))) %>%
                dplyr::pull(p)

        if(plot_file!=""){
                pdf(plot_file, width=21, height=7)
                do.call("grid.arrange", c(plots, ncol=3))
                dev.off()
        }

        var_plot <- pca_summary %>%
                dplyr::mutate(pc=factor(pc, levels=pc)) %>%
                ggplot(aes(pc,cumprop)) +
                geom_col(fill="yellow", color="yellow", alpha=0.5) +
                geom_hline(yintercept = 0.8, color="red") +
                labs(x="", y="cumulative proportion of variance") +
                theme_bw()

        plots <- c(plots,list(var_plot))

        list(pc=pc, pca_summary=pca_summary, plots=plots)
}
