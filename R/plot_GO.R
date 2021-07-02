#' GOPlot_barPlot
#'
#' @param result_GO
#' @param Clu_levels
#' @param title_name
#' @param addLine
#' @param space_width
#' @param wrap_length
#'
#' @return
#' @export
#'
#' @examples
#'
#' GO_file <- system.file("extdata", "GO.csv", package = "BioinfoPlotTK")
#' GO_info <- read.csv(GO_file)
#' p <- GOPlot_barPlot_flip(result_GO = GO_info)
#'
GOPlot_barPlot <- function(result_GO,
                           Clu_levels = NULL,
                           title_name = NULL,
                           addLine = FALSE,
                           space_width = 0.5,
                           wrap_length = 40) {

    result_GO %>%
        dplyr::distinct(.keep_all = TRUE) -> result_GO

    if (is.null(Clu_levels)){
        result_GO$Clu <- factor(result_GO$Clu, levels = unique(result_GO$Clu))
    } else {
        result_GO$Clu <- factor(result_GO$Clu, levels = Clu_levels)
    }

    result_GO %>%
        dplyr::arrange(Clu, padj) %>%
        dplyr::mutate(Description = factor(Description,
                                           levels = Description)) -> result_GO

    Clu_N <- length(table(result_GO$Clu))
    result_GO$add <- rep(seq_len(Clu_N) - 1, table(result_GO$Clu))
    result_GO$pos_x <- seq_len(nrow(result_GO)) + result_GO$add * space_width

    result_GO %>%
        ggplot(aes(x = pos_x, y = -log10(padj))) +
        geom_col(aes(fill = Clu)) +
        labs(x = " ", title = title_name) -> p

    p <- p + scale_x_continuous(breaks = result_GO$pos_x,
                                labels = function(x) str_wrap(result_GO$Description,
                                                              width = wrap_length),
                                guide = guide_axis(angle = 45))


    if (addLine) {

        line_vector <- add_groupLine(result_GO)
        tibble::tibble(x = line_vector,
                       xend = line_vector,
                       y = 0,
                       yend = max(-log10(result_GO$padj))) -> line_tb
        p <- p + geom_segment(data = line_tb,
                              aes(x = x, xend = xend,
                                  y = y, yend = yend),
                              size = 1)
    }

    return(p)

}


#' GOPlot_barPlot_flip
#'
#' @param result_GO
#' @param Clu_levels
#' @param labelInBar
#' @param title_name
#' @param space_width
#' @param wrap_length
#'
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
#'
#' GO_file <- system.file("extdata", "GO.csv", package = "BioinfoPlotTK")
#' GO_info <- read.csv(GO_file)
#' p <- GOPlot_barPlot_flip(result_GO = GO_info)
#'
#'
GOPlot_barPlot_flip <- function(result_GO,
                                Clu_levels = NULL,
                                title_name = NULL,
                                labelInBar = FALSE,
                                addLine = FALSE,
                                space_width = 0.5,
                                wrap_length = 40) {

    result_GO %>%
        dplyr::distinct(.keep_all = TRUE) -> result_GO

    if (is.null(Clu_levels)){
        result_GO$Clu <- factor(result_GO$Clu, levels = unique(result_GO$Clu))
    } else {
        result_GO$Clu <- factor(result_GO$Clu, levels = Clu_levels)
    }

    result_GO %>%
        dplyr::arrange(Clu, -padj) %>%
        dplyr::mutate(Description = factor(Description,
                                           levels = Description)) -> result_GO

    Clu_N <- length(table(result_GO$Clu))
    result_GO$add <- rep(seq_len(Clu_N) - 1, table(result_GO$Clu))
    result_GO$pos_x <- seq_len(nrow(result_GO)) + result_GO$add * space_width

    result_GO %>%
        ggplot(aes(x = pos_x, y = -log10(padj), label = Description, fill = Clu)) +
        geom_col() +
        guides(fill = guide_legend(reverse = TRUE)) +
        coord_flip() +
        labs(x = " ", title = title_name) -> p

    if (addLine) {
        p <- p + geom_vline(xintercept = add_groupLine(result_GO))
    }


    if (labelInBar){
        p <- p + geom_text(position = position_stack(vjust = 0), hjust = 0)
    } else {
        # https://stackoverflow.com/questions/21878974/wrap-long-axis-labels-via-labeller-label-wrap-in-ggplot2
        # https://stackoverflow.com/questions/5096538/customize-axis-labels
        p <- p + scale_x_continuous(breaks = result_GO$pos_x,
                                    labels = function(x) str_wrap(result_GO$Description,
                                                                  width = wrap_length))
    }


    return(p)

}




add_groupLine <- function(result_GO){
    result_GO %>%
        dplyr::group_by(Clu) %>%
        dplyr::summarise(max = max(pos_x)) %>%
        dplyr::pull(max) -> max_value

    result_GO %>%
        dplyr::group_by(Clu) %>%
        dplyr::summarise(min = min(pos_x)) %>%
        dplyr::pull(min) -> min_value

    line_vector <- vector(mode = "numeric", length(max_value) - 1)
    for (i in seq_len(length(line_vector))){
        line_vector[i] <- (max_value[i] + min_value[i + 1]) / 2
    }

    return(line_vector)
}
