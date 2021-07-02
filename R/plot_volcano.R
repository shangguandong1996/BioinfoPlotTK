#' volcanoPlot_gene
#'
#' @param result_diffGene
#' @param Fold_cutoff
#' @param padj_cutoff
#' @param gene_select
#' @param color_palette
#' @param ylim_max
#' @param max.overlaps
#'
#' @return
#' @export
#'
#' @examples
volcanoPlot_gene <- function(result_diffGene,
                             Fold_cutoff = 1,
                             padj_cutoff = 0.05,
                             gene_select,
                             color_palette = c("#669cc7", "#a6a6a6", "#d73c31"),
                             ylim_max,
                             max.overlaps = 200){

    result_diffGene %>%
        mutate(padj = case_when(
            padj == 0 ~ 1e-323,
            is.na(padj) ~ 1,
            TRUE ~ padj
        )) -> result_diffGene

    result_diffGene %>%
        mutate(color = case_when(
            abs(log2FoldChange) <= Fold_cutoff | padj >= padj_cutoff ~ "noSig",
            log2FoldChange > Fold_cutoff & padj < padj_cutoff ~ "up",
            log2FoldChange < Fold_cutoff & padj < padj_cutoff ~ "down"
        ),
        point_size = 1) %>%
        mutate(point_size = case_when(
            abs(log2FoldChange) < Fold_cutoff + 1 | padj >= padj_cutoff ~ "A",
            # if B and C location switch, there will be no C
            abs(log2FoldChange) >= Fold_cutoff + 2 &  padj < padj_cutoff ~ "C",
            abs(log2FoldChange) >= Fold_cutoff + 1 & padj < padj_cutoff ~ "B"
        )) -> data_new

    data_new %>%
        dplyr::group_by(color) %>%
        dplyr::summarise(n = n()) %>%
        pull(n) -> groupN

    # if (groupN < 3){
    #     group
    # }

    ggplot(data_new, aes(x = log2FoldChange, y = -log10(padj))) +
        ggrastr::geom_point_rast(aes(color = color, size = point_size),
                                 alpha = 0.6) +
        scale_size_manual(values=c("A" = 1,
                                   "B" = 3,
                                   "C" = 5),
                          guide = F) +
        scale_color_manual(values = c("down" = color_palette[1],
                                      "noSig" = color_palette[2],
                                      "up" = color_palette[3]),
                           labels = paste0(c("down", "nosig", "up"),
                                           "(", groupN, ")"),
                           guide = guide_legend(override.aes = list(size = 5))) +
        theme(legend.key=element_blank()) +
        geom_vline(xintercept = c(-Fold_cutoff,Fold_cutoff),
                   color="grey40",
                   linetype="longdash", lwd = 0.5) +
        geom_hline(yintercept = -log10(padj_cutoff),color="grey40",
                   linetype="longdash", lwd = 0.5) -> p

    if (missing(ylim_max)) {
        ylim_max <- as.numeric(quantile(-log10(data_new$padj),0.95))
    }

    if (!missing(gene_select)) {
        gene_select %>%
            dplyr::left_join(data_new) %>%
            dplyr::mutate(symbol = dplyr::case_when(
                is.na(symbol) ~ gene_id,
                TRUE ~ symbol
            )) -> data_select

        if (max(-log10(data_select$padj)) > ylim_max) {
            ylim_max <- max(-log10(data_select$padj)) + 5
        }

        p +
            geom_point(data = data_select,
                       alpha = 1,
                       shape = 1) +
            ggrepel::geom_text_repel(data = data_select,
                                     aes(label = symbol),
                                     size = 3, box.padding = unit(0.35, "lines"),
                                     point.padding = unit(0.3, "lines"),
                                     max.overlaps = max.overlaps) -> p

    }

    p + ylim(0, ylim_max) -> p1

    return(p1)


}
