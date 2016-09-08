library(readr,    quietly = T)
library(dplyr,    quietly = T)
library(cowplot,  quietly = T)
library(tidyr,    quietly = T)
library(magrittr, quietly = T)

chipgfp <- read_tsv("_data/results/coverage_heatmap/real_peak_summits_slop_5000bps_cov.bed",
                 col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

chipigg <- read_tsv("_data/results/coverage_heatmap/chipigg_summits_slop_5000bps_cov.bed",
                col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

inputgfp<- read_tsv("_data/results/coverage_heatmap/ip1_input_summits_slop_5000bps_cov.bed",
                       col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

inputigg <- read_tsv("_data/results/coverage_heatmap/ip2_input_summits_slop_5000bps_cov.bed",
                       col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

chipgfp$cov <- log2(chipgfp$cov + 1)
chipigg$cov <- log2(chipigg$cov + 1)

chipgfp$cov     <- (chipgfp$cov / 4127375) * 1e6
chipigg$cov     <- (chipigg$cov / 2948597) * 1e6
inputgfp$cov    <- (inputgfp$cov / 10295171) * 1e6
inputigg$cov    <- (inputigg$cov / 8330701) * 1e6

chipgfp$ratio <- log2((chipgfp$cov / inputgfp$cov) + 1)
chipigg$ratio <- log2((chipigg$cov / inputigg$cov) + 1)

get_region_order <- function(regions) {

  o <- regions %>%
    group_by(name) %>%
    summarise(avg = median(cov)) %>%
    arrange(desc(avg)) %$%
    name

  return(o)

}

plot_heatmap <- function(regions) {

  g <- regions %>%
    ggplot(aes(x = pos, y = as.factor(name))) +
    geom_tile(aes(fill=cov, width = 1.0, height = 1.0)) +
    #scale_fill_gradientn(colours = c("white", "red")) +
    scale_fill_gradientn(colours = c("blue", "black", "yellow"), limits = c(0,10), space = "rgb") +
    notheme +
    guides(fill = guide_legend(title = "Legend"))

  return(g)

}

notheme <- theme(axis.line=element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 panel.background=element_blank(),
                 panel.border=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 plot.background=element_blank())

functions <- read_tsv("_data/results/coverage_heatmap/gene_functions.tsv",
                      col_names = c("name", "func"))

# scaled <- cov %>%
#   dplyr::select(name, pos, cov) %>%
#   spread(name, cov)
#
# scaled <- scale(scaled[,2:(length(colnames(scaled)) - 1)], scale = F) %>%
#   as.data.frame()
#
# scaled$pos <- 1:1000
#
# scaled <- scaled %>%
#   gather(name, cov, -pos)
#
# order <- get_region_order(scaled)
# scaled <- transform(scaled , name = factor(name, levels = rev(order)))

# run one at a time
order <- get_region_order(chipgfp)
chipgfp <- transform(chipgfp, name = factor(name, levels = order))
chipigg <- transform(chipigg, name = factor(name, levels = order))

plot_heatmap(chipgfp)
plot_heatmap(chipigg)

