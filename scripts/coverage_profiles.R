library(readr,    quietly = T)
library(dplyr,    quietly = T)
library(cowplot,  quietly = T)
library(tidyr,    quietly = T)
library(magrittr, quietly = T)
library(stringr,  quietly = T)
library(seqinr,   quietly = T)

# import coverage data
cov_gfp <- read_tsv("_data/results/coverage_heatmap/real_peak_summits_slop_5000bps_cov.bed",
                    col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

cov_igg <- read_tsv("_data/results/coverage_heatmap/chipigg_summits_slop_5000bps_cov.bed",
                    col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

cov_input1 <- read_tsv("_data/results/coverage_heatmap/ip1_input_summits_slop_5000bps_cov.bed",
                       col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

cov_input2 <- read_tsv("_data/results/coverage_heatmap/ip2_input_summits_slop_5000bps_cov.bed",
                       col_names = c("seqnames", "start", "stop", "name", "score", "pos", "cov"))

# calculate mean of each position
cov_gfp_mean <- cov_gfp %>%
  group_by(pos) %>%
  summarize(m = mean(cov))

cov_igg_mean <- cov_igg %>%
  group_by(pos) %>%
  summarize(m = mean(cov))

cov_input1_mean <- cov_input1 %>%
  group_by(pos) %>%
  summarize(m = mean(cov))

cov_input2_mean <- cov_input2 %>%
  group_by(pos) %>%
  summarize(m = mean(cov))

# normalize by respective library sizes
lib_df <- data.frame(pos    = cov_gfp_mean$pos,
                 chipgfp    = cov_gfp_mean$m / (4127375) * 1e6,
                 chipigg    = cov_igg_mean$m / (2948597) * 1e6,
                 inputgfp   = cov_input1_mean$m / (10295171) * 1e6,
                 inputigg   = cov_input2_mean$m / (8330701) * 1e6)

lib_df$gfp_ratio <- lib_df$chipgfp / lib_df$inputgfp
lib_df$igg_ratio <- lib_df$chipigg / lib_df$inputigg

# or normalize by a scaling factor
#sca_df <- data.frame(pos    = cov_gfp_mean$pos,
#                 chipgfp    = cov_gfp_mean$m * (4127375 / 10295171),
#                 chipigg    = cov_igg_mean$m * (2948597 / 8330701),
#                 inputgfp   = cov_input1_mean$m * (4127375 / 10295171),
#                 inputigg   = cov_input2_mean$m * (2948597 / 8330701))

lib_df <- gather(lib_df, sample, cov, -pos)
#sca_df <- gather(sca_df, sample, cov, -pos)

lib_gfp <- filter(lib_df, sample == "gfp_ratio")
lib_igg <- filter(lib_df, sample == "igg_ratio")
lib <- bind_rows(lib_gfp, lib_igg)

#sca_gfp <- filter(sca_df, sample == "chipgfp" | sample == "inputgfp")
#sca_igg <- filter(sca_df, sample == "chipigg" | sample == "inputigg")

# read in fasta file
read_fasta <- function(file) {
  fasta <- as.data.frame(read.fasta(file = file,
                                    as.string = FALSE,
                                    forceDNAtolower = FALSE,
                                    seqonly = TRUE) %>%
                           str_split(""))
  colnames(fasta) <- 1:ncol(fasta)
  fasta$pos <- rownames(fasta)
  fasta <- fasta %>% gather(seq, base, -pos)
  return(fasta)
}

pwm <- fasta %>%
  group_by(as.numeric(pos)) %>%
  summarise(A = sum(base == "A")/n(), C = sum(base == "C")/n(), G = sum(base == "G")/n(), T = sum(base == "T")/n())

fasta <- read_fasta("_data/results/coverage_heatmap/real_peak_summits_slop_5000bps.fasta")
pwm <- pwm %>% mutate(GC = G+C, AT = A+T)
names(pwm)[1] <- "pos"

# calculate coverage and gc content tracks
#gfp_lib_plot <- ggplot(lib_gfp, aes(x = pos, y = cov, color = sample)) +
#  stat_smooth(method = "loess", se = F, span = 0.005, size = 1.0) +
#  geom_line(size = 1, alpha = 0.3) +
#  scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
#                     labels=c(-5000,-2500,0,2500,5000)) +
#  xlab("Position") +
#  ylab("Log2(Library Normalized Coverage + 1)") +
#  theme(legend.position="top") +
#  scale_color_manual(values = c("red", "black", "blue")) +
#  scale_y_continuous(limits=c(0.5, 1.6))

#igg_lib_plot <- ggplot(lib_igg, aes(x = pos, y = cov, color = sample)) +
#  stat_smooth(method = "loess", se = F, span = 0.005, size = 1.0) +
#  geom_line(size = 1, alpha = 0.3) +
#  scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
#                     labels=c(-5000,-2500,0,2500,5000)) +
#  xlab("Position") +
#  ylab("Log2(Library Normalized Coverage + 1)") +
#  theme(legend.position="top") +
#  scale_color_manual(values = c("blue", "black", "red")) +
#  scale_y_continuous(limits=c(0.5, 1.6))
#
lib_plot <- ggplot(lib, aes(x = pos, y = cov, color = sample)) +
  #stat_smooth(method = "loess", se = F, span = 0.001, size = 1.0) +
  geom_line(size = 1) +
  scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
                     labels=c(-5000,-2500,0,2500,5000)) +
  xlab("") +
  ylab("Fold Enrichment") +
  theme(legend.position="none") +
  scale_color_manual(values = c("red", "black")) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

gctrack <- ggplot(pwm, aes(x = pos, y = GC)) +
  stat_smooth(method = "loess", se = FALSE, span = 0.01, color="black") +
  geom_line(alpha = 0.25) +
  scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
                     labels=c(-5000,-2500,0,2500,5000)) +
  xlab("Position") +
  ylab("%GC") +
  ylim(c(0.05, 0.3))

# plot as a grid on top of one another
plot_grid(lib_plot, gctrack, nrow = 2, ncol = 1, scale = c(1, 1), rel_heights = c(2,1), align = "hv")
#g1_lib <- plot_grid(gfp_lib_plot, gctrack, nrow = 2, ncol = 1, scale = c(1, 1), rel_heights = c(2,1))
#g2_lib <- plot_grid(igg_lib_plot, gctrack, nrow = 2, ncol = 1, scale = c(1, 1), rel_heights = c(2,1))

#ggsave(filename = "library_size_gfp.png", plot = g1_lib)
#ggsave(filename = "library_size_igg.png", plot = g2_lib)


# gfp_sca_plot <- ggplot(sca_gfp, aes(x = pos, y = cov, color = sample)) +
#   stat_smooth(method = "loess", se = F, span = 0.005, size = 1.0) +
#   geom_line(size = 1, alpha = 0.3) +
#   scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
#                      labels=c(-5000,-2500,0,2500,5000)) +
#   xlab("Position") +
#   ylab("Log2(Scaled Coverage + 1)") +
#   theme(legend.position="top") +
#   scale_color_manual(values = c("red", "black")) +
#   scale_y_continuous(limits=c(1.5, 2.6))
#
# igg_sca_plot <- ggplot(sca_igg, aes(x = pos, y = cov, color = sample)) +
#   stat_smooth(method = "loess", se = F, span = 0.005, size = 1.0) +
#   geom_line(size = 1, alpha = 0.3) +
#   scale_x_continuous(breaks=c(0,2500,5000,7500,10000),
#                      labels=c(-5000,-2500,0,2500,5000)) +
#   xlab("Position") +
#   ylab("Log2(Scaled Coverage + 1)") +
#   theme(legend.position="top") +
#   scale_color_manual(values = c("blue", "black"))
#
# # plot as a grid on top of one another
# g1_sca <- plot_grid(gfp_sca_plot, gctrack, nrow = 2, ncol = 1, scale = c(1, 1), rel_heights = c(2,1))
# g2_sca <- plot_grid(igg_sca_plot, gctrack, nrow = 2, ncol = 1, scale = c(1, 1), rel_heights = c(2,1))
#
# ggsave(filename = "scaling_factor_gfp.png", plot = g1_sca)
# ggsave(filename = "scaling_factor_igg.png", plot = g2_sca)
