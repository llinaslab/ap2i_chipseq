# OPTIONS #####################################################################

options(warn=-1)

library(optparse, quietly = T, warn.conflicts = F)
library(seqinr,   quietly = T, warn.conflicts = F)
library(cowplot,  quietly = T, warn.conflicts = F)
library(readr,    quietly = T, warn.conflicts = F)
library(tidyr,    quietly = T, warn.conflicts = F)
library(dplyr,    quietly = T, warn.conflicts = F)
library(stringr,  quietly = T, warn.conflicts = F)
library(seqLogo,  quietly = T, warn.conflicts = F)

option_list <- list(
  make_option(c("-p", "--peaks_file"), dest = "peaks_file",
              help="Peaks BED file with distances to summit"),
  make_option(c("-f", "--fasta_file"), dest="fasta_file",
              help="FASTA file centered on motif"),
  make_option(c("-m", "--heatmap"), dest="heatmap", default = "motif_heatmap.eps",
              help="Output heatmap file [%default]"),
  make_option(c("-l", "--logo"), dest="logo", default = "logo.eps",
              help="Output seqLogo file [%default]"),
  make_option(c("-d", "--distances"), dest="distances", default = "distances.eps",
              help="Output distances distribution [%default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fasta_file)) {
  print_help(opt_parser)
  stop("Need to input fasta file\n", call.=FALSE)
}
if (is.null(opt$peaks_file)) {
  print_help(opt_parser)
  stop("Need to input peaks file\n", call.=FALSE)
}

# FUNCTIONS ###################################################################

read_distances <- function(file) {
  v <- read_tsv(file, col_names = F)$X7
  return(v)
}

read_fasta <- function(file) {
  fasta <- as.data.frame(read.fasta(file = file, as.string = FALSE, forceDNAtolower = FALSE, seqonly = TRUE) %>% str_split(""))
  colnames(fasta) <- 1:ncol(fasta)
  fasta$pos <- rownames(fasta)
  fasta <- fasta %>% gather(seq, base, -pos)
  return(fasta)
}

plot_heatmap <- function(fasta, out = opt$heatmap) {
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
  
  myColors <- c("#E41A1C", "#377EB8", "#F0E442", "#4DAF4A")
  names(myColors) <- c("A", "C", "G", "T")
  levels(myColors) <- c("A", "C", "G", "T")
  fill_colors <- scale_fill_manual(name = "base", values = myColors)
  outline_colors <- scale_colour_manual(name = "base", values = myColors)
  
  g <- fasta %>% ggplot(aes(x = as.numeric(pos), y = rev(as.numeric(seq)), color=base)) + 
    geom_tile(aes(fill=base)) + 
    notheme + 
    fill_colors + 
    outline_colors
  
  ggsave(plot = g, filename = out)

}

plot_logo <- function(fasta, out = opt$logo) {
  
  positions <- c(11,19)
  
  pwm <- fasta %>% 
    group_by(as.numeric(pos)) %>% 
    summarise(A = sum(base == "A")/n(), C = sum(base == "C")/n(), G = sum(base == "G")/n(), T = sum(base == "T")/n())
  
  logo <- pwm[,2:5]
  
  setEPS()
  postscript(out, width = 12)
  seqLogo(makePWM((t(logo[positions[1]:positions[2],]))), ic.scale = F)
  dev.off()
}

plot_distances <- function(dist, out = opt$distances) {
  g <- ggplot(as.data.frame(dist), aes(x = dist)) + geom_density(aes(y = ..count..), size=1.0)
  ggsave(plot = g, filename = out)
}

# MAIN ########################################################################

fasta <- read_fasta(opt$fasta_file)
distances <- read_distances(opt$peaks_file)

plot_heatmap(fasta, opt$heatmap)
plot_logo(fasta, opt$logo)
plot_distances(distances, opt$distances)
