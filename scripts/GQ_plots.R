# Libraries ====
library(data.table)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Define Functions ====

# Manhattan ==== 
# by: Pagé Goddard - https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

# QQ ====

gg_qqplot <- function(ps, ci = 0.95, snps, hlight) {
  n  <- length(ps)
  df <- data.table(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)),
    is_annotate  = 'no',
    SNP      = snps
  )
  df <- df %>% 
    mutate(is_annotate=ifelse(SNP %in% hlight, "yes", "no"))
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  ggplot(df) +
    
    geom_point(aes(expected, observed), shape = 1, alpha = 1.0, size = 2, fill="#c51f5d",color="#c51f5d") +
    geom_abline(intercept = 0, slope = 1, alpha = 1.0,color="#243447") +
    geom_line(aes(expected, cupper), linetype = 2,color="#88898c") +
    geom_line(aes(expected, clower), linetype = 2,color="#88898c") +
    geom_ribbon(aes(x=expected,ymin=clower,ymax=cupper), fill="#88898c", alpha=0.5) +
    
    scale_fill_manual(values=c('#88898c','#88898c')) +
    
    xlab(log10Pe) +
    ylab(log10Po) +
    
    geom_label_repel(
      xlim = c(0, 0.5),
       color="#243447", 
       segment.color = "#243447",
      arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
      force = 4,
      data=df[df$is_annotate=="yes",],
      aes(expected,observed,label=as.factor(x=SNP))
      
    )
  
}

# Data ====

args <- commandArgs(trailingOnly=TRUE)
#DT <- fread('hapmap.fastGWA', sep = '\t')
DT <- fread(args[1], sep = '\t')
dir <- args[2]
prefix <- args[3]
DT <- na.omit(DT)
DTp <- DT[!(DT$P>5e-8),]
ord <- DT[order(P)] 

# Variables ====

mypalette <- c("#c51f5d","#243447","#88898c") # chr color palette
mysnps <- DTp[,SNP] # snps to highlight
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

# Run Function ====
man_title = paste(prefix,"Manhattan Plot", sep =" ")
man <- gg.manhattan(DT, threshold=5e-8, hlight=mysnps, col=mypalette, ylims=c(0,10), title=man_title)

man_save = paste(dir,prefix,"_manhattan_plot.png", sep ="")
ggsave(filename = man_save, plot = man, width = 420, height = 200, units = "mm")

qq <- gg_qqplot(ps=DT$P, ci = 0.95, snps = ord$SNP, hlight = mysnps)

qq_save = paste(dir,prefix,"_qq_plot.png", sep ="")
ggsave(filename = qq_save, plot = qq, width = 420, height = 200, units = "mm")


#end