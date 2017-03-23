

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(minfi)
library(methyAnalysis)
library(sqldf)
library(qqman)

######
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Scripts/twins/data')
data_folder_2 <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder_2, '/model_data')
results_folder <- paste0(project_folder, '/Results')
cg_folder <- paste0(project_folder, '/Scripts/twins')

##########
# look at full_feat and bh results
##########
# read in full data
full_feat <- readRDS('/home/benbrew/Desktop/full_feat.rda')

# subset for reg, 48, 60, 72
feat_reg <- full_feat[full_feat$V3 == 'reg',]
feat_48 <- full_feat[full_feat$V3 == '48',]
feat_60 <- full_feat[full_feat$V3 == '60',]
feat_72 <- full_feat[full_feat$V3 == '72',]

##########
# load cg_locations
##########
# get cg_locations
cg_locations <- read.csv(paste0(cg_folder, '/cg_locations.csv'))

# rename to probe
names(cg_locations)[1] <- 'probe'

##########
# merge data
##########
reg_cg <- left_join(cg_locations,feat_reg, by = 'probe')
reg_48 <- left_join(cg_locations,feat_48, by = 'probe')
reg_60 <- left_join(cg_locations,feat_60, by = 'probe')
reg_72 <- left_join(cg_locations,feat_72, by = 'probe')

rm(list=ls(pattern="feat"))

##########
# group by chr and get counts
##########

# first remove NA 
reg_sub <- reg_cg[complete.cases(reg_cg),]

# group by chr
chr_counts <- reg_sub %>%
  group_by(seqnames) %>%
  summarise(counts= n())

# order by counts
chr_counts <- chr_counts[order(chr_counts$counts, decreasing = T),]

##########
# raw manhattan plot function
##########

getManData <- function(reg) 
{
  # remove V3
  reg$V3 <- NULL
  
  # fill NA with zero
  reg$counts[is.na(reg$counts)] <- 0
  
  ### get CHR, P, and BP columns as numeric
  reg <- reg[, c('seqnames', 'counts', 'start')]
  
  colnames(reg) <- c('CHR', 'P', 'BP')
  
  # recode CHR
  reg$CHR <- substr(reg$CHR, 4, 7)
  reg$CHR <- ifelse(grepl('X', reg$CHR), 23, 
                    ifelse(grepl('Y', reg$CHR), 24, reg$CHR))
  reg$CHR <- as.numeric(reg$CHR)
  reg$P <- jitter(reg$P, factor = 2, amount = NULL)
  
  
  return(reg)
}

reg_cg <- getManData(reg_cg)
reg_48 <- getManData(reg_48)
reg_60 <- getManData(reg_60)
reg_72 <- getManData(reg_72)


# recode CHR 


manPlot <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = adjustcolor(c("gray10", 
                                                                                          "gray60"), alpha.f = 0.6), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                     genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                     ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                                                           i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
                             ...))
  }
}


manPlot(reg_cg , 
          chr = "CHR", 
          bp = "BP",
          p = "P", 
          chrlabs = c(1:22, "X", "Y"), 
          logp = FALSE,suggestiveline = F, genomewideline = F,
          ylab = "# of times in model",
          ylim = c(0, 11),
          cex = 3,
          cex.lab=1.7, 
          cex.axis=1.7, 
          cex.main=1.7, 
          cex.sub=1.7)



# abline(a = seq(0, 10, 2))
grid(nx = 5, ny = 5, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)



##########
# generate random heatmaps
##########

#############################################################
# function to generate random proportions whose rowSums = 1 #
#############################################################
props <- function(ncol, nrow, var.names=NULL){
  if (ncol < 2) stop("ncol must be greater than 1")
  p <- function(n){
    y <- 0
    z <- sapply(seq_len(n-1), function(i) {
      x <- sample(seq(0, 1-y, by=.01), 1)
      y <<- y + x
      return(x)
    }
    )
    w <- c(z , 1-sum(z))
    return(w)
  }
  DF <- data.frame(t(replicate(nrow, p(n=ncol))))
  if (!is.null(var.names)) colnames(DF) <- var.names
  return(DF)
}
##############
# TRY IT OUT #
##############
tab <-as.matrix(props(ncol=5, nrow=5))


heatmap(tab, labRow = '', labCol = '', keep.dendro = FALSE, verbose = FALSE)
