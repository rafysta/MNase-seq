#!/usr/bin/Rscript
# calculate autocorrelation of target BED region


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--dyad"), default="NA", help="object file or dyadCov"),
  make_option(c("-i", "--in"), default="NA", help="BED file for target region, chr, start, end, name are required"),
  make_option(c("-o", "--out"), default="NA", help="output file. acf data for all target analysis, or periodicity for individual analysis"),
  make_option(c("--target_all"), default="TRUE", help="get ACF of all target (TRUE) or individual target (FALSE)"),
  make_option(c("--peaks"), default="NA", help="output file for periodicity and peak location information for all target analysis")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(IRanges)))


FILE_dyadCov <- as.character(opt["dyad"])
FILE_bed <- as.character(opt["in"])
Analysis_mode <- eval(parse(text=as.character(opt["target_all"])))
FILE_out <- as.character(opt["out"])
FILE_peaks <- as.character(opt["peaks"])


### test data
# FILE_bed <- "T:/Project/009_20190923_Eric_MNaseSeq_estimateTSS/out/2020-08-20_autocorrelation_comparison/BED/Upregulate_crf4_2_and_H3K27me.bed"
# FILE_dyadCov <- "T:/Project/009_20190923_Eric_MNaseSeq_estimateTSS/data/WT_N3753_MNase_dyadCov.rds"


O_dyadCov <- readRDS(FILE_dyadCov)
D_bed <- fread(FILE_bed, header=F)
D_bed <- D_bed[,1:4]
colnames(D_bed) <- c("chr", "start", "end", "name")

getScore <- function(i){
  score <- as.vector(unlist(O_dyadCov[[as.character(D_bed[i,"chr"])]]))
  df <- data.frame(score=score, loc=1:length(score))
  df <- df %>% filter(loc >= as.integer(D_bed[i,"start"])) %>% filter(loc <= as.integer(D_bed[i,"end"]))
  df %>% pull(score)
}

peaks<-function(series,span=100){ 
  z <- embed(series, span) 
  s <- span%/%2 
  v<- max.col(z) == 1 + s 
  result <- c(rep(FALSE,s),v) 
  result <- result[1:(length(result)-s)] 
  result 
} 

max.lag <- 1000
if(Analysis_mode){
  v <- do.call(c, pblapply(1:nrow(D_bed), getScore))
  ac <- acf(v, lag.max=max.lag, plot = FALSE)
  y <- as.vector(ac$acf)
  x <- 0:max.lag
  peaksPos <- which(peaks(y, span=100))[1:4]
  Peaks <- 1:length(peaksPos)
  fit <- lm(peaksPos ~ Peaks - 1)
  peri <- unname(round(fit$coefficients[1]))
  df <- data.frame(lag=x, acf=y)
  write.table(df, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  sink(FILE_peaks)
  cat("Periodicity =", peri, "\n")
  cat("Peaks = ", peaksPos, "\n")
  sink()
}else{
  D_peri <- do.call(rbind, pblapply(1:nrow(D_bed), function(i){
    v <- getScore(i)
    ac <- acf(v, lag.max=max.lag, plot = FALSE)
    y <- as.vector(ac$acf)
    x <- 0:max.lag
    peaksPos <- which(peaks(y, span=100))[1:4]
    Peaks <- 1:length(peaksPos)
    fit <- lm(peaksPos ~ Peaks - 1)
    peri <- unname(round(fit$coefficients[1]))
    df <- data.frame(name=as.character(D_bed[i,"name"]), periodicity=peri)
    df
  }))
  write.table(D_peri, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

