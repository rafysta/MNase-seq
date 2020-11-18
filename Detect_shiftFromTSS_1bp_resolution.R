#!/usr/bin/Rscript
# Calculate MNase shift from TSS in 1pb resolution

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="object file or dyadCov"),
  make_option(c("-o", "--out"), default="NA", help="output file shift bp for all gene"),
  make_option(c("--base_periodicity"), default="182", help="Periodicity")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(IRanges)))

### TSS location
# FILE_TSS <- "T:/Genome/data/neurospora_crassa/GeneInformation/Ensemble/TSS_list.txt"
FILE_TSS <- "/home/hidekit/Genome/data/neurospora_crassa/GeneInformation/Ensemble/TSS_list.txt"
D_tss <- fread(FILE_TSS, header=T)


FILE_dyadCov <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
NUM_base_periodicity <- as.numeric(as.character(opt["base_periodicity"]))

# ### test data
# FILE_dyadCov <- "T:/Project/009_20190923_Eric_MNaseSeq_estimateTSS/data/WT_N3753_MNase_dyadCov.rds"
# FILE_dyadCov <- "T:/Project/009_20190923_Eric_MNaseSeq_estimateTSS/data/WT_N3753_MNase_dyadCov.rds"

O_dyadCov <- readRDS(FILE_dyadCov)


peaks<-function(series,span=100){ 
  z <- embed(series, span) 
  s <- span%/%2 
  v<- max.col(z, "first") == 1 + s 
  result <- c(rep(FALSE,s),v) 
  result <- result[1:(length(result)-s)] 
  result 
} 
index <- seq(from=0, by=NUM_base_periodicity, length.out = 4)

D_shift <- NULL
for(cc in D_tss %>% distinct(chr) %>% pull(chr)){
  D_score <- data.frame(score=as.vector(unlist(O_dyadCov[[cc]])))
  D_score <- D_score %>% mutate(location = row_number())
  
  df <- do.call(rbind, pblapply(D_tss %>% filter(chr==cc) %>% pull(EnsembleID), function(gene){
    tss <- D_tss %>% filter(EnsembleID == gene) %>% pull(TSS)
    strand <- D_tss %>% filter(EnsembleID == gene) %>% pull(Strand)
    if(strand == 1){
      score <- D_score %>% filter(location >= tss - 100, location <= tss + 1000) %>% pull(score)
    }else{
      score <- D_score %>% filter(location >= tss - 1000, location <= tss + 100) %>% pull(score)
      score <- rev(score)
    }
    score.new <- sapply(1:600, function(i){
      mean(score[i + index])
    })
    shift_bp <- which(peaks(score.new, span=100))
    if(shift_bp[1] > 100){
      shift_bp <- shift_bp[1]
    }else{
      if(score.new[shift_bp[2]] > score.new[shift_bp[1]]){
        shift_bp <- shift_bp[2]
      }else{
        shift_bp <- shift_bp[1]
      }
    }
    shift_bp <- shift_bp - 100
    
    data.frame(EnsembleID=gene, chr=cc, TSS=tss, strand = strand, shift_bp=shift_bp)
  }))
  D_shift <- rbind(D_shift, df)
}

write.table(D_shift, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

