#!/usr/bin/Rscript
# Calculate MNase periodicity and shift from TSS


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("--location"), default="NA", help="location file. at least column location exist"),
  make_option(c("--template"), default="sample_XXX.txt", help="format of bigwig result filen name. XXX is used for number"),
  make_option(c("--base_periodicity"), default="182", help="Predicted periodicity"),
  make_option(c("-o", "--out"), default="NA", help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))


FILE_location <- as.character(opt["location"])
FILE_bigwig_template <- as.character(opt["template"])
NUM_base_periodicity <- as.numeric(as.character(opt["base_periodicity"]))
FILE_out <- as.character(opt["out"])

#---------------------------------------------
# test data
#---------------------------------------------
# tmp <- Sys.info()
# if(as.character(tmp["sysname"]) == "Windows"){
#   GROBAL_ROOT <- "T:/"
# }else{
#   GROBAL_ROOT <- "/home/hidekit/"
# }
# PROJECT <- paste0(GROBAL_ROOT, "Project/009_20190923_Eric_MNaseSeq_estimateTSS/")
# DIR_bigwig <- paste0(PROJECT, "out/2020-05-23_average_MNase_graph/bigwig/")
# DIR_out <- paste0(PROJECT, "out/2020-05-23_average_MNase_graph/img/")
# FILE_location <- paste0(PROJECT, "out/2020-05-23_average_MNase_graph/region/location_TSS_10bp.txt")
# FILE_bigwig_template <- paste0(DIR_bigwig, "WT_N3753_MNase_TSS_XXX_10bp.bed")
# NUM_base_periodicity <- 182
#---------------------------------------------

D_location <- fread(FILE_location, header = TRUE) %>% as.data.frame()
D_location <- D_location %>% filter(tick >= -120, tick <= 1600)


#=============================================
# load bigwig result
#=============================================
loadBigwig <- function(i){
  file <- sub("XXX", i, FILE_bigwig_template)
  df <- fread(file, header = FALSE)
  df <- df[,c(4, 7)]
  colnames(df) <- c("gene", "score")
  df <- df %>% filter(!is.na(score))
  df <- df %>% mutate(location=i)
  df
}
D_table <- do.call(rbind, pblapply(D_location %>% pull(location), loadBigwig))
D_table <- dplyr::left_join(D_table, D_location, by="location")

#=============================================
# estimating periodicity
#=============================================
#---------------------------------------------
# Make template for calculation
#---------------------------------------------
Tick_array <- D_table %>% distinct(tick) %>% arrange(tick) %>% pull(tick)
Temp_peri <- seq(120, 240)
Temp_start <- seq(-120, 120, by=10)

Peri_template <- list()
pp <- 1
for(peri in Temp_peri){
  tmp_matrix <- NULL
  for(start in Temp_start){
    tmp <- seq(start, start + 1400, by=peri)
    tmp <- tmp[1:6]
    tmp <- as.integer(tmp/10) * 10
    tmp_matrix <- rbind(tmp_matrix, Tick_array %in% tmp)
  }
  Peri_template[[pp]] <- tmp_matrix
  pp <- pp+1
}
rm(tmp_matrix)


#---------------------------------------------
# Calculate score
#---------------------------------------------
moveAvg <- function(arr, window=3){
  getAve <- function(i){
    start <- max(1, i-(window-1)/2)
    end <- min(length(arr), i+(window-1)/2)
    mean(arr[start:end])
  }
  sapply(1:length(arr), getAve)
}

calcForGene <- function(g){
  tmp <- D_table %>% filter(gene==g) %>% select(tick, score) %>% arrange(tick) %>% pull(score)
  if(length(tmp) != length(Tick_array)){
    return
  }
  
  getScore <- function(i){
    apply(t(tmp * t(Peri_template[[i]])), 1, mean, na.rm = T)
  }
  
  ### row : periodicity, column : shift from TSS
  M_score <- do.call(rbind, lapply(1:length(Temp_peri), getScore))
  
  ### Periodicity
  tmp2 <- apply(M_score, 1, max, na.rm = T)
  tmp2 <- moveAvg(tmp2)
  score_center <- tmp2[Temp_peri==NUM_base_periodicity]
  
  ##### check upstream from average
  score_upstream <- score_center
  Periodicity_upstream <- NUM_base_periodicity
  for(s in seq(NUM_base_periodicity + 1, max(Temp_peri))){
    score_new <- tmp2[Temp_peri==s]
    if(score_center >= score_new ){
      break
    }
    score_upstream <- score_new
    Periodicity_upstream <- s
  }
  
  ##### check downstream from average
  score_downstream <- score_center
  Periodicity_downstream <- NUM_base_periodicity
  for(s in seq(NUM_base_periodicity - 1, min(Temp_peri), by=-1)){
    score_new <- tmp2[Temp_peri==s]
    if(score_center >= score_new ){
      break
    }
    score_downstream <- score_new
    Periodicity_downstream <- s
  }
  
  if(score_upstream > score_downstream){
    Periodicity <- Periodicity_upstream
  }else{
    Periodicity <- Periodicity_downstream
  }
  
  ### Shift from TSS (最大で周期/2の長さだけshift可能とする)
  df <- data.frame(shift=Temp_start, score=M_score[Temp_peri==Periodicity,])
  df <- df %>% filter(shift > Periodicity/2 * -1, shift < Periodicity/2)
  max_score <- df %>% pull(score) %>% max()
  Shift_from_TSS <- df %>% filter(score == max_score) %>% pull(shift)
  if(length(Shift_from_TSS) > 1){
    Shift_from_TSS <- Shift_from_TSS[abs(Shift_from_TSS) == min(abs(Shift_from_TSS))]
  }
  
  data.frame(gene=g, Periodicity=Periodicity, Shift=Shift_from_TSS, ratio=max_score/min(M_score), stringsAsFactors = FALSE)
}
D_periodicity <- do.call(rbind, pblapply(D_table %>% distinct(gene) %>% pull(gene), calcForGene))
write.table(D_periodicity, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


