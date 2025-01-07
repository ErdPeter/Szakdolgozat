## Geno info
set.seed(1234)
mut <- function(x){
  list <- setdiff(c("A", "T", "G", "C"), x)
  return(sample(list, 1))
  
}

geno_info <- data.frame(id = paste0("G", 1:20),
                        chr = sample(as.character(1:22), 20, replace = T),
                        dist = 0,
                        pos = round(runif(20, 1, 99999999), 0),
                        A1 = sample(c("A", "T", "G", "C"), 20, replace = T))
geno_info$A2 <- unlist(lapply(geno_info$A1, mut))

saveRDS(geno_info,
        file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/geno_info.rds"),
        compress = "xz")


## Bed filek
set.seed(1234)
setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")

rm(list = ls())
gc()

geno_info <- readRDS("geno_info.rds")

fam <- data.frame(famid = NA,
                  id = paste0("S", 1:100000),
                  father = NA,
                  mother = NA,
                  sex = NA,
                  pheno = NA)

bim <- data.frame(chr = as.character(geno_info$chr),
                  id = geno_info$id,
                  dist = 0,
                  pos = geno_info$pos,
                  A1 = geno_info$A1,
                  A2 = geno_info$A2)

files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "_geno")

maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)

for(g in 1:15){
  data2 <- readRDS(files[g])
  print(paste0(g, ". beolvasas kész!"))
  
  maf_act <- maf[g]
  lambda_act <- lambda[g]
  
  for(i in 1:50){
    if(i%%10 == 0){print(i)}
    start_index <- (i - 1) * 100000 + 1
    end_index <- min(i * 100000, nrow(data2))
    
    geno <- as.matrix(data2[start_index:end_index,])
    rownames(geno) <- paste0("S", 1:100000)
    
    bed <- gaston::as.bed.matrix(geno, fam, bim)
    write.bed.matrix(bed,
                     basename = paste("c:/Tanulas/Szakdolgozat/Ero2/Data/BED/bed_maf",
                                      maf_act,
                                      "lambda",
                                      lambda_act,
                                      i,
                                      sep = "_"))
  }
}


## GWASSURVIVR
library(gwasurvivr)

setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")

rm(list = ls())
gc()

maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)

bedek <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data/BED",
                    full.names = F,
                    pattern = ".bed")
bedek <- sort(bedek)
bedek <- paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/BED/", bedek)

files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

files <- files[!grepl("geno", files)]


for(h in 1:15){
  filesorszam <- h
  data <- readRDS(paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/", files[filesorszam]))
  
  rows <- seq(0, nrow(data), by = 100000*50)
  rows1 <- rows + 1
  
  plot_result <- data.frame(p = numeric(0),
                            p5 = numeric(0),
                            p8 = numeric(0),
                            r = numeric(0),
                            point = numeric(0))
  
  for(k in 1:11){
    start <- Sys.time()
    
    gc()
    run <- data[rows1[k]:rows[k + 1], ]
    
    result <- data.frame(p = numeric(0),
                         p5 = numeric(0),
                         p8 = numeric(0),
                         r = numeric(0))
    
    bedek_tmp <- bedek[((h - 1) * 50 + 1):(min(h * 50, nrow(bedek)))]
    bedek_tmp <- bedek_tmp[c(1, 12, 23, 34, 45, 47:50, 2:11, 13:22, 24:33, 35:44, 46)]
    
    for(i in 1:50){
      start_index <- (i - 1) * 100000 + 1
      end_index <- min(i * 100000, nrow(run))
      run_seq <- run[start_index:end_index, ]
      run_seq$ID <- paste0("S", 1:100000)
      bedfile <- bedek_tmp[i]
      
      gwasurvivr::plinkCoxSurv(bed.file = bedfile,
                               covariate.file = run_seq[, c(1:2, 4:6)],
                               id.column = "ID",
                               sample.ids = run_seq$ID,
                               time.to.event = "time",
                               event = "event",
                               covariates = c("X1", "X2"),
                               print.covs = "only",
                               out.file = paste("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2/GWAS_result", maf[h], lambda[h], k, i, sep = "_"),
                               chunk.size = 20,
                               maf.filter = NULL)
      
      gwas <- data.table::fread(paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2/GWAS_result_", maf[h], "_", lambda[h], "_", k, "_", i, ".coxph"))
      result[i, "p"] <- length(gwas$PVALUE[gwas$PVALUE > 0.00005])
      result[i, "p5"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00005])
      result[i, "p8"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00000005])
      result[i, "r"] <- i
      
      print(paste0("diagram: ", h, "., pont: ", k, "., ciklus szám: ", i))
      
    }
    
    result$point <- k
    plot_result <- rbind(plot_result, result)
    
    save(result, file = "c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2_res/TEMP.RData")
    
    stop <- Sys.time()
    print(stop - start)
    
  }
  
  rm(filename)
  filename <- files[filesorszam]
  filename <- gsub(".rds", "", filename)
  filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2_res/", filename, "_gwasurvivr.rds")
  saveRDS(plot_result, file = filename, compress = "xz")
}
