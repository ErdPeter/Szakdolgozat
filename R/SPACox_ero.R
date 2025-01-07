## SPACox
library(SPACox)
library(survival)

rm(list = ls())
gc()

setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")
files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

geno_files <- files[grepl("geno", files)]
files <- files[!grepl("geno", files)]

for(h in 1:15){
  data <- readRDS(files[h])
  data2 <- readRDS(geno_files[h])
  
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
    
    for(i in 1:50){
      start_index <- (i - 1) * 100000 + 1
      end_index <- min(i * 100000, nrow(run))
      run_seq <- run[start_index:end_index, ]
      run_seq$ID <- paste0("S", 1:100000)
      
      geno <- data2[start_index:end_index, ]
      geno <- as.matrix(geno)
      rownames(geno) <- run_seq$ID
      
      run_seq$X2 <- as.factor(run_seq$X2)
      obj.null <- SPACox_Null_Model(Surv(time, event) ~ X1 + X2,
                                    data = run_seq,
                                    pIDs = run_seq$ID,
                                    gIDs = rownames(geno),
                                    length.out = 10000)
      SPACox.res <- SPACox(obj.null = obj.null, Geno.mtx = geno, G.model = "Add")
      SPACox.res <- as.data.frame(SPACox.res)
      result[i, "p"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa > 0.00005])
      result[i, "p5"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00005])
      result[i, "p8"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00000005])
      result[i, "r"] <- i
      
    }
    
    result$point <- k
    plot_result <- rbind(plot_result, result)
    
    stop <- Sys.time()
    print(stop - start)
    
  }
  
  filename <- files[h]
  filename <- gsub(".rds", "", filename)
  filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/spacox/", filename, "_spacox.rds")
  saveRDS(plot_result, file = filename, compress = "xz")
}