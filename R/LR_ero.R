## LR teszt
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
      rownames(geno) <- run_seq$ID
      
      run_seq$X2 <- as.factor(run_seq$X2)
      
      fit0 <- coxph(Surv(time, event) ~ X1 + X2,  
                    data = run_seq)
      pek <- c()
      for(g in 1:ncol(geno)){
        gcol <- colnames(geno)[g]
        fmla <- as.formula(paste0("Surv(time, event) ~ X1 + X2 + ", gcol))
        run_seq_tmp <- cbind(run_seq, geno[gcol])
        fit <- coxph(fmla,  
                     data = run_seq_tmp)
        pek[g] <- anova(fit0, fit, test = "Chisq")$`Pr(>|Chi|)`[2]
      }
      
      result[i, "p"] <- sum(pek > 0.00005)
      result[i, "p5"] <- sum(pek <= 0.00005)
      result[i, "p8"] <- sum(pek <= 0.00000005)
      result[i, "r"] <- i
      
    }
    
    result$point <- k
    plot_result <- rbind(plot_result, result)
    
    stop <- Sys.time()
    print(stop - start)
    
  }
  
  filename <- files[h]
  filename <- gsub(".rds", "", filename)
  filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/LR/", filename, "_LR.rds")
  saveRDS(plot_result, file = filename, compress = "xz")
}
