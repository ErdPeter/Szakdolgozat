library(LTFHPlus)
library(dplyr)
library(MASS)
library(future.apply)
library(tidyr)
library(survival)
library(bigstatsr)


setwd("c:/Tanulas/Szakdolgozat/4modszer/Adult")

nthreads <- 4
plan(tweak(multisession, workers = nthreads))


files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

geno_files <- files[grepl("geno", files)]
files <- files[!grepl("geno", files)]
h2 <- 0.5

for(h in 1:15){
  data <- readRDS(files[h])
  data2 <- readRDS(geno_files[h])
  
  data$time <- round(data$time, 7)
  data$aoo <- ifelse(data$event == 1, data$time, NA)
  data$role <- "g"
  gc()
  data$fid <- as.character(1:nrow(data))
  data$pid <- as.character(1:nrow(data))
  
  covariates <- as.matrix(data[, c("X1", "X2")])
  
  rows <- seq(0, nrow(data), by = 100000*50)
  rows1 <- rows + 1
  

  plot_result <- data.frame(p = numeric(0),
                            p5 = numeric(0),
                            p8 = numeric(0),
                            r = numeric(0),
                            point = numeric(0))
  
  full_p <- data.frame(full_p = numeric(0),
                       point = numeric(0))



for(k in 1:11){
  gc()
  run <- data[rows1[k]:rows[k + 1], ]
  cov <- covariates[rows1[k]:rows[k + 1], ]
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  full_p_tmp <- c()
  
  
  for(i in 1:50){
    start_index <- (i - 1) * 10000 + 1
    end_index <- min(i * 10000, nrow(run))
    run_seq <- run[start_index:end_index, ]
    cov_seq <- cov[start_index:end_index, ]
    
    run_seq$time <- round(run_seq$time, 7)
    run_seq$aoo <- round(run_seq$aoo, 7)
    
    # CIP meghatározás
    run_seq$seged <- 1
    
    cip <- run_seq %>%
      arrange(aoo) %>% 
      group_by(aoo) %>%
      summarise(event = sum(event),
                seged = sum(seged))
    
    cip$cip <- cip$event / nrow(run_seq)
    cip$cip <- cumsum(cip$cip)
    cip$cip <- round(cip$cip, 7)
    
    cip$cip <- ifelse(cip$cip == 0,
                      0,
                      ifelse(cip$cip < 2e-04, 2e-04, cip$cip))
    
    run_seq$seged <- NULL
    cip[, c("event", "seged")] <- NULL
    data_cip <- cip
    
    surv_obj <- Surv(time = run_seq$time, event = run_seq$event)
    km_fit <- survfit(surv_obj ~ 1)
    
    cip2 <- data.frame(time = summary(km_fit)$time,
                       cumhaz = summary(km_fit)$cumhaz)
    cip2$cip <- (1 - exp(-cip2$cumhaz)) / 2
    
    
    sajat_feno <- prepare_LTFHPlus_input(.tbl = run_seq,
                                         CIP = data_cip, 
                                         age_col = "time",
                                         aoo_col = "aoo",
                                         CIP_merge_columns = c("aoo"),
                                         CIP_cip_col = "cip",
                                         status_col = "event",
                                         use_fixed_case_thr = F,
                                         fam_id_col = "fid",
                                         personal_id_col = "pid",
                                         interpolation = NA,
                                         min_CIP_value = 1e-5)
    sajat_feno <- unique(sajat_feno)
    
    start <- Sys.time()
    sajat_liab <- estimate_liability(.tbl = sajat_feno,
                                     h2 = h2,
                                     pid = "pid",
                                     fam_id = "fid",
                                     role = "role",
                                     tol = 0.1,
                                     out = c(1))
    stop <- Sys.time()
    print(paste0(k, ". ciklus ", i, ". kör: ", stop - start))
    
    geno <- data2[start_index:end_index, ]
    geno <- as_FBM(geno)
    
    ########
    gwas <- big_univLinReg(X = geno,
                           y.train = sajat_liab$genetic_est,
                           ncores = 1,
                           covar.train = cov_seq)
    
    gwas$p.value <- predict(gwas, log10 = FALSE)
    full_p_tmp <- c(full_p_tmp, gwas$p.value)
    
    
    result[i, "p"] <- length(gwas$p.value[gwas$p.value > 0.00005])
    result[i, "p5"] <- length(gwas$p.value[gwas$p.value <= 0.00005])
    result[i, "p8"] <- length(gwas$p.value[gwas$p.value <= 0.00000005])
    result[i, "r"] <- i
    
    print(i)
    
  }
  
  result$point <- k
  plot_result <- rbind(plot_result, result)
  
  full_p_tmp <- data.frame(full_p = full_p_tmp)
  full_p_tmp$point <- k
  full_p <- rbind(full_p, full_p_tmp)
  

  
}
  
  filename <- files[h]
  filename <- gsub(".rds", "", filename)
  filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/ADuLT/", filename, "_ADuLT.rds")
  saveRDS(plot_result, file = filename, compress = "xz")
  saveRDS(full_p, file = paste0(filename, "_p"), compress = "xz")
  
}
