
###lasso modelling


##first data prep
#need:
#columns with different Paramtet to be modelled and different genes
#rows with different samples/ days and concentrations

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/CountsLoaded.RData")
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DaphniaAnnotation.RData")
samp_id <- rownames(cannot.s2)
conc_id_annot <- data.frame(samp_id = rownames(cannot.s2), 
                            conc = cannot.s2$conc)


day_id_annot <- data.frame(samp_id = rownames(cannot.s2),
                           days = cannot.s2$day)

rna_count <- data[,-c(1:6)]



#calculate TPM
#1.Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2.Count up all the RPK values in a sample and divide this number by 1,000,000. This is your "per million" scaling factor.
#3.Divide the RPK values by the "per million" scaling factor. This gives you TPM.

TPM <- sapply(1:dim(rna_count)[2], function(i){
  RPK <- sapply(1:dim(rna_count)[1], function(j){
    #print(j)
    RPK <- rna_count[j,i]/data$Length[i]
    #print(RPK)
    return(RPK)
  })
  s <- sum(RPK)/10**6
  TPM <- RPK/s
  return(TPM)
})

rownames(TPM) <- data$Geneid
colnames(TPM) <- colnames(rna_count)
save(TPM, file = "TPM.RData")

n_MoA <- length(list.dirs())
#load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/log2_TPM3.RData")
#log2_TPM3 <- log2_TPM3[which(colnames(log2_TPM3) %in% samp_id),]
#setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs")
setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/fits_minus2Daphnids/fits_minus2Daphnids/")
res_conc <- c(0.015,0.06,0.12,0.25,0.00,1.00,4.00)
library(readr)
library(stringr)
all_dec_energ_conduct <- data.frame()
for (i in 2:n_MoA){
  print('i')
  print(i)
  setwd(list.dirs()[i])
  
  
  parameter_mat <- matrix(nrow = length(samp_id),ncol = 2)
  rownames(parameter_mat) <- colnames(log2_TPM)
  colnames(parameter_mat) <- c("days", "conc")
  parameter_mat <- as.data.frame(parameter_mat)
  parameter_mat$days <- day_id_annot$days[match(rownames(parameter_mat), day_id_annot$samp_id)]
  parameter_mat$conc <- conc_id_annot$conc[match(rownames(parameter_mat), day_id_annot$samp_id)]
  express_data <- t(log2_TPM[,match(rownames(parameter_mat),colnames(log2_TPM))])
  parameter_mat <- cbind(parameter_mat,express_data)
 
  #parameter_mat <- parameter_mat[which(parameter_mat$conc %in% c(1.000, 4.000, 0.000)),]
  for ( j in 2:8){
    print('j')
    print(j)
    print(list.files()[j])
    results_ <- read_csv(list.files()[j], col_names = FALSE)
    
    
    colnames(results_) <- c("the structural length, cm",
                                  "energy reserve, J",
                                  "energy invested in maturity, J",
                                  "reproduction buffer, J",
                                  "internal concentration, nM",
                                  "survival probability, %", 
                                  "food density", 
                                  "maximum structural length",
                                  "cumulative embryos, #"
    )
    
    timepoint <- c(1,
                   3,
                   6,
                   8,
                   10,
                   13,
                   15,
                   17,
                   20,
                   22,
                   24,
                   27,
                   29,
                   31
    )
    conc <- rep(res_conc[j-1], 14 )
    results_ <- cbind(results_,timepoint, conc)
    all_dec_energ_conduct <- rbind(all_dec_energ_conduct, results_)
    
  }
  
  print('i')
  print(i)
  #setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs")
  setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/fits_minus2Daphnids/fits_minus2Daphnids/")
  save(all_dec_energ_conduct, file = sprintf('Param_mat%s.RData',str_remove(list.dirs()[i], './')))
  
  params_ = data.frame()
  sample_ = c()
  for (k in c(6,13,17,24)){
    
    first_id_deb <- which(all_dec_energ_conduct$timepoint == k)
    first_id_om <- which(parameter_mat$days == k)
    
    for (l in unique(parameter_mat$conc)){
      
      second_id_deb <- which(all_dec_energ_conduct[first_id_deb,]$conc == l)
      second_id_om <- which(parameter_mat[first_id_om,]$conc == l)
      sample_ <- c(sample_,rownames(parameter_mat)[first_id_om][second_id_om])
      print(rownames(parameter_mat)[first_id_om][second_id_om])
      for (t in 1:length(rownames(parameter_mat)[first_id_om][second_id_om])){
        print(all_dec_energ_conduct[first_id_deb,][second_id_deb,1:9])
        params_ <- rbind(params_,all_dec_energ_conduct[first_id_deb,][second_id_deb,1:9])
      }
    }
  }
  
  parameter_mat <- cbind(params_[match(rownames(parameter_mat),sample_),], parameter_mat)
  save(parameter_mat, file = sprintf('Mat_%s.RData',str_remove(list.dirs()[i], './')))

  all_dec_energ_conduct <- data.frame()
}


#c_0 <- c(0.0001213,0.008145,4.053e-05,1.322e-05,3.867e-05,0.1818,4.272e-05,1.322e-05,3.64e-06)
#c_T <- c(0.005718,0.7139,0.9982,0.897,1.166,2.21,0.7448, 0.897,0.03187)
c_0 <- c(6.225e-07,0.0001501,7.216e-07,3.988e-05,3.625e-05,7.614e-05,0.01192)
c_T <- c(0.001301,0.001047,0.000608,5.607e-05,5.709e-05,0.001133,1.183)

c_0 <- c(0.00013,0.0001223,0.0008062,4.053e-05,4.363e-05,4.363e-05,0.1884)
c_T <- c(0.169,0.005695,0.7281,0.9982,1.279,1.279,2.174)
stress_ls <- c()
for (i in 1:7) { #iterate over all six MoAs
  
  load(list.files()[grep('Mat',list.files())[i]])
  print(list.files()[grep('Mat',list.files())[i]])
  
  for ( j in 1:dim(parameter_mat)[1]){
    print(j)
    c_V <- parameter_mat$`internal concentration, nM`[j]
    print(c_V)
    print(c_T[i])
    print(c_0[i])
    stress <- (1/c_T[i])*max(c(0,c_V-c_0[i]))
    print(stress)
    stress_ls <-c(stress_ls,stress)
  }
  parameter_mat <- cbind(stress_ls,parameter_mat)
  save(parameter_mat, file = sprintf('New_%s_2.RData',str_remove(list.files()[grep('Mat',list.files())[i]],'.RData')))
  setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/fits_minus2Daphnids/fits_minus2Daphnids/")
  stress_ls <- c()
}

pmoas <- c('dec_feed','dec_energ','inc_all','inc_cost_growth','inc_cost_rep','inc_haz_oog','inc_maint')
setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/results concentration ranges/results concentration ranges")
for(i in 1:7){
  load(sprintf('results 0.015-0.25/%s',list.files('results 0.015-0.25/')[grep('New_Mat',list.files('results 0.015-0.25/'))][i]))
  PMOA_low <- parameter_mat
  for (j in 1:7){
    if ( i == j){
      load(sprintf('results 1&4/%s',list.files('results 1&4/')[grep('New_Mat',list.files('results 1&4/'))][i]))
      PMOA_high <- parameter_mat
      merged_pmoas <- rbind(PMOA_low,PMOA_high)
      save(merged_pmoas, file = sprintf('merged_pmoas_%s_%s.RData', pmoas[i], pmoas[j]))
    }
    
  }
}


library(glmnet)
library(gridExtra)
library(igraph)

setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs")

cols_ <- c()

for (i in 1:6) { #iterate over all six MoAs
  cols_ <- c()
  cv_out = list()
  best_lam = list()
  lasso_res <- as.data.frame(rep(0,27350))
  pred = list()
  
  
  load(list.files()[grep('Mat',list.files())[i]])
  print(list.files()[grep('Mat',list.files())[i]])
  
  for (a in 1:length(colnames(parameter_mat))){
    
    parameter_mat[,a] <- as.numeric(parameter_mat[,a])
  }
  
  for (j in 1:11){
    stop = 'no'
    print(colnames(parameter_mat)[j])
    x_vars <- parameter_mat[,-(1:11)]
    y_var <- parameter_mat[,j]
    
    lambda_seq <- 10^seq(2, -2, by = -.1)
    
    # Splitting the data into test and train
    set.seed(86)
    train = sample(1:nrow(x_vars), 5*nrow(x_vars)/6)
    x_test = (-train)
    y_test = y_var[x_test]
    
    if (var(y_var[train]) == 0){
      print('still no variance')
      stop = 'yes'
    }
    
    else if(var(y_var[train]) < 0.001 & stop == 'no' ){
      print("no variance")
      y_var <- log(parameter_mat[,j] + 0.01) 
      x_test = (-train)
      y_test = y_var[x_test]
    }
    
    if (stop == 'no' ){
  
      cols_ <- c(cols_,colnames(parameter_mat)[j])
      
      cv_output <- cv.glmnet(as.matrix(x_vars[train,]), y_var[train],
                             alpha = 1, family = 'gaussian', intercept = TRUE, lambda = lambda_seq, 
                             nfolds = 10)
      cv_out[colnames(parameter_mat)[j]] <- cv_output
      # identifying best lamda
      best_lam[colnames(parameter_mat)[j]] <- cv_output$lambda.min
      #best_lam
      
      lasso_best <- cv_output$glmnet.fit
      lasso_res <- cbind(lasso_res,lasso_best$beta[,which(as.character(lambda_seq) == as.character(cv_output$lambda.min))])
      pred_res <- predict(lasso_best, s = unname(best_lam[colnames(parameter_mat)[j]])[[1]], newx = as.matrix(x_vars[x_test,]))
      
      pred[colnames(parameter_mat)[j]] <- pred_res
      
      print(cor(pred_res, parameter_mat[x_test,j]))
      
      #print(which(lasso_best$beta != 0))
      
      if (length(which(lasso_best$beta[,which(as.character(lambda_seq) == as.character(cv_output$lambda.min))] != 0)) == 0){
        stop = 'yes'
        print('stop')
      }
    
    }

  }
  
  
  lasso_results <- lasso_res
  lasso_results <- lasso_results[,-1]
  print('finish')
  colnames(lasso_results) <- cols_
  lasso_results <- as.matrix(lasso_results)
  
  setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs")
  
  save(cv_out, file = sprintf("cv_out_%s.RData",str_remove(list.dirs()[1+i], './')))
  save(best_lam, file = sprintf("best_lam_%s.RData",str_remove(list.dirs()[1+i], './')))
  save(lasso_res, file = sprintf("lasso_res_%s.RData",str_remove(list.dirs()[1+i], './')))
  save(pred, file = sprintf("pred_%s.RData",str_remove(list.dirs()[1+i], './')))
  save(lasso_results, file = sprintf("lasso_results_%s.RData",str_remove(list.dirs()[1+i], './')))
  
  p = list()
  for (l in 1:length(colnames(lasso_results))){
    print(colnames(lasso_results)[l])
    
    energ_res <- sparseMatrix(i = c(), j = c(),dims=list(27350,27350))
    energ_res <- cbind(lasso_results[,l],energ_res )
    
    
    energ_res_add <- sparseMatrix(i = c(), j = c(),dims=list(1,27351))
    energ_res <- rbind(energ_res_add,energ_res)
    colnames(energ_res) <- c(colnames(lasso_results)[l],rownames(lasso_results))
    rownames(energ_res) <- c(colnames(lasso_results)[l],rownames(lasso_results))
    
    if ( length(which(energ_res != 0)) == 0){
      print('no predictors')
    }
    
    else {
      print('predictors')
      print(length(which(energ_res != 0)))
      graph_energ_res <- graph_from_adjacency_matrix(energ_res, mode = 'directed', weighted = TRUE)
      
      energ_res1 <- energ_res[c(1,which(rowSums(energ_res) > 0)),c(1,which(rowSums(energ_res) > 0))]
      graph_energ_res <- graph_from_adjacency_matrix(energ_res1, mode = 'directed', weighted = TRUE)
  
      V(graph_energ_res)$color <- rep('orange', length(V(graph_energ_res)$name))
      V(graph_energ_res)$color[1] <- 'red'
      V(graph_energ_res)$size <- energ_res1[,1]*round(10 / energ_res1[,1])
      V(graph_energ_res)$size[1] <- 10
      
      png(file=sprintf("net_plot_%s_%s.png",str_remove(list.dirs()[1+i], './'),colnames(lasso_results)[l]),width = 1000, height = 1000, units = "px")
      plot(graph_energ_res, edge.width=E(graph_energ_res)$weight, vertex.color = V(graph_energ_res)$color, layout=layout.lgl, edge.arrow.size = 0.5, vertex.label.color="black",  vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2)
      dev.off()
      
      p[l] <- graph_energ_res
    }
  }
  
  save(p, file = sprintf("p_%s.RData", str_remove(list.dirs()[1+i], './')))
  
}




lm_deb_energy_reserve <- lm(`energy reserve, J` ~ APZ42_000356 + APZ42_000535 + APZ42_001439 + APZ42_003672 + APZ42_004462 + APZ42_005087 + APZ42_005487 + APZ42_005794 + APZ42_006310 + APZ42_012103 + APZ42_013147 + APZ42_014675 + APZ42_018112 + APZ42_019570 + APZ42_020235 + APZ42_021071 + APZ42_022478 + APZ42_022811 + APZ42_027238 + APZ42_028278 + APZ42_032809 + APZ42_032875, data = parameter_mat)
summary(lm_deb_energy_reserve)
hist(lm_deb_energy_reserve$residuals, freq = FALSE)
plot(lm_deb_energy_reserve)
shapiro.test(lm_deb_energy_reserve$residuals)
plot(parameter_mat[,c('energy reserve, J',names(which(lasso_results[,1] != 0)))])


lm_deb_reproduction_buffer <- lm(`reproduction buffer, J` ~ APZ42_000530 + APZ42_003672 + APZ42_008549 + APZ42_010559 + APZ42_013147 + APZ42_017283 + APZ42_017894 + APZ42_017940 + APZ42_018124 + APZ42_018287 + APZ42_019693 + APZ42_022504 + APZ42_028344 + APZ42_028671 + APZ42_029556 + APZ42_032229, data = parameter_mat)
summary(lm_deb_reproduction_buffer)
hist(lm_deb_reproduction_buffer$residuals, freq = FALSE)
plot(lm_deb_reproduction_buffer)
shapiro.test(lm_deb_reproduction_buffer$residuals)
plot(parameter_mat[,c('reproduction buffer, J',names(which(lasso_results[,2] != 0)))])