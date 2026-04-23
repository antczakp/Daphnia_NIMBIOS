#### calcualte jaccard similarity for each deb parameter between the different PMOAs

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


#list.files()[grep('ranger.RData',list.files())]

df_jacc <- matrix(nrow = 7, ncol = 7)
rownames(df_jacc) <- c('Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost', 'Increased_cost_of_reproduction')
colnames(df_jacc) <- c('Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost', 'Increased_cost_of_reproduction')

for (i in 1:7){
  for (k in 1:7){
    if (i != k){
      pmoa_genes1 <- load(list.files()[grep('ranger.RData',list.files())][i])
      pmoa_genes1 <- get(pmoa_genes1)
      pmoa_genes2 <- load(list.files()[grep('ranger.RData',list.files())][k])
      pmoa_genes2 <- get(pmoa_genes2)
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
    load(list.files()[grep('New_Mat',list.files())][k+1])
    print(list.files()[grep('New_Mat',list.files())][k+1])
    load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
    print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
    lasso_stab_new_res <- lasso_stab_new_res[,-1]
    dim(lasso_stab_new_res)
    colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
    rownames(lasso_stab_new_res) <- colnames(parameter_mat)
    
    num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    regs <- names(b)[1:num[[1]]]
    genes_[[k]] <- regs

}

save(genes_, file = 'lasso_ranger_stress_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_stress_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
save(df_jacc, file = 'df_jacc_lasso_ranger_stress.RData')
  
## stress
lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

hclust(df_jacc)


save(genes_, file = 'lasso_ranger_stress_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_stress_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
save(df_jacc, file = 'df_jacc_lasso_ranger_stress.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Stress Variables Dendogram', xlab = 'jaccard distance')


### structural length
i = 2

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_structural_length_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_structural_length_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_structural_length.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Structural length Variables Dendogram', xlab = 'jaccard distance')



### energy reserve
i = 3

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_senergy_reserve_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_senergy_reserve_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_energy_reserve.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Energy reserve Variables Dendogram', xlab = 'jaccard distance')


### reproduction buffer
i = 5

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_rep_buffer_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_rep_buffer_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_rep_buffer.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Reproduction buffer Variables Dendogram', xlab = 'jaccard distance')


### internal concentration
i = 6

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_int_conc_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_int_conc_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_int_conc.RData')


hc <- hclust(as.dist(df_jacc))
plot(hc,  main = 'Internal Concentration Variables Dendogram', xlab = 'jaccard distance')


### food density
i = 8

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_food_density_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_food_density_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_food_density.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Food density Variables Dendogram', xlab = 'jaccard distance')




### maximum structural length
i = 9

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_max_structural_length_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_max_structural_length_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_max_structural_length.RData')


hc <- hclust(as.dist(1-df_jacc))
plot(hc,  main = 'Maximum structural length Variables Dendogram', xlab = 'jaccard distance')


### cumulative embryos
i = 10

lasso_nums <- c(7,9,1,2,8,3,4,5,6)
genes_ <- list()
nr_params_ <- c(4,6,8,10,12,14,16,18,20,22)
ng <- c(2,0,3,4,5,6,7,8,9)
for (k in c(1,3:9)){
  load(list.files()[grep('New_Mat',list.files())][k+1])
  print(list.files()[grep('New',list.files())][k+1])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  
  num <- as.vector(lasso_caret_res[nr_params_[i],ng[k]])
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num[[1]]]
  genes_[[k]] <- regs
  
}

save(genes_, file = 'lasso_ranger_cum_embryos_genes.RData')

df_jacc <- matrix(nrow = 9, ncol = 9)

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_cum_embryos_genes.RData")

for (i in c(1,3:9)){
  for (k in c(1,3:9)){
    if (i != k){
      pmoa_genes1 <- genes_[[i]]
      pmoa_genes2 <- genes_[[k]]
      jac_ind <- jaccard(pmoa_genes1,pmoa_genes2)
      df_jacc[i,k] <- jac_ind
      df_jacc[k,i] <- jac_ind
    }
    else {
      df_jacc[i,k] <- 0
      df_jacc[k,i] <- 0
    }
  }
}

df_jacc <- df_jacc[,-2]
df_jacc <- df_jacc[-2,]
rownames(df_jacc) <- c( 'Increased_cost_of_reproduction', 'Dec_allocation_to_soma', 'Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')
colnames(df_jacc) <- c('Increased_cost_of_reproduction', 'Dec_allocation_to_soma','Dec_energy_conductance', 'Dec_feeding', 'Hazard_during_ooegenesis', 'Inc_allocation_to_soma', 'Inc_cost_for_growth_and_repro', 'Inc_maintenance_cost')

save(df_jacc, file = 'df_jacc_lasso_ranger_cum_embryos.RData')


hc <- hclust(as.dist(1- df_jacc))
plot(hc,  main = 'Cumulative embryos Variables Dendogram', xlab = 'jaccard distance')
