#for every MoA one big network!

library(glmnet)
library(gridExtra)
library(igraph)
library(c060)


setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/results concentration ranges/results concentration ranges")

cols_ <- c()

for (i in 1:42) { #iterate over all six MoAs
  load(list.files()[grep('merged',list.files())][i])
  print(list.files()[grep('merged',list.files())][i])
  parameter_mat <- merged_pmoas
  cols_ <- c()
  lasso_stab_new_res <- as.data.frame(rep(0,length(colnames(parameter_mat))))
  lasso_stab_new_res <- Matrix(as.matrix(lasso_stab_new_res),sparse = TRUE)
  rownames(lasso_stab_new_res) = colnames(parameter_mat)
 
 #length(colnames(parameter_mat)) 
  for (a in 1:length(colnames(parameter_mat))){
    
    parameter_mat[,a] <- as.numeric(parameter_mat[,a])
  }
  count = 0
  for (j in 1:12){
    count = count +1
    stop = 'no'
    #print(colnames(parameter_mat)[j])
    rem <- c(1:12,j)
    x_vars <- parameter_mat[,-unique(rem)]
    y_var <- parameter_mat[,j]
  
    if (var(y_var) <= 0.0011 & j <= 11) {
      print('no variance')
      y_var <- log(parameter_mat[,j]+0.01)
    }
    
    if (var(y_var) <= 0.0011 | length(which(y_var != 0)) <= 4) {
      #print('still no variance')
      stop = 'yes'
      cols_ <- c(cols_,colnames(parameter_mat)[j])
      lasso_stab_new_res <- cbind(lasso_stab_new_res,rep(0, length(colnames(parameter_mat))))
    }
    
    
    if (stop == 'no' ){
      
      cols_ <- c(cols_,colnames(parameter_mat)[j])
      
      stab_path <- stabpath(y = as.matrix(y_var),
               x = as.matrix(x_vars),
               family = "gaussian", 
               alpha = 1, #lasso
               nlambda = 100, #number of lambdas calculated
               steps = 100,
               mc.cores = 4)
      
      stab_sel = stabsel(stab_path, error = 0.05)
      add_var <- rep(0,length(unique(rem)))
      names(add_var) <- colnames(parameter_mat)[unique(rem)]
      add_lasso_res <- c(stab_path$x[,stab_sel$lpos],add_var)

      lasso_stab_new_res <- cbind(lasso_stab_new_res,add_lasso_res[match(names(add_lasso_res),rownames(lasso_stab_new_res))])
      
    }
    #%in% c(11,seq(1000, 27000, by = 1000),27361)
   # if (count == 11 ){
    #  save(lasso_stab_res, file = sprintf('lasso_stab_res0.05_%s_%s.RData',i,count))
    #  print(count)
    #}
    
    
    
  }
  save(lasso_stab_new_res, file = sprintf('lasso_stab_new_res0.05_%s_%s.RData',i,count))
  print(count)
  

  
  
}


