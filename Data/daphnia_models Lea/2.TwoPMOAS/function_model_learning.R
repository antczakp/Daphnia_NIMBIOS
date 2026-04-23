


lasso_caret_learn <- function(path){
  
  library(caret)
  #need folder with aracne outputs, lasso outputs and parameter mats
  res_df <- matrix(ncol = length(list.files())/3+1, nrow = 23)
  rownames(res_df) <- c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean')
 
  
  for (k in 1:length(list.files())/3){
    
    load(list.files()[grep('New_Mat',list.files())][k])
    #print(list.files()[grep('New_Mat',list.files())][k])
    load(list.files()[grep('lasso_stab',list.files())][k])
    #print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
    lasso_stab_new_res <- lasso_stab_new_res[,-1]
    #dim(lasso_stab_new_res)
    colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
    rownames(lasso_stab_new_res) <- colnames(parameter_mat)
    #for (i in c(1:10)){
    
    if (i == 4| i == 7){
      
      print('no predictions')
      print('no predictions')
    }
    
    else {
      a <- lasso_stab_new_res[-c(1:12),i]
      b <- sort(a, decreasing = TRUE )
      r_squared <- c()
      j = 0
      
      if (i == 1 | i == 6){
        while (j < 75){
          j = j + 1
          data = as.matrix(parameter_mat[,c(colnames(parameter_mat)[i],names(b)[1:j])])
          colnames(data)[i] <- 'g'
          plsFit <- train(log(g+0.1) ~.,
                          data = data,
                          method = 'lm',
                          trControl = ctrl
          )
          
          r_squared <- c(r_squared,plsFit$results$Rsquared)
          c <- r_squared[j] - r_squared[j-1]
          
        }
        #print(r_squared[which.max(r_squared)])
        #print(which.max(r_squared))
        
      }
    
      else{
        while (j < 75){
          j = j + 1
          data = as.matrix(parameter_mat[,c(colnames(parameter_mat)[i],names(b)[1:j])])
          colnames(data)[i] <- 'g'
          plsFit <- train(g ~.,
                          data = data,
                          method = 'lm',
                          trControl = ctrl
          )
          
          r_squared <- c(r_squared,plsFit$results$Rsquared)
          c <- r_squared[j] - r_squared[j-1]
          
        }
        #print(r_squared[which.max(r_squared)])
        #print(which.max(r_squared))
        
      }
    
     
    res_df[params_[i],k - 17] <- r_squared[which.max(r_squared)]
    res_df[params_[i]+1,k - 17] <- which.max(r_squared)
  }
  
  
}


res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    inc_cost_growth_dec_energ = rep(0,23),
                    inc_cost_growth_dec_feed = rep(0,23),
                    inc_cost_growth_inc_all = rep(0,23),
                    inc_cost_growth_inc_cost_rep = rep(0,23),
                    inc_cost_growth_inc_haz_oog	= rep(0,23),
                    inc_cost_growth_inc_maint = rep(0,23)
)
params_ <- c(3,5,7,9,11,13,15,17,19,21)
lasso_nums <- c(7,1,2,3,4,5,6)
nums_ <- c(0,0,9,0,0,8)
for (k in c(2,4:10)){
  
  load(list.files()[grep('New_Mat',list.files())][k])
  print(list.files()[grep('New_Mat',list.files())][k])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  print(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  #for (i in c(1:10)){
  
  if (i == 4| i == 7){
    
    print('no predictions')
    print('no predictions')
  }
  
  else {
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    r_squared <- c()
    j = 0
    
    if (i == 1 | i == 6){
      y_var <- log(parameter_mat[,i]+0.1)
    }
    else{
      y_var <- parameter_mat[,i]
    }
    while (j < 75){
      j = j + 1
      data = as.matrix(parameter_mat[,c(colnames(parameter_mat)[i],names(b)[1:j])])
      plsFit <- train(stress_ls~.,
                      data = data,
                      method = 'lm',
                      trControl = ctrl
      )
      
      r_squared <- c(r_squared,plsFit$results$Rsquared)
      c <- r_squared[j] - r_squared[j-1]
      
    }
    print(r_squared[which.max(r_squared)])
    print(which.max(r_squared))
    
  }
  #}
  res_df[params_[i],k - 17] <- r_squared[which.max(r_squared)]
  res_df[params_[i]+1,k - 17] <- which.max(r_squared)
}
lasso_caret_inc_cost_growth2 <- res_df
for (k in 2:7){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_caret_inc_cost_growth2[params_[i],k]
  }
  lasso_caret_inc_cost_growth2[23,k] <- sum(factorized)/10
}