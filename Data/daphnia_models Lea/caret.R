install.packages("caret", dependencies = c("Depends", "Suggests"))
library(caret)
library(mlbench)


set.seed(107)
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_decreased feeding_2.RData")
load('lasso_stab_new_res0.05_1_12.RData')
lasso_stab_new_res <- lasso_stab_new_res[,-1]
dim(lasso_stab_new_res)
colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)
a <- lasso_stab_new_res[-c(1:12),1]
b <- sort(a, decreasing = TRUE )

inTrain <- createDataPartition(y = parameter_mat[,1],
                               p=.75,
                               list = FALSE)
str(inTrain)
training <- parameter_mat[inTrain,]
testing <- parameter_mat[-inTrain,]
nrow(training)
nrow(testing)
ctrl <- trainControl(method = 'boot',
                    number = 100
                    )

data = parameter_mat[,c(colnames(parameter_mat)[1],names(b)[1:66])]
plsFit <- train(stress_ls ~.,
                data = data,
                method = 'lm',
                trControl = ctrl
                )
plsFit$results$Rsquared

ctrl <- trainControl(method = 'boot',
                     number = 100
)



res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    Cost_Reproduction = rep(0,23),
                    bla = rep(0,23),
                    Decreased_energy = rep(0,23),
                    Decreased_feeding = rep(0,23),
                    Decreased_fraction = rep(0,23),
                    Hazard_oogenesis	= rep(0,23),
                    Increased_Allocatio = rep(0,23),
                    Increased_cost = rep(0,23),
                    Increased_Maintenance = rep(0,23)
)
res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    dec_all_inc_rep = rep(0,23),
                    dec_ener_inc_rep = rep(0,23),
                    dec_feed_inc_rep = rep(0,23),
                    inc_all_inc_rep	= rep(0,23),
                    inc_maint_inc_rep = rep(0,23)
)
params_ <- c(3,5,7,9,11,13,15,17,19,21)
lasso_nums <- c(9,1,2,3,4)
nums_ <- c(0,0,9,0,0,8)
for (k in 1:5){
  
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
        plsFit <- train(`cumulative embryos, #`~.,
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
  res_df[params_[i],k + 1] <- r_squared[which.max(r_squared)]
  res_df[params_[i]+1,k +1] <- which.max(r_squared)
}
lasso_caret_2Pmoas_res <- res_df
for (k in 2:6){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_caret_2Pmoas_res[params_[i],k]
  }
  lasso_caret_2Pmoas_res[23,k] <- sum(factorized)/8
}
r_sq <- c(0,0,0,0,0,0,0)
names(r_sq) <- c('dec_energ','dec_feed','inc_all','inc_cost_growth','inc_cost_rep','inc_haz_oog','inc_maint')
names_r <- c()
for (i in 1:7){ 
  lasso_df <- load(list.files()[grep('lasso_caret',list.files())][i])
  lasso_df <- get(lasso_df)
  r_sq[i] <- lasso_df[23,which.max(lasso_df[23,])]
  names_r[i] <- colnames(lasso_df)[which.max(lasso_df[23,])]
}
names(r_sq) <- names_r
as.data.frame(r_sq)
for (k in 2:8){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_ranger_res[params_[i],k] * factors_params[[k-1]][i] 
  }
  lasso_ranger_res[23,k] <- sum(factorized)/factors_params[[k-1]][11]
}

res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    Cost_Reproduction = rep(0,23),
                    bla = rep(0,23),
                    Decreased_energy = rep(0,23),
                    Decreased_feeding	= rep(0,23),
                    Decreased_fraction = rep(0,23),
                    Hazard_oogenesis = rep(0,23),
                    Increased_Allocation = rep(0,23),
                    Increased_cost = rep(0,23),
                    Increased_Maintenance = rep(0,23)
)

res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    dec_all_inc_rep = rep(0,23),
                    dec_energyinc_rep = rep(0,23),
                    dec_feedinginc_rep = rep(0,23),
                    inc_allocatioinc_rep	= rep(0,23),
                    inc_maintenanceinc_rep = rep(0,23)
)

params_ <- c(3,5,7,9,11,13,15,17,19,21)
lasso_nums <- c(9,1,2,3,4)
for (k in 1:5){
  
  load(list.files()[grep('New_Mat_',list.files())][k])
  print(list.files()[grep('New_Mat_',list.files())][k])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',lasso_nums[k]))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  for (i in c(1:3,5:6,8:10)){
    
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    #print(colnames(parameter_mat)[i])
    r_squared <- c()
    j = 1
    if (i == 1 | i == 6){
      y_var <- log(parameter_mat[,i]+0.001)
    }
    else{
      y_var <- parameter_mat[,i]
    }
    j = 1
    while (j < 75){
      j = j + 1
      rf1 <- ranger(y = y_var, x = as.matrix(parameter_mat[,names(b)[1:j]]), importance = 'impurity', num.trees = 10000)
      r_squared <- c(r_squared,rf1$r.squared)
      
    }
    
    print(max(r_squared))
    print(which.max(r_squared))
    
    
    res_df[params_[i],k+1] <- r_squared[which.max(r_squared)]
    res_df[params_[i]+1,k+1] <- which.max(r_squared)
  }
}
p <- ggplot(parameter_mat, aes(x = conc, y = stress_ls, col = as.character(days)))
p + geom_point(size = 3)
lasso_ranger_2PMOAS_res <- res_df

for (k in 2:6){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_ranger_2PMOAS_res[params_[i],k]
  }
  lasso_ranger_2PMOAS_res[23,k] <- sum(factorized)/8
}

res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    dec_ener_inc_rep = rep(0,23),
                    dec_feed_inc_rep = rep(0,23),
                    inc_all_inc_rep	= rep(0,23),
                    inc_maint_inc_rep = rep(0,23)
)
params_ <- c(3,5,7,9,11,13,15,17,19,21)
for (k in 1:4){
  
  load(list.files()[grep('New_Mat',list.files())][k])
  print(list.files()[grep('New_Mat',list.files())][k])
  load(sprintf('lasso_stab_new_res0.05_%s_12.RData',k))
  lasso_stab_new_res <- lasso_stab_new_res[,-1]
  dim(lasso_stab_new_res)
  colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
  rownames(lasso_stab_new_res) <- colnames(parameter_mat)
  for (i in c(1:3,5:6,8:10)){
    
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    #print(colnames(parameter_mat)[i])
    r_squared <- c()
    j = 1
    while (j < 75){
      j = j + 1
      rf1 <- ranger(y = parameter_mat[,i], x = as.matrix(parameter_mat[,names(b)[1:j]]), importance = 'impurity', num.trees = 10000)
      r_squared <- c(r_squared,rf1$r.squared)
      
    }
      
      print(max(r_squared))
      print(which.max(r_squared))
    
    
    res_df[params_[i],k] <- r_squared[which.max(r_squared)]
    res_df[params_[i]+1,k] <- which.max(r_squared)
  }
}

lasso_ranger_res <- res_df

for (k in 2:8){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_ranger_res[params_[i],k] * factors_params[[k-1]][i] 
  }
  lasso_ranger_res[23,k] <- sum(factorized)/factors_params[[k-1]][11]
}

params <-(stress_ls)

enrichment_fn <- function(signifGenes,pathwayGenes,genomeBackground=15208,useEASE=T,retPvalOnly=F,...){
  FP <- length(intersect(signifGenes,pathwayGenes))
  if(FP == 0)
    return(NA)
  PT <- length(pathwayGenes)
  GB <- genomeBackground
  FT <- length(signifGenes)
  if(useEASE){
    mm <- matrix(c(FP-1,PT,FT-FP,GB-PT),nrow=2,byrow=T)
  }else{
    mm <- matrix(c(FP,PT,FT-FP,GB-PT),nrow=2,byrow=T)
  }
  rownames(mm) <- c("In Pathway","Not in Pathway")
  colnames(mm) <- c("User Genes","Genome")
  #print(mm)
  enrich_genes <- intersect(signifGenes,pathwayGenes)
  fs <- fisher.test(mm,...)
  if(retPvalOnly){
    return(c(fs$p.value,enrich_genes))
  }else{
    return(fs)
  }
}


load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])
load(sprintf('lasso_stab_new_res0.05_%s_12.RData',1))
lasso_stab_new_res <- lasso_stab_new_res[,-1]
dim(lasso_stab_new_res)
colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)

all_regs <- c()
num1 <- c(7,54,66,0,63,7,0,18,63,37)
num2 <- c(17,57,46,0,4,64,0,11,28,35)
num3 <- c(17,66,72,0,9,3,0,10,63,40)
num4 <- c(11,49,15,0,74,38,0,63,52,72)
num5 <- c(16,36,51,0,69,20,0,64,65,21)
num6 <- c(28,15,15,0,65,35,0,50,48,37)
num7 <- c(69,67,50,0,71,67,0,26,74,72)
for (i in 1:10){
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num2[i]]
  all_regs <- c(all_regs,regs)
}
all_regs_unique <- unique(all_regs)



genes_inc_rep_ranger <- all_regs_unique
genes_dec_energy_ranger <- all_regs_unique
genes_dec_feed_ranger <- all_regs_unique
genes_haz_oog_ranger <- all_regs_unique
genes_inc_all_ranger <- all_regs_unique
genes_inc_cost_growth <- all_regs_unique
genes_inc_maint_ranger <- all_regs_unique

sort(table(all_regs))

names(which(sort(table(all_regs)) > 1))
sort(table(all_regs))[which(sort(table(all_regs)) > 1)]
#View(kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5])
sort(table(all_regs))[which(sort(table(all_regs)) > 1)][which(names(which(sort(table(all_regs)) > 1)) %in% kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5]$gene)]

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})


load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])
load(sprintf('lasso_stab_new_res0.05_%s_12.RData',k))
lasso_stab_new_res <- lasso_stab_new_res[,-1]
dim(lasso_stab_new_res)
colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)

num1 <- c(19,16,36,0,5,23,0,13,28,8)
num2 <- c(23,12,44,0,10,35,0,22,31,20)
num3 <- c(11,10,17,0,5,11,0,11,11,6)
num4 <- c(11,22,10,0,19,16,0,22,22,7)
num5 <- c(14,21,16,0,16,11,0,21,20,5)
num6 <- c(4,13,10,0,4,7,0,66,4,17)

for (j in 1:10){
  a <- lasso_stab_new_res[-c(1:12),j]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num6[j]]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
      print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    }
  }
}


#lm with ranger

load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])
load(sprintf('lasso_stab_new_res0.05_%s_12.RData',k))
lasso_stab_new_res <- lasso_stab_new_res[,-1]
dim(lasso_stab_new_res)
colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)

all_regs <- c()
num1 <- c(18,58,48,0,5,63,0,12,29,36)
num2 <- c(18,67,73,0,10,4,0,11,64,40)
num3 <- c(12,52,16,0,75,38,0,71,52,71)
num4 <- c(17,37,50,0,74,22,0,65,56,25)
num5 <- c(29,74,16,0,58,36,0,55,41,38)
num6 <- c(71,10,71,0,72,69,0,66,74,73)
for (i in 1:10){
  a <- lasso_stab_new_res[-c(1:12),i]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num6[i]]
  all_regs <- c(all_regs,regs)
}
all_regs_unique <- unique(all_regs)
sort(table(all_regs))

names(which(sort(table(all_regs)) > 1))
sort(table(all_regs))[which(sort(table(all_regs)) > 1)]
#View(kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5])
sort(table(all_regs))[which(sort(table(all_regs)) > 1)][which(names(which(sort(table(all_regs)) > 1)) %in% kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5]$gene)]

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})


load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])
load(sprintf('lasso_stab_new_res0.05_%s_12.RData',k))
lasso_stab_new_res <- lasso_stab_new_res[,-1]
dim(lasso_stab_new_res)
colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)

num1 <- c(23,62,22,0,5,3,0,14,28,36)
num2 <- c(18,36,74,0,69,5,0,11,71,40)
num3 <- c(14,56,16,0,74,40,0,69,54,19)
num4 <- c(19,50,54,0,71,24,0,54,27,28)
num5 <- c(14,75,38,0,7,37,0,74,51,75)
num6 <- c(73,13,51,0,72,27,0,60,48,36)

for (j in 1:10){
  a <- lasso_stab_new_res[-c(1:12),j]
  b <- sort(a, decreasing = TRUE )
  regs <- names(b)[1:num6[j]]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
      print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    }
  }
}




aracne_mat <- read.table(file = sprintf('aracne_outputs/%s',list.files(path = "C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_outputs")[k]),sep="\t", fill = TRUE, header = TRUE)

load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])

all_regs <- c()
num1 <- c(5,12,4,0,10,1,0,12,4,7)
num2 <- c(1,7,17,0,6,3,0,2,1,11)
num3 <- c(6,11,5,0,8,1,0,6,5,13)
num4 <- c(7,17,19,0,7,5,0,8,19,36)
num5 <- c(4,28,17,0,13,7,0,3,5,2)
num6 <- c(1,9,9,0,5,2,0,7,5,14)
for (i in 1:10){
  regs <- aracne_mat[which(aracne_mat[,2] == i),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == i),]$MI, decreasing = TRUE, index.return = TRUE)
  b <- regs[b$ix]
  b <- b[which(b > 12)]
  names(b) <- colnames(parameter_mat)[b]
  all_regs <- c(all_regs,names(b))
}
all_regs_unique <- unique(all_regs)
sort(table(all_regs))

names(which(sort(table(all_regs)) > 1))
sort(table(all_regs))[which(sort(table(all_regs)) > 1)]
#View(kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5])
sort(table(all_regs))[which(sort(table(all_regs)) > 1)][which(names(which(sort(table(all_regs)) > 1)) %in% kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5]$gene)]

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
  }
  return(enriched_paths)
  
})


for (j in 1:10){ 
  regs <- aracne_mat[which(aracne_mat[,2] == j),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == j),]$MI, decreasing = TRUE, index.return = TRUE)
  regs <- regs[b$ix]
  regs <- regs[which(regs > 12)]
  regs <- colnames(parameter_mat)[regs]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
    }
  }
}



aracne_mat <- read.table(file = sprintf('aracne_outputs/%s',list.files(path = "C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_outputs")[k]),sep="\t", fill = TRUE, header = TRUE)

load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])

all_regs <- c()
num1 <- c(6,6,11,0,3,13,0,7,8,1)
num2 <- c(1,7,9,0,4,10,0,4,6,3)
num3 <- c(13,4,6,0,7,2,0,11,2,3)
num4 <- c(9,13,4,0,5,5,0,1,6,3)
num5 <- c(12,13,26,0,20,7,0,2,8,1)
num6 <- c(1,6,3,0,18,3,0,8,4,7)
for (i in 1:10){
  regs <- aracne_mat[which(aracne_mat[,2] == i),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == i),]$MI, decreasing = TRUE, index.return = TRUE)
  b <- regs[b$ix]
  b <- b[which(b > 12)]
  names(b) <- colnames(parameter_mat)[b]
  all_regs <- c(all_regs,names(b))
}
all_regs_unique <- unique(all_regs)
sort(table(all_regs))

names(which(sort(table(all_regs)) > 1))
sort(table(all_regs))[which(sort(table(all_regs)) > 1)]
#View(kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5])
sort(table(all_regs))[which(sort(table(all_regs)) > 1)][which(names(which(sort(table(all_regs)) > 1)) %in% kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5]$gene)]

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})

#[1:num1[j]]
for (j in 1:10){ 
  regs <- aracne_mat[which(aracne_mat[,2] == j),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == j),]$MI, decreasing = TRUE, index.return = TRUE)
  regs <- regs[b$ix]
  regs <- regs[which(regs > 12)]
  regs <- colnames(parameter_mat)[regs[1:num1[j]]]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
      print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    }
  }
}
kannot[which(kannot$gene %in% c("APZ42_013661"   ,    "APZ42_013795")),2:3]



net_dec_energ_norm <- read.table(file = 'network_dec_energ_norm.txt',sep="\t", fill = TRUE, header = TRUE)
net_dec_energ_noTarg <- read.table(file = 'network_dec_energ_noTarg.txt',sep="\t", fill = TRUE, header = TRUE)
net_dec_energ_noTarg_p1.1 <- read.table(file = 'network_dec_energ_noTrag_p1.1.txt',sep="\t", fill = TRUE, header = TRUE)

#network_dec_feed_norm.txt

setwd("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs")
#list.files()[grep('noTa',list.files())][k]
res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
              Cost_Reproduction = rep(0,23),
              Decreased_energy = rep(0,23),
              Decreased_feeding	= rep(0,23),
              Hazard_oogenesis = rep(0,23),
              Increased_Allocation = rep(0,23),
              Increased_cost = rep(0,23),
              Increased_Maintenance = rep(0,23)
)
params_ <- c(3,5,7,9,11,13,15,17,19,21)
library(caret)
ctrl <- trainControl(method = 'boot',
                     number = 100
)
for (k in 2:8){
  
  load(list.files()[grep('New_Mat',list.files())][k])
  print(list.files()[grep('New_Mat',list.files())][k])
  aracne_mat <- read.table(file = sprintf('aracne_outputs/%s',list.files(path = "C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_outputs")[k - 1]),sep="\t", fill = TRUE, header = TRUE)
  max_r <- c()
 # for (i in 1:10){
  regs_num <- aracne_mat[which(aracne_mat[,2] == i),]$Regulator
  regs_MI <- aracne_mat[which(aracne_mat[,2] == i),]$MI
  names(regs_MI) <- colnames(parameter_mat)[regs_num]
  regs_MI <- regs_MI[which(regs_num > 12)]
  regs_num2 <- aracne_mat[which(aracne_mat[,1] == i),]$Target
  regs_MI2 <- aracne_mat[which(aracne_mat[,1] == i),]$MI
  names(regs_MI2) <- colnames(parameter_mat)[regs_num2]
  regs_MI2 <- regs_MI2[which(regs_num2 > 12)]
  regs <- c(regs_MI,regs_MI2[which(!names(regs_MI2) %in% names(regs_MI))])
  b <- sort(regs, decreasing = TRUE, index.return = TRUE)
  b <- regs[b$ix]
  mean_r <- c()
  if (length(names(b)) > 0){
    j = 0
    r_squared <- c()
    while (j < length(names(b))){
      j = j + 1
      # print(j)
      data = as.matrix(parameter_mat[,c(colnames(parameter_mat)[i],names(b)[1:j])])
      plsFit <- train(`cumulative embryos, #`~.,
                      data = data,
                      method = 'lm',
                      trControl = ctrl
      )
      
      r_squared <- c(r_squared,plsFit$results$Rsquared)
    }
    
    max_r <- c(max_r,which.max(r_squared))
    print(r_squared[which.max(r_squared)])
    print(which.max(r_squared))
    mean_r <- c(mean_r,r_squared[which.max(r_squared)])
  }
  
  else {
    #print(colnames(parameter_mat)[i])
    print('no predictors')
    print(length(regs[which(regs > 12)]))
  }
  
 # }
  print(sum(mean_r)/length(mean_r))
  res_df[params_[i],k] <- r_squared[which.max(r_squared)]
  res_df[params_[i]+1,k] <- which.max(r_squared)
}

for (k in 2:8){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- res_df[params_[i],k]
  }
  res_df[23,k] <- sum(factorized)/10
}

for (k in 2:8){
 
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- res_df[params_[i],k] * factors_params[[k-1]][i] 
  }
  res_df[23,k] <- sum(factorized)/factors_params[[k-1]][11]
}

res_df <- res_df[-24,]

aracne_caret_res <- res_df
save(aracne_caret_res, file = 'arcane_caret_res.RData')

res_df = data.frame(DEB_Param = c('MRE','SMSE','stress_ls','Variables','the structural length, cm', 'NR. Variables', 'energy reserve, J', 'NR. Variables', 'energy invested in maturity, J', 'NR. Variables', 'reproduction buffer, J' , 'NR. Variables' ,'internal concentration, nM', 'NR. Variables', 'survival probability, %', 'NR. Variables', 'food density', 'NR. Variables', 'maximum structural length', 'NR. Variables', 'cumulative embryos, #', 'NR. Variables','mean'),
                    Cost_Reproduction = rep(0,23),
                    Decreased_energy = rep(0,23),
                    Decreased_feeding	= rep(0,23),
                    Hazard_oogenesis = rep(0,23),
                    Increased_Allocation = rep(0,23),
                    Increased_cost = rep(0,23),
                    Increased_Maintenance = rep(0,23)
)
params_ <- c(3,5,7,9,11,13,15,17,19,21)
for (k in 2:8){
  
  load(list.files()[grep('New_Mat',list.files())][k])
  print(list.files()[grep('New_Mat',list.files())][k])
  aracne_mat <- read.table(file = sprintf('aracne_outputs/%s',list.files(path = "C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_outputs")[k - 1]),sep="\t", fill = TRUE, header = TRUE)
  #print(list.files()[grep('noTrag',list.files())][k])
  for (i in c(1:3,5:6,8:10)){
    regs_num <- aracne_mat[which(aracne_mat[,2] == i),]$Regulator
    regs_MI <- aracne_mat[which(aracne_mat[,2] == i),]$MI
    names(regs_MI) <- colnames(parameter_mat)[regs_num]
    regs_MI <- regs_MI[which(regs_num > 12)]
    regs_num2 <- aracne_mat[which(aracne_mat[,1] == i),]$Target
    regs_MI2 <- aracne_mat[which(aracne_mat[,1] == i),]$MI
    names(regs_MI2) <- colnames(parameter_mat)[regs_num2]
    regs_MI2 <- regs_MI2[which(regs_num2 > 12)]
    regs <- c(regs_MI,regs_MI2[which(!names(regs_MI2) %in% names(regs_MI))])
    b <- sort(regs, decreasing = TRUE, index.return = TRUE)
    regs <- regs[b$ix]
    
    if (length(regs) > 0){
      
      x_vars <- parameter_mat
      y_vars <- parameter_mat[,i]
      
      train = sample(1:nrow(x_vars), 5*nrow(x_vars)/6)
      x_test = (-train)
      y_test = y_vars[x_test]
      
      r_squared <- c()
      
      rf1 <- ranger(y = as.matrix(y_vars), x = as.matrix(x_vars[,names(regs)]), importance = 'impurity', num.trees = 10000)
      b <- sort(rf1$variable.importance, decreasing = TRUE )
      for (j in 2:length(regs)){
        rf2 <- ranger(y = parameter_mat[,i], x = as.matrix(parameter_mat[,names(regs)[1:j]]), importance = 'impurity',num.trees = 10000)
        r_squared <- c(r_squared,rf2$r.squared)
      }
      
      print(max(r_squared))
      print(which.max(r_squared))
    }
    
    else {
      #print(colnames(parameter_mat)[i])
      print('no predictors')
      print(length(regs[which(regs > 12)]))
    }
    
    res_df[params_[i],k] <- r_squared[which.max(r_squared)]
    res_df[params_[i]+1,k] <- which.max(r_squared)
  }
}


for (k in 2:8){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- res_df[params_[i],k]
  }
  res_df[23,k] <- sum(factorized)/10
}


for (k in 2:8){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- res_df[params_[i],k] * factors_params[[k-1]][i]
  }
  res_df[23,k] <- sum(factorized)/factors_params[[k-1]][11]
}

aracne_ranger_res <- res_df

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})

for (j in 1:10){ 
  regs <- aracne_mat[which(aracne_mat[,2] == j),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == j),]$MI, decreasing = TRUE, index.return = TRUE)
  regs <- regs[b$ix]
  regs <- regs[which(regs > 12)]
  regs <- colnames(parameter_mat)[regs[1:num1[j]]]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
      print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    }
  }
}

aracne_mat <- read.table(file = list.files()[grep('norm',list.files())][k],sep="\t", fill = TRUE, header = TRUE)

load(list.files()[grep('New',list.files())][k])
print(list.files()[grep('New',list.files())][k])

all_regs <- c()
num1 <- c(6,6,11,0,3,13,0,7,8,1)
num2 <- c(1,7,9,0,4,10,0,4,6,3)
num3 <- c(13,4,6,0,7,2,0,11,2,3)
num4 <- c(9,13,4,0,5,5,0,1,6,3)
num5 <- c(12,13,26,0,20,7,0,2,8,1)
num6 <- c(1,6,3,0,18,3,0,8,4,7)
for (i in 1:10){
  regs <- aracne_mat[which(aracne_mat[,2] == i),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == i),]$MI, decreasing = TRUE, index.return = TRUE)
  b <- regs[b$ix]
  b <- b[which(b > 12)]
  names(b) <- colnames(parameter_mat)[b]
  all_regs <- c(all_regs,names(b))
}
all_regs_unique <- unique(all_regs)
sort(table(all_regs))

names(which(sort(table(all_regs)) > 1))
sort(table(all_regs))[which(sort(table(all_regs)) > 1)]
#View(kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5])
sort(table(all_regs))[which(sort(table(all_regs)) > 1)][which(names(which(sort(table(all_regs)) > 1)) %in% kannot[which(rownames(kannot) %in% names(which(sort(table(all_regs)) > 1))),1:5]$gene)]

enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs_unique,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})

#[1:num1[j]]
for (j in 1:10){ 
  regs <- aracne_mat[which(aracne_mat[,2] == j),]$Regulator
  b <- sort(aracne_mat[which(aracne_mat[,2] == j),]$MI, decreasing = TRUE, index.return = TRUE)
  regs <- regs[b$ix]
  regs <- regs[which(regs > 12)]
  regs <- colnames(parameter_mat)[regs[1:num1[j]]]
  print(colnames(parameter_mat)[j])
  for ( i in 1:length(path.map.s4)){
    enrich_res <- enrichment_fn(regs,path.map.s4[[i]],retPvalOnly=T)
    if (! is.na(enrich_res) && enrich_res <= 0.05){
      print(i)
      print(enrich_res)
      enriched_paths <- c(enriched_paths,i)
      print(names(path.map.s4)[i])
      print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    }
  }
}
kannot[which(kannot$gene %in% c("APZ42_013661"   ,    "APZ42_013795")),2:3]



net_dec_energ_norm <- read.table(file = 'network_dec_energ_norm.txt',sep="\t", fill = TRUE, header = TRUE)
net_dec_energ_noTarg <- read.table(file = 'network_dec_energ_noTarg.txt',sep="\t", fill = TRUE, header = TRUE)
net_dec_energ_noTarg_p1.1 <- read.table(file = 'network_dec_energ_noTrag_p1.1.txt',sep="\t", fill = TRUE, header = TRUE)


#heatmap r
dec_ener_net <- read.table("network_dec_feed_norm.txt",sep="\t", fill = TRUE, header = TRUE)

names <- colnames(parameter_mat)
numbers <- 1:length(colnames(parameter_mat))
names_numb <- cbind(names,numbers)
names_numb[,2]<- as.numeric(names_numb[,2])

library(igraph)
adj_mat_dec_ener<- sparseMatrix(i = dec_ener_net[,1], j = dec_ener_net[,2], x = as.numeric(dec_ener_net[,3]), dims=list(27362,27362))
colnames(adj_mat_dec_ener) <- names_numb[,1]
rownames(adj_mat_dec_ener) <- names_numb[,1]
adj_mat_dec_ener_graph <- graph_from_adjacency_matrix(adj_mat_dec_ener, mode = 'directed', weighted = TRUE)


for (i in 1:12){
  regs <- dec_ener_net[which(dec_ener_net[,2] == i),]$Regulator
  regs  <- c(regs,dec_ener_net[which(dec_ener_net[,1] == i),]$Target)
  regs <- unique(regs)
  if (length(regs[which(regs > 12)]) > 0){
    lms <- lm(parameter_mat[,i] ~ as.matrix(parameter_mat[,regs[which(regs > 12)]]))
    sum_lms <- summary(lms)
    r_squared <- sum_lms$adj.r.squared
    #print(colnames(parameter_mat)[i])
    print(r_squared)
    print(length(regs[which(regs > 12)]))
  }
  else {
    #print(colnames(parameter_mat)[i])
    print('no predictors')
    print(length(regs[which(regs > 12)]))
  }
}

dec_feed_net <- dec_ener_net
all_regs <- c()
for (i in 1:10){
  regs <- dec_feed_net[which(dec_feed_net[,1] == i),]$Target
  b <- sort(dec_feed_net[which(dec_feed_net[,1] == i),]$MI, decreasing = TRUE, index.return = TRUE)
  b <- regs[b$ix]
  b <- b[which(b > 12)]
  names(b) <- colnames(parameter_mat)[b]
  all_regs <- c(all_regs,names(b))
}
all_regs <- unique(all_regs)
#all_regs <- colnames(parameter_mat)[all_regs]
signifGenes <- all_regs
pathwayGenes <- path.map.s4[[2]]


enriched_paths <- c()
res_enrich <- sapply(1:length(path.map.s4), function(i){
  #print(i)
  enrich_res <- enrichment_fn(all_regs,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
  }
  return(enriched_paths)
  
})

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_decreased feeding_2.RData")
heatmap(ab)
ab <- parameter_mat[,c("stress_ls", "the structural length, cm", "energy reserve, J", "reproduction buffer, J", "internal concentration, nM", "food density", "maximum structural length", "cumulative embryos, #", "days", "conc" ,"APZ42_001840", "APZ42_013303", "APZ42_021677")]
typeof(ab)
ab <- as.data.frame(ab)
ab <-t(ab)
?heatmap

library(corrplot)
corrplot(cor(ab), tl.cex = 0.5)
?corrplot
rownames(ab) <- sapply(1:length(ab$days), function(i){
  x <- sprintf('%s_%s', ab$days[i], i)
})
rownames(ab) <- sapply(1:length(ab$conc), function(i){
  x <- sprintf('%s_%s_%s', ab$days[i], ab$conc[i], i)
})
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
hc <- hclust(as.dist(1-cor(as.matrix(ab[order(ab$conc, ab$days),]))))
plot(hc)
hc2 <- hclust(as.dist(1-cor(t(as.matrix(ab[order(ab$conc, ab$days),])))))
plot(hc2)
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc2), scale = 'col',margins = c(12,3))


heatmap.2(as.matrix(ab[,order( ab['conc',], ab['days',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Concentration first')
heatmap.2(as.matrix(ab[,order(ab['days',],  ab['conc',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Days first')


load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_hazard during oogenesis_2.RData")

ab <- parameter_mat[,c("stress_ls", "the structural length, cm", "energy reserve, J", "reproduction buffer, J", "internal concentration, nM", "food density", "maximum structural length", "cumulative embryos, #", "days", "conc" ,"APZ42_020354", "APZ42_017398")]
ab <- as.data.frame(ab)
ab <-t(ab)
rownames(ab) <- sapply(1:length(ab$conc), function(i){
  x <- sprintf('%s_%s_%s', ab$days[i], ab$conc[i], i)
})
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
hc <- hclust(as.dist(1-cor(as.matrix(ab[order(ab$conc, ab$days),]))))
plot(hc)
hc2 <- hclust(as.dist(1-cor(t(as.matrix(ab[order(ab$conc, ab$days),])))))
plot(hc2)
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Colv = as.dendrogram(hc), Rowv = as.dendrogram(hc2), scale = 'col',margins = c(12,3))
corrplot(cor(ab), tl.cex = 0.8)

heatmap.2(as.matrix(ab[,order( ab['conc',], ab['days',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Concentration first')
heatmap.2(as.matrix(ab[,order(ab['days',],  ab['conc',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Days first')




#hippo signaling pathway
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_increased maintenance costs_2.RData")
ab <- parameter_mat[,c("stress_ls", "the structural length, cm", "energy reserve, J", "reproduction buffer, J", "internal concentration, nM", "food density", "maximum structural length", "cumulative embryos, #", "days", "conc" ,"APZ42_025278", "APZ42_026254")]
ab <- as.data.frame(ab)
ab <-t(ab)
rownames(ab) <- sapply(1:length(ab$conc), function(i){
  x <- sprintf('%s_%s_%s', ab$days[i], ab$conc[i], i)
})
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
hc <- hclust(as.dist(1-cor(as.matrix(ab[order(ab$conc, ab$days),]))))
plot(hc)
hc2 <- hclust(as.dist(1-cor(t(as.matrix(ab[order(ab$conc, ab$days),])))))
plot(hc2)
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc2), scale = 'col',margins = c(12,3))
corrplot(cor(ab), tl.cex = 0.8)

heatmap.2(as.matrix(ab[,order( ab['conc',], ab['days',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Concentration first')
heatmap.2(as.matrix(ab[,order(ab['days',],  ab['conc',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Days first')



load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_increased allocation to soma_2.RData")
ab <- parameter_mat[,c("stress_ls", "the structural length, cm", "energy reserve, J", "reproduction buffer, J", "internal concentration, nM", "food density", "maximum structural length", "cumulative embryos, #", "days", "conc" ,"APZ42_021840", "APZ42_021840")]
ab <- as.data.frame(ab)
ab <-t(ab)
rownames(ab) <- sapply(1:length(ab$conc), function(i){
  x <- sprintf('%s_%s_%s', ab$days[i], ab$conc[i], i)
})
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
hc <- hclust(as.dist(1-cor(as.matrix(ab[order(ab$conc, ab$days),]))))
plot(hc)
hc2 <- hclust(as.dist(1-cor(t(as.matrix(ab[order(ab$conc, ab$days),])))))
plot(hc2)
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc2), scale = 'col',margins = c(12,3))
corrplot(cor(ab), tl.cex = 0.8)

?heatmap.2
ab <- t(ab)

heatmap.2(as.matrix(ab[,order( ab['conc',], ab['days',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Concentration first')
heatmap.2(as.matrix(ab[,order(ab['days',],  ab['conc',])]), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Days first')

#-----------------------------------------------------------------------------
# lasso ranger 2 PMOA results
ab <- lasso_ranger_res2PMOA[params_,]
ab <- ab[,2:5]
for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- lasso_ranger_res2PMOA[params_,]$DEB_Param
ab <- as.data.frame(t(ab))
ab <- rbind(c(0,0,0,0,0,0,0,0,0,0),ab)
ab <- rbind(c(1,1,1,1,1,1,1,1,1,1),ab)



library(RColorBrewer)
coul <- brewer.pal(4, "Set2")
colors_border <- coul
library(scales)
colors_in <- alpha(coul,0.1)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart( ab,
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.8 
)
# Add a legend
legend(x=0.9, y=0.8, legend = rownames(ab[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=0.8, pt.cex=2)


#-------------------------------------------------------------
#Aracne caret results

ab <- aracne_caret_res[params_,]
ab <- ab[,2:8]
for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- aracne_caret_res[params_,]$DEB_Param
ab <- as.data.frame(t(ab))
ab <- rbind(c(0,0,0,0,0,0,0,0,0,0),ab)
ab <- rbind(c(1,1,1,1,1,1,1,1,1,1),ab)

library(RColorBrewer)
coul <- brewer.pal(7, "Set2")
colors_border <- coul
library(scales)
colors_in <- alpha(coul,0.1)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart( ab,
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=0.8, y=-0.5, legend = rownames(ab[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=0.8, pt.cex=2)



#-------------------------------------------------------------
#Aracne ranger results
#aracne_ranger_res <-aracne_ranger_res[-1,]
ab <- aracne_ranger_res[params_,]
ab <- as.matrix(ab[,2:8])

for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- aracne_ranger_res[params_,]$DEB_Param
ab <- as.data.frame(t(ab))
ab <- rbind(c(0,0,0,0,0,0,0,0,0,0),ab)
ab <- rbind(c(1,1,1,1,1,1,1,1,1,1),ab)


library(RColorBrewer)
coul <- brewer.pal(7, "Set2")
colors_border <- coul
library(scales)
colors_in <- alpha(coul,0.1)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart( ab,
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.7 
)

# Add a legend
legend(x=0.8, y=-0.5, legend = rownames(ab[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=0.8, pt.cex=2)





#-------------------------------------------------------------
#Lasso caret results

ab <- Lasso_caret_res[params_,]
ab <- as.matrix(ab[,2:8])

for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- Lasso_caret_res[params_,]$DEB
ab <- as.data.frame(t(ab))
ab <- rbind(c(0,0,0,0,0,0,0,0,0,0),ab)
ab <- rbind(c(1,1,1,1,1,1,1,1,1,1),ab)


library(RColorBrewer)
coul <- brewer.pal(7, "Set2")
colors_border <- coul
library(scales)
colors_in <- alpha(coul,0.1)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart( ab,
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.7 
)

# Add a legend
legend(x=0.8, y=-0.5, legend = rownames(ab[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=0.8, pt.cex=2)


#-------------------------------------------------------------
#Lasso ranger results
#lasso_ranger_res <- lasso_ranger_res[,-c(2,4,7,11)]
ab <- lasso_ranger_res[params_,]
ab <- as.matrix(ab[,2:8])

for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- lasso_ranger_res[params_,]$DEB_Param
ab <- as.data.frame(t(ab))
ab <- rbind(c(0,0,0,0,0,0,0,0,0,0),ab)
ab <- rbind(c(1,1,1,1,1,1,1,1,1,1),ab)


library(RColorBrewer)
coul <- brewer.pal(7, "Set2")
colors_border <- coul
library(scales)
colors_in <- alpha(coul,0.1)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
radarchart( ab,
            pcol=colors_border , pfcol=colors_in , plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=0.7 
)

# Add a legend
legend(x=0.8, y=-0.5, legend = rownames(ab[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=0.8, pt.cex=2)

library(fmsb)



## ggradar

#Lasso ranger results
#lasso_ranger_res <- lasso_ranger_res[,-c(2,4,7,11)]
params_ <- c(3,5,7,11,13,17,19,21,23)
ab <- lasso_ranger_res[params_,]
ab <- as.matrix(ab[,2:9])

for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- lasso_ranger_res[params_,]$DEB_Param
ab <- as.data.frame(t(ab))

group = rownames(ab)

ab <- cbind(group, ab)
ab
cd <- lasso_ranger_res[c(4,6,8,12,14,18,20,22),]

ps <- c()
for (i in 2:dim(cd)[2]){
  ps <- c(ps,as.vector(cd[,i]),1,cd[1,i])
}

ps <- unlist(ps)
ggradar(
  ab, 
  grid.label.size = 5,
  axis.label.size	= 5,
  base.size = 1,
  values.radar = c("0", "0.5", "1"),
  grid.min = 0, grid.mid = 0.5, grid.max = 1,
  group.line.width = 1, 
  group.point.size = ps/10,
  legend.text.size = 10,
  legend.position = 'bottom'
)



#Lasso caret results
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/Model_results/lasso_caret_res.RData")
library(ggradar)

params_ <- c(3,5,7,11,13,17,19,21,23)
ab <- lasso_caret_res[params_,]
ab <- as.matrix(ab[,2:9])

for (i in 1:length(rownames(ab))){
  ab[i,] <- as.numeric(ab[i,])
}
rownames(ab) <- lasso_caret_res[params_,]$DEB_Param
ab <- as.data.frame(t(ab))

group = rownames(ab)

ab <- cbind(group, ab)
ab
cd <- lasso_caret_res[c(4,6,8,12,14,18,20,22),]

ps <- c()
for (i in 2:dim(cd)[2]){
  ps <- c(ps,as.vector(cd[,i]),1,cd[1,i])
}

ps <- unlist(ps)
ggradar(
  ab, 
  grid.label.size = 5,
  axis.label.size	= 4,
  base.size = 1,
  values.radar = c("0", "0.5", "1"),
  grid.min = 0, grid.mid = 0.5, grid.max = 1,
  group.line.width = 1, 
  group.point.size = ps/10,
  legend.text.size = 10,
  legend.position = 'bottom'
)

params_1 <- c(3,5,7,11,13,17,19,21,23)
params_2 <- c(1,3,5,7,9)
for (i in c(1:20,22:25)){
  group_colors <- c('#ff303099','#4c709399','#7f6a7c99','#e2b9b399','#ffae1999','#60ce8099','#54f4d999','#0c472a99')
  if (i == 16){
    
    params_ <- params_2
    pars_ <- c(2,4,6,8)
  }
  
  #else{
    
    params_ <- params_1
    pars_ <- c(4,6,8,12,14,18,20,22)
  #}
  
  print(list.files()[i])
  p_df <- load(list.files()[i])
  p_df <- get(p_df)
 
  ab <- p_df[params_,]
  ab <- as.matrix(ab[,2:dim(ab)[2]])
  ab <- as.matrix(ab[,1:dim(ab)[2]])
  for (j in 1:length(rownames(ab))){
    ab[j,] <- as.numeric(ab[j,])
  }
  ab <- apply(ab, 2,            # Specify own function within apply
              function(x) as.numeric(as.character(x)))
  rownames(ab) <- p_df[params_,]$DEB_Param
  rownames(ab) <- rownames(p_df[params_,])
  ab <- as.data.frame(t(ab))
  rownames(ab) <- c('Inc reprodcution cost & decreased energy conductance','Inc reproduction cost & deccreased feeding','Inc reproduction cost & increased allocation to soma', 'Inc reprodcution cost & increased cost for growth and reproduction','Inc reprodction cost & increased hazard during oogenesis', 'Inc reproduction cost & increased maintenace cost')
  rownames(ab)[c(1,4)] <- c('Increased cost for reprodcution','decreased fraction to soma')
  group = factor(rownames(ab), levels = rownames(ab))
  
  ab <- cbind(str_wrap(group,20), ab)
  ab
  ab$`str_wrap(group, 20)`<- factor(ab$`str_wrap(group, 20)`, levels = ab$`str_wrap(group, 20)`)
  cd <- p_df[pars_,-1]
 
  ps <- c()
  for (k in 1:dim(cd)[2]){
    ps <- c(ps,as.vector(cd[,k]),1,cd[1,k])
  }
  
  ps <- unlist(ps)
  ps <- as.numeric(ps)
  for (t in 2:dim(ab)[2]){
    print(ab[i,-1])
    print(log(ab[i,-1]+1))
    ab[,t] <- log2(ab[,t]+0.01)
    ab[,t] <- rescale(as.matrix(ab[,t]),c(0,1))
  }
  
  g <- ggradar(
        ab, 
        grid.label.size = 6,
        axis.label.size	= 6,
        axis.labels = str_wrap(colnames(ab)[-1],10),
        base.size = 0.5,
        values.radar = c("0", "0.5", "1"),
        grid.min = 0, grid.mid = 0.5, grid.max = 1,
        group.line.width = 1, 
        group.point.size = log(ps),
        legend.text.size = 12,
        legend.position = 'right',
        gridline.mid.colour = '#0000006f',
        group.colours	= group_colors[1:dim(ab)[1]]
        
        
      ) + theme(legend.key.height = unit(1.5,'cm'))
  print(list.files()[i])
  ggsave(file = sprintf('Plots3/%s_.png',list.files()[i]), device = 'png',width = 10, height = 10)
}

group_colors <- c('#89898980','#ff303099','#4c709399','#7f6a7c99','#e2b9b399','#ffae1999','#60ce8099','#54f4d999','#0c472a99')
add.transparency <- function(x,perc){
  if(perc > 1 & perc <= 100){
    perc <- floor(255*(perc/100))
  }else if(perc >= 0 & perc <= 1){
    perc <- floor(255*perc)
  }else{
    stop("Error wrong percentage given.")
  }
  x <- paste(substr(x,1,7),as.hexmode(perc),sep="")
  return(x)
}


group_colors <- sapply(group_colors,function(x){
  g <- add.transparency(x,perc = 100)
  return(g)
})

vars <- c()
for (i in 1:dim(aracne_ranger_res)[2]){ 
  
  a <- sum(as.numeric(aracne_ranger_res[pars_,i]))/ length(aracne_ranger_res[pars_,i])
  vars <- c(vars, a)
  }

for (i in 1:dim(ab)[1]){
  #print(ab[i,-1])
  #print(log(ab[i,-1]+1))
  ab[i,-1] <- log2(ab[i,-1]+1)
}


