#welche genes in increased cost of reproduction?

##look at model of lasso ranger 

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_res.RData")
lasso_ranger_res[,1]
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/lasso_stab_new_res0.05_cost_rep.RData")
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/New_Mat_cost_for_reproduction.RData")

lasso_stab_new_res <- lasso_stab_new_res[,-1]
#dim(lasso_stab_new_res)

colnames(lasso_stab_new_res) <- colnames(parameter_mat)[1:12]
rownames(lasso_stab_new_res) <- colnames(parameter_mat)

rownames(lasso_ranger_res)
nr_vars <- c(4,6,8,10,12,14,16,18,20,22)
genes <- list()
all_genes <- c()
for (i in 1:10){
  print(i)
  if (lasso_ranger_res[nr_vars[i],1] > 0){
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    genes[[i]] <- names(b)[1:lasso_ranger_res[nr_vars[i],1]]
    all_genes <- c(all_genes, names(b)[1:lasso_ranger_res[nr_vars[i],1]])
  }
}

names(genes) <- rownames(lasso_ranger_res)[-c(1:2,nr_vars,23:25)]

mean_vars <- c()
for (k in 1:dim(lasso_ranger_res)[2]){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- lasso_ranger_res[nr_vars[i],k]
  }
  mean_vars[k] <- sum(factorized)/8
}

lasso_ranger_res<- rbind(lasso_ranger_res,mean_vars)

save(lasso_ranger_res, file = "lasso_ranger_res.RData")


enrichment_df <- data.frame(colnames = 'KEGGID','P_value','ratio')
for (i in 1:length(path.map.s4)){
  
  enrich_res <- enrichment_fn(all_genes,path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    new <- c(names(path.map.s4)[i],enrich_res[1],c(length(enrich_res[-1])/length(path.map.s4[[i]])))
    #enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    enrichment_df <- rbind(enrichment_df,new)
  }
  
}
colnames(enrichment_df) <- enrichment_df[1,]
enrichment_df <- enrichment_df[-1,]
dim(enrichment_df)
enrichment_df[3,]
enrichment_df[4,]
enrichment_df <- enrichment_df[1:3,]
enrichment_df  <-  enrichment_df[order(enrichment_df$P_value),]
save(enrichment_df, file = 'enrichment_df_inc_cost_rep.RData')

de_genes_inc_cost_rep <- read_excel("C:/Uni/Master/Project Modul2/de_genes_inc_cost_rep.xlsx")
#de_genes_inc_cost_rep <- de_genes_inc_cost_rep[-1,]
de_genes_inc_cost_rep$Name <-factor(de_genes_inc_cost_rep$Name, levels = de_genes_inc_cost_rep$Name[order(de_genes_inc_cost_rep$ratio)])
de_genes_inc_cost_rep$P_value <-factor(round(de_genes_inc_cost_rep$P_value,3), levels = unique(sort(round(de_genes_inc_cost_rep$P_value,3))))
p1 <- ggplot(data = de_genes_inc_cost_rep, aes(x=Name,y= ratio))+
  geom_bar(stat="identity", aes(fill = P_value))+
  scale_fill_brewer(palette = "Reds", direction = -1)+
  labs( title = "Predictive genes of Increased cost of reporduction")+
  coord_flip()
plot(p1)



enrichment_df <- data.frame(colnames = 'Param' , 'KEGGID', 'Genes','P_value','ratio','Nr. of Genes Paths')
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    print(k)
    for (i in 1:length(path.map.s4)){
      print(i)
      enrich_res <- enrichment_fn(genes[[k]],path.map.s4[[i]],retPvalOnly=T)
      if (! is.na(enrich_res) && enrich_res <= 1.0){
        
        print(enrich_res)
        genes_enrich <-  enrich_res[-1]
        genes_enrich[which(!is.na(match(genes_enrich,kannot$gene)))] <- kannot$description[which(!is.na(match(genes_enrich,kannot$gene)))]
        genes_enrich <- paste(genes_enrich, collapse = ', ')
        new <- c(names(genes)[k],names(path.map.s4)[i],genes_enrich, enrich_res[1],c(length(enrich_res[-1])/length(path.map.s4[[i]])),length(enrich_res[-1]))
        #enriched_paths <- c(enriched_paths,i)
        print(names(path.map.s4)[i])
        print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
        enrichment_df <- rbind(enrichment_df,new)
      }
    }
  }
}



gene_df <- data.frame('Parameter'= c(),
                      'Nr. of genes' = c(),
                      'KEGG paths' = c())
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    
    print(k)
    gene_df <- rbind(gene_df,rep(names(genes)[k],length(which(kannot$gene %in% genes[[k]]))),kannot[which(kannot$gene %in% genes[[k]]),2:3])

  }
}
genes_over_conc <- c()
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    
    print(k)
    genes_over_conc[k] <- length(which(genes[[k]] %in% de_genes))
    
  }
}


##look at model of aracne ranger

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_ranger_res.RData")
aracne_ranger_res[,2]
aracne_mat <- read.table(file = 'train_model/1.original/aracne_cost_reproduction.txt',sep="\t", fill = TRUE, header = TRUE)
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/New_Mat_cost_for_reproduction.RData")


nr_vars <- c(4,6,8,10,12,14,16,18,20,22)
genes <- list()
all_genes <- c()
for (i in 1:10){
  print(i)
  if (aracne_ranger_res[nr_vars[i],1] > 0){
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
    genes[[i]] <- names(b)[1:aracne_ranger_res[nr_vars[i],1]]
    all_genes <- c(all_genes,names(b)[1:aracne_ranger_res[nr_vars[i],1]])
  }
}

names(genes) <- rownames(aracne_ranger_res)[-c(1:2,nr_vars,23:26)]

mean_vars <- c()
for (k in 1:dim(aracne_ranger_res)[2]){
  
  factorized <- c()
  for (i in 1:10){
    factorized[i] <- aracne_ranger_res[nr_vars[i],k]
  }
  mean_vars[k] <- sum(as.numeric(factorized))/8
}
mean_vars[1] <- 'mean_vars'
aracne_ranger_res<- rbind(aracne_ranger_res,mean_vars)

save(aracne_ranger_res, file = "aracne_ranger_res.RData")


enrichment_df <- data.frame(colnames = 'KEGGID','P_value','ratio')
for (i in 1:length(path.map.s4)){
  
  enrich_res <- enrichment_fn(unique(all_genes),path.map.s4[[i]],retPvalOnly=T)
  if (! is.na(enrich_res) && enrich_res <= 0.05){
    print(i)
    print(enrich_res)
    new <- c(names(path.map.s4)[i],enrich_res[1],c(length(enrich_res[-1])/length(path.map.s4[[i]])))
    #enriched_paths <- c(enriched_paths,i)
    print(names(path.map.s4)[i])
    print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
    enrichment_df <- rbind(enrichment_df,new)
  }
  
}

enrichment_df <- data.frame(colnames = 'Param' , 'KEGGID','P_value','ratio')
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    print(k)
    for (i in 1:length(path.map.s4)){
      print(i)
      enrich_res <- enrichment_fn(genes[[k]],path.map.s4[[i]],retPvalOnly=T)
      if (! is.na(enrich_res) && enrich_res <= 0.9){
        
        print(enrich_res)
        new <- c(names(genes)[k],names(path.map.s4)[i],enrich_res[1],c(length(enrich_res[-1])/length(path.map.s4[[i]])))
        #enriched_paths <- c(enriched_paths,i)
        print(names(path.map.s4)[i])
        print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
        enrichment_df <- rbind(enrichment_df,new)
      }
    }
  }
}

colnames(enrichment_df) <- enrichment_df[1,]
enrichment_df <- enrichment_df[-1,]
save(enrichment_df, file = 'enrichment_df_aracne_cost_rep.RData')


Name <- df1$Name[match(enrichment_df$X.KEGGID., df1$KEGGID)]
enrichment_df_new <- cbind(enrichment_df[,1:4],Name,enrichment_df[5:8])
write.xlsx(enrichment_df_new,'cost_rep_genes_and_path5.xlsx')

enrichment_df <- data.frame(colnames = 'Param' ,'Nr. of Genes', 'Genes', 'KEGGID', 'Genes','P_value','ratio','Nr. of Genes Paths')
for (k in 1:length(genes)) {
  #print(k)
  if (length(genes[[k]]) > 1){
    print(k)
    for (i in 1:length(path.map.s4)){
      print(i)
      enrich_res <- enrichment_fn(genes[[k]],path.map.s4[[i]],retPvalOnly=T)
      if (! is.na(enrich_res) && enrich_res <= 1.0){
        
        print(enrich_res)
        genes_enrich <-  enrich_res[-1]
        genes_enrich[which(genes_enrich %in% kannot$gene)] <- kannot$description[which(kannot$gene %in% genes_enrich)]
        genes_enrich <- paste(genes_enrich, collapse = ', ')
        genes_param <- genes[[k]]
        genes_param[which(!is.na(match(genes_param,kannot$gene)))] <- kannot$description[which(!is.na(match(genes_param,kannot$gene)))]
        genes_param <- paste(genes_param, collapse = ', ')
        new <- c(names(genes)[k],length(genes[[k]]),genes_param,names(path.map.s4)[i],genes_enrich, enrich_res[1],c(length(enrich_res[-1])/length(path.map.s4[[i]])),length(enrich_res[-1]))
        #enriched_paths <- c(enriched_paths,i)
       # print(names(path.map.s4)[i])
       # print(kannot[which(kannot$gene %in% enrich_res[-1]),2:3])
        enrichment_df <- rbind(enrichment_df,new)
      }
    }
  }
}
genes_path <- c()
for (i in 1:dim(cost_rep_genes_and_path3)[1]){
  
  genes_path <- c(genes_path,paste(cost_rep_genes_and_path3$Name[i],' (',cost_rep_genes_and_path3$X.Nr..of.Genes.Paths.[i],')'))
  
}
gene_df <- data.frame()
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    
    print(k)
    gene_df <- rbind(gene_df,rep(names(genes)[k],length(which(kannot$gene %in% genes[[k]]))),kannot[which(kannot$gene %in% genes[[k]]),2:3])
    
  }
}
