#welche genes in decreased allocation to soma?

##look at model of lasso ranger 

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/lasso_ranger_res.RData")
lasso_ranger_res[,4]
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/lasso_stab_new_res0.05_dec_fract.RData")
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/New_Mat_decreased_fraction_to_soma.RData")

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
  if (lasso_ranger_res[nr_vars[i],4] > 0){
    a <- lasso_stab_new_res[-c(1:12),i]
    b <- sort(a, decreasing = TRUE )
    genes[[i]] <- names(b)[1:lasso_ranger_res[nr_vars[i],4]]
    all_genes <- c(all_genes, names(b)[1:lasso_ranger_res[nr_vars[i],4]])
  }
}

names(genes) <- rownames(lasso_ranger_res)[-c(1:2,nr_vars,23:24)]

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


de_genes_decreased_fraction_to_soma <- read_excel("C:/Uni/Master/Project Modul2/de_genes_decreased_fraction_to_soma.xlsx")
#de_genes_decreased_fraction_to_soma <- de_genes_decreased_fraction_to_soma[1,]
de_genes_decreased_fraction_to_soma$Name <-factor(de_genes_decreased_fraction_to_soma$Name, levels = de_genes_decreased_fraction_to_soma$Name[order(de_genes_decreased_fraction_to_soma$ratio)])
de_genes_decreased_fraction_to_soma$P_value <-factor(round(de_genes_decreased_fraction_to_soma$P_value,4), levels = unique(sort(round(de_genes_decreased_fraction_to_soma$P_value,4))))
p1 <- ggplot(data = de_genes_decreased_fraction_to_soma, aes(x=Name,y= ratio))+
  geom_bar(stat="identity", aes(fill = P_value))+
  scale_fill_brewer(palette = "Reds", direction = -1)+
  labs( title = "Predictive genes of Decreased_fraction to soma")+
  coord_flip()
plot(p1)


enrichment_df <- data.frame(colnames = 'Param' , 'KEGGID','P_value','ratio')
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    print(k)
    for (i in 1:length(path.map.s4)){
      print(i)
      enrich_res <- enrichment_fn(genes[[k]],path.map.s4[[i]],retPvalOnly=T)
      if (! is.na(enrich_res) && enrich_res <= 0.05){
        
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
save(enrichment_df, file = 'enrichment_df_dec_fraction_lasso.RData')


gene_df <- data.frame()
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    
    print(k)
    gene_df <- rbind(gene_df,rep(names(genes)[k],length(which(kannot$gene %in% genes[[k]]))),kannot[which(kannot$gene %in% genes[[k]]),2:3])
    
  }
}




##look at model of aracne ranger

load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/aracne_ranger_res.RData")
aracne_ranger_res[,5]
aracne_mat <- read.table(file = 'train_model/1.original/aracne_dec_fraction.txt',sep="\t", fill = TRUE, header = TRUE)
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/train_model/1.original/New_Mat_decreased_fraction_to_soma.RData")


nr_vars <- c(4,6,8,10,12,14,16,18,20,22)
genes <- list()
all_genes <- c()
for (i in 1:10){
  print(i)
  if (aracne_ranger_res[nr_vars[i],2] > 0){
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
    genes[[i]] <- names(b)[1:aracne_ranger_res[nr_vars[i],2]]
    all_genes <- c(all_genes,names(b)[1:aracne_ranger_res[nr_vars[i],5]])
  }
}

names(genes) <- aracne_ranger_res[-c(1:2,nr_vars,23:24),1]


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
enrichment_df[2,]
enrichment_df[3,]
enrichment_df <- enrichment_df[1:2,]
enrichment_df  <-  enrichment_df[order(enrichment_df$P_value),]
save(enrichment_df, file = 'enrichment_df_decreased_fraction_to_soma.RData')

de_genes_dec_fraction_to_soma_aracne <- read_excel("C:/Uni/Master/Project Modul2/de_genes_dec_fraction_to_soma_aracne.xlsx")
#de_genes_decreased_fraction_to_soma <- de_genes_decreased_fraction_to_soma[1,]
de_genes_dec_fraction_to_soma_aracne$Name <-factor(de_genes_dec_fraction_to_soma_aracne$Name, levels = de_genes_dec_fraction_to_soma_aracne$Name[order(de_genes_dec_fraction_to_soma_aracne$ratio)])
de_genes_dec_fraction_to_soma_aracne$P_value <-factor(round(de_genes_dec_fraction_to_soma_aracne$P_value,4), levels = unique(sort(round(de_genes_dec_fraction_to_soma_aracne$P_value,4))))
p1 <- ggplot(data = de_genes_dec_fraction_to_soma_aracne, aes(x=Name,y= ratio))+
  geom_bar(stat="identity", aes(fill = P_value))+
  scale_fill_brewer(palette = "Reds", direction = -1)+
  labs( title = "Predictive genes of Decreased_fraction to soma")+
  coord_flip()
plot(p1)


enrichment_df <- data.frame(colnames = 'Param' , 'KEGGID','P_value','ratio')
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    print(k)
    for (i in 1:length(path.map.s4)){
      print(i)
      enrich_res <- enrichment_fn(genes[[k]],path.map.s4[[i]],retPvalOnly=T)
      if (! is.na(enrich_res) && enrich_res <= 0.5){
        
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
save(enrichment_df, file = 'enrichment_df_aracne_dec_fraction.RData')


gene_df <- data.frame()
for (k in 1:length(genes)) {
  if (length(genes[[k]]) > 1){
    
    print(k)
    gene_df <- rbind(gene_df,rep(names(genes)[k],length(which(kannot$gene %in% genes[[k]]))),kannot[which(kannot$gene %in% genes[[k]]),2:3])
    
  }
}
