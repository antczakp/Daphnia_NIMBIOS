
load("C:/Uni/Master/Project Modul2/Daphnia_NIMBIOS/Daphnia_NIMBIOS/DEB Modelling Outputs/model outputs/New_Mat_cost_for_reproduction.RData")
genes <- c("APZ42_014481" ,"APZ42_015046", "APZ42_028914",'days','conc')
length(genes)
length(unique(genes))

ab <- parameter_mat[,unique(genes)]
typeof(ab)
ab <- as.data.frame(ab)
#ab <-t(ab)
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
heatmap(as.matrix(t(ab[order(ab$conc, ab$days),])), Rowv = order(ab$days) ,scale = 'row', margins = c(12,3))
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Rowv = order(ab$conc) ,scale = 'col', margins = c(12,3))
hc <- hclust(as.dist(1-cor(as.matrix(ab[order(ab$conc, ab$days),]))))
plot(hc)
hc2 <- hclust(as.dist(1-cor(t(as.matrix(ab[order(ab$conc, ab$days),])))))
plot(hc2)
heatmap(as.matrix(ab[order(ab$conc, ab$days),]), Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc2), scale = 'col',margins = c(12,3))


ab$days <- factor(ab$days, levels = c('6','13','17','24'))

for (i in 1:dim(ab)[2]){
  
  ab[,i] <- as.numeric(ab[,i])
  
}
ab$conc <- log(ab$conc + 0.001)
colnames(ab)[5] <- 'log concentration'
heatmap.2(as.matrix(t(ab[order(ab$`log concentration`, ab$days ),c(1:3)])), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = as.character(parameter_mat$conc[order(parameter_mat$conc,parameter_mat$days)]), cexRow = 1.5, cexCol = 1.5, xlab = 'concentration nM [Hg]')
heatmap.2(t(as.matrix(ab[order(ab$days, ab$conc),])), Colv = FALSE, col = greenred(100), scale = 'row', margins = c(6,12), trace="none", labCol = FALSE, main = 'Samples sorted by Days first')


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


