library(stringr)
library(ggplot2)
library(ggradar)
params_1 <- c(3,5,7,11,13,17,19,21,23)
params_2 <- c(1,3,5,7,9)
for (i in c(1:20,22:25)){
  group_colors <- c('#ff303099','#4c709399','#7f6a7c99','#e2b9b399','#ffae1999','#60ce8099','#54f4d999','#0c472a99')
  if (i == 16){
    
    params_ <- params_2
    pars_ <- c(2,4,6,8)
  }
  
  else{
  
    params_ <- params_1
    pars_ <- c(4,6,8,12,14,18,20,22)
  
  }
  print(list.files()[i])
  p_df <- load(list.files()[i])
  p_df <- get(p_df)
  
  ab <- p_df[params_,]
  ab <- as.matrix(ab[,2:dim(ab)[2]])
  
  for (j in 1:length(rownames(ab))){
    ab[j,] <- as.numeric(ab[j,])
  }
  ab <- apply(ab, 2,          
              function(x) as.numeric(as.character(x)))
  rownames(ab) <- p_df[params_,]$DEB_Param
  ab <- as.data.frame(t(ab))
  group = factor(rownames(ab), levels = rownames(ab))
  
  ab <- cbind(group, ab)
  ab$group<- factor(ab$group, levels = ab$group)
  cd <- p_df[pars_,-1]
  cd <- as.data.frame(cd)
  ps <- c()
  for (k in 1:dim(cd)[2]){
    ps <- c(ps,as.vector(cd[,k]),1,cd[1,k])
  }
  
  ps <- unlist(ps)
  ps <- as.numeric(ps)
  
  #log scaled
  
  #for (t in 2:dim(ab)[2]){
  #  print(ab[i,-1])
  # print(log(ab[i,-1]+1))
  #ab[,t] <- log2(ab[,t]+0.01)
  #ab[,t] <- rescale(as.matrix(ab[,t]),c(0,1))
  #}
  
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
    
  )
  print(list.files()[i])
  ggsave(file = sprintf('SPider_plots/%s_.png',list.files()[i]), device = 'png',width = 10, height = 10)
}