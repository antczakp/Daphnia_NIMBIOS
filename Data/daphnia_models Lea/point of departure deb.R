##identify point of deprature


library(ggplot2)


p1 <- ggplot(data = parameter_mat[which(parameter_mat$days == 24),], aes(y=stress_ls,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p1)

p2 <- ggplot(data = parameter_mat[which(parameter_mat$days == 24),], aes(y=`the structural length, cm`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p2)

p3 <- ggplot(data = parameter_mat, aes(y=`energy reserve, J`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p3)

p4 <- ggplot(data = parameter_mat, aes(y=`energy invested in maturity, J`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p4)

p5 <- ggplot(data = parameter_mat[which(parameter_mat$days == 24),], aes(y=`reproduction buffer, J`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p5)


p7 <- ggplot(data = parameter_mat, aes(y=`internal concentration, nM`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p7)

p8<- ggplot(data = parameter_mat, aes(y=`survival probability, %`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p8)

p9 <- ggplot(data = parameter_mat, aes(y=`food density`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p9)

p10 <- ggplot(data = parameter_mat[which(parameter_mat$days == 24),], aes(y=`maximum structural length`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p10)

p11 <- ggplot(data = parameter_mat[which(parameter_mat$days == 24),], aes(y=`cumulative embryos, #`,x= conc))+
  geom_point(aes(col = as.factor(days)))
plot(p11)



log_scaled_conc <- rescale(log2(all_dec_energ_conduct$conc+0.001), c(0,1))
all_dec_energ_conduct <- cbind(all_dec_energ_conduct,log_scaled_conc)

p1 <- ggplot(data = all_dec_energ_conduct, aes(y=stress_ls,x= log_scaled_conc))+
  geom_point(aes(col = timepoint)) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position="none")+
  xlab(str_wrap('log scaled concentration nM [HG]',30 ))+
  ylab('stress')
plot(p1)

p2 <- ggplot(data = all_dec_energ_conduct, aes(y=`the structural length, cm`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position="none")+
  xlab(str_wrap('log scaled concentration nM [HG]',30 ))
plot(p2)

p3 <- ggplot(data = all_dec_energ_conduct, aes(y=`energy reserve, J`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p3)

p4 <- ggplot(data = all_dec_energ_conduct, aes(y=`energy invested in maturity, J`,x= log_scaled_conc))+
  geom_point(aes(col =timepoint))
plot(p4)

p5 <- ggplot(data = all_dec_energ_conduct, aes(y=`reproduction buffer, J`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p5)


p7 <- ggplot(data = all_dec_energ_conduct, aes(y=`internal concentration, nM`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p7)

p8<- ggplot(data = all_dec_energ_conduct, aes(y=`survival probability, %`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p8)

p9 <- ggplot(data = all_dec_energ_conduct, aes(y=`food density`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p9)

p10 <- ggplot(data = all_dec_energ_conduct, aes(y=`maximum structural length`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))
plot(p10)

p11 <- ggplot(data = all_dec_energ_conduct[which(all_dec_energ_conduct$timepoint <= 24),1:13], aes(y=`cumulative embryos, #`,x= log_scaled_conc))+
  geom_point(aes(col = timepoint))+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position="none")+
  xlab(str_wrap('log scaled concentration nM [HG]',30 ))
plot(p11)




p1 <- ggplot(data = all_dec_energ_conduct[], aes(y=stress_ls,x= conc))+
  geom_point(aes(col = timepoint))
plot(p1)

p2 <- ggplot(data = all_dec_energ_conduct, aes(y=`the structural length, cm`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p2)

p3 <- ggplot(data = all_dec_energ_conduct, aes(y=`energy reserve, J`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p3)

p4 <- ggplot(data = all_dec_energ_conduct, aes(y=`energy invested in maturity, J`,x= conc))+
  geom_point(aes(col =timepoint))
plot(p4)

p5 <- ggplot(data = all_dec_energ_conduct, aes(y=`reproduction buffer, J`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p5)


p7 <- ggplot(data = all_dec_energ_conduct, aes(y=`internal concentration, nM`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p7)

p8<- ggplot(data = all_dec_energ_conduct, aes(y=`survival probability, %`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p8)

p9 <- ggplot(data = all_dec_energ_conduct, aes(y=`food density`,x= conc))+
  geom_point(aes(col = timepoint))
plot(p9)

p10 <- ggplot(data = all_dec_energ_conduct, aes(y=`maximum structural length`,x= log(conc)))+
  geom_point(aes(col = timepoint))
plot(p10)

p11 <- ggplot(data = all_dec_energ_conduct[which(all_dec_energ_conduct$timepoint <= 24),], aes(y=`cumulative embryos, #`,x= log(conc)))+
  geom_point(aes(col = timepoint))
plot(p11)
