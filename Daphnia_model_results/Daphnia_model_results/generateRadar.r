rannot <- c(
	"Cost_Reproduction"="Increased cost for reproduction",
	"Decreased_energy"="Decreased energy conductance",
	"Decreased_feeding"="Decreased feeding",
	"Decreased_fraction"="Decreased fraction to soma",
	"Hazard_oogenesis"="Hazard during oogenesis",
	"Increased_Allocation"="Increased allocation to soma",
	"Increased_cost"="Increased cost for growth\nand reproduction",
	"Increased_Maintenance"="Increased maintenance cost"
)

cannot <- c(
	"group",
	"stress_ls"="stress_ls",
	"the structural length, cm"="the structural length, cm",
	"energy reserve, J"="energy reserve, J",
	"energy invested in maturity, J"="energy invested in maturity, J",
	"reproduction buffer, J"="reproduction buffer, J",
	"internal concentration, nM"="internal concentration, nM",
	"survival probability, %"="survival probability, %",
	"food density"="food density",
	"maximum structural length"="maximum structural length",
	"cumulative embryos, #"="cumulative embryos, #",
	"mean"="mean"
)

generateRadar <- function(var,var.sel=seq(3,22,by=2),rem.zero=T,row.annot=rannot,col.annot=cannot,wrap.col=10,wrap.row=wrap.col,legend.position="right",group.line.width=2,group.point.size=4,values.radar=c(0,0.5,1),...){
	require(ggradar)
	require(scales)
	require(stringr)
	require(ggplot2)
	vard <- var[c(var.sel,nrow(var)),-1]
	rownames(vard) <- var[c(var.sel,nrow(var)),1]
	if(rem.zero)
		vard <- vard[which(vard[,1] != 0),]
	vard <- t(vard)
	vard.df <- data.frame(vard)
	if(!is.null(col.annot))
		colnames(vard.df) <- str_wrap(col.annot[colnames(vard)],wrap.col)
	if(!is.null(row.annot)){
		vard.df$group <- str_wrap(row.annot[rownames(vard)],wrap.row)
	}else{
		vard.df$group <- rownames(vard.df)
	}
	vard.df <- vard.df[,c("group",colnames(vard.df)[-length(colnames(vard.df))])]
	gpl <- ggradar(vard.df,values.radar=values.radar,group.line.width=group.line.width,group.point.size=group.point.size,legend.position=legend.position,...)+
	theme(
		plot.margin=margin(0.5,0,0.5,0.5,unit="cm"),
		legend.margin=margin(0,0,0,0,unit="cm"),
		legend.box.margin=margin(-5,0.5,-5,0.5,unit="cm")
	)+
	coord_cartesian(clip = "off")
	return(gpl)
}