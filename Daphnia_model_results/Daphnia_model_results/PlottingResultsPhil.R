g1 <- ggradar(dts[-c(9,10),],values.radar=c(0,0.5,1),group.line.width=2,group.point.size=4,legend.position="right")+
	theme(
		plot.margin=margin(0.5,0,0.5,0,unit="cm"),
		legend.margin=margin(0,0,0,0,unit="cm"),
		legend.box.margin=margin(-5,0,-5,-1,unit="cm")
	)+
	coord_cartesian(clip = "off")



+guides(fill=guide_legend(nrow=4,byrow=TRUE))

dts[,1] <- c(	"Increased cost\nfor reproduction",
			"Decreased energy\nconductance",
			"Decreased\nfeeding",
			"Decreased\nfraction to soma",
			"Hazard during\noogenesis",
			"Increased\nallocation to soma",
			"Increased cost\nfor growth\nand reproduction",
			"Increased\nmaintenance cost",
			"Max conductance",
			"Min conductance")

colnames(dts) <- c(	"group",
				"stress_ls",
				"the structural\nlength, cm",
				"energy\nreserve, J",
				"reproduction\nbuffer, J",
				"internal\nconcentration, nM",
				"food\ndensity",
				"maximum\nstructural\nlength",
				"cumulative\nembryos, #",
				"mean")
