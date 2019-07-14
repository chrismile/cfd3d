library("ggplot2")
library("reshape2")
require(ggplot2)
require(ggthemes)

setwd("/home/stefan/Uni/CFD_Lab/cfd3d/tests/")

#function for RuntimeGraph
printGraph <- function(fileName, subtitle){
  data <- read.csv(fileName, check.names=FALSE, header=TRUE, sep=',')
  data <- melt(data, id="solver")
  
  p <- ggplot(data=data, mapping=aes(x=reorder(solver, -solver), y=value, fill=variable))
  p <- p + scale_fill_manual(values=c("#333333","#E37222","#0065BD", "#A2AD00"))
  p <- p + coord_flip()
  p <- p + guides(fill = guide_legend(reverse=TRUE))
  p <- p + labs(title = subtitle)
  p <- p + labs(x = "Solver", y = "Time (s)", fill = "Order")
  p <- p + geom_bar(stat="identity", position = "dodge")
  p <- p + theme_tufte()
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + theme(plot.subtitle = element_text(hjust = 0.5))
  p
}

p <- printGraph("results/rayleigh_benard_convection_8-2-1_runtime.csv", "Runtime Rayleigh-Benard-Convection")
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1.svg', p, width=7, height=4, units='in')
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1.pdf', p, width=7, height=4, units='in')