library("ggplot2")
library("reshape2")
require(ggplot2)
require(ggthemes)

setwd("/home/stefan/Uni/CFD_Lab/cfd3d/tests/")

#function for RuntimeGraph
printGraph <- function(fileName, subtitle, x_axis_label, fill_legend){
  data <- read.csv(fileName, check.names=FALSE, header=TRUE, sep=',')
  levels(data$group)
  data
  p <- ggplot(data=data, mapping=aes(x=reorder(platform, -time), y=time, fill=as.factor(precision)))
  p <- p + scale_fill_manual(values=c("#E37222","#0065BD", "#A2AD00", "#333333"))
  p <- p + guides(fill = guide_legend(reverse=TRUE))
  p <- p + coord_flip()
  p <- p + labs(title = "Runtime Rayleigh Benard Convection")
  p <- p + labs(x = x_axis_label, y = "time (s)", fill = fill_legend)
  p <- p + geom_bar(stat="identity", position="dodge")
  p <- p + theme_tufte()
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + theme(plot.subtitle = element_text(hjust = 0.5))
  p
}

p <- printGraph("results/rayleigh_benard_convection_8-2-1_runtime_solver.csv", "Solvers", "Parallelisation", "Solver")
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1_solver.svg', p, width=7, height=4, units='in')
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1_solver.pdf', p, width=7, height=4, units='in')

p <- printGraph("results/rayleigh_benard_convection_8-2-1_runtime_platforms.csv", "Platforms", "Platform", "Precision")
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1_platforms.svg', p, width=7, height=4, units='in')
ggsave('graphs/runtime_rayleigh_benard_convection_8-2-1_platforms.pdf', p, width=7, height=4, units='in')