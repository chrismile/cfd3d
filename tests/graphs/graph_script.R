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
  p <- p + labs(subtitle = "(80x20x10)")
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


    #Speedup Graph
    data <- read.csv("results/speedup_OpenMP_double.csv", check.names=FALSE, header=TRUE, sep=',')
    
    for(i in 1:nrow(data)) 
    {
      data[i, "speedup"] <- data[data$threads == 1 & data$group == data[i,"group"], "time"]/data[i, "time"]
    }
    levels(data$group)
    #data$group <- factor(data$group, levels = rev(levels(data$group)))
    data
    
    p <- ggplot(data=data, mapping=aes(x=threads, y=speedup, linetype=group))
    p <- p + scale_colour_manual(values=c("#333333","#E37222","#0065BD", "#A2AD00"))
    p <- p + labs(title = "Speedup OpenMP")
    p <- p + labs(subtitle = "Rayleigh Benard Convection (80x20x10)")
    p <- p + labs(x = "Threads", y = "Speedup")
    p <- p + theme_linedraw()
    p <- p + theme(legend.title = element_blank())
    p <- p + guides(linetype = guide_legend(reverse=TRUE))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_continuous(breaks=c(1,2,4,8,16,28))
    p <- p + scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + theme(plot.subtitle = element_text(hjust = 0.5))
    p
    
    ggsave('graphs/speedup_OpenMP_double.svg', p, width=7, height=4, units='in')
    ggsave('graphs/speedup_OpenMP_double.pdf', p, width=7, height=4, units='in')
    
    #Speedup Graph MPI
    data <- read.csv("results/speedup_MPI.csv", check.names=FALSE, header=TRUE, sep=',')
    
    for(i in 1:nrow(data)) 
    {
      data[i, "speedup"] <- data[data$threads == 1, "time"]/data[i, "time"]
    }
    levels(data$group)
    data
    
    p <- ggplot(data=data, mapping=aes(x=threads, y=speedup,  linetype=group))
    p <- p + scale_colour_manual(values=c("#333333","#E37222","#0065BD", "#A2AD00"))
    p <- p + labs(title = "Speedup MPI")
    p <- p + labs(subtitle = "Rayleigh Benard Convection (80x20x10)")
    p <- p + labs(x = "Processes", y = "Speedup")
    p <- p + theme_linedraw()
    p <- p + theme(legend.title = element_blank())
    p <- p + guides(linetype = guide_legend(reverse=TRUE))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + scale_x_continuous(breaks=c(1,27,56))
    p <- p + scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20))
    p <- p + theme_tufte()
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + theme(plot.subtitle = element_text(hjust = 0.5))
    p
    
    ggsave('graphs/speedup_MPI.svg', p, width=7, height=4, units='in')
    ggsave('graphs/speedup_MPI.pdf', p, width=7, height=4, units='in')
