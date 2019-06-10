library(ggplot2)
library(dunn.test) #https://cran.r-project.org/web/packages/dunn.test/dunn.test.pdf

#Figure 2A
fig.2A.data<-read.table(file="Figure2A.distances.by.label.txt", header=TRUE) 
p<-ggplot(fig.2A.data, aes(x=fig.2A.data$Label, y=fig.2A.data$Distance))
p +geom_boxplot(fatten=1) + geom_jitter(shape=16, size=2, position=position_jitter(0.2))+ theme_bw()+ labs(x = NULL, y = expression(paste("Distance (", theta,"YC)")))+ scale_x_discrete(labels=c("Between_pairs" = "Between pairs", "Within_pairs" = "Within pairs", "Within_subject" = "Within subjects"))+ylim(0,1.2)
ggsave("fig.2A.eps", width = 5, height = 5)

kruskal.test(fig.2A.data$Distance ~ fig.2A.data$Label)
dunn.test(fig.2A.data$Distance, g=fig.2A.data$Label, method="bonferroni", kw=TRUE, label=TRUE, wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

#Figure 2B
del.mode<-read.table(file="pair.daughter_delivery_mode.txt", header=TRUE)
average.within.pair.distances<-read.table(file="average.within.pair.distances.plus.txt", header=TRUE)
del.mode.average.within.pair.distances<-merge(del.mode,average.within.pair.distances,by.x=c("Pair"), by.y=c("Pair"))
del.mode.average.within.pair.distances$Daughter_delivery_mode<-as.factor(del.mode.average.within.pair.distances$Daughter_delivery_mode)
p<-ggplot(del.mode.average.within.pair.distances, aes(x=reorder(del.mode.average.within.pair.distances$Daughter_delivery_mode_label, del.mode.average.within.pair.distances$Distance), y=del.mode.average.within.pair.distances$Distance))
p +geom_boxplot(fatten=1) + geom_jitter(shape=16, size=2, position=position_jitter(0.2))+ theme_bw()+ labs(x = NULL, y = expression(paste("Distance (", theta,"YC)")))+ylim(0,1.2)
ggsave("fig.2B.eps", width = 5, height = 5)

wilcox.test(del.mode.average.within.pair.distances$Distance ~ del.mode.average.within.pair.distances$Daughter_delivery_mode_label)
