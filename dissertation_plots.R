overlaps <- read.csv("Branco-work/overlaps.csv")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Checking for overlap with full length L1Hs sequences. 
blen100 <- overlaps[overlaps$blen == 1,]
blen95 <- overlaps[overlaps$blen == 0.95,]
blen0 <- overlaps[overlaps$blen == 0,]

newcol = rainbow(4)[(blen100$macs_q)]
plot(blen100$min_macs_peak_intensity, blen100$full.length.percent, col=newcol)

library(ggplot2)
plot_blen100 <- ggplot(blen100, aes(min_macs_peak_intensity, full.length.percent, colour=factor(macs_q) )) + 
  geom_point() + 
  labs(y="% of Peaks at 6KB LINEs", x="Min Peak Intensity") + 
  ggtitle("Filtered by 100% Matches") + 
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.045)) + 
  theme_minimal()

plot_blen95 <- ggplot(blen95, aes(min_macs_peak_intensity, full.length.percent, colour=factor(macs_q) )) +
  geom_point() + labs(y="% of Peaks at 6KB LINEs", x="Min Peak Intensity") +
  ggtitle("Filtered by 95% Matches") + 
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.045))+ 
  theme_minimal()

plot_blen0 <-ggplot(blen0, aes(min_macs_peak_intensity, full.length.percent, colour=factor(macs_q) )) + geom_point() + 
  labs(y="% of Peaks at 6KB LINEs", x="Min Peak Intensity") + 
  ggtitle("Minimap2 short-read preset alignment") + 
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.045))+ 
  theme_minimal()

multiplot(plot_blen0, plot_blen95, plot_blen100, cols=3)
#plot(overlaps95nona$min_macs_peak_intensity, overlaps95nona$overlapping.full.length.L1HS)
#plot(overlaps$min_macs_peak_intensity, overlaps$full.length.percent)
install.packages("tidyverse")
library(forcats)
library(dplyr)

#Genome Composition Plot 
genome_comp_NA <- read.csv("Branco-work/genome_composition.csv")
genome_comp_all <- na.omit(genome_comp_NA)
genome_comp <- genome_comp_all[genome_comp_all$Type != "Unmasked",]

genome_comp %>%
  mutate(name = fct_reorder(genome_comp$Type, genome_comp$Bp)) %>%
  ggplot( aes(x=name, y=genome_comp$Bp)) +
  geom_bar(stat="identity", alpha=.6, width=.4, fill="green4") +
  ggtitle("Annotated Genome Composition") +
  coord_flip() +
  xlab("Types of Sequence") +
  ylab("Genome coverage in base pairs") +
  theme_minimal()


shuffled_best <- read.csv("Branco-work/random_shuffle.csv")

  
shuffled3 <- ggplot(shuffled_best, aes(x= reorder(file, X_L1), y=X_L1, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random Specificity LINE-1") +
  labs(y="percent of peaks overlapping LINE-1", x="Model Input Type")+
  geom_text(aes(label=X_L1), color="black", size = 3, vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

  
shuffled1 <- ggplot(shuffled_best, aes(x= reorder(file, X_L1hs), y=X_L1hs, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random Specificity L1hs") +
  labs(y="percent of peaks overlapping L1Hs", x="Model Input Type")+
  geom_text(aes(label=X_L1hs), color="black", size = 3 ,vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

  
shuffled2 <- ggplot(shuffled_best,aes(x= reorder(file, X_FL_L1Hs), y=X_FL_L1Hs, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random Specificity Full Length L1hs") +
  labs(y="percent of peaks overlapping Full Length L1Hs", x="Model Input Type")+
  geom_text(aes(label=X_FL_L1Hs), color="black", size = 3, vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()


sensitivity1 <- ggplot(shuffled_best, aes(x= reorder(file, X_L1), y=L1_sens, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random LINE-1 Sensitivity") +
  labs(y="percent of LINE-1 found", x="Model Input Type")+
  geom_text(aes(label=L1_sens), color="black", size = 3, vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()


sensitivity2 <- ggplot(shuffled_best, aes(x= reorder(file, X_L1hs), y=L1Hs_sens, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random L1hs Sensitivity") +
  labs(y="percent of L1Hs found", x="Model Input Type")+
  geom_text(aes(label=L1Hs_sens), color="black", size = 3 ,vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()


sensitivity3 <- ggplot(shuffled_best,aes(x= reorder(file, X_FL_L1Hs), y=FL_L1Hs_sens, fill=status)) +
  geom_bar(stat = "identity", width= 0.8, position = position_dodge()) + 
  ggtitle("Model vs Random Full Length L1hs Sensitivity") +
  labs(y="percent of Full Length L1Hs found", x="Model Input Type")+
  geom_text(aes(label=FL_L1Hs_sens), color="black", size = 3, vjust=0, position = position_dodge(0.8)) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

multiplot(shuffled3, sensitivity1, shuffled1, sensitivity2, shuffled2, sensitivity3,cols=3)

multi_seq_scores <- shuffled_best[shuffled_best$file == "All L1Hs",]


allseq1 <- ggplot(multi_seq_scores, aes(x= "All LINE-1", y=L1_sens, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle("Multiple Sequence Model Performance") +
  labs(y="Sensitivity", x=" ") +
#  geom_text(aes(label="random"), color="black", size = 3 ,vjust=0.5) + 
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()
allseq1

allseq2 <- ggplot(multi_seq_scores, aes(x= "L1Hs", y=L1Hs_sens, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle(" ") +
  labs(y=" ", x=" ") +
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

allseq3 <- ggplot(multi_seq_scores, aes(x= "Full length L1Hs", y=FL_L1Hs_sens, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle(" ") +
  labs(y=" ", x=" ") +
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

allseq4 <- ggplot(multi_seq_scores, aes(x= "All LINE-1", y=X_L1, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle("Multiple Sequence Model Performance") +
  labs(y="Specificity", x=" ") +
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

allseq5 <- ggplot(multi_seq_scores, aes(x= "L1Hs", y=X_L1hs, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle(" ") +
  labs(y=" ", x=" ") +
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

allseq6 <- ggplot(multi_seq_scores, aes(x= "Full length L1Hs", y=X_FL_L1Hs, fill=staus)) +
  geom_bar(stat = "identity", width= 0.2, fill = "palegreen") + 
  ggtitle(" ") +
  labs(y=" ", x=" ") +
  scale_fill_brewer(palette="Greens")+
  coord_cartesian(ylim = c(0, 100))+
  theme_minimal()

multiplot(allseq1, allseq4, allseq2,allseq5, allseq3, allseq6, cols = 3)
