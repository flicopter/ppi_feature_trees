library(ggplot2)
library(plyr)
library(reshape2)
library(scales)

setwd("~/git/ppi_feature_trees/Data/Results/ismb2015/")
data = read.csv("11receptors_rocs_bedrocs.csv")

auc_roc.mean  = data[, grep("^Receptor|^AUC.*mean$",  colnames(data))]
auc_roc.stdev = data[, grep("^Receptor|^AUC.*stdev$", colnames(data))]

bedroc.mean  = data[, grep("^Receptor|^BED.*mean$",  colnames(data))]
bedroc.stdev = data[, grep("^Receptor|^BED.*stdev$", colnames(data))]

methods = c("ExtraTrees", "SVM", "Z-score")
colnames(auc_roc.mean)[2:4]  = methods
colnames(auc_roc.stdev)[2:4] = methods
colnames(bedroc.mean)[2:4]   = methods
colnames(bedroc.stdev)[2:4]  = methods

auc_roc = merge(
  melt(auc_roc.mean,  id.vars="Receptor", variable.name="Method", value.name = "Mean"),
  melt(auc_roc.stdev, id.vars="Receptor", variable.name="Method", value.name = "StDev")
)

bedroc = merge(
  melt(bedroc.mean,  id.vars="Receptor", variable.name="Method", value.name = "Mean"),
  melt(bedroc.stdev, id.vars="Receptor", variable.name="Method", value.name = "StDev")
)

#### Plotting ####
p <- ggplot(auc_roc, aes(x=Receptor, y=Mean, fill=Method)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean-StDev, ymax=Mean+StDev, color=Method),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("AUC-ROC")
p
phoriz = p + coord_flip()

ggsave(filename="auc_roc.png", plot=p, width=11, height=5)
ggsave(filename="auc_roc.horiz.png", plot=phoriz, width=5, height=6)

## BEDROC
p <- ggplot(bedroc, aes(x=Receptor, y=Mean, fill=Method)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean*0.5, ymax=Mean+StDev, color=Method),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("BEDROC (alpha = 20)")
p
phoriz = p + coord_flip()

ggsave(filename="bedroc.png", plot=p, width=11, height=5)
ggsave(filename="bedroc.horiz.png", plot=phoriz, width=5, height=6)

## both at the same time
D = rbind(cbind(auc_roc, Measure="AUC-ROC"), cbind(bedroc, Measure="BEDROC"))
p.all <- ggplot(D, aes(x=Receptor, y=Mean, fill=Method)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean*0.5, ymax=Mean+StDev, color=Method),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ylab("Score") +
  facet_grid(.~Measure) +
  coord_flip()

p.all
ggsave(filename="aucroc_and_bedroc.png", plot=p.all, width=8, height=5)

