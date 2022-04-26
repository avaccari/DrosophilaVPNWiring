library(ggplot2)
source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function

allLC_plots = readRDS("output/LC_kmeans/allLC_plots.rds")

for (lc in names(allLC_plots)) {
  multiplot(allLC_plots[[lc]])
}

