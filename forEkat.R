# For Ekat
install.packages("grImport")
library(grDevices)
library(grImport)
library(igraph)
library(tcltk)
devtools::install_github("ParkerICI/panorama")

library(panorama)
panorama()

path = "/Volumes/Inflammation_research/People/Dimitri/clustered_single_samples_100/unsupervised_graph/100_unsupervised.graphml"

g = read_graph(path, format = "graphml")

tkplot(g, canvas.width = 450, canvas.height = 450)

PostScriptTrace("Rplot.eps", "new")
r = readPicture("new")

plot.igraph(g)
