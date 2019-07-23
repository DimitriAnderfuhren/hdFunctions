
library(plotly)
library(microbenchmark)
# quickDensityPlot
#install.packages("microbenchmark")

nObs = 32
n = 1000000
obs = 1:nObs
y = matrix(nrow = n, ncol = nObs)
for(i in obs){
  y[,i] = rnorm(n = n, mean = runif(1,-1,1), sd = 1)
}



microbenchmark::microbenchmark(
  ds = apply(y,2,density)
)

ds = apply(y,2,density)

plot(ds[[1]])
for(i in 2:length(ds)){
  lines(ds[[i]])
}

ys = sapply(ds, function(x){x[["y"]]})

p = plot_ly(x = ds[[1]]$x,y = ds[[10]]$y, type = "scatter", mode = "line")
for(i in 2:length(ds)){
  p = add_trace(p,x = ds[[i]]$x, y = ds[[i]]$y)
}
p

