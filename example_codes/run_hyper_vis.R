library(PPCI)
library(hyperclusvis)


#hyperbolice OV example

## load Fisher's Iris data set
data(iris)
d <- iris[, 1:4]
dat1 <- as.matrix(d)

set.seed(1000000)
hDat <- hyperbolize(dat1, 3, preprocess="kmeans", scale="std")
plot_disk(hDat,3)


## load wine
library(pgmm)
data("wine")
d <- wine[, -1]
dat1 <- as.matrix(d)
seed <- 1000000
set.seed(seed)
## apply the hyperbolic mapper and extract the returned information
hDat <- hyperbolize(dat1, 3, preprocess="pgmm", scale="std",use_input_order =F)
plot_disk(hDat,3)


###EV example#####
library(HDclassif)
#yale
Dp <- PrepareData(yale$x)
set.seed(42)
hddc_r<- hddc(Dp,10)
km_r <-kmeans(Dp,10)
yale_hddc=hddc_r$class
yale_km=km_r$cluster
cluster_performance(yale_hddc,yale$c)
cluster_performance(yale_km,yale$c)



set.seed(42)
origin_yale=hyperbolize(yale$x,preprocess =yale$c,scale='norm')
origin_yale$points=origin_yale$points[-1,]
plot_disk(origin_yale,10)

set.seed(42)
out_yale_true2 <- hyperbolizeEV(yale$x, numk=10,
                                           preprocess=yale$c,
                                           scale="norm",
                                           local_model="kmeansWrap",
                                           kg_per_class=10,
                                           spread_scale =0.6)
plot_disk(out_yale_true2,10)


set.seed(42)
out_yale_hddc2 <- hyperbolizeEV(yale$x, numk=10,
                                           preprocess=yale_hddc,
                                           scale="norm",
                                           local_model="kmeansWrap",
                                           kg_per_class=10,
                                           spread_scale =0.6)
plot_disk(out_yale_hddc2,10)

set.seed(42)
out_yale_km2<- hyperbolizeEV(yale$x, numk=10,
                                        preprocess=yale_km,
                                        scale="norm",
                                        local_model="kmeansWrap",
                                        kg_per_class=10,
                                        spread_scale =0.6)
plot_disk(out_yale_km2,10)



#digits
d_idx=which(colSums(optidigits$x)==0)
digit_data=optidigits$x[,-d_idx]
Dp <- PrepareData(digit_data)

set.seed(42)
hddc_r<- hddc(Dp,10)
km_r <-kmeans(Dp,10)
digit_hddc=hddc_r$class
digit_km=km_r$cluster
cluster_performance(digit_hddc,optidigits$c)
cluster_performance(digit_km,optidigits$c)

set.seed(42)
origin_digit=hyperbolize(digit_data,preprocess =optidigits$c+1,scale='norm')
origin_digit$points=origin_digit$points[-1,]
plot_disk(origin_digit,10)

set.seed(42)
out_digit_true2 <- hyperbolizeEV(digit_data, numk=10,
                                            preprocess=optidigits$c+1,
                                            scale="norm",
                                            local_model="kmeansWrap",
                                            kg_per_class=10,
                                            spread_scale =0.5)
plot_disk(out_digit_true2,10)



set.seed(42)
out_digit_hddc <- hyperbolizeEV(digit_data, numk=10,
                                           preprocess=digit_hddc,
                                           scale="norm",
                                           local_model="kmeansWrap",
                                           kg_per_class=10,
                                           spread_scale =0.5)
plot_disk(out_digit_hddc,10)


set.seed(42)
out_digit_km <- hyperbolizeEV(digit_data, numk=10,
                                         preprocess=digit_km,
                                         scale="norm",
                                         local_model="kmeansWrap",
                                         kg_per_class=10,
                                         spread_scale = 0.5)
plot_disk(out_digit_km,10)


#simulated

library(clusterGeneration)
## number of clusters
nK <- 7

##
seed <- (1000000 + 7)
set.seed(seed)

## generate random data sets
datasets <- clusterGeneration::genRandomClust(
  numClust = nK,
  sepVal = 0.1,
  numNonNoisy = 100,
  numReplicate = 3,
  fileName = file.path(tempdir(), "temp")
)

data  <- datasets$datList[[1]]
label <- datasets$memList[[1]]
## apply the hyperbolic mapper and extract the returned information
seed <- (1000000 + 7)
set.seed(seed)
hDat <- hyperbolize(data, nK, preprocess="hddc", scale="norm")
hDat2 <- hyperbolizeEV(data, nK, preprocess="hddc", scale="norm",local_model="kmeansWrap",kg_per_class=10)
plot_disk(hDat2,7)
plot_disk(hDat,7)

