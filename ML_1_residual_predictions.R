#
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
# Author: Guido Cervone (cervone@psu.edu) and Fangcao Xu (xfangcao@psu.edu)
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

library(gtools)
library(RColorBrewer)
library(plotly)
source("Z_global_variables.R")
source("Z_MODTRAN_functions.R")
source("Z_plot_functions.R")
source("Z_general_functions.R")


### This script is for the predicted matrix csv files with the local time

# read predicted files
# 15330
predictfiles.trans <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/TransmissionCSV",full.names=TRUE))
predictfiles.down <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/DownwellingCSV",full.names=TRUE))
predictfiles.grnd <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/GroundReflectedCSV",full.names=TRUE))
predictfiles.up1 <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/PathThermalEmissionCSV",full.names=TRUE))
predictfiles.up2 <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/PathThermalScatteringCSV",full.names=TRUE))
predictfiles.surf <- mixedsort(list.files("~/disk10TB/DARPA/predictedCSV/MidLatitude/SurfaceEmissionCSV",full.names=TRUE))

# plot predicted and original matrixplot
originalfilesfolder <- "~/disk10TB/DARPA/MatrixCSV/MidLatitude"
originalfiles.trans <- find.originalfiles(predictfiles.trans,originalfilesfolder)
originalfiles.down <- find.originalfiles(predictfiles.down,originalfilesfolder)
originalfiles.grnd <- find.originalfiles(predictfiles.grnd,originalfilesfolder)
originalfiles.up1 <- find.originalfiles(predictfiles.up1,originalfilesfolder)
originalfiles.up2<- find.originalfiles(predictfiles.up2,originalfilesfolder)
originalfiles.surf<- find.originalfiles(predictfiles.surf,originalfilesfolder)
# find the original total radiance files, all componemts correspond to the same associated total input files 
totalInputfiles <- find.originaltotalInput(predictfiles.down, originalfilesfolder)

# Create a list where the first component are the predicted tables, and the second component are the original tables
predorigfiles.trans.read <- read.predorigmatrix(predictfiles.trans, originalfiles.trans)
predorigfiles.grnd.read <- read.predorigmatrix(predictfiles.grnd, originalfiles.grnd)
predorigfiles.surf.read <- read.predorigmatrix(predictfiles.surf, originalfiles.surf)
predorigfiles.down.read <- read.predorigmatrix(predictfiles.down, originalfiles.down)
predorigfiles.up1.read <- read.predorigmatrix(predictfiles.up1, originalfiles.up1)
predorigfiles.up2.read <- read.predorigmatrix(predictfiles.up2, originalfiles.up2)
# read the original total radiance files
totalInputfiles.read <- read.originalmatrix(totalInputfiles)

# predicted files stored in disk10TB are not complete and same. #15530 total files 
vars.trans <-return.index(predictfiles.trans, 2:4, names = c("Day","Time","Reflectivity"))
vars.grnd <-return.index(predorigfiles.grnd, 2:4, names = c("Day","Time","Reflectivity"))
vars.surf<-return.index(predorigfiles.surf, 2:4, names = c("Day","Time","Reflectivity"))
vars.down<-return.index(predorigfiles.down, 2:4, names = c("Day","Time","Reflectivity"))
vars.up1 <-return.index(predorigfiles.up1, 2:4, names = c("Day","Time","Reflectivity"))
vars.up2 <-return.index(predorigfiles.up2, 2:4, names = c("Day","Time","Reflectivity"))
vars.total <- return.index(totalInputfiles, 2:4, names = c("Day","Time","Reflectivity"))

# # Plot everything together
# comparisonMatrixPlot(predorigfiles.down.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/GroundReflected/")
# comparisonMatrixPlot(predorigfiles.up1.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/PathThermalEmission/")
# comparisonMatrixPlot(predorigfiles.up2.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/PathThermalScattering/")
# comparisonMatrixPlot(predorigfiles.surf.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/SurfaceEmission/")

# Let's plot radiance as a function of wavelength, angles, and Reflectivity
day <- 107
time <- 14  # here the time is the local time

# prepare the sliced input data based on the assigned variables
valid.index <- which(vars.trans$Day == day & vars.trans$Time == time)
#valid.index <- which(vars$Day == day & vars$Time == time)
input.trans <- list(predorigfiles.trans.read[[1]][valid.index],predorigfiles.trans.read[[2]][valid.index])
valid.index <- which(vars.grnd$Day == day & vars.grnd$Time == time)
input.grnd <- list(predorigfiles.grnd.read[[1]][valid.index],predorigfiles.grnd.read[[2]][valid.index])
valid.index <- which(vars.surf$Day == day & vars.surf$Time == time)
input.surf <- list(predorigfiles.surf.read[[1]][valid.index],predorigfiles.surf.read[[2]][valid.index])
valid.index <- which(vars.down$Day == day & vars.down$Time == time)
input.down <- list(predorigfiles.down.read[[1]][valid.index],predorigfiles.down.read[[2]][valid.index])
valid.index <- which(vars.up1$Day == day & vars.up1$Time == time)
input.up1 <- list(predorigfiles.up1.read[[1]][valid.index],predorigfiles.up1.read[[2]][valid.index])
valid.index <- which(vars.up2$Day == day & vars.up2$Time == time)
input.up2 <- list(predorigfiles.up2.read[[1]][valid.index],predorigfiles.up2.read[[2]][valid.index])
valid.index <- which(vars.total$Day == day & vars.total$Time == time)
input.total <-  totalInputfiles.read[valid.index]

# get the minimum and maximum residual over four components 
input.trans.diff <- mapply(function(x,y) diff.matrices(x, y), input.trans[[1]], input.trans[[2]], SIMPLIFY = FALSE)
input.grnd.diff <- mapply(function(x,y) diff.matrices(x, y), input.down[[1]], input.down[[2]], SIMPLIFY = FALSE)
input.surf.diff <- mapply(function(x,y) diff.matrices(x, y), input.surf[[1]], input.surf[[2]], SIMPLIFY = FALSE)
input.down.diff <- mapply(function(x,y) diff.matrices(x, y), input.down[[1]], input.down[[2]], SIMPLIFY = FALSE)
input.up1.diff <- mapply(function(x,y) diff.matrices(x, y), input.up1[[1]], input.up1[[2]], SIMPLIFY = FALSE)
input.up2.diff <- mapply(function(x,y) diff.matrices(x, y), input.up2[[1]], input.up2[[2]], SIMPLIFY = FALSE)
input.diff.minmax <- minmax.crosstables(c(unlist(input.down.diff, recursive = FALSE), unlist(input.up1.diff, recursive = FALSE),unlist(input.up2.diff, recursive = FALSE)))

# transmission
comparisonMatrixPlot(input.trans, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/Transmission/",residual.opt = "T",residual.square.opt = F, zname="Transmission", zunit=NULL)
# ground reflected component
comparisonMatrixPlot(input.grnd, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/GroundReflected/",residual.opt = "T",residual.square.opt = F)
# surface component
comparisonMatrixPlot(input.surf, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/SurfaceEmission/",residual.opt = "T", residual.square.opt = F)
# downwelling component
comparisonMatrixPlot(input.down, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/Downwelling/",residual.opt = "T",residual.square.opt = F, zlim.residual = input.diff.minmax)
# path thermal emission component
comparisonMatrixPlot(input.up1, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/PathThermalEmission/", residual.opt = "T", residual.square.opt = F, zlim.residual = input.diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(input.up2, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/PathThermalScattering/",residual.opt = "T", residual.square.opt = F,  zlim.residual = input.diff.minmax)
# total radiance component
mapply(function(x,y) MODTRAN.matplot(as.matrix(x), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/June30th/TotalRadiance/", y, ".jpg",sep = "")), input.total, names(input.total))


### Check one variable by fixing all other variables for the radiance residual
### average plots for down componet
marray.original <- statistic.prep(input.down[[1]])
marray.predict <- statistic.prep(input.down[[2]])
marray.residual<- statistic.prep(input.down, check.residual=T)
# average radiance of different wavelength
output.avg <- "/amethyst/s0/fbx5002/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/averageMatrixPlots/"
matplot.average(marray.predict, avg = "wavelength",nlevels = 10,plotContour = T, plotKey = T, filename=paste(output.avg,"predictAvgWavelength_down.jpg",sep = ""))
matplot.average(marray.original, avg = "wavelength",nlevels = 10,plotContour = T, plotKey = T,filename=paste(output.avg,"originalAvgWavelength_down.jpg",sep = ""))
matplot.average(marray.residual, avg = "wavelength",nlevels = 10,plotContour = T, plotKey = T, filename=paste(output.avg,"residualAvgWavelength_down.jpg",sep = ""))
#
matplot.average(marray.predict, avg = "angle",nlevels = 10,plotContour = T, plotKey = T, filename=paste(output.avg,"originalAvgAngle_down.jpg",sep = ""))
matplot.average(marray.original, avg = "angle",nlevels = 10,plotContour = T, plotKey = T,filename=paste(output.avg,"predictAvgAngle_down.jpg",sep = ""))
matplot.average(marray.residual, avg = "angle",nlevels = 10,plotContour =F, plotKey = T,filename=paste(output.avg,"residualAvgAngle_down.jpg",sep = ""))
#
matplot.average(marray.predict, avg = "reflectivity",nlevels = 10,plotContour = F, plotKey = T,filename=paste(output.avg,"originalAvgReflectivity_down.jpg",sep = ""))
matplot.average(marray.original, avg = "reflectivity",nlevels = 10,plotContour = T, plotKey = T,filename=paste(output.avg,"predictAvgReflectivity_down.jpg",sep = ""))
matplot.average(marray.residual, avg = "reflectivity",nlevels = 10,plotContour = F, plotKey = T,filename=paste(output.avg,"residualAvgReflectivity_down.jpg",sep = ""))


### Fix other variables to check the distribution
# check the surf component
day.check.surf <- residual.prep(predorigfiles.surf.read, var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60))
time.check.surf <- residual.prep(predorigfiles.surf.read, var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60))
reflectivity.check.surf <- residual.prep(predorigfiles.surf.read, var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60))
angle.check.surf <- residual.prep(predorigfiles.surf.read, var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10))
matplot.residual.check(as.matrix(day.check.surf), var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60), 
                 filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/daycheck_14_10_60_surf.jpg"), zlim=c(-5.5e-5, 5.5e-5),nlevels=10)
matplot.residual.check(as.matrix(time.check.surf), var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/timecheck_107_10_60_surf.jpg", nlevels=10)
matplot.residual.check(as.matrix(reflectivity.check.surf), var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/reflectivitycheck_107_14_60_surf.jpg", nlevels=10)
matplot.residual.check(as.matrix(angle.check.surf), var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10), 
                 filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/anglecheck_107_14_10_surf.jpg", nlevels=10)
# check the down component
day.check <- residual.prep(predorigfiles.down.read, var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60))
time.check <- residual.prep(predorigfiles.down.read, var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60))
reflectivity.check <- residual.prep(predorigfiles.down.read, var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60))
angle.check <- residual.prep(predorigfiles.down.read, var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10))
matplot.residual.check(as.matrix(day.check), compDown=TRUE, var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60), 
                 filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/daycheck_14_10_60_down.jpg"), nlevels=10)
matplot.residual.check(as.matrix(time.check), compDown=TRUE, var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/timecheck_107_10_60_down.jpg", nlevels=10)
matplot.residual.check(as.matrix(reflectivity.check), compDown=TRUE, var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/reflectivitycheck_107_14_60_down.jpg", nlevels=10)
matplot.residual.check(as.matrix(angle.check), compDown=TRUE, var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10), 
                 filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/anglecheck_107_14_10_down.jpg", nlevels=10)
# check the up1 component
day.check <- residual.prep(predorigfiles.up1.read, var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60))
time.check <- residual.prep(predorigfiles.up1.read, var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60))
reflectivity.check <- residual.prep(predorigfiles.up1.read, var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60))
angle.check <- residual.prep(predorigfiles.up1.read, var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10))
matplot.residual.check(as.matrix(day.check), var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60), 
                 filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/daycheck_14_10_60_up1.jpg"), zlim=c(-3.1e-6, 3.1e-6), nlevels=10)
matplot.residual.check(as.matrix(time.check), var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/timecheck_107_10_60_up1.jpg", nlevels=10)
matplot.residual.check(as.matrix(reflectivity.check), var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/reflectivitycheck_107_14_60_up1.jpg", nlevels=10)
matplot.residual.check(as.matrix(angle.check), var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10), 
                 filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/anglecheck_107_14_10_up1.jpg", nlevels=10)
# check the up2 component 
day.check <- residual.prep(predorigfiles.up2.read, var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60))
time.check <- residual.prep(predorigfiles.up2.read, var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60))
reflectivity.check <- residual.prep(predorigfiles.up2.read, var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60))
angle.check <- residual.prep(predorigfiles.up2.read, var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10))
matplot.residual.check(as.matrix(day.check), var="day", fix.var= c("time"=14,"reflectivity"=10,"angle"=60), 
                 filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/daycheck_14_10_60_up2.jpg"), zlim=c(-6.8e-8, 6.8e-8), nlevels=10)
matplot.residual.check(as.matrix(time.check), var="time", fix.var= c("day"=107,"reflectivity"=10,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/timecheck_107_10_60_up2.jpg", nlevels=10)
matplot.residual.check(as.matrix(reflectivity.check), var="reflectivity", fix.var= c("day"=107,"time"=14,"angle"=60), 
                 filename= "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/reflectivitycheck_107_14_60_up2.jpg", nlevels=10)
matplot.residual.check(as.matrix(angle.check), var="angle", fix.var= c("day"=107,"time"=14,"reflectivity"=10), 
                 filename="~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/residuals/anglecheck_107_14_10_up2.jpg", nlevels=10)
# # batch plot all residuals
batch.matplot.residual.check(predorigfiles.down.read, comp ="down", check="day")
# batch.matplot.residual.check(predorigfiles.down.read, comp ="down", check="angle")
# batch.matplot.residual.check(predorigfiles.down.read, comp ="down", check="reflectivity")
# batch.matplot.residual.check(predorigfiles.down.read, comp ="down", check="time")
# batch.matplot.residual.check(predorigfiles.up1.read, comp ="up1", check="day")
# batch.matplot.residual.check(predorigfiles.up1.read, comp ="up1", check="angle")
# batch.matplot.residual.check(predorigfiles.up1.read, comp ="up1", check="time")
# batch.matplot.residual.check(predorigfiles.up1.read, comp ="up1", check="reflectivity")
# batch.matplot.residual.check(predorigfiles.up2.read, comp ="up2", check="day")
# batch.matplot.residual.check(predorigfiles.up2.read, comp ="up2", check="angle")
# batch.matplot.residual.check(predorigfiles.up2.read, comp ="up2", check="time")
# batch.matplot.residual.check(predorigfiles.up2.read, comp ="up2", check="reflectivity")
# batch.matplot.residual.check(predorigfiles.surf.read, comp ="surf", check="day")
# batch.matplot.residual.check(predorigfiles.surf.read, comp ="surf", check="angle")
# batch.matplot.residual.check(predorigfiles.surf.read, comp ="surf", check="time")
# batch.matplot.residual.check(predorigfiles.surf.read, comp ="surf", check="reflectivity")


