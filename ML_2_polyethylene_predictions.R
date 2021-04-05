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

############### Polyethylene Surface; Day: 107; Time: 14
temp.scan <-  mixedsort(list.files('~/disk10TB/DARPA/MODTRANSimulated/Stage1_7.5_12um/Polyethylene',pattern="*scan.csv",full.names=TRUE)) # geometric elevation angle, day, time, reflectivity (i.e., 13*365*6*7)
files.scan <- MODTRAN.readfiles.scan(temp.scan)
table.scan <-  MODTRAN.combineGeo(files.scan, eles)
table.scan<- MODTRAN.calculate.downwelling(table.scan)
# generate matrix CSV
mat.total_rad <- MODTRAN.mat(table.scan,"total_rad")
mat.up1 <- MODTRAN.mat(table.scan,"path_emiss")
mat.up2 <- MODTRAN.mat(table.scan,"path_thermal_scat")
mat.grnd <- MODTRAN.mat(table.scan,"grnd_rflt")
mat.down <- MODTRAN.mat(table.scan,"downwelling")
mat.surf <- MODTRAN.mat(table.scan,"surface_emiss")
mat.trans <- MODTRAN.mat(table.scan,"path_trans")
# mat.emiss <- MODTRAN.mat(table.scan,"surface_emiss")

outputfolder <- '~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Polyethylene'
write.csv(mat.total_rad, file = paste(outputfolder,"/Radiance_107_14_total.csv",sep = ""))
write.csv(mat.trans, file = paste(outputfolder,"/Radiance_107_14_trans.csv",sep = ""))
write.csv(mat.up1, file = paste(outputfolder,"/Radiance_107_14_up1.csv",sep = ""))
write.csv(mat.up2, file = paste(outputfolder,"/Radiance_107_14_up2.csv",sep = ""))
write.csv(mat.down, file = paste(outputfolder,"/Radiance_107_14_down.csv",sep = ""))
write.csv(mat.grnd, file = paste(outputfolder,"/Radiance_107_14_grnd.csv",sep = ""))
write.csv(mat.surf, file = paste(outputfolder,"/Radiance_107_14_surf.csv",sep = ""))
# predicted files
predictfile.trans <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_trans.csv"
predictfile.surf <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_surf.csv"
predictfile.grnd <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_grnd.csv"
predictfile.down <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_down.csv"
predictfile.up1 <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_up1.csv"
predictfile.up2 <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/predictedCSV/Radiance_107_14_up2.csv"
# original files
originalfile.trans<- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_trans.csv"
originalfile.surf<- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_surf.csv"
originalfile.grnd <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_grnd.csv"
originalfile.down <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_down.csv"
originalfile.up1 <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_up1.csv"
originalfile.up2<- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_up2.csv"
totalInputfile <- "~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/Polyethylene/Radiance_107_14_total.csv"
# read predicted and original
predorigfile.trans.read <- read.predorigmatrix(predictfile.trans, originalfile.trans)
predorigfile.surf.read <- read.predorigmatrix(predictfile.surf, originalfile.surf)
predorigfile.grnd.read <- read.predorigmatrix(predictfile.grnd, originalfile.grnd)
predorigfile.down.read <- read.predorigmatrix(predictfile.down, originalfile.down)
predorigfile.up1.read <- read.predorigmatrix(predictfile.up1, originalfile.up1)
predorigfile.up2.read <- read.predorigmatrix(predictfile.up2, originalfile.up2)
totalInputfile.read <- read.originalmatrix(totalInputfile)
# calculate zlim
grnd.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.grnd.read[[1]], predorigfile.grnd.read[[2]], SIMPLIFY = FALSE)
down.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.down.read[[1]], predorigfile.down.read[[2]], SIMPLIFY = FALSE)
up1.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.up1.read[[1]], predorigfile.up1.read[[2]], SIMPLIFY = FALSE)
up2.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.up2.read[[1]], predorigfile.up2.read[[2]], SIMPLIFY = FALSE)
surf.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.surf.read[[1]], predorigfile.surf.read[[2]], SIMPLIFY = FALSE)
diff.minmax <- minmax.crosstables(c(grnd.diff, up1.diff, up2.diff))
#diff.minmax <- c(-2.189695e-05,1.428151e-05)
# generate plots
comparisonMatrixPlot(predorigfile.trans.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/",residual.opt = "T",residual.square.opt = F,zname="Transmission", zunit=NULL)
comparisonMatrixPlot(predorigfile.grnd.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
comparisonMatrixPlot(predorigfile.down.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/",residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax)
comparisonMatrixPlot(predorigfile.surf.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/",residual.opt = "T", residual.square.opt = F)
# path thermal emission component
comparisonMatrixPlot(predorigfile.up1.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/", residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
# path thermal scattering component
comparisonMatrixPlot(predorigfile.up2.read, "~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/",residual.opt = "T", residual.square.opt = F, zlim.residual = diff.minmax)
MODTRAN.matplot(as.matrix(totalInputfile.read[[1]]), filename=paste("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/Polyethylene/", names(totalInputfile.read), ".jpg",sep = ""))


### Different temperatures, 300, 310, 330, 340, 350, 360, emissivity 0.9
folder.scan <- list.files('~/disk10TB/DARPA/SimulatedData/DifferentTemp',full.names=TRUE)
temp.scan <- lapply(folder.scan, function(x) list.files(x,pattern="*scan.csv",full.names=TRUE))
temp.scan <- lapply(temp.scan, mixedsort) # geometric elevation angle, day, time, reflectivity (i.e., 13*365*6*7)
files.scan <-lapply(temp.scan, MODTRAN.readfiles.scan)
table.scan <-  lapply(files.scan, function(x) MODTRAN.combineGeo.scan(x, eles)[[1]])
table.scan<-  lapply(table.scan, MODTRAN.calculate.downwelling)
# generate matrix CSV
mat.total_rad <- lapply(table.scan, function(x) MODTRAN.mat(x,"total_rad"))
mat.up1 <- lapply(table.scan, function(x) MODTRAN.mat(x,"path_emiss")) 
mat.up2 <- lapply(table.scan, function(x) MODTRAN.mat(x,"path_thermal_scat"))
mat.grnd <- lapply(table.scan, function(x) MODTRAN.mat(x,"grnd_rflt"))
mat.down <- lapply(table.scan, function(x) MODTRAN.mat(x,"downwelling"))
mat.surf <- lapply(table.scan, function(x) MODTRAN.mat(x,"surface_emiss"))
outputfolder <- list.files('~/disk10TB/DARPA/DifferentTemp/originalCSV',full.names=TRUE)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","total",".csv",sep = "")), mat.total_rad, outputfolder)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","up1",".csv",sep = "")), mat.up1, outputfolder)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","up2",".csv",sep = "")), mat.up2, outputfolder)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","down",".csv",sep = "")), mat.down, outputfolder)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","grnd",".csv",sep = "")), mat.grnd, outputfolder)
mapply(function(x, y) write.csv(x, file = paste(y,"/Radiance_107_14","_","surf",".csv",sep = "")), mat.surf, outputfolder)
# predicted files
predictfile.grnd <-list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/predictedCSV", pattern = "*grnd.csv", recursive=TRUE, full.names=TRUE)
predictfile.down <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/predictedCSV", pattern = "*down.csv", recursive=TRUE, full.names=TRUE)
predictfile.up1 <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/predictedCSV", pattern = "*up1.csv", recursive=TRUE, full.names=TRUE)
predictfile.up2 <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/predictedCSV", pattern = "*up2.csv", recursive=TRUE, full.names=TRUE)
predictfile.surf <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/predictedCSV", pattern = "*surf.csv", recursive=TRUE, full.names=TRUE)
# original files
originalfile.grnd <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*grnd.csv", recursive=TRUE, full.names=TRUE)
originalfile.down <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*down.csv", recursive=TRUE, full.names=TRUE)
originalfile.up1 <- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*up1.csv", recursive=TRUE, full.names=TRUE)
originalfile.up2<- list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*up2.csv", recursive=TRUE, full.names=TRUE)
originalfile.surf<-list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*surf.csv", recursive=TRUE, full.names=TRUE)
totalInputfile <-list.files("~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentTemp/originalCSV", pattern = "*total.csv", recursive=TRUE, full.names=TRUE)
# read predicted and original
predorigfile.grnd.read <- read.predorigmatrix(predictfile.grnd, originalfile.grnd) # for six different temperatures 
predorigfile.down.read <- read.predorigmatrix(predictfile.down, originalfile.down)
predorigfile.up1.read <- read.predorigmatrix(predictfile.up1, originalfile.up1)
predorigfile.up2.read <- read.predorigmatrix(predictfile.up2, originalfile.up2)
predorigfile.surf.read <- read.predorigmatrix(predictfile.surf, originalfile.surf)
totalInputfile.read <- read.originalmatrix(totalInputfile)
# calculate zlim
grnd.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.grnd.read[[1]], predorigfile.grnd.read[[2]], SIMPLIFY = FALSE)
down.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.down.read[[1]], predorigfile.down.read[[2]], SIMPLIFY = FALSE)
up1.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.up1.read[[1]], predorigfile.up1.read[[2]], SIMPLIFY = FALSE)
up2.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.up2.read[[1]], predorigfile.up2.read[[2]], SIMPLIFY = FALSE)
surf.diff <- mapply(function(x,y) diff.matrices(x, y), predorigfile.surf.read[[1]], predorigfile.surf.read[[2]], SIMPLIFY = FALSE)
diff.minmax <- minmax.crosstables(c(unlist(down.diff, recursive = FALSE), unlist(up1.diff, recursive = FALSE),unlist(up2.diff, recursive = FALSE)))
diff.minmax.grnd <- minmax.crosstables(unlist(grnd.diff, recursive = FALSE))
minmax.grnd <- minmax.crosstables(unlist(predorigfile.grnd.read, recursive = FALSE))
diff.minmax.surf <- minmax.crosstables(unlist(surf.diff, recursive = FALSE))
# generate plots
plotoutput <- list.files("~/geolab_storage_V3/data/DARPA/Plots/matrixPlotsML/DifferentTemp",full.names=TRUE)
plotoutput <- paste(plotoutput, "/",sep="")
# lapply(predorigfile.grnd.read, '[', 2)[[1]]==predorigfile.grnd.read[[1]][[2]] extract the 2nd element with its name of each list in a nest list
lapply(1:length(plotoutput), function(i) comparisonMatrixPlot(lapply(predorigfile.grnd.read, '[', i), plotoutput[i], residual.opt = "T",residual.square.opt = F, zlim= minmax.grnd, zlim.residual = diff.minmax.grnd))
lapply(1:length(plotoutput), function(i) comparisonMatrixPlot(lapply(predorigfile.down.read, '[', i), plotoutput[i], residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax))
lapply(1:length(plotoutput), function(i) comparisonMatrixPlot(lapply(predorigfile.up1.read, '[', i), plotoutput[i], residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax))
lapply(1:length(plotoutput), function(i) comparisonMatrixPlot(lapply(predorigfile.up2.read, '[', i), plotoutput[i], residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax))
lapply(1:length(plotoutput), function(i) comparisonMatrixPlot(lapply(predorigfile.surf.read, '[', i), plotoutput[i], residual.opt = "T",residual.square.opt = F, zlim.residual = diff.minmax.surf))
lapply(1:length(plotoutput), function(i) MODTRAN.matplot(as.matrix(totalInputfile.read[[i]]), filename=paste(plotoutput[i], names(totalInputfile.read)[i], ".jpg",sep = "")))

