library(MetaboMate) # to compute the PCA
library(hastaLaVista) # to visualize the data using hastaLaVista
library(car) # to compute the ellipses

get.idx=function(range=c(1,5), ppm){
  range=sort(range, decreasing=T);
  which(ppm<=range[1] & ppm>=range[2])}

# X <- "your matrix of data"
X <- bariatricRat.binned.4$binned.data

# x_axis <- "your x-axis"
x_axis <- bariatricRat.binned.4$binned.ppm

# ID <- "a vector with unique ids"
ID <- seq(1, nrow(X))

# group <- country # or the metadata you wish to use for coloring
group <- bariatricRat$metadata$Class

# metadata <- "a data.frame with your metadata"
metadata <- bariatricRat$metadata

# define the blocks
rangeList <- list(c(0.5, 1.000), 
                  c(1.0005, 1.5), 
                  c(1.5005, 2), 
                  c(2.0005, 2.5), 
                  c(2.5005, 3), 
                  c(3.0005, 3.5), 
                  c(3.5005, 4.0), 
                  c(6.5005, 7),
                  c(7.0005, 7.5),
                  c(7.5005, 8),
                  c(8.0005, 8.5),
                  c(8.5005, 9),
                  c(9.0005, 9.49))

# ONLY IF YOU HAVE VERY LARGE FILES
# if your data are large you may save each spectrum as individual file to avoid
# loading them all in the browser. Modify and uncomment the following lines to fit
# you system. The location of the data should be within the package folder 
# visu/data/json folder you can find out the path using the ".libPathts()" command.

# # create json spectra
# path = "/home/jul/R/x86_64-redhat-linux-gnu-library/3.6/visualizeR/visu/data/json/"
# for (i in seq_along(ID)) {
#   print(i)
#   filename = paste0(ID[i],".json")
#   saveJSON(list("y" = as.numeric(t(X[i,]))), path, filename)
# }
# metadata['JcampUrl'] <- paste0('http://127.0.0.1:5474/data/json/', filename, '.json')
# END OF LARGE FILE SECTION

# this line create a visualization object
v <- new("visualization")

# the name given to the data file stored in the visualizer package folder /visu/data
v@data <- "multiblocking.data.json"

# the name of the view file stored in the visualizer package folder /visu/view
# choose this view if data are small and if the dataMatrix at line 50 is not commented
v@view <- "modelExplorer_1_1.view.json" 

# choose the view below instead, if data are large and original spectra are retrieved individually 
# See line 38-41 and comment line 50
#v@view <- "modelExplorer_1_0.view.json" 

#======================================================================================
# no modification should be necessary below this line
#======================================================================================

color <- sapply(group, function(x) getColor2(as.character(x)))

d = list()
c <- data.frame(ID = ID,
                group = group,
                color = color,
                "_highlight" = seq_along(group) - 1,
                dataMatrix = I(unname(X)), # if data are large comment this line (you need to provide 
                # metadata['JcampUrl'] in order to be able to see your original data instead)
                metadata = I(metadata),
                check.names = FALSE
)

d <- appendData(data = d, variableName = "data", variable = c, type = "table")
remove(c)
d <- appendData(data = d, variableName = "xAxis", variable = x_axis, type = "table")

block <- list()
model <- list()
            
it <- 1
for (ran in rangeList) {
idx <- get.idx(range(ran), as.numeric(x_axis))
pcaName <- paste0("model",rangeList[it])

mod <- MetaboMate::pca(X[, idx], pc = 3, method = 'nipals')

ellipse <- dataEllipse(mod@t[,1], mod@t[,2], levels=0.80, draw = FALSE)

ellipseChart <- data.frame("x" = ellipse[,1],
                           "y" = ellipse[,2],
                           "color" = rep('black', length(ellipse[,1])))

c <- list(name = pcaName,
    scores = mod@t,
    expVar = mod@R2,
    loadings = cov(mod@t, X[, idx]),
    loadingsXaxis = x_axis[idx],
    loadingsColor = abs(cor(mod@t, X[, idx])),
    ellipses = ellipseChart
)

model[[it]] <- c
remove(c)
block[[it]] <- X[, idx]
it <- it + 1
} ## end of loop

## we compute the pca for the whole spectra

pcaName <- paste0("model"," full spectra")

mod <- MetaboMate::pca(X, pc = 3, method = 'nipals')

ellipse <- dataEllipse(mod@t[,1], mod@t[,2], levels=0.80, draw = FALSE)

ellipseChart <- data.frame("x" = ellipse[,1],
                           "y" = ellipse[,2],
                           "color" = rep('black', length(ellipse[,1])))

c <- list(name = pcaName,
          scores = mod@t,
          expVar = mod@R2,
          loadings = cov(mod@t, X),
          loadingsXaxis = x_axis,
          loadingsColor = abs(cor(mod@t, X)),
          ellipses = ellipseChart
)

model[[it]] <- c
remove(c)

## we compute the pca of the superscores

it <- it + 1

pcaName <- paste0("model"," super score")

newLength <- length(model)-1
newX <- matrix(unlist(lapply(model[1:newLength], function(x) x$scores)), nrow(X), 3 * newLength)
mod <- MetaboMate::pca(newX, pc = 4, method = 'nipals')

ellipse <- dataEllipse(mod@t[,1], mod@t[,2], levels=0.80, draw = FALSE)

ellipseChart <- data.frame("x" = ellipse[,1],
                           "y" = ellipse[,2],
                           "color" = rep('black', length(ellipse[,1])))

c <- list(name = pcaName,
          scores = mod@t,
          expVar = mod@R2,
          loadings = cov(mod@t, X),
          loadingsXaxis = x_axis,
          loadingsColor = abs(cor(mod@t, X)),
          ellipses = ellipseChart
)

model[[it]] <- c
remove(c)

d[['model']] <- model

push(v, type="data", d)
print(v) # print out the link 
visualize(v) # start the server and point your browser to the link
