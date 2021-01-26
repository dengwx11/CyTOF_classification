#https://biosurf.org/cytof_data_scientist.html

require(ncdfFlow)
require(cydar)
require(flowCore)

files = list.files(path = './fcs_data',pattern = "*\\.fcs", full.names = TRUE)
fcs <- read.FCS(filename=files[1], transformation=FALSE)
exprs <- fcs@exprs
fcs@parameters@data
markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
colnames(exprs)[which(!is.na(markers))] <- markers[which(!is.na(markers))]