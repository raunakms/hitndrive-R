### Load Libraries --------------------------------------------------------------
library("bigmemory")

### Load MutMatrix --------------------------------------------------------------
cat(paste(Sys.time()), "LOADING ABERRATION DATA ....", "\n", sep="\t")
	dat.abr <- read.table(file.abr, sep="\t", header=T, stringsAsFactors=F)
cat(paste(Sys.time()), "LOADED ....", "\n", sep="\t")

### Load OutlierMatrix ----------------------------------------------------------
cat(paste(Sys.time()), "LOADING OUTLIER DATA ....", "\n", sep="\t")
	dat.outlier <- read.table(file.outlier, sep="\t", header=T, stringsAsFactors=F)
cat(paste(Sys.time()), "LOADED ....", "\n", sep="\t")

### Load Influence Graph ---------------------------------------------------------
cat(paste(Sys.time()), "LOADING INFLUENCE MATRIX ....", "\n", sep="\t")
	nodes <- read.table(file.nodes, header=F, as.is="V1")$V1
	influence.graph <- read.big.matrix(filename=file.network, sep = "\t", header = TRUE, 
						col.names = nodes, has.row.names=TRUE, ignore.row.names=FALSE,
						type = "double", skip = 0)
	influence.graph <- as.matrix(influence.graph)		
cat(paste(Sys.time()), "LOADED ....", "\n", sep="\t")

### BipartiteConstruction ---------------------------------------------------------
cat(paste(Sys.time()), "BI-PARTITE GRAPH CONSTRUCTION STARTED ....", "\n", sep="\t")
	source(file.path(dir.script, "get_bipartite_graph_parallel.R"))
	BiPartiteGraph <- getBiPartiteGraph(
					dat.abr = dat.abr, 
					dat.outlier = dat.outlier, 
					influence.graph = influence.graph)
cat(paste(Sys.time()), "BI-PARTITE GRAPH CONSTRUCTION FINISHED ....", "\n", sep="\t")

### Output: BipartiteGraph ---------------------------------------------------------
cat(paste(Sys.time()), "WRITING FILE ....", "\n", sep="\t")				
	write.table(BiPartiteGraph, file=file.bip, sep="\t", col.names=T, row.names=F, quote=F)		
cat(paste(Sys.time()), "FILE CREATED ....", "\n", sep="\t")

### Output: Compression -------------------------------------------------------------
cat(paste(Sys.time()), "FILE COMPRESSING STARTED....", "\n", sep="\t")		
	cmd <- paste("gzip", file.bip, sep=" ")
	system(cmd)
cat(paste(Sys.time()), "FILE COMPRESSING COMPLETED....", "\n", sep="\t")	
cat(paste(Sys.time()), "FILE GENERATED:", file.bip, "\n", sep="\t")				
