### Load Libraries --------------------------------------------------------------
library("stringr")
library("foreach")
library("doParallel")

getBiPartiteGraph <- function(dat.abr, dat.outlier, influence.graph){
### Get GeneList --------------------------------------------------------------
	genelist.inf <- colnames(influence.graph)
	genelist.abr <- unique(dat.abr$Gene)
	genelist.out <- unique(dat.outlier$Gene)

### Match Genelist in MutMatrix (or OutlierMatrix) with GeneList of InfGraph ----
	genes.abr.inf <- intersect(genelist.abr, genelist.inf)
	genes.out.inf <- intersect(genelist.out, genelist.inf)

### Create Subset of MutMatrix and OutlierMatrix --------------------------------
	dat.abr <- subset(dat.abr, dat.abr$Gene %in% genes.abr.inf)	
	dat.out <- subset(dat.outlier, dat.outlier$Gene %in% genes.out.inf)
	
	upatients <- intersect(unique(dat.abr$SampleID), unique(dat.out$SampleID))
	dat.abr <- subset(dat.abr, dat.abr$SampleID %in% upatients)
	dat.out <- subset(dat.out, dat.out$SampleID %in% upatients)

### Combine SampleID;OutlierGene for labelling --------------------------------	
	dat.out$OutlierEvent <- paste(dat.out$SampleID, dat.out$Gene, sep=";")

### BIPARTITE GRAPH ------------------------------------------------------------		
	cat(paste(Sys.time()), "COMPUTING BIPARTITE GRAPH ....", "\n", sep="\t")
	
	#no_cores <- 30
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
	
	cat(paste(Sys.time()), "NUMBER OF CORES IN USE: ", no_cores, "\n", sep="\t")
	
	bip.lpar <- foreach(i = 1:nrow(dat.out), .combine=rbind, .inorder=FALSE, .packages="stringr") %dopar%{
					events.out <- str_split(dat.out$OutlierEvent[i], ";")[[1]]
					id <- events.out[1]
					outgene <- events.out[2]
					abr.genelist <- subset(dat.abr, dat.abr$SampleID == id)$Gene
					inf.val <- influence.graph[abr.genelist, outgene]
					
					if(any(colnames(dat.out) == "Weight")){
						list.dat <- data.frame(MutGene=abr.genelist,
											OutlierEvent=dat.out$OutlierEvent[i],
											Influence=inf.val,
											OutlierWt=dat.out$Weight[i])
					}else{
						list.dat <- data.frame(MutGene=abr.genelist,
											OutlierEvent=dat.out$OutlierEvent[i],
											Influence=inf.val)
					}
	}
	stopCluster(cl)
	rownames(bip.lpar) <- c(1:nrow(bip.lpar))

	bip.lpar <- subset(bip.lpar, bip.lpar$Influence != 0)
	
	cat(paste(Sys.time()), "BIPARTITE GRAPH CREATED ....", "\n", sep="\t")				
	return(bip.lpar)
}
