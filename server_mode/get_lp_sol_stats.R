### LOAD LIBRARIES --------------------------------------------------------------------------
library("stringr")

get_lp_sol_stats <- function(dir.sol.parse, dir.batch, name.batch, file.bipartite){

### List of files in the directory ---------------------------------------------------------- 	
	filenames <- list.files(dir.sol.parse, pattern="*.txt", full.names=TRUE)
	name <- unlist(lapply(str_split(filenames, "/"), function(x) x[length(x)]))
	f <- as.numeric(unlist(lapply(lapply(str_split(name, "_"), function(x) x[length(x)]), function(x) str_extract(x, "\\d{1,15}"))))

	val.x <- list()
	val.y <- list()
	val.z <- list()

	par$TotalAbr <- rep(0, nrow(par)) #TOTAL ABERRANT GENES
	par$TotalEdg <- rep(0, nrow(par)) #TOTAL OUTLIERS
	par$TotalOut <- rep(0, nrow(par)) #TOTAL EDGES BETWEEN ABR-OUT
	
	par$x <- rep(0, nrow(par))
	par$y <- rep(0, nrow(par))
	par$z <- rep(0, nrow(par))
	par$DriverGene <- rep("", nrow(par))


### PARSE SOL FILES ----------------------------------------------------------------------------	
	row.num <- 1
	for(i in 1:length(f)){
		file.sol <- filenames[i]
	
		### LOAD SOL FILES ----------------------------------------------------------------------
		dat <- read.delim(file.sol, header=F, stringsAsFactors=F)
		colnames(dat) <- c("Solution","Class","Name","Value")
	
		dat$Value <- round(dat$Value, 1)
		TotalSol <- unique(dat$Solution)
		#TotalSol <- 1
	
		for(j in 1:length(TotalSol)){
			dat.temp <- dat[which(dat$Solution == j),]
	
			dat.Edges <- subset(dat.temp, dat.temp$Class == "EDGE")[-c(1,2)]
			dat.Zeta <- subset(dat.temp, dat.temp$Class == "ZETA")[-c(1,2)]
			dat.Node <- subset(dat.temp, dat.temp$Class == "NODE")[-c(1,2)]
				
			### Create Matrix
			if(i == 1){
				dat.x <- matrix(0, nrow=length(dat.Node$Name), ncol=length(par$alpha[f]), dimnames=list(dat.Node$Name, paste(par$gamma[f], par$alpha[f], sep="_")))
				dat.y <- matrix(0, nrow=length(dat.Edges$Name), ncol=length(par$alpha[f]), dimnames=list(dat.Edges$Name, paste(par$gamma[f], par$alpha[f], sep="_")))
				dat.z <- matrix(0, nrow=length(dat.Zeta$Name), ncol=length(par$alpha[f]), dimnames=list(dat.Zeta$Name, paste(par$gamma[f], par$alpha[f], sep="_")))
			}
			
			DriverGene <- paste(dat.Node$Name[which(dat.Node$Value == 1)], collapse=",")
		
			val.x[[row.num]] <- dat.Node$Name[which(dat.Node$Value == 1)]
			val.y[[row.num]] <- dat.Edges$Name[which(dat.Edges$Value == 1)]
			val.z[[row.num]] <- dat.Zeta$Name[which(dat.Zeta$Value == 1)]
				
			par$TotalAbr[f[i]] <- nrow(dat.Node)
			par$TotalEdg[f[i]] <- nrow(dat.Edges)
			par$TotalOut[f[i]] <- nrow(dat.Zeta)
			
			par$x[f[i]] <- length(val.x[[row.num]])
			par$y[f[i]] <- length(val.y[[row.num]])
			par$z[f[i]] <- length(val.z[[row.num]])
			par$DriverGene[f[i]] <- DriverGene
			#par$Outliers[f[i]] <- paste(val.z[[row.num]],collapse=",")
			#par$Edges[f[i]] <- paste(val.y[[row.num]],collapse=",")
		
			row.num <- row.num  + 1
		}
	
		cat(paste(Sys.time()), "PARSED", i, "OF", length(f),"\n", sep="\t")
	}	
	colnames(par) <- c("gamma","alpha","beta","iter","TotalAbr","TotalEdg","TotalOut","x","y","z","DriverGene")
	
	
	### Fill Binary Data Matrix for X, Y and Z
	for(ctr in 1:length(par$alpha[f])){
		index.x <- which(rownames(dat.x) %in% val.x[[ctr]])
		index.y <- which(rownames(dat.y) %in% val.y[[ctr]])
		index.z <- which(rownames(dat.z) %in% val.z[[ctr]])
		
		dat.x[index.x,ctr] <- 1
		dat.y[index.y,ctr] <- 1
		dat.z[index.z,ctr] <- 1
	}
	
	dat.x <- dat.x[,order(colnames(dat.x),decreasing=F)]
	dat.y <- dat.y[,order(colnames(dat.y),decreasing=F)]
	dat.z <- dat.z[,order(colnames(dat.z),decreasing=F)]
	
	
	### Write X,Y,Z matrix files
	file.x <- file.path(dir.batch, paste(name.batch, "_drivers.txt", sep=""))
	file.y <- file.path(dir.batch, paste(name.batch, "_edges.txt", sep=""))
	file.z <- file.path(dir.batch, paste(name.batch, "_outliers.txt", sep=""))
	
	write.table(dat.x, file.x, sep="\t", row.names=T, col.names=NA, quote=F)
	write.table(dat.y, file.y, sep="\t", row.names=T, col.names=NA, quote=F)
	write.table(dat.z, file.z, sep="\t", row.names=T, col.names=NA, quote=F)

	par <- par[!(par$x == 0),]
	rownames(par) <- c(1:nrow(par))
	
	gene.sp <- str_split(par$DriverGene, ",")

	#func.x <- function(x){
		#val.driver.cgc.name <- intersect(dat.cgc$Symbol, x)
		#val.driver.cgc.num <- length(val.driver.cgc.name)
	
		#val.driver.cosmic.name <- intersect(dat.cosmic$COSMIC_GENE_NAME, x)
		#val.driver.cosmic.num <- length(val.driver.cosmic.name)
	
		#val.cgc.ratio <- round(val.driver.cgc.num/length(x),2)
		#val.cosmic.ratio <- round(val.driver.cosmic.num/length(x),2)
	
		#result <- c(val.driver.cgc.num, val.cgc.ratio, val.driver.cosmic.num, val.cosmic.ratio)
		#return(result)
	#}
	
	#driver.cgc.list <- lapply(gene.sp, func.x)
	#driver.cgc.list <- do.call(rbind.data.frame, driver.cgc.list)
	#colnames(driver.cgc.list) <- c("CGC.DriverNum","CGC.DriverRatio","COSMIC.DriverNum","COSMIC.DriverRatio")

	#par$CGC.DriverNum <- driver.cgc.list$CGC.DriverNum
	#par$CGC.DriverRatio <- driver.cgc.list$CGC.DriverRatio
	#par$COSMIC.DriverNum <- driver.cgc.list$COSMIC.DriverNum
	#par$COSMIC.DriverRatio <- driver.cgc.list$COSMIC.DriverRatio

	par <- par[order(par$gamma, decreasing=F),]
	
	#Write Stats file
	file.stats <- file.path(dir.batch, paste(name.batch, "_ilpsol_stats.txt", sep=""))
	write.table(par, file.stats, sep="\t", row.names=F, col.names=T, quote=F)
	
	source(file.path(dir.script, "get_patient_specific_driver.R"))
}
