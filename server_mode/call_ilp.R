### LOAD LIBRARIES --------------------------------------------------------
library("stringr")
library("foreach")
library("doParallel")

#### FUNCTION TO COMPUTE TOTAL INFLUENCE FOR EACH OUTLIER -----------------
func.inf_out <- function(BiPartiteGraph){
	Unique.Outliers <- sort(unique(BiPartiteGraph$OutlierEvent), decreasing=FALSE)

	cat(paste(Sys.time()), "TOTAL ITERATION REQUIRED :", length(Unique.Outliers),"\n", sep="\t")
	cat(paste(Sys.time()), "ITERATING ... PLEASE WAIT ...","\n", sep="\t")
	
	#no_cores <- 30
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
	
	cat(paste(Sys.time()), "NUMBER OF CORES IN USE: ", no_cores, "\n", sep="\t")
	
	lpar <- foreach(ctr = 1:length(Unique.Outliers), .inorder=FALSE, .combine=rbind) %dopar% {
				inf.out <- BiPartiteGraph$Influence[which(BiPartiteGraph$OutlierEvent == Unique.Outliers[ctr])]		
				#sumInfluence <- format(sum(inf.out), scientific = NA)
				lpar <- data.frame(Outlier=Unique.Outliers[ctr], sumInfluence=sum(inf.out))
	}
	stopCluster(cl)
	
	lpar$Outlier <- as.character(lpar$Outlier)
	lpar$sumInfluence <- as.numeric(lpar$sumInfluence)
	
	return(lpar)
}


#### FUNCTION TO CALL ILP FORMULATION -------------------------------------
call.ilp <- function(file.bipartite, par, beta, dir.input, name.batch, file.driverWt, seed.genes=NA){
	
#### LOAD BIPARTITE GRAPH -------------------------------------------------
	cat(paste(Sys.time()), "LOADING BIPARTITE GRAPH IN MEMORY ...","\n", sep="\t")
	bip <- read.delim(file.bipartite, header=TRUE, stringsAsFactors=FALSE)
	
	bip$MutGene <- gsub("-","_", bip$MutGene)             # Convert DASH "-" symbol to UNDERSCORE "_" symbol
	bip$OutlierEvent <- gsub("-","_", bip$OutlierEvent)   # Convert DASH "-" symbol to UNDERSCORE "_" symbol
  
	cat(paste(Sys.time()), "LOADED ...","\n", sep="\t")
	cat("\n")
  
#### COMPUTE TOTAL INFLUENCE ----------------------------------------------
	cat(paste(Sys.time()), "COMPUTING TOTAL INFLUENCE ON EACH OUTLIER ...","\n", sep="\t")
	influence.outlier <- func.inf_out(bip)
	cat(paste(Sys.time()), "COMPUTED TOTAL INFLUENCE ON EACH OUTLIER ...","\n", sep="\t")
	cat("\n")

### PARSE BIPARTITE GRAPH --------------------------------------------------
	cat(paste(Sys.time()), "PARSING BIPARTITE GRAPH ...","\n", sep="\t")
	y1 <- which(colnames(bip) == "MutGene")
	y2 <- which(colnames(bip) == "OutlierEvent")
	z1 <- which(colnames(bip) == "OutlierEvent")
	
	edge <- apply(bip, 1, function(x) paste("Y",paste(x[y1],x[y2], sep=";"), sep=";"))
	zeta <- apply(bip, 1, function(x) paste("Z",x[z1],sep=";"))

	bip$Edge <- edge
	bip$Zeta <- zeta
  
	patient.sep <- str_split(bip$OutlierEvent, ";")
	patient.sep <- do.call(rbind.data.frame, patient.sep)
	colnames(patient.sep) <- c("Patient","OutlierGene")
  
	bip$Patient <- patient.sep$Patient
	bip$Patient <- as.character(bip$Patient)
	
	cat(paste(Sys.time()), "PARSING BIPARTITE GRAPH : COMPLETE...","\n", sep="\t")
	cat("\n")

### CALL FUNCTION MODULE-I -----------------------------------------------------
	file.driverWt <- file.driverWt
	list.module_i <- func.module_i(bip, beta, file.driverWt, seed.genes)

### CALL FUNCTION MODULE-I -----------------------------------------------------
	list.module_ii <- func.module_ii(bip, par, influence.outlier)

	cat(paste(Sys.time()), "MERGE ALL THE CONSTRAINTS START: STEP 7 OF 8  ...","\n", sep="\t")
	
### MERGE ALL THE CONSTRAINTS ---------------------------------------------------
	oper.A <- "MINIMIZE"
	oper.B <- "SUBJECT TO"
	oper.C <- "BINARY"
	#oper.D <- "GENERAL"
	oper.E <- "END"

	opt.func <- list.module_i$opt.func
	constraint_2 <- list.module_i$constraint_2
	if(any(colnames(bip) == "OutlierWt")){
		constraint_3_beta_lambda <- list.module_i$constraint_3_beta_lambda
	}
	constraint_4 <- list.module_i$constraint_4
	constraint_5 <- list.module_i$constraint_5
	constraint_6 <- list.module_i$constraint_6
	
	if(!is.na(seed.genes[1])){
		constraint_seed <- list.module_i$constraint_seed
	}
	
	list.constraint <- list()
	for(i in 1:nrow(par)){
		constraint_1 <- list.module_ii$constraint_1[[i]]
		constraint_3 <- list.module_ii$constraint_3[[i]]

			if(any(colnames(bip) == "OutlierWt")){
				if(is.na(seed.genes[1])){
					list.constraint[[i]] <- list(oper.A, opt.func,
						oper.B, constraint_1, constraint_2, constraint_3, constraint_3_beta_lambda,
						oper.C, constraint_4, constraint_5, constraint_6,
						oper.E)
				} else{
					list.constraint[[i]] <- list(oper.A, opt.func,
						oper.B, constraint_1, constraint_2, constraint_3, constraint_3_beta_lambda, constraint_seed,
						oper.C, constraint_4, constraint_5, constraint_6,
						oper.E)
				}
			} else{
				list.constraint[[i]] <- list(oper.A, opt.func,
					oper.B, constraint_1, constraint_2, constraint_3,
					oper.C, constraint_4, constraint_5, constraint_6,
					oper.E)
			}

		cat(paste(Sys.time()), "ITERATION COMPLETE: ", i, "OF", nrow(par), "\n", sep="\t")
	}
	cat(paste(Sys.time()), "STEP 8 OF 9 COMPLETE ...","\n", sep="\t")

#### GENERATE ILP FILES ------------------------------------------------------------
	cat(paste(Sys.time()), "WRITING ILP FILES: STEP 9 OF 9 ...","\n", sep="\t")
	for(i in 1:nrow(par)){
		file.output <- file.path(dir.input, paste(name.batch,"_0",i,".lp",sep=""))

		file.out <- file(file.output, open="a")
		for(j in 1:length(list.constraint[[i]])){
			write.table(list.constraint[[i]][[j]], file.out, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
		close(file.out)
	
		cmd <- paste("gzip", file.output, sep=" ")
		system(cmd)
	
		cat("ILP FILE GENERATED :", i, "OF", nrow(par), paste(file.output,".gz", sep=""), "\n", sep=" ")
	}
	cat(paste(Sys.time()), "STEP 9 OF 9 COMPLETE ...","\n", sep="\t")
	cat(paste(Sys.time()), "END OF LINEAR PROGRAMING FORMULATION ...","\n", sep="\t")
	cat("\n")
	#return(list.constraint)
}

#### FUNCTION TO COMPUTE TOTAL INFLUENCE FOR EACH OUTLIER
#func.inf_out <- function(BiPartiteGraph)
#{
#	Unique.Outliers <- sort(unique(BiPartiteGraph$OutlierEvent), decreasing=FALSE)
#	Influence.Outlier <- data.frame(row.names=NULL)
#	
#	cat(paste(Sys.time()), "TOTAL ITERATION REQUIRED :", length(Unique.Outliers),"\n", sep="\t")
#	cat(paste(Sys.time()), "ITERATING ... PLEASE WAIT ...","\n", sep="\t")
#	for(i in 1:length(Unique.Outliers))
#	{
#		inf.out <- BiPartiteGraph$Influence[which(BiPartiteGraph$OutlierEvent == Unique.Outliers[i])]
#		Influence.Outlier[i,1] <- Unique.Outliers[i]
#		Influence.Outlier[i,2] <- format(sum(inf.out),scientific = NA)
#		cat("COMPLETE :", i, "OF", length(Unique.Outliers), "\n", sep="\t")
#	}
#	colnames(Influence.Outlier) <- c("Outlier","sumInfluence")
#	
#	#min.inf <- min(Influence.Outlier$sumInfluence)
#	
#	return(Influence.Outlier)
#}
