######## HIT'nDRIVE ILP FORMULATION ########
## AUTHOR: RAUNAK SHRESTHA (sraunak@gmail.com)
##
## OPTIMIZATION FUNCTION :
## 	MINIMIZE : sum_i(x_i)
##
## OPTIMIZATION FUNCTION : DRIVER WEIGHT VERSION.
## 	MINIMIZE : sum_i(s_i * x_i)
##
## SUBJECT TO:
## 	CONSTRAINT-1    : for each eta_j: sum_i(w_ij * y_ij) >= g * eta_j * sum_j(w_ij)
## 	CONSTRAINT-1_wt : for each eta_j: sum_i(w_ij * y_ij) >= g * lambda_j * eta_j * sum_j(w_ij)
## 	CONSTRAINT-2    : for each x_i  : y_ij = x_i 
## 	CONSTRAINT-3    : sum_j(eta_j) >= a * N
## 	CONSTRAINT-4    : x_i = [0,1]
## 	CONSTRAINT-5    : eta_j = [0,1]
## 	CONSTRAINT-6    : y_ij = [0,1]
##
## HERE, x_i		: aberrant genes
##       eta_j		: outlier genes
##       y_ij		: edges connecting x_i and eta_j in the bipartite graph
##		 w_ij		: influence falue of a aberrant gene over an outlier gene
##       lambda_j	: weight of an outlier gene
##		 s_i 		: Score or Weight (or p-value) of the potential driver gene (x_i). Lower scores/weight are stronger.
##       N			: total number of patient-specific outlier genes
##       a			: alpha (fraction of the total outliers to be covered)
##       g			: gamma (fraction of the total influence that a selected driver gene must at leaset have ouver the outlier)
###############################################

######### FUNCTION TO GENERATE ILP FORMULATION #########
#### FUNCTION MODULE-I
#### MODULE-I THOSE ILP FUNC. THAT DOES NOT REQUIRE GAMMA OR ALPHA ITERATIONS STEPS
func.module_i <- function(BiPartiteGraph, beta, file.driverWt, seed.genes=NA){
	BiPartiteGraph <- BiPartiteGraph
	beta <- beta
	
	cat(paste(Sys.time()), "COMPUTING LINEAR PROGRAMING FORMULATION ...","\n", sep="\t")
  
  		# CALL OPTIMIZATION FUNCTION ----
  		if(!is.na(file.driverWt)){
  			dat.nodeWt <- read.delim(file.driverWt, header=TRUE, stringsAsFactors=FALSE)
  			opt.func <- func.optimization.nodeWt(BiPartiteGraph, dat.nodeWt)
  		} else{
  			opt.func <- func.optimization(BiPartiteGraph)		
  		}
		
	
	cat(paste(Sys.time()), "STEP 1 OF 9: CALL FUNCTION: OPTIMIZATION  COMPLETE ...","\n", sep="\t")

		# CALL FUNCTION: CONSTRAINT-2 ----
		constraint_2 <- func.constraint_2(BiPartiteGraph)

	cat(paste(Sys.time()), "STEP 2 OF 9: CALL FUNCTION: CONSTRAINT-2  COMPLETE ...","\n", sep="\t")

		# CALL FUNCTION: CONSTRAINT-3 WITH OUTLIER WEIGHT ----
		if(any(colnames(BiPartiteGraph) == "OutlierWt")){
			constraint_3_beta_lambda <- func.constraint_3_beta_lambda(BiPartiteGraph, beta)
		}

	cat(paste(Sys.time()), "STEP 3 OF 9: CALL FUNCTION: CONSTRAINT-3  COMPLETE ...","\n", sep="\t")
	
		# CALL FUNCTION: CONSTRAINT-4 ----
		constraint_4 <- func.constraint_4(BiPartiteGraph)

	cat(paste(Sys.time()), "STEP 4 OF 9: CALL FUNCTION: CONSTRAINT-4 COMPLETE ...","\n", sep="\t")

		# CALL FUNCTION: CONSTRAINT-5 ----
		constraint_5 <- func.constraint_5(BiPartiteGraph)

	cat(paste(Sys.time()), "STEP 5 OF 9: CALL FUNCTION: CONSTRAINT-5 COMPLETE ...","\n", sep="\t")
	
		# CALL FUNCTION: CONSTRAINT-6 ----
		constraint_6 <- func.constraint_6(BiPartiteGraph)

	cat(paste(Sys.time()), "STEP 6 OF 9: CALL FUNCTION: CONSTRAINT-6 COMPLETE ...","\n", sep="\t")

	if(!is.na(seed.genes[1])){
		constraint_seed <- func.constraint_seed(seed.genes)
	}
	
	list.module_i <- list()
	list.module_i$opt.func <- opt.func
	list.module_i$constraint_2 <- constraint_2
	if(any(colnames(BiPartiteGraph) == "OutlierWt")){
		list.module_i$constraint_3_beta_lambda <- constraint_3_beta_lambda
	}
	list.module_i$constraint_4 <- constraint_4
	list.module_i$constraint_5 <- constraint_5
	list.module_i$constraint_6 <- constraint_6
	
	if(!is.na(seed.genes[1])){
		list.module_i$constraint_seed <- constraint_seed
	}
	
	return(list.module_i)
}

#### FUNCTION MODULE-II
#### MODULE-II THOSE ILP FUNC. THAT REQUIRES GAMMA OR ALPHA ITERATIONS STEPS
func.module_ii <- function(BiPartiteGraph, par, influence.outlier){
	BiPartiteGraph <- BiPartiteGraph
	Influence.Outlier <- influence.outlier

	cat(paste(Sys.time()), "STARTING LOOP FOR PARAMETER (gamma,alpha) SIMULATION: STEP 7 OF 9  ...","\n", sep="\t")

	list.module_ii <- list()

	for(i in 1:nrow(par)){
		cat(paste(Sys.time()), "ITERATION START: ", i, "OF", nrow(par), "\n", sep="\t") 
		
		gamma <- par$gamma[i]
		alpha <- par$alpha[i]
		
		if(any(colnames(BiPartiteGraph) == "OutlierWt")){
			constraint_1 <- func.constraint_1_wt(BiPartiteGraph, gamma, Influence.Outlier)
		}else{
			constraint_1 <- func.constraint_1(BiPartiteGraph, gamma, Influence.Outlier)
		}
		
		constraint_3 <- func.constraint_3(BiPartiteGraph, alpha)

		list.module_ii$constraint_1[[i]] <- constraint_1
		list.module_ii$constraint_3[[i]] <- constraint_3

		cat(paste(Sys.time()), "ITERATION COMPLETE: ", i, "OF", nrow(par), "\n", sep="\t")    
	}
	cat(paste(Sys.time()), "STEP 7 OF 9 COMPLETE ...","\n", sep="\t")

	return(list.module_ii)
}
########################################################

#### OPTIMIZATION FUNCTION ---------------------------------------------------------------------------------
func.optimization <- function(BiPartiteGraph){
	AbGenes <- sort(unique(BiPartiteGraph$MutGene), decreasing=F)
	AbGenes <- paste("@", AbGenes, sep="")     # Place @ symbol before every genename
	
	q <- length(AbGenes) %/% 10
	r <- length(AbGenes) %% 10
	
	if(r == 0){
		n <- q
	} else {
		n <- q + 1
	}
	
	opt.func <- rep(NA, n)
	a <- 0
	for(i in 1:n){	
		start <- a + 1
		stop <- a + 10
		
		if(stop > length(AbGenes)){
			if(i >= 2){
				opt.func[i] <- paste("+", paste(AbGenes[start:length(AbGenes)], collapse=" + "), sep=" ")
			} else {
				opt.func[i] <- paste(AbGenes[start:length(AbGenes)], collapse=" + ", sep=" ")
			}	
		} else {
			if(i >= 2){
				opt.func[i] <- paste("+", paste(AbGenes[start:stop], collapse=" + "), sep=" ")
			} else {
				opt.func[i] <- paste(AbGenes[start:stop], collapse=" + ", sep=" ")
			}
			a <-  stop
		}
	}
	return(opt.func)
}	

#### OPTIMIZATION FUNCTION: NODE WEIGHT -------------------------------------------------------
func.optimization.nodeWt <- function(BiPartiteGraph, dat.nodeWt){

	# GET ALTERED GENES ---
	AbGenes <- sort(unique(BiPartiteGraph$MutGene), decreasing=F)
	
	# MATCH GENES IN df.nodeWt ---
	dat.nodeWt <- subset(dat.nodeWt, dat.nodeWt$Gene %in% AbGenes)
	dat.nodeWt <- dat.nodeWt[match(AbGenes, dat.nodeWt$Gene),]

	nodeWt <- dat.nodeWt[,2]
	nodeWt[which(is.na(nodeWt))] <- 0
	nodeWt[which(nodeWt == 1)] <- 0

	# Place @ symbol before every genename ---
	AbGenes <- paste("@", AbGenes, sep="")     
	
	# GET QUOTIENT & REMAINDER ---
	q <- length(AbGenes) %/% 10
	r <- length(AbGenes) %% 10
	
	if(r == 0){
		n <- q
	} else {
		n <- q + 1
	}
	
	# PREPARE DATA ---
	opt.func <- rep(NA, n)
	a <- 0
	for(i in 1:n){	
		start <- a + 1
		stop <- a + 10
		
		if(stop > length(AbGenes)){
			if(i >= 2){
				opt.func[i] <- paste("+", paste(paste(nodeWt[start:length(AbGenes)], AbGenes[start:length(AbGenes)], sep=" x "), collapse=" + "), sep=" ")
			} else {
				opt.func[i] <- paste(paste(nodeWt[start:length(AbGenes)], AbGenes[start:length(AbGenes)], sep=" x "), collapse=" + ", sep=" ")
			}	
		} else {
			if(i >= 2){
				opt.func[i] <- paste("+", paste(paste(nodeWt[start:length(AbGenes)], AbGenes[start:length(AbGenes)], sep=" x "), collapse=" + "), sep=" ")
			} else {
				opt.func[i] <- paste(paste(nodeWt[start:length(AbGenes)], AbGenes[start:length(AbGenes)], sep=" x "), collapse=" + ", sep=" ")
			}
			a <-  stop
		}
	}
	return(opt.func)
}	

#### CONSTRAINT-1 ---------------------------------------------------------------------------------------------------------------
func.constraint_1 <- function(BiPartiteGraph, gamma, Influence.Outlier){
	g <- gamma
	Unique.Outliers <- sort(unique(BiPartiteGraph$OutlierEvent), decreasing=FALSE)

	#no_cores <- 30
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
		
	lpar <- foreach(i = 1:length(Unique.Outliers), .combine=c, .inorder=FALSE) %dopar% {	
				outlier <- Unique.Outliers[i]
				Sum.Inf <- Influence.Outlier$sumInfluence[which(Influence.Outlier$Outlier == outlier)]
				Sum.Inf.frac <- g * as.numeric(Sum.Inf) #### Cover g fraction
		
				bip.sub <- BiPartiteGraph[which(BiPartiteGraph$OutlierEvent == outlier),]
		
				z.func <- paste("-", Sum.Inf.frac, bip.sub$Zeta[1], sep=" ")
				Edg.Inf <- paste(bip.sub$Influence, bip.sub$Edge, sep=" ")
		
				q <- length(Edg.Inf) %/% 4
				r <- length(Edg.Inf) %% 4

				if(r == 0){
					n <- q
				} else {
					n <- q + 1
				}	
	
				dat.constraint_1 <- rep(NA, n)
				a <- 0
				for(j in 1:n){
					start <- a + 1
					stop <- a + 4
		
					if(stop > length(Edg.Inf)){
						if(j >= 2){
							dat.constraint_1[j] <- paste("+", paste(Edg.Inf[start:length(Edg.Inf)], collapse=" + "), sep=" ")
						} else {
							dat.constraint_1[j] <- paste(Edg.Inf[start:length(Edg.Inf)], collapse=" + ", sep=" ")
						}	
					} else {
						if(j >= 2){
							dat.constraint_1[j] <- paste("+", paste(Edg.Inf[start:stop], collapse=" + "), sep=" ")
						} else {
							dat.constraint_1[j] <- paste(Edg.Inf[start:stop], collapse=" + ", sep=" ")
						}
						a <-  stop
					}
				}
		
				dat.constraint_1[n] <- paste(dat.constraint_1[n],z.func, ">= 0")
				lpar <- dat.constraint_1
			}
	stopCluster(cl)
	return(lpar)
}


#### CONSTRAINT-1_WT
func.constraint_1_wt <- function(BiPartiteGraph, gamma, Influence.Outlier){
	g <- gamma
	Unique.Outliers <- sort(unique(BiPartiteGraph$OutlierEvent), decreasing=FALSE)
	
	#no_cores <- 30
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
		
	lpar <- foreach(i = 1:length(Unique.Outliers), .combine=c, .inorder=FALSE) %dopar% {	
				outlier <- Unique.Outliers[i]
				Sum.Inf <- Influence.Outlier$sumInfluence[which(Influence.Outlier$Outlier == outlier)]
				bip.sub <- subset(BiPartiteGraph, BiPartiteGraph$OutlierEvent == outlier)
				out.wt <- unique(bip.sub$OutlierWt)
				Sum.Inf.frac <- g * out.wt * Sum.Inf #### Cover g fraction

				z.func <- paste("-", Sum.Inf.frac, unique(bip.sub$Zeta), sep=" ")
				Edg.Inf <- paste(bip.sub$Influence, bip.sub$Edge, sep=" ")
		
				q <- length(Edg.Inf) %/% 4
				r <- length(Edg.Inf) %% 4

				if(r == 0){
					n <- q
				} else {
					n <- q + 1
				}	
	
				dat.constraint_1 <- rep(NA, n)
				a <- 0
				for(j in 1:n){
					start <- a + 1
					stop <- a + 4
		
					if(stop > length(Edg.Inf)){
						if(j >= 2){
							dat.constraint_1[j] <- paste("+", paste(Edg.Inf[start:length(Edg.Inf)], collapse=" + "), sep=" ")
						} else {
							dat.constraint_1[j] <- paste(Edg.Inf[start:length(Edg.Inf)], collapse=" + ", sep=" ")
						}	
					} else {
						if(j >= 2){
							dat.constraint_1[j] <- paste("+", paste(Edg.Inf[start:stop], collapse=" + "), sep=" ")
						} else {
							dat.constraint_1[j] <- paste(Edg.Inf[start:stop], collapse=" + ", sep=" ")
						}
						a <-  stop
					}
				}
		
				dat.constraint_1[n] <- paste(dat.constraint_1[n],z.func, ">= 0")
				lpar <- dat.constraint_1
			}
	stopCluster(cl)
	return(lpar)
}


#### CONSTRAINT-2
func.constraint_2 <- function(BiPartiteGraph){
	y_ij <- which(colnames(BiPartiteGraph) == "Edge")
	x_i <- which(colnames(BiPartiteGraph) == "MutGene")
	constraint_2 <- apply(BiPartiteGraph, 1, function(x) paste(paste(x[y_ij],paste("@", x[x_i], sep=""), sep=" - "), "= 0", sep=" "))
	return(constraint_2)
}


#### CONSTRAINT-3
func.constraint_3 <- function(BiPartiteGraph, alpha){
	Unique.zeta <- sort(unique(BiPartiteGraph$Zeta), decreasing=F)
	a.func <- prod(alpha, length(unique(BiPartiteGraph$OutlierEvent)))
	
	q <- length(Unique.zeta) %/% 10
	r <- length(Unique.zeta) %% 10
	
	if(r == 0){
		n <- q
	} else {
		n <- q + 1
	}
	
	constraint_3 <- rep(NA, n)
	a <- 0
	for(i in 1:n){	
		start <- a + 1
		stop <- a + 10
		
		if(stop > length(Unique.zeta)){
			if(i >= 2){
				constraint_3[i] <- paste("+", paste(Unique.zeta[start:length(Unique.zeta)], collapse=" + "), sep=" ")
			} else {
				constraint_3[i] <- paste(Unique.zeta[start:length(Unique.zeta)], collapse=" + ", sep=" ")
			}	
		} else {
			if(i >= 2){
				constraint_3[i] <- paste("+", paste(Unique.zeta[start:stop], collapse=" + "), sep=" ")
			} else {
				constraint_3[i] <- paste(Unique.zeta[start:stop], collapse=" + ", sep=" ")
			}
			a <-  stop
		}
	}
	constraint_3[n] <- paste(constraint_3[n], a.func, sep=" >= ")
	return(constraint_3)
}

#### CONSTRAINT-3 modified #### MODIFIED (same fraction of outlier covered in every patient)
func.constraint_3_modified <- function(BiPartiteGraph, beta){
	pat <- as.character(unique(BiPartiteGraph$Patient))
	
	constraint_3_modified <- list()
	
	for(ctr in 1:length(pat)){
		bip.temp <- subset(BiPartiteGraph, as.character(BiPartiteGraph$Patient) == pat[ctr])
		
		Unique.zeta <- sort(unique(bip.temp$Zeta), decreasing=F)
		b.func <- prod(beta, length(unique(bip.temp$OutlierEvent)))
		
		q <- length(Unique.zeta) %/% 10
		r <- length(Unique.zeta) %% 10
	
		if(r == 0){
			n <- q
		} else {
			n <- q + 1
		}
	
		constraint_3 <- rep(NA, n)
		a <- 0
		for(i in 1:n)
		{	
			start <- a + 1
			stop <- a + 10
		
			if(stop > length(Unique.zeta)){
				if(i >= 2){
					constraint_3[i] <- paste("+", paste(Unique.zeta[start:length(Unique.zeta)], collapse=" + "), sep=" ")
				} else {
					constraint_3[i] <- paste(Unique.zeta[start:length(Unique.zeta)], collapse=" + ", sep=" ")
				}	
			} else {
				if(i >= 2){
					constraint_3[i] <- paste("+", paste(Unique.zeta[start:stop], collapse=" + "), sep=" ")
				} else {
					constraint_3[i] <- paste(Unique.zeta[start:stop], collapse=" + ", sep=" ")
				}
				a <-  stop
			}
		}
		constraint_3[n] <- paste(constraint_3[n], b.func, sep=" >= ") ## original

		constraint_3_modified[[ctr]] <- constraint_3
	}
	constraint_3_modified <- unlist(constraint_3_modified)
	
	return(constraint_3_modified)
}


#### CONSTRAINT-3 BETA_LAMBDA #### 
func.constraint_3_beta_lambda_discard <- function(BiPartiteGraph, beta){
	pat <- unique(BiPartiteGraph$Patient)
	
	constraint_3_beta_lambda <- list()
	
	for(ctr in 1:length(pat)){
		bip.temp <- subset(BiPartiteGraph, BiPartiteGraph$Patient == pat[ctr])
		
		z_j <- which(colnames(bip.temp) == "Zeta")
		l_j <- which(colnames(bip.temp) == "OutlierWt")
		dat.outwt <- bip.temp[,c(z_j,l_j)]
		dat.outwt <- dat.outwt[!duplicated(dat.outwt),]
		dat.outwt <- dat.outwt[order(dat.outwt$OutlierWt, decreasing=T),]
		
		select.num <- round(beta * nrow(dat.outwt) , 0)
		select.zeta <- dat.outwt$Zeta[c(1:select.num)]
		
		q <- length(select.zeta) %/% 10
		r <- length(select.zeta) %% 10
	
		if(r == 0){
			n <- q
		} else {
			n <- q + 1
		}
	
		constraint_3 <- rep(NA, n)
		a <- 0
		for(i in 1:n){	
			start <- a + 1
			stop <- a + 10
		
			if(stop > length(select.zeta)){
				if(i >= 2){
					constraint_3[i] <- paste("+", paste(select.zeta[start:length(select.zeta)], collapse=" + "), sep=" ")
				} else {
					constraint_3[i] <- paste(select.zeta[start:length(select.zeta)], collapse=" + ", sep=" ")
				}	
			} else {
				if(i >= 2){
					constraint_3[i] <- paste("+", paste(select.zeta[start:stop], collapse=" + "), sep=" ")
				} else {
					constraint_3[i] <- paste(select.zeta[start:stop], collapse=" + ", sep=" ")
				}
				a <-  stop
			}
		}
		constraint_3[n] <- paste(constraint_3[n], length(select.zeta), sep=" >= ") ## original

		constraint_3_beta_lambda[[ctr]] <- constraint_3
	}
	constraint_3_beta_lambda <- unlist(constraint_3_beta_lambda)
	
	return(constraint_3_beta_lambda)
}

#### CONSTRAINT-3 BETA_LAMBDA #### 
func.constraint_3_beta_lambda <- function(BiPartiteGraph, beta){
	pat <- unique(BiPartiteGraph$Patient)

	#no_cores <- 30
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
		
	lpar <- foreach(ctr = 1:length(pat), .combine=c, .inorder=FALSE) %dopar% {
				bip.temp <- subset(BiPartiteGraph, BiPartiteGraph$Patient == pat[ctr])
		
				z_j <- which(colnames(bip.temp) == "Zeta")
				l_j <- which(colnames(bip.temp) == "OutlierWt")
				dat.outwt <- bip.temp[,c(z_j,l_j)]
				dat.outwt <- dat.outwt[!duplicated(dat.outwt),]
				dat.outwt <- dat.outwt[order(dat.outwt$OutlierWt, decreasing=T),]
		
				select.num <- round(beta * nrow(dat.outwt), 0)
				select.zeta <- dat.outwt$Zeta[c(1:select.num)]
		
				lpar <- paste(select.zeta, "1", sep=" = ")
	}
	stopCluster(cl)
	
	return(lpar)
}

#### CONSTRAINT-4
func.constraint_4 <- function(BiPartiteGraph){
	AbGenes <- sort(unique(BiPartiteGraph$MutGene), decreasing=F)
	constraint_4 <- paste("@", AbGenes, sep="")
	return(constraint_4)
}


#### CONSTRAINT-5
func.constraint_5 <- function(BiPartiteGraph){
	constraint_5 <- sort(unique(BiPartiteGraph$Zeta), decreasing=F)
	return(constraint_5)
}

#### CONSTRAINT-6
func.constraint_6 <- function(BiPartiteGraph){
	Unique.Edge <- sort(unique(BiPartiteGraph$Edge), decreasing=F)
	
	constraint_6 <- Unique.Edge
	return(constraint_6)
}

#### CONSTRAINT SEED ##FOR TESTING PERPOSE ONLY
func.constraint_seed <- function(seed.genes){
	seed.genes <- gsub("-","_", seed.genes)             # Convert DASH "-" symbol to UNDERSCORE "_" symbol
	constraint_seed <- paste(paste("@", seed.genes, sep=""), " = 1", sep="")
	return(constraint_seed)
}
