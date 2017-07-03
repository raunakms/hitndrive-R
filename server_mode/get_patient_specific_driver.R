cat(paste(Sys.time()), "PATIENT SPECIFIC DRIVERS: COMPUTING ... ","\n", sep="\t")

library("stringr")

file.stats <- file.path(dir.batch, paste(name.batch, "_ilpsol_stats.txt", sep=""))
file.x <- file.path(dir.batch, paste(name.batch, "_drivers.txt", sep=""))
file.y <- file.path(dir.batch, paste(name.batch, "_edges.txt", sep=""))
file.z <- file.path(dir.batch, paste(name.batch, "_outliers.txt", sep=""))
	
###Load Data File
dat <-  read.table(file.stats, sep="\t", header=T, as.is="DriverGene")

### Load X,Y,Z data matrix
dat.x <- read.delim(file.x, header=T, row.names=1)
dat.y <- read.delim(file.y, header=T, row.names=1)
dat.z <- read.delim(file.z, header=T, row.names=1)

colnames(dat.x) <- paste(dat$gamma, dat$alpha, sep="_")
colnames(dat.y) <- paste(dat$gamma, dat$alpha, sep="_")
colnames(dat.z) <- paste(dat$gamma, dat$alpha, sep="_")

y.sp <- str_split(rownames(dat.y),";")
y.sp  <- do.call(rbind.data.frame, y.sp)
colnames(y.sp) <- c("Y","AbrGene", "Patient","OutlierGene")
y.sp <- y.sp[,-1]
dat.y <- cbind(dat.y, y.sp)

z.sp <- str_split(rownames(dat.z),";")
z.sp  <- do.call(rbind.data.frame, z.sp)
colnames(z.sp) <- c("Z","Patient","OutlierGene")
z.sp <- z.sp[,-1]
dat.z <- cbind(dat.z, z.sp)

### Load bipartiteartite Graph
bipartite <- read.delim(file.bipartite, header=T, as.is=c("MutGene","OutlierEvent"))

out.event <- str_split(bipartite$OutlierEvent, ";")
out.event <- do.call(rbind.data.frame, out.event)
colnames(out.event) <- c("Patient","OutlierGene")
bipartite$Patient <- out.event$Patient
bipartite$OutlierGene <- out.event$OutlierGene

samples <- as.character(unique(bipartite$Patient))

#### Compute Total Alteration per samples
total.abr <- data.frame(Patient=samples)
for(i in 1:length(samples)){
  bipartite.temp1 <- subset(bipartite, bipartite$Patient == samples[i])
  total.abr$TotalAbr[i] <- length(unique((bipartite.temp1$MutGene)))
}  


#### Compute Drivers Per samples per alpha
dat.list <- list()

count <- 1
for(i in 1:nrow(dat)){
	g <- dat$gamma[i]
	a <- dat$alpha[i]
	
	name.col <- paste(g, a, sep="_")
	sol.geneset <- dat$DriverGene[i]
	sol.genes <- str_split(sol.geneset, ",")[[1]]
  
	sol.outliers <- rownames(dat.z)[which(dat.z[,name.col] == 1)]
	sol.patients <- unique(as.character(dat.z$Patient[which(dat.z[,name.col] == 1)]))
  
  #z.sol <- do.call(rbind.data.frame, str_split(sol.outliers, ";"))
  #colnames(z.sol) <- c("Z","Patient","Outlier")
  #z.sol <- z.sol[,-1]
  
  for(j in 1:length(sol.patients)){
    sol.pat <- sol.patients[j]
    dat.temp1 <- subset(bipartite, bipartite$Patient == sol.pat)
    dat.temp2 <- subset(dat.temp1, dat.temp1$MutGene %in% sol.genes)
    driver.genes <- sort(unique(dat.temp2$MutGene), decreasing=F)
    
    dat.list[[count]] <- data.frame(gamma=dat$gamma[i],
									alpha=dat$alpha[i],
									beta=dat$beta[i],
                                    patient=sol.pat,
                                    driver.genes=paste(driver.genes, collapse=","),
                                    driver.num=length(driver.genes))
    count <- count+1
  }  
}  
dat.list <- do.call(rbind.data.frame, dat.list)

file.patient.driver.genes <- file.path(dir.batch, paste(name.batch, "_patientdrivers.txt", sep=""))
write.table(dat.list, file.patient.driver.genes, sep="\t", row.names=F, col.names=T, quote=F)

### Compute Drivers Per sample per alpha
driver.num <- matrix(0, nrow=length(samples), ncol=ncol(dat.x), dimnames=list(samples, colnames(dat.x)))
driver.val <- matrix("", nrow=length(samples), ncol=ncol(dat.x), dimnames=list(samples, colnames(dat.x)))

for(i in 1:nrow(dat.list)){
	name.col <- paste(dat.list$gamma[i], dat.list$alpha[i], sep="_")
	driver.num[as.character(dat.list$patient)[i],name.col] <- dat.list$driver.num[i]
	driver.val[as.character(dat.list$patient)[i],name.col] <- as.character(dat.list$driver.genes)[i]
}

driver.ratio <- round(driver.num/total.abr$TotalAbr,4)

file.driver.val <- file.path(dir.batch, paste(name.batch, "_driver_value.txt", sep=""))
file.driver.num <- file.path(dir.batch, paste(name.batch, "_driver_number.txt", sep=""))
file.driver.ratio <- file.path(dir.batch, paste(name.batch, "_driver_ratio.txt", sep=""))

write.table(driver.val,file.driver.val, sep="\t", row.names=T, col.names=NA, quote=F)
write.table(driver.num,file.driver.num, sep="\t", row.names=T, col.names=NA, quote=F)
write.table(driver.ratio, file.driver.ratio, sep="\t", row.names=T, col.names=NA, quote=F)


### Create Outlier Number matrix
outlier.num <- matrix(0, nrow=length(samples), ncol=ncol(dat.x), dimnames=list(samples, colnames(dat.x)))
outlier.ratio <- matrix(0, nrow=length(samples), ncol=ncol(dat.x), dimnames=list(samples, colnames(dat.x)))

for(i in 1:length(colnames(dat.x))){
  for(j in 1:length(samples)){
    pat <- samples[j]
    dat.temp <- subset(dat.z, as.character(dat.z$Patient) == pat)
    outlier.num[j,i] <- length(which(dat.temp[,i] == 1))
    outlier.ratio[j,i] <- length(which(dat.temp[,i] == 1))/nrow(dat.temp)
  }  
}  

file.outlier.num <- file.path(dir.batch, paste(name.batch, "_outlier_number.txt", sep=""))
file.outlier.ratio <- file.path(dir.batch, paste(name.batch, "_outlier_ratio.txt", sep=""))

write.table(outlier.num,file.outlier.num, sep="\t", row.names=T, col.names=NA, quote=F)
write.table(outlier.ratio, file.outlier.ratio, sep="\t", row.names=T, col.names=NA, quote=F)

cat(paste(Sys.time()), "PATIENT SPECIFIC DRIVERS: COMPLETE ... ","\n", sep="\t")
