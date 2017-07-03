get.ilp.bash.script <- function(dir.input, dir.output, par, name.batch, file.bash.script)
{
	cmd_all <- list()

	for(i in 1:nrow(par))
	{
		filename.ilp <- paste(name.batch,"_0",i,".lp.gz",sep="")
		filename.sol <- paste(name.batch,"_0",i,".sol",sep="")
	
		file.ilp <- file.path(dir.input,filename.ilp)
		file.sol <- file.path(dir.output,filename.sol)
	
		cmd_read <- paste("read",file.ilp, sep=" ")
		cmd_write <- paste("write",file.sol, sep=" ")
	
		cmd_file <- c(cmd_read, "optimize", cmd_write)
	
		cmd_all[[i]] <- cmd_file
	}

	cmd_all <- unlist(cmd_all)
	cmd_all <- c(cmd_all, "quit")

	write.table(cmd_all, file.bash.script, sep="\t", quote=F, row.names=F, col.names=F)

	cat("CPLEX BATCH SCRIPT GENERATED :", file.bash.script, "\n", sep="\t")
}