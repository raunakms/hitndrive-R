run.cplex <- function(file.bash.script){
	cat(paste(Sys.time()), "STARTING CPLEX TO SOLVE ILP ...","PLEASE WAIT ...","\n", sep="\t")
		#path.cplex <- "/Data/Raunak/CPLEX/cplex/bin/x86_sles10_4.1/cplex"
		path.cplex <- "/home/stas/Software/CPLEX/cplex/bin/x86-64_linux/cplex"
	
		cmd <- paste(path.cplex, "<", file.bash.script, ">", file.cplex.log, sep=" ")
		system(cmd)
	cat(paste(Sys.time()), "TERMINATING CPLEX ...","\n", sep="\t")
}
