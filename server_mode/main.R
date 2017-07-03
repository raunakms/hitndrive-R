cat("# ------------- HIT'nDRIVE -------------- #", "\n", sep="\t")
cat("# CONSTRUCT BI-PARTITE GRAPH :", task.create_bipartite, "\n", sep="\t")
cat("# GENERATE ILP FORMULATION   :", task.generate_ilp, "\n", sep="\t")
cat("# GENERATE CPLEX BASH SCRIPT :", task.generate_bash, "\n", sep="\t")
cat("# RUN CPLEX                  :", task.run_cplex, "\n", sep="\t")
cat("# PARSE CPLEX SOLUTION FILE  :", task.prase_sol, "\n", sep="\t")
cat("# GENERATE STATISTICS        :", task.generate_stat, "\n", sep="\t")
cat("# --------------------------------------- #", "\n", sep="\t")
cat("\n")

### TASK CREATE DIRECTORY STRUCTURE
source(file.path(dir.script, "get_dir.R"))

#### CALL BI-PARTITE
if(task.create_bipartite == "T"){
	cat(paste(Sys.time()), "STARTING TASK: CREATE BI-PARTITE GRAPH ...","\n", sep="\t")
		source(file.path(dir.script, "call_bipartite_graph.R"))
	cat(paste(Sys.time()), "ENDING TASK: CREATE BI-PARTITE GRAPH ...","\n", sep="\t")	
}

if(task.generate_ilp == "T" | task.generate_bash == "T" | task.run_cplex == "T" | task.prase_sol == "T" | task.generate_stat == "T"){
	### TASK CREATE DIRECTORY STRUCTURE
	#if(task.create_dir == "T")
	#{
	#	source(file.path(dir.script, "get_dir.R"))
	#}
	
	### LOAD REQUIRED FILES
	cat(paste(Sys.time()), "REQUIRED FUNCTIONS : LOADING ...","\n", sep="\t")
		source(file.path(dir.script, "call_ilp.R"))
		source(file.path(dir.script, "get_ilp.R"))
		source(file.path(dir.script, "get_ilp_bash_script.R"))
		source(file.path(dir.script, "run_cplex.R"))
		source(file.path(dir.script, "get_lp_sol_stats.R"))
	cat(paste(Sys.time()), "REQUIRED FUNCTIONS : LOADED ...","\n", sep="\t")
	cat("\n")

	file.par <- file.path(dir.batch, paste(name.batch, "_par.txt", sep=""))
	file.bash.script <- file.path(dir.batch, paste(name.batch, "_batch_lp.sh", sep=""))
	file.cplex.log <- file.path(dir.batch, paste(name.batch, "_batch_lp.log", sep=""))
	
	### CREATE PARAMETER FILES
	cat(paste(Sys.time()), "PARAMETER FILE : CREATING ...","\n", sep="\t")
		par <- expand.grid(gamma, alpha, beta)
		colnames(par) <- c("gamma","alpha","beta")
		par <- par[order(par$gamma, par$beta, decreasing=F),]
		par <- par[which(apply(par, 1, function(x) ifelse(x[2] >= x[3], 1, 0)) == 1),]
		
		if(nrow(par) == 0){
			stop(paste(Sys.time(), "We expect alpha >= beta. Please correct simulation values.", sep="\t"))
		}
		
		par$iter <- c(1:nrow(par))
		print(par)
		write.table(par, file.par, sep="\t", row.names=F, col.names=T, quote=F)
	cat(paste(Sys.time()), "PARAMETER FILE : CREATED ...","\n", sep="\t")
	cat("\n")

	cat(paste(Sys.time()), "DIRECTORY DATASET :", dir.dataset, "\n", sep=" ")
	cat(paste(Sys.time()), "DIRECTORY BATCH:", dir.batch, "\n", sep=" ")
	cat(paste(Sys.time()), "FILE PARAMETER:", file.par, "\n", sep=" ")
	cat(paste(Sys.time()), "FILE BI-PARTITE:", file.bipartite, "\n", sep=" ")
	cat("\n")

	#### CALL ILP
	if(task.generate_ilp == "T"){
		if(!exists("seed.genes")){
			seed.genes <- NA
		}

		list.constraint <- call.ilp(file.bipartite, par, beta, dir.input, name.batch, file.driverWt, seed.genes)
	}

	#### GENERATE lp BASH SCRIPT
	if(task.generate_bash == "T"){
		get.ilp.bash.script(dir.input, dir.output, par, name.batch, file.bash.script)
	}

	#### RUN CPLEX
	if(task.run_cplex == "T"){
		run.cplex(file.bash.script)
	}

	### PARSE LP SOL FILES
	if(task.prase_sol == "T"){
		cmd <- paste("perl", file.path(dir.script, "parse_lp_sol.pl"), dir.output, dir.sol.parse, name.batch, sep=" ")
		system(cmd)
	}

	### COMPUTE SOL STATS
	if(task.generate_stat == "T"){
		get_lp_sol_stats(dir.sol.parse, dir.batch, name.batch, file.bipartite)
	}
}
