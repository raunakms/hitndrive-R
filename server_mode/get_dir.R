dir.main <- file.path(dir.dataset, "hitndrive_files")
dir.batch <- file.path(dir.main,name.batch)
dir.input <- file.path(dir.batch, "ilp_input")
dir.output <- file.path(dir.batch,"ilp_output")
dir.sol.parse <- file.path(dir.batch, "ilp_sol_parse")

file.bip <- file.path(dir.batch, paste(name.batch, "_BiPartiteGraph", ".txt", sep=""))
file.bipartite <- file.path(dir.batch, paste(name.batch, "_BiPartiteGraph", ".txt.gz", sep=""))

cat(paste(Sys.time()), "DIRECTORY STRUCTURE : CREATING ...","\n", sep="\t")
	#### CREATE DIRECTORIES ####
	dir.create(dir.main, showWarnings=FALSE)
	dir.create(dir.batch, showWarnings=FALSE)
	dir.create(dir.input, showWarnings=FALSE)
	dir.create(dir.output, showWarnings=FALSE)
	dir.create(dir.sol.parse, showWarnings=FALSE)
cat(paste(Sys.time()), "DIRECTORY STRUCTURE : CREATED ...","\n", sep="\t")

