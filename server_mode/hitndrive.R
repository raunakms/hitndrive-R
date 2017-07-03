##########################################################
# SETTINGS #
###########################################################
dir.work <- "/Data/Raunak/HITnDRIVE"
dir.network <- file.path(dir.work, "InteractionNetwork")
dir.script <- file.path(dir.work,"Scripts/R/HITnDRIVE")
dir.dataset <- file.path(dir.work, "TCGAdata/GBM")
dir.analysis <- file.path(dir.dataset, "Analysis")
dir.db <- file.path(dir.work, "DriverGeneDB")

name.batch <- "TCGA_GBM_testb"

file.abr <- file.path(dir.analysis,"TCGA_GBM_AbrGenes.txt")
file.outlier <- file.path(dir.analysis,"TCGA_GBM_GE_outlier.txt")
file.network <- file.path(dir.network,"influenceMatrix_20000_influences.txt.gz")

file.cgc <- file.path(dir.db, "cancer_gene_census.tsv")
file.cosmic <- file.path(dir.db, "CosmicHGNC_v66_160713.tsv")

###########################################################

#### SET PARAMETERS ####
## alpha = fraction of the total outliers required to be covered
## gamma = fraction of the influence that the selected driver gene must at least have over the outlier
alpha <- c(0.8,0.9)
#alpha <- seq(0.1,1, by=0.1)
gamma <- 0.7
beta <- 0.2

########### HIT'nDRIVE TASK ###########
#                                     #
# TASK-1 : CONSTRUCT BI-PARTITE GRAPH #
# TASK-2 : CREATE DIRECTORY STRUCTURE #
# TASK-3 : GENERATE ILP FORMULATION   #
# TASK-4 : GENERATE CPLEX BASH SCRIPT #
# TASK-5 : RUN CPLEX                  #
# TASK-6 : PARSE CPLEX SOLUTION FILE  #
# TASK-7 : GENERATE STATISTICS        #
#                                     #
########### HIT'nDRIVE TASK ###########

# ENTER "T" FOR TRUE AND "F" FOR FALSE FOR EACH TASK
task.create_bipartite <- "T"
#task.create_dir <- "F"
task.generate_ilp <- "F"
task.generate_bash <- "F"
task.run_cplex <- "F"
task.prase_sol <- "F"
task.generate_stat <- "F"

##### START HIT'nDRIVE ####
source(file.path(dir.script, "main.R"))