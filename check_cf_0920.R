############input RDS file generated in each cycle;input the serial number of the RDS file (the number followed by '.RDS');#############
##########four files would be generated (direct file from RDS file;raw crossfeeding file; crossfeeding file with annotation of compounds ID;output and input of each individual)######
library(lattice)
library(ReacTran)
library(rootSolve)
library(deSolve)
library(shape)
library(sybil)
library(Matrix)
library(BacArena)
library(parallel)
library(plyr)
option_list = list(
  optparse::make_option(c("-o", "--outfile_name"),                 type="character",      default='opt_ex',       help="part of output file name"),
  optparse::make_option(c("-i","--infile"),type="character",help="input file should be RDS file"),
  optparse::make_option(c("-t", "--iter"),                  type="character",         help="iteration"));

opt_parser = optparse::OptionParser(option_list=option_list, add_help_option=FALSE);
opt = optparse::parse_args(opt_parser);
infile     = opt$infile
outfile_name=opt$outfile_name
iter=opt$iter
# ############################ Download RDS file and run on MacOS ############################
getwd <- getwd()
getwd
simulation_loop <- readRDS(paste(getwd,'/',infile,'.RDS',sep = '')) ####BacArena_community_0806_none_400grids_mineral_sw_4.RDS

# 10.1. exchangeslist: A list of containing exchanges per time step (unit: mmol/(h x g_dw)?).
# By calculating the sum of each column, you will get the net flux of the compound in the current iteration.
tp = 2
df_exchangeslist <- as.data.frame(simulation_loop[[1]]@exchangeslist[[tp]]) #A list of containing exchanges per time step.
write.csv(df_exchangeslist, file = paste(getwd,'/',outfile_name,"_",iter,'.csv',sep = ''), row.names = F)
## This output file includes all metabolite exchanges of each individual "cell" in the area.

## Then we need to use the output file as the input file in this step, to get the metabolite exchange for each species.
# St1. Run "python3 /Get_network_graph_Willis.py -indir -infile"
cmd0="module load python/3.7.3"
system(command = cmd0)
cmd1 = "python /srv/scratch/z5245780/scripts_and_shell_command/Bacarena/crossfeeding_check_0919.py"
mydir = paste(getwd,sep = '')
infile = paste('/',outfile_name,"_",iter,'.csv',sep = '')
system(command = paste(cmd1," -mydir ", mydir," -infile ", infile, sep = ''))