
library(argparse)

create_np_model_dir <- function(basedir, model_name, overwrite = FALSE){
  if(dir.exists(basedir)){
    model_dir <- file.path(basedir, model_name)
    if(dir.exists(model_dir) & !overwrite){
      stop('Model directory already exists. See help.')
    } else if(dir.exists(model_dir) & overwrite) {
      message('Overwriting model files in ', model_dir, '...')
      message('Note, this does not _delete_ any other files!')
    } else {
      dir.create(model_dir)  
    }
  } else {
    stop('Directory does not exist: ', basedir)
  }
  return(model_dir)
}

make_model_header <- function(y, x, within_prefix = 'WCEN_', between_prefix = 'GCEN_'){
  
  x_b_name <- paste0(between_prefix, x)
  x_w_name <- paste0(within_prefix, x)
  
  header <- '
library(nlme)
library(clubSandwich)

processVoxel <- function(v) {
  BRAIN <- voxeldat[,v]
  NOVAL <- 999
  retvals <- c()
  retnames <- c()'
  
  if(x == 'BRAIN'){
    header <- paste0(header,'
  brain_df <- data.frame(BRAIN, idnum) #idnum is defined by neuropoint
  
  # Compute Within-person Mean of brain data
  win_mean <- aggregate(brain_df, by = list(idnum), FUN = mean)[,1:2]
  colnames(win_mean) <- c("idnum", "WIN_BRAIN_MEAN")
  brain_df <- merge(brain_df, win_mean, by="idnum")

  grandmean <- mean(ag$win_mean)

  # This is a bit of an extra step... Tara had this here so I\'m leaving it
  brain_df[, "', x_b_name,'"] <- brain_df$WIN_BRAIN_MEAN - grandmean
  brain_df[, "', x_w_name,'"] <- BRAIN - brain_df$WIN_BRAIN_MEAN
  
  # attach to these new variables
  ', x_b_name,' <- brain_df[, "', x_b_name,'"]
  ', x_w_name,' <- brain_df[, "', x_w_name,'"]
') 
  }


  
  
  # Compute Grand-Mean Centered Variable
  grandmean <- mean(ag$WCEN_BRAIN_MEAN)
  dat$GCEN_BRAIN <- dat$WCEN_BRAIN_MEAN - grandmean
  dat$WCEN_BRAIN <- BRAIN - dat$WCEN_BRAIN_MEAN
  # attach to these new variables
  GCEN_BRAIN=dat$GCEN_BRAIN
  WCEN_BRAIN=dat$WCEN_BRAIN
  retvals <- c()
  retnames <- c()
  e <- try( DEP_BETWEEN_WITHIN <-lme(PHQ9_TOT ~ TIMECENTER + GCEN_BRAIN + GCEN_BRAIN * TIMECENTER +  WCEN_BRAIN + WCEN_BRAIN * TIMECENTER,random=~1|idnum,  method="ML", na.action=na.omit) )

}
make_lme_model_syntax <- function(y, x, within_prefix = 'WCEN_', between_prefix = 'GCEN_'){
  
}






parser <- ArgumentParser(description='Create a neuropoint model directory to be run')
parser$add_argument('base_dir', type="character",
                    help='This is where the model directory will be created')
parser$add_argument('model_name', type="character", help='Name of the model')
parser$add_argument('--IV', type="character",
                    help = 'IV variable name. Use "BRAIN" to specify BOLD contrast as the IV',
                    required = TRUE)
parser$add_argument('--DV', type="character",
                    help = 'DV variable name. Use "BRAIN" to specify BOLD contrast as the DV',
                    required = TRUE)
parser$add_argument('--overwrite', action='store_true', help = 'Overwrite existing model files?')
parser$add_argument('--win_pre', type="character", help='Within-person IV name prefix', default = 'WCEN_')
parser$add_argument('--bw_pre', type="character", help='Between-person IV name prefix', default = 'GCEN_')
parser$print_help()

args <- parser$parse_args(c('~/', 'newmodelthing', '--IV', 'EPISODICTOT', '--DV', 'BRAIN'))

if(! 'BRAIN' %in% c(args$DV, args$IV)){
  stop('Either --DV or --IV must be "BRAIN"... else why are you using neuropoint?')
}

model_dir <- create_np_model_dir(args$base_dir, args$model_name, overwrite = args$overwrite)

