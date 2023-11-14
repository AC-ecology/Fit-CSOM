# BirdNET outputs a csv for each hour (1200 3-sec clips) of recording 
# Raw data is in site folders - contains 1 file per hour recording, starting at 
# morning, with no lunch recordings = 4 files per day 

# FILES - Nrow = no. 3 secs in 1 hour (20) x no. species included in e-bird checklist (110)
# Ncol = start time (secs), sp name, confidence score, 

# For this example practice fit, I am using 5 sites (N = 5), 1 species (Wren),
# and 6720 clips per site (see below)
## Use data for 7 days only - whole data is 2400 at dawn, 2400 at dusk - 
# From this, select every 5th clip (480 clips at dawn, 480 at dusk) 
# equivalent to sampling 12 min of every hour (discont.)

## NOTE - in future, species specific subsample can be drawn knowing ecology of the species
### i.e. dawn or dusk, repetitiveness of call (thinner subsample) and commonness of call


# For CSOM - need a sample of 3-sec clips per day 
x <- list()
x$NSITES <- 6       # Number of sites
x$NFILES <- 6720     # Number of 3-sec clips per site

x$true_score <- matrix(nrow = x$NSITES, ncol = x$NFILES, data = NA) # vector of shape NSITES x NFILES containing a continuous score for each clip - the real-numbered score (logit of the 0-1 classifier score) produced by the ML classifier for each clip.

#####################################
#  Format data to get x$true_score  #
#####################################

# EXPLORING OUTSIDE OF FUNCTION BELOW

# Format for each site separately 
site <- "TM50"

wd <- setwd(paste("C:/Users/ajpc1/Desktop/FIT CSOM TENTSMUIR/Data/Raw",
                  "/", site, sep = ""))

hour_files <- list.files(wd)  # list of files in site folder

hour_dat <- read.csv(hour_files[1]) # first hour

hour_dat <- hour_dat[hour_dat$Common.name == "Eurasian Wren", ]

range(hour_dat$Confidence) # very low confidence - 0.01 median, 0.11 mean
# keeping 1 out of every 5 risks keeping only very low conf dets 
# surely we just care if a wren was present in the file?!

# TWO OPTIONS - keep every nth row, or take a stratified sample based on conf scores

# take a stratifed sample of 20% based on conf scores instead of time
library(splitstackshape)

# round confidence to 2dp
hour_dat$Confidence <- round(hour_dat$Confidence, 2)

set.seed(1)
strat_dat <- stratified(hour_dat, group = "Confidence", 0.2, replace = FALSE)
range(strat_dat$Confidence)

# how about taking every 5th sample?
new_dat <- hour_dat[seq(1, nrow(hour_dat), 5), ]
range(new_dat$Confidence)


# FUNCTION TO RUN 

# GENERATE FOR ONE SITE - we want to process the first 7 x 4 = 28 files for just 1 week of dets

# Function to generate x_dat to fit CSOM from hour length BN output - run for each
# site to generate a vector of scores 
# INPUTS:
##  Path - directory that houses raw csv files, organised into site folders
##  Site - name of site folder 
##  Species - focal species of exact format as BN output
##  N - to retain every Nth clip from the 1200 3-sec clips, default to keep 20%
##  days - days of recordings to process, where each day has 4 recordings 
scores <- NULL

site_xdat <- function(path = "C:/Users/ajpc1/Desktop/FIT CSOM TENTSMUIR/Data/Raw", 
                      site, species = "Eurasian Wren", N = 5, days = 7) {
  
  wd <- setwd(paste(path,
                    "/", site, sep = ""))
  
  hour_files <- list.files(wd) # list all files in the site folder
  
  n_hour_files <- 4 * days # specify how many hour files are to be processed 

  for(i in 1:n_hour_files) {
    
    # get first hour recording dat
    hour_dat <- read.csv(hour_files[i])
    
    # get full hour data for a single focal species 
    hour_dat <- hour_dat[hour_dat$Common.name == paste(species), ]
    # if(nrow(hour_dat) = 1200) stop("File does not contain 3600 3-sec clips")
    
    # take every Nth sample
    new_dat <- hour_dat[seq(1, nrow(hour_dat), N), ]
    if(nrow(new_dat) != (1200/N)) stop("Has not retained every Nth 3-sec clip")
    
    # retain the scores from new_dat 
    scores <- c(scores, new_dat$Confidence)
  }
  return(scores)
}

# Run in loop across sites, and store each site element in a list x$true_score

x$true_score <- matrix(nrow = x$NSITES, ncol = x$NFILES, data = NA)

sitenames <- list.files("C:/Users/ajpc1/Desktop/FIT CSOM TENTSMUIR/Data/Raw")

scores <- NULL

for(j in 1:5) {
  site <- sitenames[j]
  xdat <- site_xdat(site = site)
  x$true_score[j , ] <- xdat
}

x$true_score[1 , ] <- site_xdat(site = "TM50")
x$true_score[2 , ] <- site_xdat(site = "TM51")
x$true_score[3 , ] <- site_xdat(site = "TM52")
x$true_score[4 , ] <- site_xdat(site = "TM53")
x$true_score[5 , ] <- site_xdat(site = "TM54")
x$true_score[6 , ] <- site_xdat(site = "TM55")
# TO INVESTIGATE - Get error first time run each above line for any site other than TM50 then it runs fine 2nd time?

hist(x$true_score)
hist(nimble::logit(x$true_score))

# Add rest of data elements to list x
x$z_data <- rep(NA, x$NSITES) # vector of length NSITES containing either NA, 0, or 1 for each site) - the “naive” Z information for each site, where NA = sites that have not had a positive-annotated clip at them yet, and 1 = sites that have had a positive-annotated clip at them yet (i.e. a 1 for sites that are “confirmed occupied by annotator”)

x$annotation_all <- matrix(nrow = x$NSITES, ncol = x$NFILES, data = NA) # vector of shape NSITES x NFILES containing either NA, 0, or 1 for each clip) - the annotation status of each clip. NA = unannotated, 0 = annotated negative, 1 = annotated positive.
dim(x$annotation_all)

x$constraint <- 1

View(x)

# Save x as .Rdata object
setwd("C:/Users/ajpc1/Desktop/FIT CSOM TENTSMUIR/Data")

saveRDS(x, "x.rds")
