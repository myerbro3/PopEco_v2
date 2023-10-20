##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

#install relevant packages
install.packages( 'Rtools' )
install.packages( "nmixgof" ) #for evaluating N-mixture models
#load relevant packages
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #runs N-mixture models
library( MuMIn ) #calculates pseudo-R^2
library( AICcmodavg) #gof tests (Duarte et al. 2018)
library( nmixgof ) #more gof tests (Knape et al. 2018)
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Week_8/Data/", sep = "" )
# load cleaned data
closeddf <- read.csv( file = paste( datadir, "closed_counts.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (counts):
umf <- unmarkedFramePCount( 
  y = as.matrix( closeddf[ ,c("count.j1", "count.j2", "count.j3")]),
  # Define predictors at the site level:
  siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
  # Define predictors at the survey level as a list:
  obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")],
                  time = closeddf[ ,c('time.j1', 'time.j2', 'time.j3') ],
                  time2 = cbind( (closeddf$time.j1)^2, (closeddf$time.j2)^2, 
                                 (closeddf$time.j3)^2 ) ) ) 

# View summary of unmarked dataframe:
summary( umf )
#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
#now for time:
timesc <- as.vector(scale( obsCovs(umf)[2] ))
#replace with scaled values:
obsCovs(umf)[2] <- timesc
time2sc <- as.vector(scale( obsCovs(umf)[3] ))
#replace with scaled values:
obsCovs(umf)[3] <- time2sc

#recheck
summary( umf )
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. We start with a full model:
fm.closed <- pcount( ~ 1 + obsv + time + time2
                     ~ 1 + sagebrush + cheatgrass, 
                     #Define the maximum possible abundance
                     #during the primary occasion
                     K = 40,
                     data = umf, 
                     mixture="NB")
# Note that we start with the observation submodel #
#We then define the ecological submodel 
# You should try alternative K values to make sure that your model 
#isn't sensitive to the value you provided. 

# View model results:
fm.closed

# Estimate confidence intervals:
confint( fm.closed, type = "state" )
# coefficients for detection submodel:
confint( fm.closed, type = 'det' )
#

#
#What is the mean abundance from our model?
exp(coef(fm.closed[1]))


#############end full model ###########
##########################################################################
# Model fit and evaluation -----------------------------------------------
# We start with goodness of fit (GoF) outlined by Duarte et al. 2018 #
# Ecological modelling 374:51-59 and available via AICmodavg package #
# The Nmix.gof.test relies on a Pearson chi-square to assess the fit of #
# N-mixture models. The approach uses bootstrapping to estimate the p values #
# The test also estimates a c-hat measure of overdispersion, as the  #
# observed test statistic divided by the mean of the simulated test statistics #

# Let's compute observed chi-square, assess significance, and estimate c-hat
gof.boot <- Nmix.gof.test( fm.closed, nsim = 100, 
                           print.table = TRUE )
#view
gof.boot

# We run the null model
fm.null <- pcount( ~ 1 ~ 1,
                   K = 40, data = umf )
#view
fm.null
# Now build the list with the two models of interest:
rms <- fitList( 'full' = fm.closed,
                'null' = fm.null )
# Then use model selection function from unmarked, defining which is the null:
unmarked::modSel(rms, nullmod = "null" )

#
# Now use gof checks outlined in Knape et al. 2018 MEE 9:2102-2114
# We start by estimating overdispersion metrics 
chat( fm.closed, type = 'marginal' )
chat( fm.closed, type = 'site-sum' )
chat( fm.closed, type = 'observation' )

# Plot rq residuals against untransformed numeric predictors. This may help
# detect problems with the functional form assumed in a model
residcov( fm.closed )


# Plot residuals against fitted values. Site-sum randomized quantile residuals
# are used for site covariates while marginal residuals are used for
# observation covariates. 
residfit( fm.closed, type = 'site-sum' )
# Plot the observation model residuals
residfit( fm.closed, type = 'observation' )

# Nw plot rq plots of randomized residuals against standard normal quantiles. #
# Under a good fit, residuals should be close to the identity line. 
residqq( fm.closed, type = 'site-sum' )
residqq( fm.closed, type = 'observation' )



############# REST OF ASSIGNMENT CONTINUED IN BM_CountAnalysis #############