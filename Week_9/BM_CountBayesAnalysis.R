#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for a single year of point count   ##
#  observations for Piute ground squirrels at the NCA and run a      ##
## closed population N-mixture analysis. The model is hierarchical    #
#  with : (1) an ecological submodel linking abundance to             #
## environmental predictors at each site; (2) an observation submodel #
## linking detection probability to relevant predictors.             ##
##  We rerun previous single season model from unmarked in JAGS Bayesian.##
# Abundance is expected to be higher in sites with more sagebrush     #
# and lower in those with more cheatgrass.                            #                                        #
# Detection may be related to observer effects and to time of day     #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse ) 
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI )
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Week_9/Data/", sep = "" )
# load cleaned data
closeddf <- read.csv( file = paste( datadir, "closed_counts.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# extract broad parameters of interest
#number of sites:
I <- length( closeddf$o.sites )
#number of visits
J = 3

#Create a new dataframe to scale predictors starting with site level predictors #
XI <- closeddf %>% 
  select( cheatgrass, sagebrush ) %>% 
  mutate( cheatgrass = scale( cheatgrass),
          sagebrush = scale( sagebrush) )

#for sitexvisit predictor we create separate objects for each 
#create index of columns of interest
tidx <- grep("time", colnames( closeddf), value = FALSE)
oidx <- grep( "obs", colnames(closeddf), value = FALSE)
#also for our observations
yidx <- grep( "count.j", colnames( closeddf ), value = FALSE)
#now use those indices to automatically select correct columns and scale
time_sc <- scale( closeddf[ ,tidx] )
#we want a quadratic term for time also
time2_sc <- scale( (closeddf[ ,tidx])^2 )
# for observers, we will include them as a random intercept instead
# of a fixed effect so we convert them to a number to represent index
obvs <- as.matrix( closeddf[,oidx] )
obvs <- as.numeric( as.factor(obvs) )
#turn back into matrix
obvs <- structure( obvs, dim = dim( closeddf[,oidx] ), class = "matrix" )
#number of observers 
S <- max(obvs, na.rm = TRUE)

#For homework you need to run only one of the three model versions below
# List your choice here including some reasoning for your choice:
# Answer: ORIGINALLY I CHOSE M3 - I'M STILL NEW TO BAYESIAN SO IT APPEARED THE 
# SIMPLEST FOR ME TO UNDERSTAND. I GOT AN ERROR "NO OBJECT WIN.DATA" 
# IT LOOKS LIKE WE DEFINE WIN.DATA IN M1 - NOT SURE IF WIN.DATA IS THE SAME
# FOR ALL MODELS BUT M1 RAN SUCCESSFULLY SO I WILL PROCEED WITH M1. 
#
 
##########################################################################
####### m1 single season N-mixture abundance model    #
# ecological predictors: cheatgrass and sagebrush #
# detection predictors: time and observer as random intercept #
############################################################################
############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #random intercept for observer
      for( s in 1:S ){
        eps.det[s] ~ dnorm( 0, pres.det ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.det <- 1/ ( sigma.det * sigma.det )
      #sigma prior specified as a student t half-normal:
      sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
      
      #priors for detection coefficients:
      #define as a slightly informative prior
      for( a in 1:A){
        alpha[ a ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #prior for abundance model intercept 
      int.lam ~ dgamma( 0.01, 0.01 ) T(-10, 10 )
      
    
    # ecological model of abundance
    for( i in 1:I ){
      #relative abundance modelled as a Poisson distribution 
      N[ i ] ~ dpois( lambda[ i ] )
      
      #linking abundance rate to predictor
      log( lambda[ i ] ) <- #intercept for abundance 
                          int.lam + 
                 # fixed effects of sagebrush and cheatgrass
                        inprod( beta, XI[ i, ]  )
    }
      #observation model
     for( i in 1:I ){  
      for( j in 1:J ){
        #model for probability of detection
        logit( p[i,j] ) <- #intercept for detection 
                          int.det + 
                  #random intercept for observer effect
                  eps.det[ obvs[i,j] ] +
                  #quadratic effect of time of day
                  alpha[1] * time[i,j] +
                  alpha[2] * time2[i,j]
                  
        #observed counts distributed as a Binomial:
        y_obs[ i, j ] ~ dbin( p[i,j], N[i] )  
                  
        #Model evaluation: We calculate Chi-squared discrepancy
        
        #start with expected abundance
        eval[i,j] <- p[i,j] * N[i]
        #compare vs observed counts
        E[i,j] <- pow( ( y_obs[i,j] - eval[i,j] ), 2 ) /
                  ( eval[i,j] + 0.001 )
        # Generate replicate data and compute fit stats
        #expected counts
        y_hat[ i, j ] ~ dbin( p[i,j], N[i] )
        #compare vs expected counts
        E.new[i,j] <- pow( ( y_hat[i,j] - eval[i,j] ), 2 ) /
                  ( eval[i,j] + 0.001 )
      
    } #close J
    } #close I
        
    #derived estimates of model fit
    fit <- sum( E[,] )
    fit.new <- sum( E.new[,] )
    
    } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################     
modelname <- "m1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'eps.det' #random intercepts in detection
              , 'sigma.det' #error for random intercept
              , 'beta' #abundance coefficients
              , 'p' #estimate of detection probability
              , 'y_hat' #predicted observations
              , 'N' #estimates of abundance
              , 'fit' #estimate of fit for observed data
              , 'fit.new' #estimate of fit for predicted data
)

#initial values defined as max counts
Nst <- apply( closeddf[ ,yidx], 1, max, na.rm = TRUE )
#replace 0s with 1
Nst[which(Nst== 0 )] <- 1

#how many ecological predictors that are fixed effects
B <- dim(XI)[2]
#how many detection predictors that are fixed effects
A <- 2
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           N = Nst ) }

#define data that will go in the model
str( win.data <- list( y_obs = as.matrix( closeddf[ ,yidx] ),
                       #number of sites, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, S = S,
                       #site level habitat predictors
                       XI = XI,
                       #observation predictors:
                       time = time_sc,
                       time2 = time2_sc,
                       obvs = obvs
) )

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                n.chains = 5, n.thin = 10, n.burnin = 20000,
                iter.increment = 10000, max.iter = 500000, 
                Rhat.limit = 1.05,
                save.all.iter = FALSE, parallel = TRUE ) 

#view results 
summary(m1)
plot(m1)

#chat
hist( m1$sims.list$fit / m1$sims.list$fit.new )
#mean chat
mean( m1$mean$fit ) / mean( m1$mean$fit.new )
#Bayesian pvalue
plot( m1$sims.list$fit, m1$sims.list$fit.new )
###### end m1 ########




############################################################################
################## Save your data and workspace ###################
# Save workspace:
save.image( "CountBayesResults.RData" )

######### End of saving section ##################################

############# END OF SCRIPT #####################################