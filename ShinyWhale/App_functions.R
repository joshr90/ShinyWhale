whale_data <- function(Cap_H, month){
  #'
  #'---
  #'import data provided by the North Atlantic right whale Consortium
  whales  <- Cap_H 
  #'check that all unwanted information from the spreadsheet is removed
  str(whales)
  #'---
  #'remove duplicate SightingIDs
  whales.2<-whales[!duplicated(whales[, -1]), ]
  head(whales.2) ##Check that the duplicates have been removed
  #'---
  #'remove rows where SightingEGNO is NA
  whales.3<-whales.2[complete.cases(whales.2[ ,2]),]
  head(whales.3)
  #'---
  #'Set 1 Nov as the start of a year, so Nov or Dec (Year)==Year+1
  #'code logic: Year, if SightingMonth is greater than or equal to 11,
  #'then: Year is SightingYear+1; else: Year is SightingYear
  whales.4<-mutate(whales.3, Year=ifelse(SightingMonth>= month, yes=SightingYear+1, no=SightingYear))
  head(whales.4) ##Check the years are correct
  #'---
  #'Creates a column to indicate if the whale has been recovered dead
  whales.5 <- whales.4 %>% mutate(Dead = case_when(str_detect(Behaviors, "DEAD") ~ "DEAD",
                                                   TRUE ~ "ALIVE"))
  #'Concatenates the dead column by EGNO and Year in order to know that individual died when the data is collapsed 
  whales.5 <- whales.5 %>% 
    group_by(SightingEGNo, Year) %>% 
    mutate(DEAD = paste0(Dead, collapse = ",")) 
  #'---
  #'remove duplicate years when EGNO is the same
  whales.6<-distinct(whales.5,SightingEGNo,Year,.keep_all = TRUE)
  head(whales.6)

  #'---
  #'converts the DEAD column into a single value based upon if the individual was seen alive or dead
  whales.7 <- whales.6 %>% mutate(Dead = case_when(str_detect(DEAD, "DEAD") ~ "DEAD",
                                                   TRUE ~ "ALIVE"))
  #Sort the data by EGNo and SightingYear
  whales.8<-arrange(whales.7, SightingEGNo, Year)
  whales.8 <- whales.8 %>% mutate(Dead = case_when(str_detect(Dead, "DEAD") ~ 2,
                                               TRUE ~ 1))
  whales.9 <- whales.8 %>%
    group_by(SightingEGNo) %>%
    filter(!(Dead == 2 & n() > 1))
  return(whales.9)
}

ch <- function(whales.9) {
  mr <- acast(whales.9, SightingEGNo~Year, value.var="Dead")
  #'replace any NA's or 0's in the data with 3's, indicating the individual was neither seen nor recovered
  mr[is.na(mr)] = 3
  mr[mr == 0] <- 3
  #'Convert the data to numeric
  datachr <- as.numeric(mr)
  #'d is designed to find the dimensions of the data, number of individuals and the number of years, then we convert the numeric data into a matrix, 
  #'with the number of columns set to the number of years in the data set
  d <- dim(mr)
  datachr <- matrix(datachr, ncol = d[2])
  first_year <- rep(3, d[1])
  datachr <- cbind(first_year, datachr)
  return(datachr)
}


aug <- function(ch, da){
  d <- dim(ch)
  #'next we add pseudo-individuals for data augmentation. This is to account for individuals who may be in the population but never seen,
  #'for the NARW we add 300 additional individuals that can enter the population at any time over the study period, though this number can be changed depending on the length of the study.
  d1 <- d[2]  # this creates the right amount of columns for the number of years in the study period, plus the first year
  aug <- matrix(rep(3), nrow = da, ncol = d1) # this creates a matrix of 3's (not seen nor recovered) for the study period, with 300 individuals
  #'creates a csv of the data for record
  datachr <- as.matrix(aug)
  
  return(datachr)
}

jags_data <- function(ch) {
#'create the capture history matrix from the condensed NARW data set
#'CH - capture history
#'ch - the name of the data set created in "NARWDataInitialPrep.R"
#'ncol - number of years in the data 
CH <- ch
#'---
#'Compute date of first capture
get.first <- function(x) min(which(x!=3)) #a function to find the occasion that each individual was first seen 
get.last <- function(x) max(which(x!=3))  #a function to find the last occasion each individual was seen 
f <- apply(ch, 1, get.first) #creates a vector indicating the first column for each individual that is not equal to 3 (neither seen or recovered)
l <- apply(ch, 1, get.last) #creates a vector indicating the last column for each individual that is not equal to 3 (neither seen or recovered)
#'Bundle data - combines the data into a list that can be used by JAGS


#'Initial values
#'In order to run the model, we must provide sound initial values for the true latent state z. 
#'The difficulty is that observed states do not always correspond to the true latent state. 
#'For example, in our observation the state 3 refers to an individuals who has been reported as dead, 
#'while the true state 3 refers to an individuals that is alive, but outside the study area. 
#'Consequently, we cannot use the observed states as the initial values for the true state in JAGS. 
#'The important things are i) that the observations correspond to the true state (i.e. all "2" are converted into "4"), ii) that states after a "4" are all "5", 
#'iii) that all non-observations between "2's" become "2", and iv) that all remaining originally "3" after the first observation become "1".


z.data <- CH                    #duplicate the capture history to convert for the true state matrix
z.data[z.data==3] <- NA         #converts individuals who were not seen to NA's
z.data[z.data==2] <- 4          #changes observation 2 (recovered dead)  to state 4 (recovered dead) 
z.data[z.data==1] <- 2          #changes observation 1 (seen alive) to state 2 (alive in the study area)

for (i in 1:dim(z.data)[1]){
  z.data[i,f[i]:l[i]-1] <- 2       #for each individual, this loop makes all states between the first observation and last observation a 2 (alive in the study area)
  }
for (i in 1:dim(z.data)[1]){
  z.data[i,1:f[i]-1] <- 1       #for each individual, makes the state from the first column, until the individual was first observes a 1 (not yet entered)
}

NAs <- is.na(z.data)            #finds the occurrences of NA's in the matrix


z.data[NAs==T] <- 3             #converts NA's to 3's
z.data[z.data==3] <- 5          #converts 3 (not seen) to state 5 (dead) 

for (i in 1:dim(z.data)[1]){
z.data[i,1] <- NA                #converts the first time step to NA, as this state is defined in the model, and the model will not work if there is a number assigned to this occasion.
}
#'Next we need to define the initial values used in the model, for these we use uniform vague priors between 0-1, as as we are calculating them for each time step - 1, we need to provide priors for each.
#'for each mean.# below, set the first value of the runif as the number of columns in your datach matrix - 1, and the other two values represent the range of the priors (0-1). z.data is the true state matrix.
return(z.data)
}


inits <- function(n_yr, z.data){list(mean.s = runif(n_yr, 0, 1), mean.f = runif(n_yr, 0, 1), mean.p = runif(n_yr, 0, 1), mean.r = runif(n_yr, 0, 1), z = z.data)}  

