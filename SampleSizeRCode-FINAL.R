
#'########################################################
#'           BIOL 5337P – MSc. RESEARCH PROJECT.         #
#'                    JUMA MINYA © 2022                  #
#'########################################################


#' code to simulate age/sex distribution based on 2014 elephant aerial census + 
#' age distribution data from Jones et, al (2018)


#Set working directory
setwd('C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification')

# Fixing the codes to run the same way every time
set.seed(104)

#Loading the library to manipulate and plot the data
library (dplyr)# Manipulating the data
library (ggplot2) #Plotting the data
library (ggpubr) #Arranging the graphs into panels
library(grid)# Manipulating the labels


# Age classes distribution (years): 0-4, 5-9, 10-14, 15-19, 20-24, 25-39, 40+
fageclass <- c('0','5','10','15','20','25','40') # Female group age classes
mageclass <- c('15','20','25','40')# Male group age classes

########################################################
#####           1. SERENGETI POPULATION            #####
########################################################

# Loading the Serengeti elephant survey data and assume similar trends with other populations
sere.grp <- read.csv('ele2014.csv') 

#Preview of the data
head(sere.grp)

# Finding the group sizes
fem <- length(sere.grp$Corrected[sere.grp$Species=='ELF'])
bul <- length(sere.grp$Corrected[sere.grp$Species=='ELB'])

# Calculating the proportions of the sexes relative to the population
fem.prop <- fem/(fem+bul) #for females
bul.prop <- bul/(fem+bul)#for males

# Fit the distribution of the groups
# females
nbinom <- fitdistrplus::fitdist(data=sere.grp$Corrected[sere.grp$Species=='ELF'],distr = 'nbinom')
summary(nbinom)

# males
nbinom <- fitdistrplus::fitdist(data=sere.grp$Corrected[sere.grp$Species=='ELB'],distr = 'nbinom')
summary(nbinom)

# looking at distribution of the groups
hist(sere.grp$Corrected[sere.grp$Species=='ELF'],probability = T,breaks = 10,col='cornflowerblue',ylim=c(0,0.2))
hist(sere.grp$Corrected[sere.grp$Species=='ELB'],probability = T,breaks = 10,add=T,col='darkolivegreen')

# Adopting the age classes (years): 0–4 5–9 10–14 15–19 20–24 25–39 40+
Sfaead <- c(72,30,24,12,15,43,17) ### female age Distribution - from Jones data in Jones et al. 2018
faeadS <- Sfaead/sum(Sfaead) #make female age distribution proportional to total
Smaead <- c(71,30,18,6,7,17,2) ### male age distribution - from Jones data data
maeadS <- Smaead/sum(Smaead) #make male age proportional to total

# input parameters
S.popsize <- 6100 # approx pop size of Serengeti elephants (reported 6087)
samplesize <- c(0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.9) * S.popsize # create a sequence sample sizes ranging from 1% to 90% of the total pop
nsim <- 1000 # number of times to repeat the simulation for each combination of values

# function to calc where to stop sampling (when running out)
stopfun <- function(vTR,element) {
  val = cumsum(vTR)
  return(length(val[val<=element])+1)# return the next element in the length
}
  
  fin <- data.frame()# create an empty dataframe where to write in simulations
  
  # start loop across sample sizes
  for(k in samplesize){
    
    # loop across iterations
    for(j in 1:nsim){  
      
      # create empty dataframe
      df <- data.frame()
      
      ##### step 1 -- Sex and group size sampling
      # sample sex of groups
      sex <- rbinom(n = k,size=1,prob = fem.prop) # 0=male group 1=female group
      
      # sample group sizes for each sex
      grp <- rep(NA,k)
      grp[sex==0] <- rnbinom(n = length(sex[sex==0]),size = 2.247386,mu = 2.339081) #Values from the distribution summary of Serengeti male group
      grp[sex==1] <- rnbinom(n = length(sex[sex==1]),size = 1.703228,mu = 13.517563) #Values from the distribution summary of Serengeti female group
      
      # remove any groups with size=0
      sex <- sex[grp>0]
      grp <- grp[grp>0]
      
      # reduce group size vector
      sex <- sex[1:stopfun(grp,k)]
      grp <- grp[1:stopfun(grp,k)]
      
      # change sex of groups to be M and F
      sex[sex==0] <- 'M'
      sex[sex==1] <- 'F'

      ### step 2 randomly sample age distribution
      for(i in 1:length(grp)){
        sex.grp <- rep(sex[i],grp[i])
        
        if(sex[i]=='F'){
          agedist <- sample(x = fageclass,size = grp[i],replace=T,prob=faeadS)
          sex.grp[agedist %in% c("0","5","10")] <- sample(c('M','F'),
                                                          size = sum(agedist %in% c("0","5","10")),
                                                          prob = c(0.5,0.5),replace=T) #assuming equal sex ratio 
        }
        if(sex[i]=='M'){
          
          agedist <- sample(x = mageclass,size = grp[i],replace=T,prob=maeadS[4:7]) 
        }
        
        # compile all data into a single dataframe
        df <- rbind(df,
                    data.frame(agedist=agedist,
                               sex=sex.grp))
        
      } # end i loop - sampling agedist
      
      ## compare samples to 'truth' for k samples
      df <- df[1:k,]
      Females <- df[df$sex=='F',] 
      Females$agedist <- factor(Females$agedist,levels=fageclass)
      
      Males <- df[df$sex=='M',] 
      Males$agedist <- factor(Males$agedist,levels=fageclass) # To incorporate all age categories
      
      S.obs.age.distf <- table(Females$agedist)/nrow(Females)
      S.obs.age.distm <- table(Males$agedist)/nrow(Males)
      
      #check if the ages for females and males are in the correct order
      if(!all(names(table(Females$agedist)) == fageclass)){
        S.obs.age.distf <- S.obs.age.distf[match(fageclass,names(table(Females$agedist)))]
      }
      
      if(!all(names(table(Males$agedist)) == fageclass)){
        S.obs.age.distm <- S.obs.age.distm[match(fageclass,names(table(Males$agedist)))]
      }
      
      # calculate the proportional difference (Bias)
      bfS <- mean((S.obs.age.distf - faeadS)/faeadS,na.rm=T) #females
      bmS <- mean((S.obs.age.distm - maeadS)/maeadS,na.rm=T) #males
      
      # calculate the root mean square error (Precision)-
      rmsefS <- sqrt(mean((S.obs.age.distf - faeadS)^2,na.rm=T)) #females
      rmsemS <- sqrt(mean((S.obs.age.distm - maeadS)^2,na.rm=T)) #males
      
      # extract the information on age and sex and combine it in a dataframe fin
      fin <- rbind(fin,
                   data.frame(           
                     age0f = S.obs.age.distf[[1]],
                     age5f = S.obs.age.distf[[2]],
                     age10f = S.obs.age.distf[[3]],
                     age15f = S.obs.age.distf[[4]],
                     age20f = S.obs.age.distf[[5]],
                     age25f = S.obs.age.distf[[6]],
                     age40f = S.obs.age.distf[[7]],
                     age0m = S.obs.age.distm[[1]],
                     age5m = S.obs.age.distm[[2]],
                     age10m = S.obs.age.distm[[3]],
                     age15m = S.obs.age.distm[[4]],
                     age20m = S.obs.age.distm[[5]],
                     age25m = S.obs.age.distm[[6]],
                     age40m = S.obs.age.distm[[7]],
                     bfS = bfS,
                     bmS = bmS,
                     rmsefS = rmsefS,
                     rmsemS = rmsemS,
                     samplesize=k, nsim=j))                 
      
    } # end j loop  - # iterations
  }  # end k loop - sample size, i.e. proportion of the population
  
#Save the resultant data frame
write.csv(fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/S.Finalsimul.csv", row.names = TRUE)

# read simulated data
S.Finalsimul <- read.csv('S.Finalsimul.csv')

#Grouping and summarizing the data
S.sum.fin<- S.Finalsimul %>% 
  group_by(samplesize) %>% 
  summarise(bfS_mean = mean(bfS,na.rm=T), bmS_mean = mean(bmS,na.rm=T),
            rmsefS_mean = mean(rmsefS,na.rm=T),rmsemS_mean = mean(rmsemS,na.rm=T),
            bfS_sd = sd(bfS,na.rm=T),bmS_sd = sd(bmS,na.rm=T),
            rmsefS_sd = sd(rmsefS,na.rm=T),rmsemS_sd = sd(rmsemS,na.rm=T))

#Save the resultant data frame
write.csv(S.sum.fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/sum.SFinsim.csv", row.names = TRUE)

# read the summary of the final simulated data
sum.SFinsim <- read.csv('sum.SFinsim.csv')

#Combined plot
#Variation of RMSE Sample size
Serengeti <- ggplot(sum.SFinsim,)+ 
  geom_line(aes(x=samplesize,y=rmsemS_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=rmsefS_mean,linetype = "Females"))+
  labs(x="Sample size", y="RMSE")+
  theme_classic() #white background no grid lines

#Variation of Bias with sample size
SerengetiB <- ggplot(sum.SFinsim)+ 
  geom_line(aes(x=samplesize,y=bmS_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=bfS_mean,linetype = "Females"))+
  labs(x="Sample size", y="Bias")+
  theme_classic() #white background no grid lines


########################################################
#####           2. RUAHA POPULATION                #####
########################################################

# read Ruaha elephant survey data from Jones study
Ruaha <- read.csv('Ruaha.csv') 

#Preview of the population
head(Ruaha)

# look at distribution of group sizes
hist(log(Ruaha$Total))

# For Ruaha population
R.faead <- c(60,24,25,9,16,46,4) ### Female age distribution for Ruaha ecosystem
faeadR <- R.faead/sum(R.faead) #make proportional of the female to total population
R.maead <- c(60,23,23,13,14,10,2) ### Male age distribution for Ruaha ecosystem
maeadR <- R.maead/sum(R.maead) #make proportional of the male to total population

# input parameters
R.popsize <- 14400 # approx pop size of Ruaha elephants (reported 14384)
samplesize <- c(0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.9) * R.popsize # create a sequence sample sizes ranging from 1% to 90% of the total pop
nsim <- 1000 # number of times to repeat the simulation for each combination of values

# function to calc where to stop sampling
stopfun <- function(vTR,element) {
  val = cumsum(vTR)
  return(length(val[val<=element])+1)# return the next element in the length
}

fin <- data.frame()

# start loop across sample sizes
for(k in samplesize){
  
  # loop across iterations
  for(j in 1:nsim){  
    
    # create empty dataframe
    df <- data.frame()
    
    ##### step 1 -- Sex and group size sampling
    # sample sex of groups
    sex <- rbinom(n = k,size=1,prob = fem.prop) # 0=male group 1=female group
    
    # sample group sizes for each sex
    grp <- rep(NA,k)
    grp[sex==0] <- rnbinom(n = length(sex[sex==0]),size = 2.247386,mu = 2.339081) 
    grp[sex==1] <- rnbinom(n = length(sex[sex==1]),size = 1.703228,mu = 13.517563) 
    
    # remove any groups with size=0
    sex <- sex[grp>0]
    grp <- grp[grp>0]
    
    # reduce group size vector
    sex <- sex[1:stopfun(grp,k)]
    grp <- grp[1:stopfun(grp,k)]
    
    # change sex of groups to be M and F
    sex[sex==0] <- 'M'
    sex[sex==1] <- 'F'
    
    ### step 2 randomly sample age distribution
    for(i in 1:length(grp)){
      sex.grp <- rep(sex[i],grp[i])
      
      if(sex[i]=='F'){
        agedist <- sample(x = fageclass,size = grp[i],replace=T,prob=faeadR)
        sex.grp[agedist %in% c("0","5","10")] <- sample(c('M','F'),
                                                        size = sum(agedist %in% c("0","5","10")),
                                                        prob = c(0.5,0.5),replace=T)  
      }
      if(sex[i]=='M'){
        
        agedist <- sample(x = mageclass,size = grp[i],replace=T,prob=maeadR[4:7]) 
      }
      
      # compile all data into a single dataframe
      df <- rbind(df,
                  data.frame(agedist=agedist,
                             sex=sex.grp))
      
    } # end i loop - sampling agedist
    
    ## compare samples to 'truth' for k samples
    df <- df[1:k,]
    Females <- df[df$sex=='F',] 
    Females$agedist <- factor(Females$agedist,levels=fageclass)
    
    Males <- df[df$sex=='M',] 
    Males$agedist <- factor(Males$agedist,levels=fageclass) # To incorporate all age categories
    
    R.obs.age.distf <- table(Females$agedist)/nrow(Females)
    R.obs.age.distm <- table(Males$agedist)/nrow(Males)
    
    #check if the ages for females and males are in the correct order
    if(!all(names(table(Females$agedist)) == fageclass)){
      R.obs.age.distf <- R.obs.age.distf[match(fageclass,names(table(Females$agedist)))]
    }
    
    if(!all(names(table(Males$agedist)) == fageclass)){
      R.obs.age.distm <- R.obs.age.distm[match(fageclass,names(table(Males$agedist)))]
    }
    
    # calculate the proportional difference (Bias)
    bfR <- mean((R.obs.age.distf - faeadR)/faeadR,na.rm=T)
    bmR <- mean((R.obs.age.distm - maeadR)/maeadR,na.rm=T)
    
    # calculate the root mean square error (Precision)-
    rmsefR <- sqrt(mean((R.obs.age.distf - faeadR)^2))
    rmsemR <- sqrt(mean((R.obs.age.distm - maeadR)^2))
    
    fin <- rbind(fin,
                 data.frame(           
                   age0f = R.obs.age.distf[[1]],
                   age5f = R.obs.age.distf[[2]],
                   age10f = R.obs.age.distf[[3]],
                   age15f = R.obs.age.distf[[4]],
                   age20f = R.obs.age.distf[[5]],
                   age25f = R.obs.age.distf[[6]],
                   age40f = R.obs.age.distf[[7]],
                   age0m = R.obs.age.distm[[1]],
                   age5m = R.obs.age.distm[[2]],
                   age10m = R.obs.age.distm[[3]],
                   age15m = R.obs.age.distm[[4]],
                   age20m = R.obs.age.distm[[5]],
                   age25m = R.obs.age.distm[[6]],
                   age40m = R.obs.age.distm[[7]],
                   bfR = bfR,
                   bmR = bmR,
                   rmsefR = rmsefR,
                   rmsemR = rmsemR,
                   samplesize=k, nsim=j))                 
    
  } # end j loop  - # iterations
}  # end k loop - sample size, i.e. proportion of the population

#Save the resultant data frame
write.csv(fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/R.Finalsimul.csv", row.names = TRUE)

# read simulated data
R.Finalsimul <- read.csv('R.Finalsimul.csv')

#Grouping and summarizing the data
sum.fin <- R.Finalsimul %>% 
  group_by(samplesize) %>% 
  summarise(bfR_mean = mean(bfR), bmR_mean = mean(bmR),
            rmsefR_mean = mean(rmsefR),rmsemR_mean = mean(rmsemR),
            bfR_sd = sd(bfR),bmR_sd = sd(bmR),
            rmsefR_sd = sd(rmsefR),rmsemR_sd = sd(rmsemR))

#Save the resultant data frame
write.csv(sum.fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/sum.RFinsim.csv", row.names = TRUE)

# read the summary of the final simulated data
sum.RFinsim <- read.csv('sum.RFinsim.csv')

#Combined plot
#Variation of RMSE with Sample size
Ruaha <- ggplot(sum.RFinsim,)+
  geom_line(aes(x=samplesize,y=rmsemR_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=rmsefR_mean,linetype = "Females"))+
  labs(x="Sample size", y="RMSE")+ 
  theme_classic() #white background no grid lines

#variation of Bias with sample size
RuahaB <- ggplot(sum.RFinsim,)+
  geom_line(aes(x=samplesize,y=bmR_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=bfR_mean,linetype = "Females"))+
  labs(x="Sample size", y="Bias")+ 
  theme_classic() #white background no grid lines


########################################################
#####           3. UGALLA POPULATION               #####
########################################################


## read Ugalla elephant survey data
Ugalla <- read.csv('Ugalla.csv') 

#preview the data
head(Ugalla)

# look at distribution of group sizes
hist(log(Ugalla$Total))


#Age sex distribution and proportional to the population
U.faead <- c(19,8,17,15,22,25,1) ### Female age distribution for Ugalla
faeadU <- U.faead/sum(U.faead) #make proportional to total
U.maead <- c(19,7,10,8,1,1,0.00001) ### Male age distribution for Ugalla
maeadU <- U.maead/sum(U.maead) #make proportional to total

# input parameters
U.popsize <- 700 # approx pop size of Ugalla elephants (reported 659)
samplesize <- c(0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.9) * U.popsize # create a sequence sample sizes ranging from 1% to 90% of the total pop
nsim <- 1000 # number of times to repeat the simulation for each combination of values

# function to calc where to stop sampling
stopfun <- function(vTR,element) {
  val = cumsum(vTR)
  return(length(val[val<=element])+1)
}

fin <- data.frame()

# start loop across sample sizes
for(k in samplesize){
  
  # loop across iterations
  for(j in 1:nsim){  
    
    # create empty dataframe
    df <- data.frame()
    
    ##### step 1 -- Sex and group size sampling
    # sample sex of groups
    sex <- rbinom(n = k,size=1,prob = fem.prop) # 0=male group 1=female group
    
    # sample group sizes for each sex
    grp <- rep(NA,k)
    grp[sex==0] <- rnbinom(n = length(sex[sex==0]),size = 2.247386,mu = 2.339081) 
    grp[sex==1] <- rnbinom(n = length(sex[sex==1]),size = 1.703228,mu = 13.517563) 
    
    # remove any groups with size=0
    sex <- sex[grp>0]
    grp <- grp[grp>0]
    
    # reduce group size vector
    sex <- sex[1:stopfun(grp,k)]
    grp <- grp[1:stopfun(grp,k)]
    
    # change sex of groups to be M and F
    sex[sex==0] <- 'M'
    sex[sex==1] <- 'F'
    
    ### step 2 randomly sample age distribution
    for(i in 1:length(grp)){
      sex.grp <- rep(sex[i],grp[i])
      
      if(sex[i]=='F'){
        agedist <- sample(x = fageclass,size = grp[i],replace=T,prob=faeadU)
        sex.grp[agedist %in% c("0","5","10")] <- sample(c('M','F'),
                                                        size = sum(agedist %in% c("0","5","10")),
                                                        prob = c(0.5,0.5),replace=T)  
      }
      if(sex[i]=='M'){
        
        agedist <- sample(x = mageclass,size = grp[i],replace=T,prob=maeadU[4:7]) 
      }
      
      # compile all data into a single dataframe
      df <- rbind(df,
                  data.frame(agedist=agedist,
                             sex=sex.grp))
      
    } # end i loop - sampling agedist
    
    ## compare samples to 'truth' for k samples
    df <- df[1:k,]
    Females <- df[df$sex=='F',] 
    Females$agedist <- factor(Females$agedist,levels=fageclass)
    
    Males <- df[df$sex=='M',] 
    Males$agedist <- factor(Males$agedist,levels=fageclass) # To incorporate all age categories
    
    U.obs.age.distf <- table(Females$agedist)/nrow(Females)
    U.obs.age.distm <- table(Males$agedist)/nrow(Males)
    
    #check if the ages for females and males are in the correct order
    if(!all(names(table(Females$agedist)) == fageclass)){
      U.obs.age.distf <- U.obs.age.distf[match(fageclass,names(table(Females$agedist)))]
    }
    
    if(!all(names(table(Males$agedist)) == fageclass)){
      U.obs.age.distm <- U.obs.age.distm[match(fageclass,names(table(Males$agedist)))]
    }
    
    # calculate the proportional difference (Bias)
    bfU <- mean((U.obs.age.distf - faeadU)/faeadU)
    bmU <- mean((U.obs.age.distm - maeadU)/maeadU)
    
    # calculate the root mean square error (Precision)-
    rmsefU <- sqrt(mean((U.obs.age.distf - faeadU)^2))
    rmsemU <- sqrt(mean((U.obs.age.distm - maeadU)^2))
    
    fin <- rbind(fin,
                 data.frame(           
                   age0f = U.obs.age.distf[[1]],
                   age5f = U.obs.age.distf[[2]],
                   age10f = U.obs.age.distf[[3]],
                   age15f = U.obs.age.distf[[4]],
                   age20f = U.obs.age.distf[[5]],
                   age25f = U.obs.age.distf[[6]],
                   age40f = U.obs.age.distf[[7]],
                   age0m = U.obs.age.distm[[1]],
                   age5m = U.obs.age.distm[[2]],
                   age10m = U.obs.age.distm[[3]],
                   age15m = U.obs.age.distm[[4]],
                   age20m = U.obs.age.distm[[5]],
                   age25m = U.obs.age.distm[[6]],
                   age40m = U.obs.age.distm[[7]],
                   bfU = bfU,
                   bmU = bmU,
                   rmsefU = rmsefU,
                   rmsemU = rmsemU,
                   samplesize=k, nsim=j))                 
    
  } # end j loop  - # iterations
}  # end k loop - sample size, i.e. proportion of the population

#Save the resultant data frame
write.csv(fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/U.Finalsimul.csv", row.names = TRUE)

# read simulated data
U.Finalsimul <- read.csv('U.Finalsimul.csv')

# Grouping and summarizing the data
U.sum.fin<- U.Finalsimul %>% 
  group_by(samplesize) %>% 
  summarise(bfU_mean = mean(bfU,na.rm=T), bmU_mean = mean(bmU,na.rm=T),
            rmsefU_mean = mean(rmsefU,na.rm=T),rmsemU_mean = mean(rmsemU,na.rm=T),
            bfU_sd = sd(bfU,na.rm=T),bmU_sd = sd(bmU,na.rm=T),
            rmsefU_sd = sd(rmsefU,na.rm=T),rmsemU_sd = sd(rmsemU,na.rm=T))

#Save the resultant data frame
write.csv(U.sum.fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/sum.UFinsim.csv", row.names = TRUE)

# read the summary of the final simulated data
sum.UFinsim <- read.csv('sum.UFinsim.csv')

#Combined plot
#variation of RMSE with sample size
Ugalla <- ggplot(sum.UFinsim,)+
  geom_line(aes(x=samplesize,y=rmsemU_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=rmsefU_mean,linetype = "Females"))+
  labs(x="Sample size", y="RMSE")+
  theme_classic() #white background no grid lines

#variation of Bias with sample size
UgallaB <- ggplot(sum.UFinsim,)+
  geom_line(aes(x=samplesize,y=bmU_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=bfU_mean,linetype = "Females"))+
  labs(x="Sample size", y="Bias")+
  theme_classic() #white background no grid lines


########################################################
#####           4. KATAVI POPULATION               #####
########################################################

# read Katavi elephant survey data 
Katavi <- read.csv('Katavi.csv') 

head(Katavi)

# look at distribution of group sizes
hist(log(Katavi$Total))

# For Katavi population
K.faead <- c(58,41,44,14,25,55,6) ### Female age distribution for Katavi ecosystem
faeadK <- K.faead/sum(K.faead) #make proportional of the female to total population
K.maead <- c(58,40,23,10,25,14,0.00001) ### Male age distribution for Katavi ecosystem
maeadK <- K.maead/sum(K.maead) #make proportional of the male to total population

# input parameters
K.popsize <- 3200 # approx pop size of Katavi elephants (reported 3128)
samplesize <- c(0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,0.7,0.9) * K.popsize # create a sequence sample sizes ranging from 1% to 90% of the total pop
nsim <- 1000 # number of times to repeat the simulation for each combination of values

# function to calc where to stop sampling
stopfun <- function(vTR,element) {
  val = cumsum(vTR)
  return(length(val[val<=element])+1)
}

fin <- data.frame()

# start loop across sample sizes
for(k in samplesize){
  
  # loop across iterations
  for(j in 1:nsim){  
    
    # create empty dataframe
    df <- data.frame()
    
    ##### step 1 -- Sex and group size sampling
    # sample sex of groups
    sex <- rbinom(n = k,size=1,prob = fem.prop) # 0=male group 1=female group
    
    # sample group sizes for each sex
    grp <- rep(NA,k)
    grp[sex==0] <- rnbinom(n = length(sex[sex==0]),size = 2.247386,mu = 2.339081) 
    grp[sex==1] <- rnbinom(n = length(sex[sex==1]),size = 1.703228,mu = 13.517563) 
    
    # remove any groups with size=0
    sex <- sex[grp>0]
    grp <- grp[grp>0]
    
    # reduce group size vector
    sex <- sex[1:stopfun(grp,k)]
    grp <- grp[1:stopfun(grp,k)]
    
    # change sex of groups to be M and F
    sex[sex==0] <- 'M'
    sex[sex==1] <- 'F'
    
    ### step 2 randomly sample age distribution
    for(i in 1:length(grp)){
      sex.grp <- rep(sex[i],grp[i])
      
      if(sex[i]=='F'){
        agedist <- sample(x = fageclass,size = grp[i],replace=T,prob=faeadK)
        sex.grp[agedist %in% c("0","5","10")] <- sample(c('M','F'),
                                                        size = sum(agedist %in% c("0","5","10")),
                                                        prob = c(0.5,0.5),replace=T)  
      }
      if(sex[i]=='M'){
        
        agedist <- sample(x = mageclass,size = grp[i],replace=T,prob=maeadK[4:7]) 
      }
      
      # compile all data into a single dataframe
      df <- rbind(df,
                  data.frame(agedist=agedist,
                             sex=sex.grp))
      
    } # end i loop - sampling agedist
    
    ## compare samples to 'truth' for k samples
    df <- df[1:k,]
    Females <- df[df$sex=='F',] 
    Females$agedist <- factor(Females$agedist,levels=fageclass)
    
    Males <- df[df$sex=='M',] 
    Males$agedist <- factor(Males$agedist,levels=fageclass) # To incorporate all age categories
    
    K.obs.age.distf <- table(Females$agedist)/nrow(Females)
    K.obs.age.distm <- table(Males$agedist)/nrow(Males)
    
    #check if the ages for females and males are in the correct order
    if(!all(names(table(Females$agedist)) == fageclass)){
      K.obs.age.distf <- K.obs.age.distf[match(fageclass,names(table(Females$agedist)))]
    }
    
    if(!all(names(table(Males$agedist)) == fageclass)){
      K.obs.age.distm <- K.obs.age.distm[match(fageclass,names(table(Males$agedist)))]
    }
    
    # calculate the proportional difference (Bias)
    bfK <- mean((K.obs.age.distf - faeadK)/faeadK)
    bmK <- mean((K.obs.age.distm - maeadK)/maeadK)
    
    # calculate the root mean square error (Precision)-
    rmsefK <- sqrt(mean((K.obs.age.distf - faeadK)^2))
    rmsemK <- sqrt(mean((K.obs.age.distm - maeadK)^2))
    
    fin <- rbind(fin,
                 data.frame(           
                   age0f = K.obs.age.distf[[1]],
                   age5f = K.obs.age.distf[[2]],
                   age10f = K.obs.age.distf[[3]],
                   age15f = K.obs.age.distf[[4]],
                   age20f = K.obs.age.distf[[5]],
                   age25f = K.obs.age.distf[[6]],
                   age40f = K.obs.age.distf[[7]],
                   age0m = K.obs.age.distm[[1]],
                   age5m = K.obs.age.distm[[2]],
                   age10m = K.obs.age.distm[[3]],
                   age15m = K.obs.age.distm[[4]],
                   age20m = K.obs.age.distm[[5]],
                   age25m = K.obs.age.distm[[6]],
                   age40m = K.obs.age.distm[[7]],
                   bfK = bfK,
                   bmK = bmK,
                   rmsefK = rmsefK,
                   rmsemK = rmsemK,
                   samplesize=k, nsim=j))                 
    
  } # end j loop  - # iterations
}  # end k loop - sample size, i.e. proportion of the population


#Save the resultant data frame
write.csv(fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/K.Finalsimul.csv", row.names = TRUE)

# read simulated data
K.Finalsimul <- read.csv('K.Finalsimul.csv')

# Grouping and summarizing the data
K.sum.fin<- K.Finalsimul %>% 
  group_by(samplesize) %>% 
  summarise(bfK_mean = mean(bfK,na.rm=T), bmK_mean = mean(bmK,na.rm=T),
            rmsefK_mean = mean(rmsefK,na.rm=T),rmsemK_mean = mean(rmsemK,na.rm=T),
            bfK_sd = sd(bfK,na.rm=T),bmK_sd = sd(bmK,na.rm=T),
            rmsefK_sd = sd(rmsefK,na.rm=T),rmsemK_sd = sd(rmsemK,na.rm=T))

#Save the resultant data frame
write.csv(K.sum.fin,"C:/Users/Juma/OneDrive - University of Glasgow/Documents/MSc/Research Project/Project Morrison/Power Analysis/classification/sum.KFinsim.csv", row.names = TRUE)

# read the summary of the final simulated data
sum.KFinsim <- read.csv('sum.KFinsim.csv')

#Combined plot
#variation of RMSE with sample size
Katavi <- ggplot(sum.KFinsim,)+
  geom_line(aes(x=samplesize,y=rmsemK_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=rmsefK_mean,linetype = "Females"))+
  labs(x="Sample size", y="RMSE")+
  theme_classic() #white background no grid lines

#variation of Bias with sample size
KataviB <- ggplot(sum.KFinsim,)+
  geom_line(aes(x=samplesize,y=bmK_mean,linetype = "Males"))+
  geom_line(aes(x=samplesize,y=bfK_mean,linetype = "Females"))+
  labs(x="Sample size", y="Bias")+
  theme_classic() #white background no grid lines

#Arrange the graphs into panels
#For RMSE against Sample size
figureRSME<- ggarrange(Ugalla + rremove("ylab") + rremove("xlab") + xlim(0, 650),# remove axis labels and set x scale from the plot
          Katavi + rremove("ylab") + rremove("xlab") + xlim(0, 3000),
          Serengeti+ rremove("ylab") + rremove("xlab")+ xlim(0, 7000),
          Ruaha + rremove("ylab") + rremove("xlab")+ xlim(0, 15000), 
          labels = c('A','B','C','D'),ncol = 2, 
          nrow = 2,common.legend = TRUE,legend="top")
annotate_figure(figureRSME, left = textGrob("RMSE", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Sample size", gp = gpar(cex = 1)))

#For Bias against Sample size
figure<- ggarrange(UgallaB + rremove("ylab") + rremove("xlab")+ xlim(0, 650),
                   KataviB + rremove("ylab") + rremove("xlab")+ xlim(0, 3000),
                   SerengetiB + rremove("ylab") + rremove("xlab")+ xlim(0, 7000),
                   RuahaB + rremove("ylab") + rremove("xlab")+ xlim(0, 15000), 
                   labels = c('A','B','C','D'),ncol = 2, 
                   nrow = 2,common.legend = TRUE,legend="top")

annotate_figure(figure, left = textGrob("Bias", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Sample size", gp = gpar(cex = 1)))
