## devtools::install_github("alicebalard/parasiteLoad", ref = "v2.1")
library(parasiteLoad)
library(optimx)
library(fitdistrplus)
library(parallel)

library(ggplot2)

# for reproducibility
set.seed(8)

>######################################
### Load data from the 4 field studies
### For Sage et al. and Moulia et al.: extract negbin parameters for simulation of bigger study

# Sage et al. 1986
# no easily interpretable individual data, so I took farm level data
SAGdata <- read.csv("data4fieldStudies/Sage1986.csv")
SAGdata$HI <- (SAGdata$Hioriginal +8) /16
SAGdata$load <- SAGdata$mean_abundance

descdist(SAGdata$load)
hist(SAGdata$load)
fitSAG <- fitdist(data = round(SAGdata$load), distr = "nbinom")
fitSAG

# Moulia et al. 1993
MOUdata <- read.csv("data4fieldStudies/moulia1993.csv") # N=120
MOUdata$HI <- 0.5
MOUdata$HI[MOUdata$genotype == "Mmm"] <- 1
MOUdata$HI[MOUdata$genotype == "Mmd"] <- 0

hist(MOUdata$load)
descdist(MOUdata$load)
fitMOU <- fitdist(data = round(MOUdata$load), distr = "nbinom")
fitMOU

# Baird et al. 2  012
BAIdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/WATWMdata.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BAIdata <- BAIdata[!is.na(BAIdata$HI) & !is.na(BAIdata$Sex) & !is.na(BAIdata$Aspiculuris.Syphacia),]

# Balard et al. 2020
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex) & !is.na(BALdata$Aspiculuris_Syphacia),]

######################################
### Simulate small & big datasets for each study 
### Shuffle the the nematode load to make it random along HI

### -> to do from there: BS

make_data_shuffle <- function(size=100, mu = 90, theta = 0.5, real.data=FALSE){
    if(isFALSE(real.data)) {
        df <- data.frame(HI = sample(seq(0,1, 0.001), size = size),
                         load = rnegbin(n = 100, mu = mu,   theta = theta),
                         Sex = as.factor(sample(rep(c("F", "M"), 1000), size = size)))
        return(df)
    }
    if(isTRUE(real.data)) {
        ## Sage et al: simulate N=100, simulate N=600
        SAG100_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 100),
                                load = rnegbin(n = 100, mu = fitSAG$estimate[["mu"]],
                                               theta = fitSAG$estimate[["size"]]),
                                Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 100)))
        SAG600_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 600),
                                load = rnegbin(n = 600, mu = fitSAG$estimate[["mu"]],
                                               theta = fitSAG$estimate[["size"]]),
                                Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 600)))
        ## Moulia et al: simulate N=100 & N=600
        MOU100_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 100),
                                load = rnegbin(n = 100, mu = fitMOU$estimate[["mu"]],
                                               theta = fitMOU$estimate[["size"]]),
                                Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 100)))
        MOU600_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 600),
                                load = rnegbin(n = 600, mu = fitMOU$estimate[["mu"]],
                                               theta = fitMOU$estimate[["size"]]),
                                Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 600)))
        ## Baird et al: subsample N=100 & N=600
        BAI100_df <- data.frame(HI = sample(BAIdata$HI, size = 100), 
                                load = sample(BAIdata$Aspiculuris.Syphacia, size = 100),
                                Sex = as.factor(sample(factor(BAIdata$Sex), size = 100)))
        BAI600_df <- data.frame(HI = sample(BAIdata$HI, size = 600), 
                                load = sample(BAIdata$Aspiculuris.Syphacia, size = 600),
                                Sex = as.factor(sample(factor(BAIdata$Sex), size = 600)))
        ## Balard et al: subsample N=100 & N=600 (585 for BAL)
        BAL100_df <- data.frame(HI = sample(BALdata$HI, size = 100), 
                                load = sample(BALdata$Aspiculuris_Syphacia, size = 100),
                                Sex = as.factor(sample(factor(BALdata$Sex), size = 100)))
        BAL585_df <- data.frame(HI = sample(BALdata$HI, size = 585), 
                                load = sample(BALdata$Aspiculuris_Syphacia, size = 585),
                                Sex = as.factor(sample(factor(BALdata$Sex), size = 585)))
        return(list(SAG100_df = SAG100_df, SAG600_df = SAG600_df, 
                    MOU100_df = MOU100_df, MOU600_df = MOU600_df, 
                    BAI100_df = BAI100_df, BAI600_df = BAI600_df, 
                    BAL100_df = BAL100_df, BAL585_df = BAL585_df))
    }
}

######################################
# Apply the three test on each df:
runTests <- function(df, which=c(1, 2, 3)){
    results <- list()
    if(1%in%which) {
        ## Test 1: chi2 wormy VS non wormy, cut at 250, hybrids 12.5% < HI < 87.5%
        ## 250 is the ~ 0.9 quantile, better doing 0.8 for more power
        chi2 <- chisq.test(df$load>quantile(df$load,0.8),
                           df$HI<0.125|df$HI>0.875)
        results["chi2"]  <- chi2$p.value
    }
    ## Test 2: kruskal-Wallis test on load, 20% < HI < 60%
    if(2%in%which) {
        df$genotype <- "Hybrid"
        df$genotype[df$HI < 0.2] <- "Mmd"
        df$genotype[df$HI > 0.6] <- "Mmm"
        KW <- kruskal.test(load ~ genotype, data = df)
        results["KW"] <- KW$p.value
    }
    ## Test 3: maximum likelihood estimation and LRT, H0 (no diff in sex and both sides)
    if(3%in%which) {
        ML <- parasiteLoad::analyse(data = df, response = "load", model = "negbin",
                                    group = "Sex", hybridIndex = "HI")
        results["ML"] <- ML$H0$Gtest$pvalue
    }
    return(results)
}



mybootstrap <- function(numBS){
    mclapply(1:numBS, function(i) {
        d <- make_data_shuffle()
        lapply(d, runTests, c(1,2))
    }, mc.cores=64)
}

numBS=100

try <- mybootstrap(numBS)


## resultsDF <-  do.call(rbind, lapply(try, unlist))

resultsDF <- reshape::melt(try)

tapply(resultsDF$value, list(resultsDF$L3, resultsDF$L2) ,
       function(x) length(x[x<0.05]))

head(resultsDF)

pheatmap(resultsDF)

write.csv(x = resultsDF,
          file = "/home/alice/Desktop/Git/Article3_review/output/resultsDF.csv",
          quote = "F")


## 240 F2 with HI->rbinom(240, 40, 0.5)/40 + 48 HI = 0 and 48 HI = 1
