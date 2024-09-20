if (!require(pacman)) install.packages("pacman")
suppressPackageStartupMessages(pacman::p_load(
  tidyverse,
  ggplot2,
  zoo, 
  table1,
  minpack.lm,
  broom,
  nlmixr2,
  pander,
  cowplot,
  kableExtra,
  tidyvpc,
  shinythemes
))


##############################################################
##############################################################

#Dataset formatting - PK
moxi_pk_source <- readRDS("Moxi_NHP_PKdata.RDS")
moxi_pd_source <- readRDS("Moxi_NHP_PDdata.RDS")

colnames(moxi_pk_source) <- toupper(colnames(moxi_pk_source))
moxi_pk <- moxi_pk_source %>%
  mutate(
    ID = as.integer(ID),
    GENDER = as.character(GENDER)
    #GENDER = factor("GENDER", levels=c(0,1), labels=c("Male","Female"))
  )
moxi_pk$SEX <- factor(moxi_pk$GENDER, levels=c("Male","Female"), labels=c(0,1))
moxi_pk$TIME <- 24*(moxi_pk$DAY-1) + moxi_pk$HOUR
moxi_pk$EVID <- 0
moxi_pk$AMT <- NA
moxi_pk$DV <- moxi_pk$CONCENTRATION


moxi_pk <- moxi_pk %>% 
  select(c("ID", "TIME", "DAY", "HOUR", "DV", "AMT", "EVID", "CONCENTRATION", "SEX", "TREATMENT"))

moxi_pk_doses <- moxi_pk %>% distinct(ID,DAY,.keep_all = T)
moxi_pk_doses <- moxi_pk_doses %>% 
  mutate(
    TIME = (DAY-1)*24,
    HOUR = 0,
    DV = NA,
    AMT = 80000,
    EVID = 1,
    CONCENTRATION = NA
  )

moxi_pk <- rbind(moxi_pk,moxi_pk_doses) %>%
  arrange(AMT) %>% 
  arrange(DAY) %>% 
  arrange(TIME) %>% 
  arrange(ID)

moxi_pk$MDV <- ifelse(moxi_pk$DV==0 | is.na(moxi_pk$DV),1,0)

moxi_pk$OCC1 <- ifelse(moxi_pk$DAY==3, 1, 0)
moxi_pk$OCC2 <- ifelse(moxi_pk$DAY==10, 1, 0)
moxi_pk$OCC3 <- ifelse(moxi_pk$DAY==15, 1, 0)
moxi_pk$OCC4 <- ifelse(moxi_pk$DAY==17, 1, 0)


moxi_pk_formatted <- moxi_pk[, c("ID","TIME","DV","AMT","EVID","MDV","OCC1","OCC2","OCC3","OCC4","SEX","TREATMENT","DAY","HOUR")]

moxi_pk_formatted$DV <- ifelse(moxi_pk_formatted$DV==0 & moxi_pk_formatted$DAY!=15, .01, moxi_pk_formatted$DV)

moxi_weights <- read_csv("Weight_Data.csv")
moxi_weights$WT <- moxi_weights$AVG
moxi_weights <- moxi_weights %>% select("ID","WT")

moxi_pk_formatted <- merge(moxi_pk_formatted,moxi_weights,by="ID")
moxi_pk_formatted$AMT <- moxi_pk_formatted$AMT * moxi_pk_formatted$WT
moxi_pk_formatted$WTALOM <- (moxi_pk_formatted$WT/3.5)
moxi_pk_formatted$WTMED <- moxi_pk_formatted$WT-3.5
#moxi_pk_formatted$WTLOG <- log(moxi_pk_formatted$WT)

moxi_pk_formatted <- moxi_pk_formatted %>%
  arrange(AMT) %>% 
  arrange(DAY) %>% 
  arrange(TIME) %>% 
  arrange(ID)


moxi_pk_formatted$MDV <- ifelse(moxi_pk_formatted$DV != 0 & moxi_pk_formatted$TIME==336 & moxi_pk_formatted$EVID==0, 1, moxi_pk_formatted$MDV)
moxi_pk_formatted$DV <- ifelse(moxi_pk_formatted$DV == 0 & moxi_pk_formatted$TIME==336, NA, moxi_pk_formatted$DV)
moxi_pk_formatted$TIME <- ifelse(moxi_pk_formatted$EVID==0 & moxi_pk_formatted$TIME==384, moxi_pk_formatted$TIME-.001, moxi_pk_formatted$TIME)


remove_vomit <- function(inData, later_occ=F){
  early_occ_vomit <- c(1010, 1021, 1503, 1510, 1518, 1520, 1522, 1524)
  later_occ_vomit <- c(1501, 1510, 1521)
  
  data <- inData %>% filter((!inData$ID %in% later_occ_vomit) & (inData$DAY>=15 | !inData$ID %in% early_occ_vomit) )
}


get_pk_dataset <- function(later_occ=F,remove_vomit=F){
  data <- (moxi_pk_formatted)
  if (later_occ){
    data <-  data %>% filter(OCC3==1 | OCC4==1)
  }
  if (remove_vomit){
    data <- remove_vomit(data)
  }
  return(data)
}


##############################################################
# Dataset Formatting - PD



##############################################################
##############################################################
# Models

one.compartment.proportional <- function(){
  #Weight median = 3.5 kg
  ini({
    tka <- log(.8); label("Ka")
    tcl <- log(1); label("Cl")
    tv <- log(30); label("V")
    
    #covka <- c(-0.3,0.01, 0.3)    
    covka <- -.1
    allcl <- .75 # Unfixed OBJF = 12714.77, est. value = 1.01 | Fixed OBJF = 12721.18
    allv <- fix(1) # Fixed OBJF = 12714.61 | Unfixed OBJF = 12714.92, est. value = 0.957
    
    #etas
    eta.ka ~ .1
    
    eta.cl + eta.v ~ c(.1, .01, .1)
    
    #eta.cl ~ .05
    #eta.v ~ .1
    
    prop.sd <- 0.1
    #add.sd <- 0.1
  })
  
  model({
    #ka <- exp(tka + eta.ka + log(1+covka*(WTMED))) #Linear - 12719.71
    #ka <- exp(tka + eta.ka + covka*log(WTALOM)) #Power - 12714.45
    ka <- exp(tka + eta.ka + (covka*WTMED))  #Exponential - 12714.61 - Best since it's also good outside of range and BSV of Ka went from ~60% to 54% (5% drop)
    #ka <- exp(tka + eta.ka)
    
    cl <- exp(tcl + eta.cl + allcl*log(WTALOM))# * exp(allcl*log(WT/3.5))
    v <- exp(tv + eta.v + allv*log(WTALOM))
    
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v) * center
    
    cp <- center / v
    cp ~ prop(prop.sd) #+ add(add.sd) # Proportional OBJF = 12714.77 | Combined OBJF = 12714.4 (add.sd = 0.0846)
  })
}

two.compartment.proportional <- function(){
  ini({
    tka <- log(.8); label("Ka")
    tcl <- log(5); label("Cl")
    tq <- log(1); label("Q")
    
    tv1 <- log(30); label("V1")
    tv2 <- log(10); label("V2")
    
    allcl <- .75
    allv1 <- fix(1)
    
    allq <- fix(.75)
    allv2 <- fix(1)
    
    #lagtime <- log(.1)
    #etas
    eta.ka ~ .1
    eta.cl ~ .05
    eta.v1 ~ .1
    #eta.cl + eta.v ~ c(.1, .01, .1)
    eta.q ~ 0.05
    
    #eta.v1 ~ .1
    eta.v2 ~ fixed(0)
    
    prop.sd <- 0.1
    #add.sd <- 0.1
  })
  
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + allcl*log(WTALOM))
    q <- exp(tq + eta.q + allq*log(WTALOM))
    
    v1 <- exp(tv1 + eta.v1 + allv1*log(WTALOM))
    v2 <- exp(tv2 + eta.v2 + allv2*log(WTALOM))
    
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v1) * center + peripheral * (q / v2)
    d/dt(peripheral) <- center * (q / v1) - peripheral * (q / v2)
    #lag(center) = exp(lagtime)
    
    cp <- center / v1
    cp ~ prop(prop.sd)# + add(add.sd) #add = additive error, prop=proportional error
  })
}

##############################################################
##############################################################
# VPC Plotting

get_vpc <- function(fit_used, scale_y_log10=T, scale_x_log10=F, use_tad=F, stratify_by_sex=F, pred_corrected=F){
  vpc_sim <- vpcSim(object = fit_used[fit_used$EVID==0,])
  obs_data <- fit_used
  
  if (use_tad){
    vpc <- observed(obs_data, x=tad, y=DV) %>%
      simulated(vpc_sim, x=tad, y=sim) %>% 
      binning(bin = "jenks", nbins=44)
  }
  else{
    vpc <- observed(obs_data, x=TIME, y=DV) %>%
      simulated(vpc_sim, x=TIME, y=sim) %>% 
      binning(bin = "jenks", nbins=44)
  }
  
  if (stratify_by_sex){
    vpc <- vpc %>% stratify(~SEX)
  }
  if (pred_corrected){
    vpc <- vpc %>% predcorrect(pred=PRED)
  }
  
  vpc <- vpc %>% vpcstats()
  plot(vpc) + scale_y_log10(limits = c(.01,100000))
  
}



one.compartment.additive <- function(){
  ini({
    tka <- log(.7); label("Ka")
    tcl <- log(.72); label("Cl")
    tv <- log(4.9); label("V")
    
    allocl <- .75
    allov <- 1
    
    eta.ka ~ fixed(0)
    eta.cl ~ 0.3
    eta.v ~ 0.3
    
    add.sd <- 0.3
  })
  
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + allocl * log(WTALOM))
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v) * center
    
    cp <- center / v
    cp ~ prop(add.sd) #add = additive error, prop=proportional error
  })
}
#fit.one.add <- nlmixr2(one.compartment.additive, dataset_name,  est="focei", foceiControl(print=0))










