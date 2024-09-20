#TO DO
# Functions:
#  1: Box plots of etas faceted by sex (separate plots per eta) - not included in dataset
#  2: Scatterplots of etas vs. covariates (p value <5% = correlation) - FIGURE OUT HOW TO SEE CORRELATION
#  3: make each point distinct (weight/id) - HELP
#  4: print fit$objdf, fit$parFixed, fit$omega, and graph_goodness_of_fit/other plots

#Done:
#  1: make 1, 2, 3 compartment models for warfarin
#  2: try proportional and combined error (add(add.sd), prop(prop.sd), and add+prop)

# 9 total reports, keep same template for each model (type of error and # of compartments)




##############################################################
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
  ggpubr
))

#source("Compartment_Models.R")
##############################################################

num_patients <- 8


##############################################################

lightblack<-'#212121'
darkblack<-'#0c0c0c'

lightblue<-'#00447c'
darkblue<-'#002756'

lightcyan<-'#0b7b9e'
darkcyan<-'#035b7a'

lightbeige<-'#b09b72'
darkbeige<-'#938160'

lightwhite<-'#ededed'
darkwhite<-'#cccccc'

lightred<-'#540000'
darkred<-'#3a0101'

lightgreen<-'#004749'
darkgreen<-'#013b3f'

#Basic functions
##############################################################
graph_scatterplot <- function(inData, xAxis, yAxis, xLab='', yLab='', to_facet=NULL, facet_size=0, page=1, point=FALSE, line=FALSE, trend_line=NULL, hline=NULL, abline=FALSE, log10_y=FALSE, log10_x=FALSE){
  xAxisCol <- inData[[xAxis]]
  yAxisCol <- inData[[yAxis]]
  
  if (facet_size!=0){
    #plot <- ggplot(data=inData %>% filter(ID<5),#to_facet>=((page-1)*facet_size) & to_facet<(page*facet_size)
    #aes(x=xAxis,y=yAxis)) + facet_wrap(to_facet)
    
    plot <- ggplot(data=inData %>% filter(ID<5), 
                   aes(x=xAxisCol, y=yAxisCol, group = ID))
  }
  else{
    plot <- ggplot(data=inData ,aes(x=xAxisCol,y=yAxisCol))
  }
  
  if (!is.null(hline)){
    plot <- plot + geom_hline(yintercept=hline, color=darkbeige, linewidth=1, linetype=5)
  }
  if (abline){
    plot <- plot + geom_abline(color=lightblack, linewidth=1)
  }
  if (!is.null(trend_line)){
    plot <- plot + geom_smooth(method=trend_line, linewidth=1, color=darkbeige, fill=darkwhite)
  }
  if (point){
    plot <- plot + geom_point(color=darkgreen, size=2, alpha=.9)
  }
  if (line){
    plot <- plot + geom_line(color=darkcyan, linewidth=1)
  }
  if (log10_y){
    plot <- plot + scale_y_log10()
  }
  if (log10_x){
    plot <- plot + scale_x_log10()
  }
  if (!is.null(to_facet)){
    plot <- plot + facet_wrap(to_facet)
  }
  plot <- plot + theme_bw() + labs(x=xLab, y=yLab)
  
  plot
}


graphBoxplot <- function(inData, xAxis, yAxis, xAxisOverride, xLab='', yLab='', to_facet=NULL, facet_size=0, page=1, outlier_shape=1){
  if (is.null(xAxisOverride)){
    xAxisCol <- inData[[xAxis]]
  }
  else {
    xAxisCol <- xAxisOverride
  }
  yAxisCol <- inData[[yAxis]]
  
  if (facet_size!=0){
    plot <- ggplot(data=inData %>% filter(ID<5), 
                   aes(x=xAxisCol, y=yAxisCol, group = ID))
  }
  else{
    plot <- ggplot(data=inData ,aes(x=xAxisCol,y=yAxisCol))
  }
  
  plot <- plot + geom_boxplot(color = darkblack, fill=lightbeige, alpha=.7, linewidth=.6, 
                              outliers=FALSE,
                              outlier.color=darkgreen, outlier.size = 2, outlier.alpha=1, outlier.shape = outlier_shape) + 
    geom_jitter(width=.2,size=2,color=darkgreen)
  
  
  
  if (!is.null(to_facet)){
    plot <- plot + facet_wrap(to_facet)
  }
  plot <- plot + theme_bw() + labs(x=xLab, y=yLab)
  
  plot
}



#Scatterplots
##############################################################

##USE plot_grid

#DV vs CPRED + identity line + trend line
graph_dv_vs_cpred <- function(inData, scale_log10=F){
  min_lim <- min(c(min(na.omit(inData)$CPRED), min(na.omit(inData)$DV)))
  max_lim <- max(c(max(na.omit(inData)$CPRED), max(na.omit(inData)$DV)))
  
  if (scale_log10){
    graph_scatterplot(inData, xAxis='CPRED', yAxis='DV', xLab='CPRED', yLab='DV', abline=TRUE, point=TRUE, trend_line = 'lm', log10_y=TRUE, log10_x=TRUE) +
      scale_x_log10(limits=c(min_lim,max_lim)) +
      scale_y_log10(limits=c(min_lim,max_lim)) +
      coord_fixed(ratio=1)
  }
  else{
    graph_scatterplot(inData, xAxis='CPRED', yAxis='DV', xLab='CPRED', yLab='DV', abline=TRUE, point=TRUE, trend_line = 'lm', log10_y=TRUE, log10_x=TRUE) +
      scale_x_continuous(limits=c(min_lim,max_lim)) +
      scale_y_continuous(limits=c(min_lim,max_lim)) +
      coord_fixed(ratio=1)
  }
}
#graph_dv_vs_cpred(fit_used)

#DV vs IPRED + identity line + trend line
graph_dv_vs_ipred <- function(inData, scale_log10=F){
  min_lim <- min(c(min(na.omit(inData)$IPRED), min(na.omit(inData)$DV)))
  max_lim <- max(c(max(na.omit(inData)$IPRED), max(na.omit(inData)$DV)))
  
  if (scale_log10){
    graph_scatterplot(inData, xAxis='IPRED', yAxis='DV', xLab='IPRED', yLab='DV', abline=TRUE, point=TRUE, trend_line = 'lm', log10_y=TRUE, log10_x=TRUE) +
      scale_x_log10(limits=c(min_lim,max_lim)) +
      scale_y_log10(limits=c(min_lim,max_lim)) +
      theme(aspect.ratio=1)
  }
  else{
    graph_scatterplot(inData, xAxis='IPRED', yAxis='DV', xLab='IPRED', yLab='DV', abline=TRUE, point=TRUE, trend_line = 'lm', log10_y=TRUE, log10_x=TRUE) +
      scale_x_continuous(limits=c(min_lim,max_lim)) +
      scale_y_continuous(limits=c(min_lim,max_lim)) +
      theme(aspect.ratio=1)
  }
}
#graph_dv_vs_ipred(fit_used)

#CWRES vs CPRED + trend line + horizontal lines at -2, 0 and 2
graph_cwres_vs_cpred <- function(inData, scale_x_log10=F, trendline_type='lm'){ #trendlines = lm or loess
  plot <- graph_scatterplot(inData, 
                            xAxis='CPRED', yAxis='CWRES', 
                            xLab='CPRED', yLab='CWRES', 
                            point=TRUE, 
                            hline=c(-2, 0, 2),
                            log10_x = scale_x_log10)
  
  if (trendline_type=='both'){
    plot$layers <- c(
      geom_smooth(method='loess', linewidth=1, color=darkbeige, alpha=0, linetype="dashed"),
      geom_smooth(method='lm', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
    
  }
  else if (trendline_type=='lm') {
    plot$layers <- c(
      geom_smooth(method='lm', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  else if (trendline_type=='loess') {
    plot$layers <- c(
      geom_smooth(method='loess', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  else {
    plot$layers <- c(
      geom_smooth(method=trendline_type, linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  plot
}
#graph_cwres_vs_cpred(fit_used, F, 'lm')

#IWRES(y) vs IPRED(x)
graph_iwres_vs_ipred <- function(inData, scale_x_log10=F, trendline_type='lm'){
  plot <- graph_scatterplot(inData, 
                            xAxis='IPRED', yAxis='IWRES', 
                            xLab='IPRED', yLab='IWRES', 
                            point=TRUE, 
                            hline=c(-2, 0, 2),
                            log10_x = scale_x_log10)
  #coord_fixed(ratio=1)
  
  if (trendline_type=='both'){
    plot$layers <- c(
      geom_smooth(method='loess', linewidth=1, color=darkbeige, alpha=0, linetype="dashed"),
      geom_smooth(method='lm', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
    
  }
  else if (trendline_type=='lm') {
    plot$layers <- c(
      geom_smooth(method='lm', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  else if (trendline_type=='loess') {
    plot$layers <- c(
      geom_smooth(method='loess', linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  else {
    plot$layers <- c(
      geom_smooth(method=trendline_type, linewidth=1, color=darkbeige, fill=darkwhite),
      plot$layers
    )
  }
  plot
}
#graph_iwres_vs_ipred(fit_used, F, 'both')

#CWRES vs TIME + trend line + horizontal lines at -2, 0 and 2
graph_cwres_vs_time <- function(inData, scale_x_log10=F){
  graph_scatterplot(inData, 
                    xAxis='TIME', yAxis='CWRES', 
                    trend_line = 'loess', 
                    xLab='TIME', yLab='CWRES',  
                    point=TRUE, 
                    hline=c(-2, 0, 2), 
                    log10_x=scale_x_log10)
  #coord_fixed(ratio=1)
}
#graph_cwres_vs_time(fit_used)

#IWRES vs TIME + trend line + horizontal lines at -2, 0 and 2
graph_iwres_vs_time <- function(inData, scale_x_log10=F){
  graph_scatterplot(inData, 
                    xAxis='TIME', yAxis='IWRES', 
                    trend_line = 'loess', 
                    xLab='TIME', yLab='IWRES',  
                    point=TRUE, 
                    hline=c(-2, 0, 2), 
                    log10_x=scale_x_log10)
  #coord_fixed(ratio=1)
}
#graph_iwres_vs_time(fit_used)

#CWRES(y) vs TAD(x)
graph_cwres_vs_tad <- function(inData, scale_x_log10=F){
  graph_scatterplot(inData, 
                    xAxis='tad', yAxis='CWRES', 
                    trend_line = 'loess', 
                    xLab='TAD', yLab='CWRES',  
                    point=TRUE, 
                    hline=c(-2, 0, 2), 
                    log10_x=scale_x_log10) 
  #coord_fixed(ratio=1)
}


#CWRES(y) vs TAD(x)
graph_iwres_vs_tad <- function(inData, scale_x_log10=F){
  graph_scatterplot(inData, 
                    xAxis='tad', yAxis='IWRES', 
                    trend_line = 'loess', 
                    xLab='TAD', yLab='IWRES',  
                    point=TRUE, 
                    hline=c(-2, 0, 2), 
                    log10_x=scale_x_log10)
  #coord_fixed(ratio=1)
}

#Individual plots DV vs time + IPRED line
graph_dv_vs_time <- function(inData, pageSize=9, page=0){
  inData<- inData %>% filter(as.numeric(ID)>page*pageSize & as.numeric(ID)<=page*pageSize+pageSize)
  graph_scatterplot(inData, xAxis='TIME', to_facet=~ID, yAxis='DV', xLab='TIME', yLab='DV', point=TRUE, log10_y = TRUE) +
    geom_line(aes(y=IPRED), show.legend=TRUE, color=darkgreen, linewidth=1, linetype=1) +
    geom_line(aes(y=PRED), color=lightred, linewidth=1, linetype=5) 
}

#Individual plots DV vs tad + IPRED line
graph_dv_vs_tad <- function(inData, pageSize=9, page=0){
  inData<- inData %>% filter(as.numeric(ID)>page*pageSize & as.numeric(ID)<=page*pageSize+pageSize)
  graph_scatterplot(inData, xAxis='tad', to_facet=~ID, yAxis='DV', xLab='TAD', yLab='DV', point=TRUE, log10_y = TRUE) +
    geom_line(aes(y=IPRED), show.legend=TRUE, color=darkgreen, linewidth=1, linetype=1) +
    geom_line(aes(y=PRED), color=lightred, linewidth=1, linetype=5) 
}

graph_goodness_of_fit <- function(inData){
  p1 <- graph_dv_vs_cpred(inData)
  p2 <- graph_cwres_vs_cpred(inData)
  p3 <- graph_cwres_vs_time(inData)
  p4 <- graph_cwres_vs_tad(inData)
  
  p5 <- graph_dv_vs_ipred(inData)
  p6 <- graph_iwres_vs_ipred(inData)
  p7 <- graph_iwres_vs_time(inData)
  p8 <- graph_iwres_vs_tad(inData)
  
  plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow=2, ncol=4, scale=.985)
}


#fit.F <- theo_model(1,"add")

#graph_goodness_of_fit(tempFit)



#Covariate Scatterplots
##############################################################

#Correlation between each eta and age, weight, sex
#dataset only shows weight
graph_wt_correlation <- function(inData){
  inData <- inData %>% distinct(ID, .keep_all = TRUE)
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  
  plots_list <- list()
  
  current_num <- 1
  
  for (eta_col in eta_columns){
    col_data <- inData[[eta_col]]
    p1 <- graph_scatterplot(inData, 
                            xAxis='WT', yAxis=eta_col, 
                            trend_line = 'lm', 
                            xLab='weight', yLab=eta_col,  
                            point=TRUE, 
                            log10_x=FALSE) + 
      stat_cor(method="pearson")
    plots_list[[length(plots_list) + 1]] = p1
  }
  
  plot_grid(plotlist=plots_list, ncol=3)
}
#graph_wt_correlation(fit_used)

graph_age_correlation <- function(inData){
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  
  plots_list <- list()
  
  current_num <- 1
  
  for (eta_col in eta_columns){
    col_data <- inData[[eta_col]]
    p1 <- graph_scatterplot(inData, 
                            xAxis='AGE', yAxis=eta_col, 
                            trend_line = 'lm', 
                            xLab='age', yLab=eta_col,  
                            point=TRUE, 
                            log10_x=FALSE)
    plots_list[[length(plots_list) + 1]] = p1
  }
  
  plot_grid(plotlist=plots_list, ncol=3)
}

graph_sex_correlation <- function(inData, num_col=3){
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  eta_table <- as.data.frame(inData) %>% 
    select(c("ID",eta_columns,"SEX")) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    mutate(SEX_LABEL = ifelse(SEX == 1, "male", "female"))
  
  plots_list <- list()
  
  for (eta_col in eta_columns){
    #eta_table <- stack(select(as.data.frame(inData) %>% distinct(ID, .keep_all = TRUE), eta_col))
    
    p1 <- graphBoxplot(eta_table, xAxisOverride=inData$SEX_LABEL, xAxis = "SEX_LABEL", yAxis=eta_col, xLab=eta_col, yLab='value') +
      stat_compare_means(method="t.test", ref.group="male")
    #facet_wrap(~SEX_LABEL) + 
    #stat_compare_means()
    
    plots_list[[length(plots_list) + 1]] <- p1
  }
  
  plot_grid(plotlist=plots_list, ncol=num_col)
}

#graph_sex_correlation(fit_used)

#Barplots/histogram
##############################################################
#2*n^(1/3)
graph_histogram <- function(inData, xAxis, xLab='', yLab='Density', numBins=ceiling(2*(num_patients^(1/3))), standard_deviation){
  #Example standard deviation: sqrt(inData$omega[sd_row,sd_col])
  #lim = max(xAxis)*1.1
  ggplot(data=inData, 
         aes(x=.data[[xAxis]], 
             after_stat(density))) + 
    #scale_x_continuous(limits=c(-lim,lim)) +
    labs(x=xLab, y=yLab) +
    geom_histogram(bins=numBins, color=darkgreen, fill=lightgreen, linewidth=1.5, alpha=.8) +
    geom_function(inherit.aes=FALSE, 
                  fun=dnorm, args=list(mean=0, sd=standard_deviation), 
                  color=darkbeige, linewidth=1) +
    geom_density(trim=TRUE, color=darkred, linewidth=1) +
    geom_point(aes(y=max(density(.data[[xAxis]])$y)/20), color=darkred, alpha=.7, size=2) +
    theme_bw()
}

residual_histogram_summary <- function(inData){
  p1 <- residual_histogram(inData,"CWRES")
  p2 <- residual_histogram(inData,"IWRES")
  p3 <- residual_histogram(inData,"RES")
  
  plot_grid(p1,p2,p3, ncol=3)
}


residual_histogram <- function(inData, type){
  inData <- na.omit(inData)
  num_bins <- ceiling(2*(nrow(inData)^(1/3)))
  
  if (type == "CWRES"){
    cwres_boundaries <- max(c(max(inData$CWRES), abs(min(inData$CWRES))))
    
    plot <- graph_histogram(inData=inData, xAxis='CWRES', numBins = num_bins, standard_deviation=sd(inData$CWRES), xLab='CWRES') +
      scale_x_continuous(limits=c(-cwres_boundaries,cwres_boundaries)) +
      theme(aspect.ratio=1)
  }
  else if (type == "IWRES"){
    iwres_boundaries <- max(c(max(inData$IWRES), abs(min(inData$IWRES))))
    
    plot <- graph_histogram(inData=inData, xAxis='IWRES', numBins = num_bins, standard_deviation=sd(inData$IWRES), xLab='IWRES') +
      scale_x_continuous(limits=c(-iwres_boundaries,iwres_boundaries))  +
      theme(aspect.ratio=1)
  }
  else {
    res_boundaries <- max(c(max(inData$RES), abs(min(inData$RES))))
    
    plot <- graph_histogram(inData=inData, xAxis='RES', numBins = num_bins, standard_deviation=sd(inData$RES), xLab='RES') +
      scale_x_continuous(limits=c(-res_boundaries,res_boundaries))  +
      theme(aspect.ratio=1)
  }
  plot
}
#residual_histogram_summary(fit_used)

eta_histograms <- function(inData, num_cols=3){
  num_bins <- ceiling(2*(length(unique(inData$ID))^(1/3)))
  
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  
  plots_list <- list()
  
  current_num <- 1
  
  for (eta_col in eta_columns){
    col_data <- inData[[eta_col]]
    
    limit <- max(c(max(inData[[eta_col]]), abs(min(inData[[eta_col]]))))
    
    p1 <- graph_histogram(inData=inData, 
                          xAxis=eta_col, 
                          standard_deviation=sqrt(inData$omega[eta_col,eta_col]), 
                          xLab=eta_col,
                          numBins = num_bins) + 
      scale_x_continuous(limits=c(-limit,limit)) +
      theme(aspect.ratio=1)
    plots_list[[length(plots_list) + 1]] = p1
    
    current_num <- current_num + 1
  }
  
  plot_grid(plotlist=plots_list, ncol=num_cols)
  
}
#eta_histograms(fit_used)

#eta_histograms(tempFit)
#eta_histograms(fit.F)


#Boxplots
##############################################################
#Problem: fit.F does not have column for sex

#tempFit <- theo_model(2, "add")
#eta_columns <- colnames(tempFit)[grep("eta", colnames(tempFit))]


eta_boxplots <- function(inData){
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  eta_table <- as.data.frame(inData) %>% 
    select(c("ID",eta_columns,"SEX")) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    mutate(SEX_LABEL = ifelse(SEX == 1, "female", "male"))
  
  plots_list <- list()
  
  for (eta_col in eta_columns){
    #eta_table <- stack(select(as.data.frame(inData) %>% distinct(ID, .keep_all = TRUE), eta_col))
    limit <- max(c(max(eta_table[[eta_col]]), abs(min(eta_table[[eta_col]]))))
    
    p1 <- graphBoxplot(eta_table, xAxisOverride=eta_col, yAxis=eta_col, xLab=eta_col, yLab='value') +
      scale_y_continuous(limits=c(-limit,limit))
    
    plots_list[[length(plots_list) + 1]] <- p1
  }
  
  plot_grid(plotlist=plots_list, ncol=3)
}
#eta_boxplots(fit_used)

num_etas <- function(inData){
  length(colnames(inData)[grep("eta", colnames(inData))])
}
#eta_boxplots(fit.F)

##############################################################


# Eta Correlation Plots

eta_correlations <- function(inData){
  inData <- inData %>% distinct(ID, .keep_all = TRUE)
  
  plots_list <- list()
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  
  for (col1 in eta_columns){
    for (col2 in eta_columns){
      plot <- graph_scatterplot(inData, 
                                xAxis=col1, yAxis=col2, 
                                trend_line = 'lm', 
                                xLab=col1, yLab=col2,  
                                point=TRUE, 
                                log10_x=FALSE) + 
        stat_cor(method="pearson") +
        theme(aspect.ratio=1)
        
      plots_list[[length(plots_list) + 1]] <- plot
    }
  }
  
  plot_grid(plotlist=plots_list)
}

eta_correlations_list <- function(inData){
  inData <- inData %>% distinct(ID, .keep_all = TRUE)
  
  plots_list <- list()
  eta_columns <- colnames(inData)[grep("eta", colnames(inData))]
  
  for (col1 in eta_columns){
    for (col2 in eta_columns){
      plot <- graph_scatterplot(inData, 
                                xAxis=col1, yAxis=col2, 
                                trend_line = 'lm', 
                                xLab=col1, yLab=col2,  
                                point=TRUE, 
                                log10_x=FALSE) + 
        stat_cor(method="pearson") +
        theme(aspect.ratio=1)
      
      plots_list[[length(plots_list) + 1]] <- plot
    }
  }
  
  plots_list
}

