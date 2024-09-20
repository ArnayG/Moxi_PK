library(shiny)
library(rxode2)
library(bslib)

source("Model_Report_Functions.R", chdir=TRUE)
source("Moxi_PKPD.R", chdir=TRUE)

final_model <- function(){
  ini({
    tka <- fix(log(0.7630404)); label("Ka")
    tcl <- fix(log(3.2363909)); label("Cl")
    tv <- fix(log(33.9763030)); label("V")
    
    covka <- fix(-0.2138061)
    allcl <- fix(1.0050231)
    allv <- fix(1.0000000)
    
    eta.ka ~ fix(0.2676694)
    eta.cl + eta.v ~ fix(c(0.04400711, 0.02083013, 0.01501074))
    
    prop.sd <- fix(0.2756416)
  })
  
  model({
    WTMED <- 3.5-3.5
    WTALOM <- 3.5/3.5
    
    ka <- exp(tka + eta.ka + (covka*WTMED))
    
    cl <- exp(tcl + eta.cl + allcl*log(WTALOM))
    v <- exp(tv + eta.v + allv*log(WTALOM))
    
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v) * center
    
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

final_model_unfixed <- function(){
  ini({
    tka <- log(0.7630404); label("Ka")
    tcl <- log(3.2363909); label("Cl")
    tv <- log(33.9763030); label("V")
    
    covka <- -0.2138061
    allcl <- 1.0050231
    allv <- fix(1.0000000)
    
    eta.ka ~ 0.2676694
    eta.cl + eta.v ~ c(0.04400711, 0.02083013, 0.01501074)
    
    prop.sd <- 0.2756416
  })
  
  model({
    WTMED <- 3.5-3.5
    WTALOM <- 3.5/3.5
    
    ka <- exp(tka + eta.ka + (covka*WTMED))
    
    cl <- exp(tcl + eta.cl + allcl*log(WTALOM))
    v <- exp(tv + eta.v + allv*log(WTALOM))
    
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v) * center
    
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}


final_model <- rxode(final_model)
final_model <- final_model$simulationModel


final_omegas <- lotri::lotri({
  eta.ka ~ 0.2676694
  eta.cl + eta.v ~ c(0.04400711, 0.02083013, 0.01501074)
})


make_pk_sim <- function(wt, dose, num_doses=1, interval=24, final_time_override=0){
  final_time <- (num_doses+2)*interval
  if (final_time_override!=0){
    final_time <- final_time_override
  }
  
  ev  <- et(amountUnits="mg", timeUnits="hours") %>% 
    et(amt=dose*wt*1000, cmt="depot", addl=num_doses-1, ii=interval) %>%
    et(0:final_time) %>%
    et(id=1:100)
  
  set.seed(10)
  rxSetSeed(10)
  
  cov_df <- data.frame(id=1:100,WT=rnorm(100,wt,0.8797681))
  cov_df$WTMED <- cov_df$WT-3.5
  cov_df$WTALOM <- cov_df$WT/3.5
  
  
  sim <- rxSolve(final_model, ev, 
                 omega=final_omegas, 
                 #params = params_wt, 
                 nSub=100,
                 iCov=cov_df)
  
  return(sim)
}


make_pk_plot <- function(wt, dose, num_doses=1, interval=24, final_time_override=0){
  #dose = mg/kg
  sim <- make_pk_sim(wt, dose, num_doses, interval, final_time_override)
  
  return(
    ggplot(data=sim,aes(x=time,y=cp, group=id))+
      geom_line(alpha=.2, linewidth=1) +
      scale_y_log10(limits=c(1,100000))
    )
}


pk_simulation_given_data <- NULL

make_pk_sim_given_data <- function(inDataset){
  if (! 'WT' %in% colnames(inDataset)){
    inDataset$WT <- 3.5
  }
  
  if (! 'SEX' %in% colnames(inDataset)){
    inDataset$SEX <- NA
  }
  
  inDataset$WTMED <- inDataset$WT - 3.5
  inDataset$WTALOM <- inDataset$WT / 3.5
  
  pk_fit <- nlmixr2(final_model_unfixed, data = inDataset,  est="posthoc", table=tableControl(keep = c('SEX', 'WT')))
  
  fit_params <- tibble(
    tka = pk_fit$parFixedDf$Estimate[[1]],
    tcl = pk_fit$parFixedDf$Estimate[[2]],
    tv = pk_fit$parFixedDf$Estimate[[3]],
    covka = pk_fit$parFixedDf$Estimate[[4]],
    allcl = pk_fit$parFixedDf$Estimate[[5]],
    allv = pk_fit$parFixedDf$Estimate[[6]],
    prop.sd = pk_fit$parFixedDf$Estimate[[7]],
    WTMED = inDataset$WTMED[[1]],
    WTALOM = inDataset$WTALOM[[1]],
    eta.ka = pk_fit$eta.ka[[1]],
    eta.cl = pk_fit$eta.cl[[1]],
    eta.v = pk_fit$eta.v[[1]]
  )

  
  ev  <- et(amountUnits="mg", timeUnits="hours")
  
  for (n in 1:nrow(inDataset %>% filter(EVID==1))){
    row <- inDataset %>% filter(EVID==1)
    ev <- ev %>% 
      et(amt=row$AMT[[n]], cmt="depot", time=row$TIME[[n]])
      #et(amt=row$AMT[[n]], cmt="depot", time=row$TIME[[n]])
  }
  
  last_dose <- max((inDataset %>% filter(EVID==1))$TIME)
  
  ev <- ev %>% 
    et(seq(0,(last_dose+48),(1/60)))#ceiling(max(inDataset$TIME)/100)*100
  
  
  
  sim <- rxSolve(one.compartment.proportional, ev, 
                 params=fit_params)
  return(sim)
  
  
  
}

make_pk_plot_given_data <- function(inDataset){
  pk_simulation_given_data <<- make_pk_sim_given_data(inDataset)
  return(
    ggplot(data=pk_simulation_given_data,aes(x=time,y=cp))+
      geom_line(linewidth=1) +
      scale_y_log10() + 
      geom_point(color="red", size=2, data=inDataset, aes(x=TIME,y=DV))
  )
}


make_pk_plot_threshold <- function(wt, dose, num_doses, interval, max_threshold){
  make_pk_plot(wt, dose, num_doses=num_doses, interval=interval) + geom_hline(yintercept = max_threshold, linetype="dashed")
}


#p4_avg_weight, p4_dosing_interval, p4_dose_amount, p4_auc_window
make_pk_plot_auc <- function(wt, dose, interval, windows){
  num_doses <- ceiling((14*24)/interval)
  make_pk_plot(wt, dose, num_doses, interval, final_time_override=(14+1)*24)
}

get_auc <- function(inData){
  sum <- 0
  
  for (i in 1:nrow(inData)){
    if (i == 1 | i ==nrow(inData)){
      sum = sum + inData$`median(cp)`[[i]]
    }
    else {
      sum = sum + 2*inData$`median(cp)`[[i]]
    }
  }
  
  return(.5 * sum * 1) #*abs(inData$time[[0]]-inData$time[[1]])
}

make_auc_table <- function(wt, dose, interval, windows){
  windows <- windows
  num_doses <- ceiling((14*24)/interval)
  sim <- make_pk_sim(wt, dose, num_doses, interval, final_time_override=0)
  
  sim_medians <- as.data.frame(sim %>% group_by(time) %>% summarize(median(cp)))
  
  auc_df <- data.frame()
  
  for (i in 1:ceiling((14*24)/windows)){
    temp_auc <- get_auc(sim_medians %>% filter(time<i*windows & time>=(i-1)*windows))
    temp_interval <- paste("AUC(",(i-1)*windows,",",(i)*windows,")",sep="")
    auc_df <- rbind(auc_df, data.frame(Interval=temp_interval, Value=temp_auc))
    #auc_df$Interval[nrow(auc_df),] <- 1
    #auc_df$Value[i] <- get_auc(sim %>% filter(time<i*windows & time>=(i-1)*windows))
      
    #   get_auc(
    #   sim %>% filter(
    #     time < i*windows &
    #     time >= (i-1)*windows
    #   )
    # )
  }
  print(auc_df)
  return(auc_df)
}





ui <- page_navbar(
  theme = bs_theme(version = 5, bootswatch = "sandstone"),
  bg = "#335d88",
  title = "Moxi NHP PK Simulator",
  
  nav_panel(
    title="PK Simulation by Weight",
    sidebarLayout(
      sidebarPanel(
        # Input: Slider for the number of bins ----
        sliderInput(
          inputId = "p1_dose",
          label = "Dose (mg/kg)",
          min = 0,
          max = 200,
          value = 80,
          step = 10
        ),
        sliderInput(
          inputId = "p1_female_wt",
          label = "Female Average Weight",
          min = 0,
          max = 10,
          value = 3.15,
          step = .05
        ),
        sliderInput(
          inputId = "p1_male_wt",
          label = "Male Average Weight",
          min = 0,
          max = 10,
          value = 4.45,
          step = .05
        ),
        submitButton(text = "Run Model", icon = NULL, width = NULL)
      ),
      mainPanel(
        card(
          card_header("Female PK"),
          plotOutput(outputId = "p1_pkPlotFemale")
        ),
        card(
          card_header("Male PK"),
          plotOutput(outputId = "p1_pkPlotMale")
        )
      )
    ),
    
  ),
  nav_panel(
    title="PK Simulation Given Sampling",
    sidebarLayout(
      sidebarPanel(
        fileInput("p2_pkData", "Upload Sampling Data"),
        downloadButton("p2_simulationDownload", "Download Simulation")
      ),
      mainPanel(
        card(
          card_header("Simulation Graph"),
          plotOutput(outputId = "p2_pkPlotDataSimulation")
        )
      )
    )
  ),
  nav_panel(
    title="PK Simulation With Threshold",
    sidebarLayout(
      sidebarPanel(
        numericInput(
          inputId = "p3_threshold_weight",
          label = "Average Weight",
          min = 0,
          max = 20,
          value = 3.5
        ),
        numericInput(
          inputId = "p3_max_concentration_threshold",
          label = "Maximum Threshold",
          min = 0,
          max = 100000,
          value = 1800
        ),
        numericInput(
          inputId = "p3_num_doses",
          label = "Number of Doses",
          min = 1,
          value = 2
        ),
        numericInput(
          inputId = "p3_dosing_interval",
          label = "Dosing Interval (h)",
          value = 24
        ),
        numericInput(
          inputId = "p3_dose_amt",
          label = "Dose Amount (mg/kg)",
          min = 0,
          max = 500,
          value = 20
        ),
        submitButton(text = "Run Model", icon = NULL, width = NULL)
      ),
      mainPanel(
        card(
          card_header(""),
          plotOutput(outputId = "p3_pkPlotThreshold")
        )
      )
    )
  ), 
  nav_panel(
    title="PK AUC Simulation",
    sidebarLayout(
      sidebarPanel(
        numericInput(
          inputId = "p4_avg_weight",
          label = "Average Weight",
          min = 0,
          max = 20,
          value = 3.5
        ),
        numericInput(
          inputId = "p4_dosing_interval",
          label = "Dosing Interval",
          min = 0,
          value=24
        ),
        numericInput(
          inputId = "p4_dose_amount",
          label = "Dosing Amount",
          min = 0,
          value=20
        ),
        radioButtons(
          inputId = "p4_auc_window",
          label = "AUC Window Length (h)",
          choices=c(
            "2 Hours" = 2,
            "3 Hours" = 3,
            "4 Hours" = 4,
            "6 Hours" = 6,
            "12 Hours" = 12,
            "24 Hours" = 24
          )
        ),
        
        # numericInput(
        #   inputId = "p4_auc_window",
        #   label = "AUC Window Length",
        #   min = 0,
        #   value=2
        # ),
        submitButton(text = "Run Model", icon = NULL, width = NULL),
        downloadButton("p4_aucDownload", "Download AUC Values")
      ),
      mainPanel(
        card(
          card_header(""),
          plotOutput(outputId = "p4_pkPlotAUC")
        ),
        card(
          DT::dataTableOutput("p4_aucTable")
        )
      )
    )
  )
)

  
  
  

server <- function(input, output) {
  output$p1_pkPlotFemale <- renderPlot({
    make_pk_plot(input$p1_female_wt, input$p1_dose)
  })
  
  output$p1_pkPlotMale <- renderPlot({
    make_pk_plot(input$p1_male_wt, input$p1_dose)
  })
  
  output$p2_pkPlotDataSimulation <- renderPlot({
    print(input$p2_pkData$datapath)
    
    make_pk_plot_given_data(read_csv(input$p2_pkData$datapath))
  })
  
  output$p2_simulationDownload <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_PK_simulation", ".csv", sep="")
    },
    content = function(file) {
      write.csv(pk_simulation_given_data, file)
    }
  )
  
  output$p3_pkPlotThreshold <- renderPlot({
    make_pk_plot_threshold(input$p3_threshold_weight, input$p3_dose_amt, input$p3_num_doses, input$p3_dosing_interval, input$p3_max_concentration_threshold)
  })
  
  
  #p4_avg_weight, p4_dosing_interval, p4_dose_amount, p4_auc_window
  output$p4_pkPlotAUC <- renderPlot({
    make_pk_plot_auc(input$p4_avg_weight, input$p4_dose_amount, input$p4_dosing_interval, as.numeric(input$p4_auc_window))
  })
  
  output$p4_aucTable <- DT::renderDataTable(DT::datatable({
    make_auc_table(input$p4_avg_weight, input$p4_dose_amount, input$p4_dosing_interval, as.numeric(input$p4_auc_window))
  }))
  
  output$p4_aucDownload <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_PK_AUC_simulation", ".csv", sep="")
    },
    content = function(file) {
      write.csv(make_auc_table(input$p4_avg_weight, input$p4_dose_amount, input$p4_dosing_interval, as.numeric(input$p4_auc_window)), file)
    }
  )
  
  output$pkPlotOld <- renderPlot({
    print("hi1")
    
    final_model <- function(){
      ini({
        tka <- fix(log(0.7630404)); label("Ka")
        tcl <- fix(log(3.2363909)); label("Cl")
        tv <- fix(log(33.9763030)); label("V")
        
        covka <- fix(-0.2138061)
        allcl <- fix(1.0050231)
        allv <- fix(1.0000000)
        
        eta.ka ~ fix(0.2676694)
        eta.cl + eta.v ~ fix(c(0.04400711, 0.02083013, 0.01501074))
        
        prop.sd <- fix(0.2756416)
      })
      
      model({
        WTMED <- 3.5-3.5
        WTALOM <- 3.5/3.5
        
        ka <- exp(tka + eta.ka + (covka*WTMED))
        
        cl <- exp(tcl + eta.cl + allcl*log(WTALOM))
        v <- exp(tv + eta.v + allv*log(WTALOM))
        
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - (cl / v) * center
        
        cp <- center / v
        cp ~ prop(prop.sd)
      })
    }
    
    final_model <- rxode(final_model)
    final_model <- final_model$simulationModel
    
    final_omegas <- lotri::lotri({
      eta.ka ~ 0.2676694
      eta.cl + eta.v ~ c(0.04400711, 0.02083013, 0.01501074)
    })
    
    ev  <- et(amountUnits="mg", timeUnits="hours") %>% 
      et(amt=400000, cmt="depot") %>%
      et(0:100) %>%
      et(id=1)
    
    
    
    params_wt <- tibble(
      WTMED <- 3.5-3.5,
      WTALOM <- 3.5/3.5
    )
    
    final_model <- final_model %>% 
      ini(tka = (log(1))) %>% 
      ini(tcl = (log(3))) %>% 
      ini(tv = (log(30)))
    
    
    sim <- rxSolve(final_model, ev, omega=final_omegas, params = params_wt, nSub=100)
    ggplot(data=sim,aes(x=time,y=cp, group=sim.id))+
      geom_line(alpha=.2, linewidth=1) +
      scale_y_log10()
    
    
    
    #print("Hi")
    #plot(sim$time,sim$cp)
    #model_file_name <- paste("../",model_name,".RData",sep="")
    #load(model_file_name)
    
    #x    <- faithful$waiting
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #inData<- as.data.frame(fit_used %>% filter(ID==as.data.frame(id_list)$ID[[input$id]]))
    #graph_scatterplot(inData, xAxis='TIME', yAxis='DV', xLab='TIME', yLab='DV', point=TRUE, log10_y = TRUE) +
    #  geom_line(aes(y=IPRED), show.legend=TRUE, color=darkgreen, linewidth=1, linetype=1) +
    #  geom_line(aes(y=PRED), color=lightred, linewidth=1, linetype=5) 
    
    #graph_dv_vs_time(fit_used,5,as.data.frame(id_list)$ID[[input$id]])
    
  })
  
}


shinyApp(ui = ui, server = server)
