#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

requiredPackages <- c("unix","dplyr","reshape2","MASS","stringr","R2jags","MCMCvis","ggplot2","shiny","shinybusy")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(requiredPackages)

source('App_functions.R', local = TRUE)

if (interactive()) {
  
  ui <- fluidPage(
    
    titlePanel("ShinyWhale"),
    
    navlistPanel(
      "Getting Started",
      tabPanel("Before you start",
               helpText("This app is designed to provide simplified tool for estimating abundance and demographic information (survival, recapture, site fidelity and dead recovery)
                  using a Bayesian multi-state mark-recapture-recovery framework. With the use of this app, we aim to provide managers, NGO's, government organisations and 
                  conservationists a tool to assess the changes in population trends, without the need to understand the complete workings of Bayesian mark-recapture techniques,
                  or strong programming skills. More information on this model can be found in the paper is available in the Github ripository.
                  
                  This app was designed for North Atlantic right whales, and as such, to be use with data provided by the North Atlantic Right Whale Consortium, 
                  and uses their data structure to conduct the model.
                  The columns required for this app are as below, please ensure that the column headings match"),
               strong("Please ensure your data is set out as below:"),
               tableOutput("test.data"),
               strong("Before you start you must download JAGS (Just Another Gibbs Sampler)."),
               uiOutput("jags"),
               helpText("")),
      
      "Data",
      tabPanel("Import Data",
               helpText("File must be a csv. and in the structure shown in the Getting Started - input data tab."),
               fileInput("file1", "Choose CSV File",
                         accept = c(
                           "text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")
               ),
               dataTableOutput("df"),
      ),
      tabPanel("Capture Histories",
               helpText("Select the first year you want to model from."),
               numericInput("start.year", 
                            label = "Start Year:",
                            min = 1977, max = 3000, value = 1990, step = 1),
               helpText("Select the last year you want to model, or the last year of complete data."),
               numericInput("end.year", 
                            label = "End Year:",
                            min = 1980, max = 3000, value = 2000, step = 1),
               helpText("Select the month to be used as the start of the year"),
               selectInput("month", "Starting month", choices = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")),
               helpText("Enter the number of additional individuals you want to include for data augmentation. 
                         This is the number of individuals who could enter the population, but have never been observed."),
               numericInput("data.aug",
                            label = "Number of individuals for data augmentation",
                            min = 0, max = 1000, value = 300, step = 1),
               actionButton("button",
                            strong("Create capture histories")),
               dataTableOutput("whale_CH")
               
      ),
      "Model",
      tabPanel("Model settings",
               helpText("Select the settings for the MCMC (Markov Chain Monte Carlo) simulations to run the model"),
               helpText("The number of draws (iterations) per chain"),
               numericInput("ni", 
                            label = "Number of iterations:",
                            min = 2000, max = 100000, value = 6000, step = 100),
               helpText("How many iterations to keep (1 keeps every iteration, 10 keeps every 10, etc)"),
               numericInput("nt", 
                            label = "Amount of thinning:",
                            min = 1, max = 100, value = 1, step = 1),
               helpText("The number of iterations discarded before the model to converge 
                       (needs to be less than Number of iterations)"),
               numericInput("nb", 
                            label = "Length of burn-in:",
                            min = 500, max = 80000, value = 3000, step = 100),
               helpText("The number of times the model is run with unique starting values"),
               numericInput("nc", 
                            label = "Number of chains:",
                            min = 3, max = 10, value = 3, step = 1),
               actionButton("button2",
                            strong("Run Model")),
               verbatimTextOutput("mod"),
               add_busy_bar(color = "red", height = "8px")
               
      ),
      "Output",
      tabPanel("Abundance",
               splitLayout(
                 dataTableOutput("Ntable"),
                 plotOutput("Nplot")
               ),
               downloadButton("downloadData_abund", "Download the data"),
               downloadButton("save_abund", "save plot")
      ),
      tabPanel("Survival",
               splitLayout(
                 dataTableOutput("Stable"),
                 plotOutput("Splot")
               ),
               downloadButton("downloadData_surv", "Download the data"),
               downloadButton("save_surv", "save plot")
      ),
      tabPanel("Recapture",
               splitLayout(
                 dataTableOutput("Ptable"),
                 plotOutput("Pplot")
               ),
               downloadButton("downloadData_recp", "Download the data"),
               downloadButton("save_recp", "save plot")
      ),
      tabPanel("Site Fidelity",
               splitLayout(
                 dataTableOutput("Ftable"),
                 plotOutput("Fplot")
               ),
               downloadButton("downloadData_SF", "Download the data"),
               downloadButton("save_SF", "save plot")
      ),
      tabPanel("Dead Recovery",
               splitLayout(
                 dataTableOutput("Rtable"),
                 plotOutput("Rplot")
               ),
               downloadButton("downloadData_Drec", "Download the data"),
               downloadButton("save_Drec", "save plot")
      )
    )
  )
  
  server <- function(input, output) {
    
    output$test.data <- renderTable({
      tab <- data.frame("SightingEGNo" = c("###", "###", "###", "###", "..."),
                        "SightingYear" = c("year seen", "year seen", "year seen", "year seen", "..."),
                        "SightingMonth" = c("1-12", "1-12", "1-12", "1-12", "..."),
                        "SightingDay" = c("1-31", "1-31", "1-31", "1-31", "..."), 
                        "Behaviors" = c("###", "###", "###", "###", "..."))
      
      return(tab)
    })
    
    url <- a("JAGS", href="https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/")
    output$jags <- renderUI({
      tagList("available here:", url)
    })
    #inFile: is to import the csv of the sightings data used
    #dataf: is the output of this function
    inFile <- reactive({
      if(is.null(input$file1)) {
        return("")
      }
      read.csv(file = input$file1$datapath)
    })
    output$df <- renderDataTable({
      req(inFile())
      
      dataf <- inFile()
      return(dataf)
    })      
    
    Month <- reactive({
      switch(input$month,
             "January" = 1,
             "February" = 2,
             "March" = 3,
             "April" = 4,
             "May" = 5,
             "June" = 6,
             "July" = 7,
             "August" = 8,
             "September" = 9,
             "October" = 10,
             "November" = 11,
             "December" = 12)
    })
    
    #CH: takes the dataset that has been entered and assigns it to a vector to be used by other functions
    CH <- reactive({
      if(is.null(input$button)) {
        return("")
      }
      wh_ch <- whale_data(Cap_H = inFile(), month = Month())
      return(wh_ch)
    })
    
    #capt_hist: creates the capture histories for the individuals in the dataset
    capt_hist <- reactive({
      if(is.null(input$button)) {
        return("")
      }
      req(CH())
      ch <- CH()
      ch <- subset(ch, ch$Year >= input$start.year & ch$Year <= input$end.year)
      ch <- ch(ch)
      return(ch)
    })
    
    Aug <-  eventReactive(input$button, {
      data.aug <- aug(capt_hist(), input$data.aug)
      return(data.aug)
    })   
    
    Aug_data <-  eventReactive(input$button, {
      req(Aug())
      req(capt_hist())
      full_capt_hist <- rbind(capt_hist(), Aug())
      return(full_capt_hist)
    })  
    
    Aug_z.data <- eventReactive(input$button, {
      req(Aug())
      full_z_hist <- Aug()
      full_z_hist[full_z_hist == 3] <- 1
      for (i in 1:dim(full_z_hist)[1]){
        full_z_hist[i,1] <- NA                #converts the first time step to NA, as this state is defined in the model, and the model will not work if there is a number assigned to this occasion.
      }
      return(full_z_hist)
    })  
    
    #whale_CH: prints the capture histories for the dataset of the selected time periods and the number of data augmented individuals
    output$whale_CH <- renderDataTable({
      
      return(Aug_data())
    })
    
    #NYears: calculates the number of years within the capture histories to use in the initial values function. 
    
    model <- eventReactive(input$button2, {
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = (end.year - start.year) + 1
      jags.data = list(y = Aug_data(), n.occasions = dim(Aug_data())[2], nind = dim(Aug_data())[1])
      iv = jags_data(capt_hist())
      iv_2 = rbind(iv, Aug_z.data())
      inits = function(){list(mean.s = runif(years, 0, 1), mean.f = runif(years, 0, 1), mean.p = runif(years, 0, 1), mean.r = runif(years, 0, 1), z = iv_2)}  
      parameters = c("mean.s", "mean.f", "mean.r", "mean.p", "N")
      MRR_NARW = jags(jags.data, inits, parameters, "lifedead.jags", n.chains = input$nc, n.thin = input$nt, n.iter = input$ni, n.burnin = input$nb)
      return(MRR_NARW)
    })
    update_busy_bar(0) # update with the desire value [0-100], 100 hide the bar
    
    
    modo <- reactive({
      if(is.null(input$button2)) {
        return("") 
      }
      req(model())
      modout <- model()
      out = MCMCsummary(modout, digit = 3,
                        params = c("mean.s", "mean.f", "mean.r", "mean.p", "N"))
      out = as.data.frame(out)
      outp <- ifelse (max(out$Rhat) > 1.2,
                      print("Rerun model with more iterations for better convergence"),
                      
                      print("Model complete"))
      
      return(outp)  
    })
    
    output$mod <- renderPrint({
      return(modo())
    })
    
    Abund <- reactive({
      req(model())
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = seq.int(start.year, end.year, by = 1)
      modout <- model()
      N = MCMCsummary(modout, digit = 3,
                      params =  "N")
      N <- as.data.frame(N)
      N <- cbind(years, N)
      out = MCMCsummary(modout, digit = 3)
      out = as.data.frame(out)
      
      if (max(out$Rhat) > 1.2){
        shiny::showNotification("Re-run model with more iterations", type = "error")
        NULL
      } else {
        return(N) }
    })
    
    output$Ntable <- renderDataTable({
      return(Abund())
    })
    
    Nplot <- reactive({
      req(model())
      N <- Abund()
      plot <- ggplot(N, aes(x= years, y= mean)) + geom_point(aes(x=years, y= mean), size = 3)  + 
        geom_errorbar(aes(ymin= `2.5%`, 
                          ymax=`97.5%`), width=0.2, size = 1) +
        scale_y_continuous(name = "Abundance")+
        scale_x_continuous(name = "Year")+ 
        theme_classic()
      
      return(plot)
    })
    
    output$Nplot <- renderPlot({
      return(Nplot())
    })
    
    
    Surv <- reactive({
      req(model())
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = seq.int(start.year, end.year, by = 1)
      modout <- model()
      S = MCMCsummary(modout, digit = 3,
                      params =  "mean.s")
      S <- as.data.frame(S)
      S <- cbind(years, S)
      out = MCMCsummary(modout, digit = 3)
      out = as.data.frame(out)
      
      if (max(out$Rhat) > 1.2){
        shiny::showNotification("Re-run model with more iterations", type = "error")
        NULL
      } else {
        return(S) }
    })
    
    output$Stable <- renderDataTable({
      return(Surv())
    })
    
    Splot <- reactive({
      req(model())
      S <- Surv()
      plot <- ggplot(S, aes(x= years, y= mean)) + geom_point(aes(x=years, y= mean), size = 3)  + 
        geom_errorbar(aes(ymin= `2.5%`, 
                          ymax=`97.5%`), width=0.2, size = 1) +
        scale_y_continuous(name = "Survival probability")+
        scale_x_continuous(name = "Year")+ 
        theme_classic() 
      return(plot)
    })
    
    output$Splot <- renderPlot({
      return(Splot())
    })
    
    Rec <- reactive({
      req(model())
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = seq.int(start.year, end.year, by = 1)
      modout <- model()
      P = MCMCsummary(modout, digit = 3,
                      params =  "mean.p")
      P <- as.data.frame(P)
      P <- cbind(years, P)
      out = MCMCsummary(modout, digit = 3)
      out = as.data.frame(out)
      
      if (max(out$Rhat) > 1.2){
        shiny::showNotification("Re-run model with more iterations", type = "error")
        NULL
      } else {
        return(P) }
    })
    
    output$Ptable <- renderDataTable({
      return(Rec())
    })
    
    Pplot <- reactive({
      req(model())
      P <- Rec()
      plot <- ggplot(P, aes(x= years, y= mean)) + geom_point(aes(x=years, y= mean), size = 3)  + 
        geom_errorbar(aes(ymin= `2.5%`, 
                          ymax=`97.5%`), width=0.2, size = 1) +
        scale_y_continuous(name = "Recapture probability")+
        scale_x_continuous(name = "Year")+ 
        theme_classic()
      
      return(plot)
    })
    
    output$Pplot <- renderPlot({
      return(Pplot())
    })
    
    Fid <- reactive({
      req(model())
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = seq.int(start.year, end.year, by = 1)
      modout <- model()
      Fi = MCMCsummary(modout, digit = 3,
                       params =  "mean.f")
      Fi <- as.data.frame(Fi)
      Fi <- cbind(years, Fi)
      out = MCMCsummary(modout, digit = 3)
      out = as.data.frame(out)
      
      if (max(out$Rhat) > 1.2){
        shiny::showNotification("Re-run model with more iterations", type = "error")
        NULL
      } else {
        return(Fi) }
    })
    
    output$Ftable <- renderDataTable({
      return(Fid())
    })
    
    Fplot <- reactive({
      req(model())
      Fi <- Fid()
      plot <- ggplot(Fi, aes(x= years, y= mean)) + geom_point(aes(x=years, y= mean), size = 3)  + 
        geom_errorbar(aes(ymin= `2.5%`, 
                          ymax=`97.5%`), width=0.2, size = 1) +
        scale_y_continuous(name = "Site fidelity probability")+
        scale_x_continuous(name = "Year")+ 
        theme_classic() 
      
      return(plot)
    })
    
    output$Fplot <- renderPlot({
      return(Fplot())
    })
    
    DRec <- reactive({
      req(model())
      start.year = as.numeric(input$start.year)
      end.year = as.numeric(input$end.year)
      years = seq.int(start.year, end.year, by = 1)
      modout <- model()
      R = MCMCsummary(modout, digit = 3,
                      params =  "mean.r")
      R <- as.data.frame(R)
      R <- cbind(years, R)
      out = MCMCsummary(modout, digit = 3)
      out = as.data.frame(out)
      
      if (max(out$Rhat) > 1.2){
        shiny::showNotification("Re-run model with more iterations", type = "error")
        NULL
      } else {
        return(R) }
    })
    
    output$Rtable <- renderDataTable({
      return(DRec())
    })
    
    Rplot <- reactive({
      req(model())
      R <- DRec()
      plot <- ggplot(R, aes(x= years, y= mean)) + geom_point(aes(x=years, y= mean), size = 3)  + 
        geom_errorbar(aes(ymin= `2.5%`, 
                          ymax=`97.5%`), width=0.2, size = 1) +
        scale_y_continuous(name = "Dead recovery probability")+
        scale_x_continuous(name = "Year")+ 
        theme_classic() 
      
      return(plot)
    })
    
    output$Rplot <- renderPlot({
     return(Rplot())
    })
    #Abundance
    # Downloadable csv of selected dataset ----
    output$downloadData_abund <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Abundance_data_extraction", ".csv", sep = "")
        
      },
      content = function(file) {
        write.csv(Abund(), file, row.names = FALSE)
      }
    )
    
    output$save_abund <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Abundance_plot", ".tiff", sep = "")
        
      },
      content = function(file) {
        ggsave(Nplot(), filename = file)
        
      })
    
    #Survival
    output$downloadData_surv <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Survival_data_extraction", ".csv", sep = "")
        
      },
      content = function(file) {
        write.csv(Surv(), file, row.names = FALSE)
      }
    )
    
    output$save_surv <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Survival_plot", ".tiff", sep = "")
        
      },
      content = function(file) {
        ggsave(Splot(), filename = file)
        
      })
    
    #Recapture
    output$downloadData_recp <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Recapture_data_extraction", ".csv", sep = "")
        
      },
      content = function(file) {
        write.csv(Rec(), file, row.names = FALSE)
      }
    )
    
    output$save_recp <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Recapture_plot", ".tiff", sep = "")
        
      },
      content = function(file) {
        ggsave(Pplot(), filename = file)
        
      })
    
    #Site fidelity
    output$downloadData_SF <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Site_fidelity_data_extraction", ".csv", sep = "")
        
      },
      content = function(file) {
        write.csv(Fid(), file, row.names = FALSE)
      }
    )
    
    output$save_SF <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Site_fidelity_plot", ".tiff", sep = "")
        
      },
      content = function(file) {
        ggsave(Fplot(), filename = file)
        
      })
    
    #Dead recovery
    output$downloadData_Drec <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Dead_recovery_data_extraction", ".csv", sep = "")
        
      },
      content = function(file) {
        write.csv(DRec(), file, row.names = FALSE)
      }
    )
    
    
    output$save_Drec <- downloadHandler(
      filename = function() {
        ## paste(input$dataset, ".csv", sep = "")
        paste("Dead_recovery_plot", ".tiff", sep = "")
        
      },
      content = function(file) {
        ggsave(DRec(), filename = file)
        
      })
    
  }
  
  shinyApp(ui, server)
}

