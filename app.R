#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(PlanBsmooth)

#--------------------------------------------------------------------------------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#--------------------------------------------------------------------------------
############################################################################################
ApplyFSD <- function(surveys,          # matrix with Year in first column and surveys in remaining columns
                     npoint = 5,       # number of points in regression for slope (default = 5)
                     termyr = NA,      # terminal year to use when computing slope (default = last year)
                     Kp     = 0.75,    # gain for first derivative (default = 0.75)
                     Kd     = 0.50     # gain for second derivative (default = 0.50)
){
  
  # select data to use
  if(is.na(termyr)) termyr <- max(surveys$Year, na.rm=TRUE)
  dat.use <- filter(surveys, Year <= termyr) %>%
    drop_na()  # removes years with missing index values
  nyears <- max(dat.use$Year) - min(dat.use$Year) + 1 
  
  # get average slope from surveys
  nsurveys <- length(surveys[1,]) - 1
  sdat <- filter(dat.use, Year >= (termyr - npoint + 1)) 
  sdat1 <- filter(dat.use, Year >= (termyr - npoint), Year <= (termyr - 1))
  slopeterm <- rep(NA, nsurveys)  # slope for terminal year period for each survey
  slopeterm1 <- rep(NA, nsurveys) # slope for period prior to terminal year for each survey
  for (i in 1:nsurveys){
    mylm <- lm(log(sdat[,(i+1)]) ~ sdat[,1])
    mylm1 <- lm(log(sdat1[,(i+1)]) ~ sdat1[,1])
    slopeterm[i] <- as.numeric(coefficients(mylm)[2])
    slopeterm1[i] <- as.numeric(coefficients(mylm1)[2])
  }
  avgslope <- mean(slopeterm, na.rm=TRUE)
  deltaslope <- avgslope - mean(slopeterm1, na.rm=TRUE)
  
  # compute multiplier
  multiplier <- exp(Kp * avgslope + Kd * deltaslope)
  
  # make results list
  res <- list()
  res$avgslope <- avgslope
  res$deltaslope <- deltaslope
  res$multiplier <- multiplier
  
  return(res)
}
############################################################################################
# some constants for now (make them variable later)
R0 <- 1000
M <- 0.2
nages <- 20
WAA <- 0.000001 * (100*(1 - exp(-0.2 * seq(1,nages)))) ^ 3
maturity <- c(0,0.2,0.5,0.8,rep(1, nages-4))
selx <- c(0.1,0.2,0.5,0.8,0.9,rep(1, nages-5))
nbaseyears <- 35
minimumF <- 0.01

# calculate F.table for use in MSY reference points
nsteps <- 2001
F.table <- data.frame(Fval = double(),
                      ypr = double(),
                      spr = double()  )
for (istep in 1:nsteps){
  Fval <- (istep - 1)/1000
  yprval <- 0.0
  sprval <- 0.0
  Nval <- 1.0
  for (i in 1:nages){
    Zval <- Fval * selx[i] + M
    yprval <- yprval + Nval * WAA[i] * Fval * selx[i] * (1 - exp(-Zval)) / Zval
    sprval <- sprval + Nval * WAA[i] * maturity[i]
    Nval <- Nval * exp(-Zval)
  }
  Nval <- Nval / (1 - exp(-(Fval * selx[nages] + M)))
  yprval <- yprval + Nval * WAA[nages] * Fval * selx[nages] * (1 - exp(-Zval)) / Zval
  sprval <- sprval + Nval * WAA[nages] * maturity[nages]
  F.row <- data.frame(Fval = Fval,
                      ypr = yprval,
                      spr = sprval)
  F.table <- rbind(F.table, F.row)
}
F.table
spr0 <- filter(F.table, Fval == 0)$spr 

#--------------------------------------------------------------------------------
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Simulation Test PlanBsmooth (red) and FSD (blue)"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(6,
                 sliderInput("steepness",
                             label = "Steepness",
                             min = 0.21,
                             max = 1.0,
                             value = 0.7),
                 
                 sliderInput("Fyr1",
                             label = "F in first year",
                             min = 0,
                             max = 1.0,
                             value = 0.1),  
                 
                 sliderInput("Fyr20",
                             label = "F in year 20",
                             min = 0,
                             max = 1.0,
                             value = 0.6),  
                 
                 sliderInput("Fyrbreak",
                             label = "Year change",
                             min = 21,
                             max = 35,
                             step = 1,
                             value = 30),
                 
                 sliderInput("Ffinal",
                             label = "Final F",
                             min = 0,
                             max = 1.0,
                             value = 0.25),
                 
                 sliderInput("numprjyrs",
                             label = "# Projection Years",
                             min = 3,
                             max = 50,
                             step = 1,
                             value = 10)
                 ),
          column(6,
                 sliderInput("sigmaF",
                             label = "sigmaF",
                             min = 0,
                             max = 1,
                             value = 0.2),  
                 
                 sliderInput("sigmaR",
                             label = "sigmaR",
                             min = 0,
                             max = 1,
                             value = 0.6),
                 
                 sliderInput("sigmaS1",
                             label = "sigmaS1",
                             min = 0,
                             max = 1,
                             value = 0.2),
                 
                 sliderInput("sigmaS2",
                             label = "sigmaS2",
                             min = 0,
                             max = 1,
                             value = 0.4),
                 
                 selectInput("catch_advice_mult",
                             label="Catch Advice Multiplier",
                             choices = list("1.0"   = 1,
                                            "0.75"  = 2),
                             selected = 1),
                 
                 selectInput("assess_frequency",
                             label = "Assessment Frequency",
                             choices = list("Every Year" = 1,
                                            "Alternating Years" = 2,
                                            "Every Third Year" = 3),
                             selected = 1)
                 
                 )
        )
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("myPlots")
      )
   )
)

#--------------------------------------------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) {

   output$myPlots <- renderPlot({

     # colors
     PB_color = "tomato"
     FSD_color = "blue"
     
     # get total number of years
     nprojyears <- input$numprjyrs
     ntotyears <- nbaseyears + nprojyears
     
     # create the deviate streams
     F_devs <- rnorm(nbaseyears)
     N1_devs <- rnorm(nages)
     R_devs <- rnorm(ntotyears)
     S1_devs <- rnorm(ntotyears)
     S2_devs <- rnorm(ntotyears)
       
     # Stock-Recruitment Plot     
     x <- seq(0, R0*spr0, length.out = 1000)
     sr_alpha <- 4 * input$steepness * R0 / (5 * input$steepness - 1)
     sr_beta  <- R0 * spr0 * (1 - input$steepness) / (5 * input$steepness - 1)
     y <- sr_alpha * x / (sr_beta + x)
     srdf <- data.frame(SSB = x, 
                        Recruitment = y)
     plot1 <- ggplot(srdf, aes(x=SSB, y=Recruitment)) +
       geom_line() +
       theme_bw()

     # get reference points
     F.table.1 <- mutate(F.table, SSB = sr_alpha * spr - sr_beta)
     F.table.2 <- mutate(F.table.1, Rval = sr_alpha * SSB / (sr_beta + SSB))
     F.table.full <- mutate(F.table.2, Yield = ypr * Rval)
     ref.pts <- filter(F.table.full, Yield == max(Yield))

     # Fmultiplier during base years plot
     Fmult_base <- rep(NA, nbaseyears)
     Fmult_base[1:20] <- seq(input$Fyr1, input$Fyr20, length.out = 20)
     Fmult_base[21:input$Fyrbreak] <- seq(input$Fyr20, input$Ffinal, length.out = input$Fyrbreak - 20)
     Fmult_base[input$Fyrbreak:35] <- rep(input$Ffinal, 35 - input$Fyrbreak + 1)
     Fmult_applied <- Fmult_base * exp(F_devs * input$sigmaF - 0.5 * input$sigmaF * input$sigmaF)
     FAA <- matrix(NA, nrow=ntotyears, ncol=nages)
     FAA[1:nbaseyears,] <- outer(Fmult_applied, selx)
     ZAA <- FAA + M
     Fdf <- data.frame(Year = 1:nbaseyears,
                       Fmult = Fmult_applied)
     plot2 <- ggplot(Fdf, aes(x=Year, y=Fmult)) +
       geom_point() +
       ylim(0, NA) +
       theme_bw()
     
     # fill NAA matrix for base period
     N1 <- rep(1, nages)
     for (i in 2:nages){
       N1[i] <- N1[i-1] * exp(-ZAA[1,i-1])
     }
     N1[nages] <- N1[nages] / (1 - exp(-ZAA[1,nages]))
     
     NAA <- matrix(NA, nrow=ntotyears, ncol=nages)
     R1 <- filter(F.table.full, Fval == input$Fyr1)
     NAA[1,] <- R1$Rval * N1 * exp(N1_devs * input$sigmaR - 0.5 * input$sigmaR * input$sigmaR)
     NAA[1,nages] <- R1$Rval * N1[nages] 
     SSB <- rep(NA, ntotyears)
     SSB[1] <- sum(NAA[1,] * maturity * WAA)
     
     for (j in 2:nbaseyears){
       pred_R <- sr_alpha * SSB[j-1] / (sr_beta + SSB[j-1])
       NAA[j,1] <- pred_R * exp(R_devs[j] * input$sigmaR - 0.5 * input$sigmaR * input$sigmaR)
       for (i in 2:nages){
         NAA[j,i] <- NAA[j-1,i-1] * exp(-ZAA[j-1,i-1])
       }
       NAA[j,nages] <- NAA[j,nages] + NAA[j-1,nages] * exp(-ZAA[j-1,nages])
       SSB[j] <- sum(NAA[j,] * maturity * WAA)
     }
     
     # add SR points to the SR curve plot
     sr_obs <- data.frame(SSB = SSB[1:(nbaseyears-1)],
                          Recruitment = NAA[2:nbaseyears, 1])
     plot1 <- plot1 + geom_point(data=sr_obs, aes(x=SSB, y=Recruitment))
       
     # calculate catch in numbers and yield in weight
     CAA <- NAA * FAA * (1 - exp(-ZAA)) / ZAA
     YAA <- CAA * outer(rep(1, ntotyears), WAA)
     Yield <- apply(YAA, 1, sum)
     
     Ydf <- data.frame(Year = 1:nbaseyears,
                       Yield = Yield[1:nbaseyears])
     plot3 <- ggplot(Ydf, aes(x=Year, y=Yield)) +
       geom_point() +
       ylim(0, NA) +
       theme_bw()
     
     # make the surveys
     popB <- apply(NAA * outer(rep(1, ntotyears), WAA), 1, sum)
     survey1 <-  popB * exp(S1_devs * input$sigmaS1 - 0.5 * input$sigmaS1 * input$sigmaS1)
     survey2 <-  popB * exp(S2_devs * input$sigmaS2 - 0.5 * input$sigmaS2 * input$sigmaS2)

     surveys <- data.frame(Year = 1:nbaseyears,
                           survey1 = survey1[1:nbaseyears],
                           survey2 = survey2[1:nbaseyears])

     avgS <- data.frame(Year = 1:nbaseyears,
                        avg = 0.5 * (survey1[1:nbaseyears] + survey2[1:nbaseyears]))
     
     # duplicate starting point for FSD
     NAA_FSD <- NAA
     ZAA_FSD <- ZAA
     SSB_FSD <- SSB
     Fmult_applied_FSD <- Fmult_applied
     Yield_FSD <- Yield

     #----------------------------------------------------------------
     # loop through projection years
     PBmult <- rep(NA, ntotyears)
     catch_advice <- rep(NA, ntotyears)
     FSDmult <- rep(NA, ntotyears)
     catch_advice_FSD <- rep(NA, ntotyears)
     
     for(iyear in (nbaseyears + 1):ntotyears){
       
       # determine whether this is an assessment year or not
       runassess <- FALSE
       iyeartest <- iyear - nbaseyears - 1
       aaa <- input$assess_frequency
       if((iyeartest %% as.numeric(input$assess_frequency)) == 0) runassess <- TRUE
       
       if(runassess == TRUE){
         # get PlanBsmooth multiplier
         PBres <- ApplyPlanBsmooth(dat           = avgS,
                                   od            = ".\\",
                                   my.title      = "PlanBsmooth fit",
                                   terminal.year = iyear - 2,
                                   nyears        = 33,
                                   loess.span    = NA,
                                   saveplots     = FALSE,
                                   showplots     = FALSE)
         PBmult[iyear] <- PBres$multiplier
         if(is.na(PBmult[iyear])) PBmult[iyear] <- 1 # is this the correct default for PlanBsmooth crashing?
         
         # Calculate recent mean catch
         rec_mean_catch <- mean(Yield[(iyear-4):(iyear-2)], na.rm=TRUE)
         
         # Determine catch advice
         catch_advice_multiplier <- 1.0
         if(input$catch_advice_mult == 2) catch_advice_multiplier <- 0.75
         catch_advice[iyear] <- PBmult[iyear] * rec_mean_catch * catch_advice_multiplier
         
         # get FSD multiplier and catch advice
         FSD <- ApplyFSD(surveys = surveys,
                         npoint  = 5,
                         termyr  = NA,
                         Kp      = 0.75,
                         Kd      = 0.50)
         FSDmult[iyear] <- FSD$multiplier
         catch_advice_FSD[iyear] <- FSDmult[iyear] * Yield[(iyear-2)]
       }
       
       if(runassess == FALSE){
         catch_advice[iyear] <- catch_advice[iyear - 1]
         catch_advice_FSD[iyear] <- catch_advice_FSD[iyear - 1]
       }

#--------------------------------------------------------------
       # bring forward population to this year for PlanBsmooth 
       pred_R <- sr_alpha * SSB[iyear-1] / (sr_beta + SSB[iyear-1])
       NAA[iyear,1] <- pred_R * exp(R_devs[iyear] * input$sigmaR - 0.5 * input$sigmaR * input$sigmaR)
       for (i in 2:nages){
         NAA[iyear,i] <- NAA[iyear-1,i-1] * exp(-ZAA[iyear-1,i-1])
       }
       NAA[iyear,nages] <- NAA[iyear,nages] + NAA[iyear-1,nages] * exp(-ZAA[iyear-1,nages])
       SSB[iyear] <- sum(NAA[iyear,] * maturity * WAA)
       
       # determine F for this year
       Nvec <- NAA[iyear,]
       ctarget <- catch_advice[iyear]
       if(ctarget <= 1e-6 ){
         Fmult <- 0
       }
       if(ctarget > 1e-6) {
         Fmax <- 2
         Fvec <- Fmax * selx
         Zvec <- Fvec + M
         y <- sum(Nvec * WAA * Fvec * (1-exp(-Zvec)) / Zvec)
         if(y <= ctarget){
           Fmult <- Fmax
         }
         if (y >ctarget){
           a = 0
           b = Fmax
           for (itry in 1:25){
             Fval <- (a + b) / 2
             Fvec <- Fval * selx
             Zvec <- Fvec + M
             y <- sum(Nvec * WAA * Fvec * (1-exp(-Zvec)) / Zvec)
             if (y <= ctarget){
               a <- Fval
             }
             if (y > ctarget){
               b <- Fval  
             }
           }    
           Fmult <- (a + b) / 2
         }
       }
       Fmult_applied[iyear] <- max(Fmult, minimumF)
       Faa <- Fmult_applied[iyear] * selx
       ZAA[iyear,] <- Faa + M
       
       # calculate catch (should be equal to catch_advice)
       Caa <- NAA[iyear,] * Faa * (1 - exp(-ZAA[iyear,])) / ZAA[iyear,]
       Yaa <- Caa * WAA
       Yield[iyear] <- sum(Yaa)

       # make surveys
       popB <- sum(NAA[iyear,] * WAA)
       survey1[iyear] <-  popB * exp(S1_devs[iyear] * input$sigmaS1 - 0.5 * input$sigmaS1 * input$sigmaS1)
       survey2[iyear] <-  popB * exp(S2_devs[iyear] * input$sigmaS2 - 0.5 * input$sigmaS2 * input$sigmaS2)
       avgS_year <- data.frame(Year = iyear,
                               avg = 0.5 * (survey1[iyear] + survey2[iyear]))
       avgS <- rbind(avgS, avgS_year)
       
#--------------------------------------------------------------
       # bring forward population to this year for FSD 
       pred_R <- sr_alpha * SSB_FSD[iyear-1] / (sr_beta + SSB_FSD[iyear-1])
       NAA_FSD[iyear,1] <- pred_R * exp(R_devs[iyear] * input$sigmaR - 0.5 * input$sigmaR * input$sigmaR)
       for (i in 2:nages){
         NAA_FSD[iyear,i] <- NAA_FSD[iyear-1,i-1] * exp(-ZAA_FSD[iyear-1,i-1])
       }
       NAA_FSD[iyear,nages] <- NAA_FSD[iyear,nages] + 
         NAA_FSD[iyear-1,nages] * exp(-ZAA_FSD[iyear-1,nages])
       SSB_FSD[iyear] <- sum(NAA_FSD[iyear,] * maturity * WAA)
       
       # determine F for this year
       Nvec <- NAA_FSD[iyear,]
       ctarget <- catch_advice_FSD[iyear]
       if(ctarget <= 1e-6 ){
         Fmult <- 0
       }
       if(ctarget > 1e-6) {
         Fmax <- 2
         Fvec <- Fmax * selx
         Zvec <- Fvec + M
         y <- sum(Nvec * WAA * Fvec * (1-exp(-Zvec)) / Zvec)
         if(y <= ctarget){
           Fmult <- Fmax
         }
         if (y >ctarget){
           a = 0
           b = Fmax
           for (itry in 1:25){
             Fval <- (a + b) / 2
             Fvec <- Fval * selx
             Zvec <- Fvec + M
             y <- sum(Nvec * WAA * Fvec * (1-exp(-Zvec)) / Zvec)
             if (y <= ctarget){
               a <- Fval
             }
             if (y > ctarget){
               b <- Fval  
             }
           }    
           Fmult <- (a + b) / 2
         }
       }
       Fmult_applied_FSD[iyear] <- max(Fmult, minimumF)
       Faa <- Fmult_applied_FSD[iyear] * selx
       ZAA_FSD[iyear,] <- Faa + M
       
       # calculate catch (should be equal to catch_advice)
       Caa <- NAA_FSD[iyear,] * Faa * (1 - exp(-ZAA_FSD[iyear,])) / ZAA_FSD[iyear,]
       Yaa <- Caa * WAA
       Yield_FSD[iyear] <- sum(Yaa)
       
       # make surveys
       popB <- sum(NAA_FSD[iyear,] * WAA)
       survey1[iyear] <-  popB * exp(S1_devs[iyear] * input$sigmaS1 - 0.5 * input$sigmaS1 * input$sigmaS1)
       survey2[iyear] <-  popB * exp(S2_devs[iyear] * input$sigmaS2 - 0.5 * input$sigmaS2 * input$sigmaS2)
       surveys_year <- data.frame(Year = iyear,
                                  survey1 = survey1[iyear],
                                  survey2 = survey2[iyear])
       surveys <- rbind(surveys, surveys_year)

     } # end of projection loop
     #----------------------------------------------------------------
     
     # summarize results and add projection values to plots
     
     # add SR points to the SR curve plot
     sr_obs_proj <- data.frame(SSB = SSB[nbaseyears:(ntotyears - 1)],
                               Recruitment = NAA[(nbaseyears + 1):ntotyears, 1])
     sr_obs_proj_FSD <- data.frame(SSB = SSB_FSD[nbaseyears:(ntotyears - 1)],
                                   Recruitment = NAA_FSD[(nbaseyears + 1):ntotyears, 1])
     plot1 <- plot1 + geom_point(data=sr_obs_proj, aes(x=SSB, y=Recruitment), color=PB_color)
     plot1 <- plot1 + geom_point(data=sr_obs_proj_FSD, aes(x=SSB, y=Recruitment), color=FSD_color)

     # add Fmult projections
     Fdf_proj <- data.frame(Year = (nbaseyears + 1):ntotyears,
                       Fmult = Fmult_applied[(nbaseyears + 1):ntotyears])
     Fdf_proj_FSD <- data.frame(Year = (nbaseyears + 1):ntotyears,
                                Fmult = Fmult_applied_FSD[(nbaseyears + 1):ntotyears])
     plot2 <- plot2 + geom_point(data=Fdf_proj, aes(x=Year, y=Fmult), color=PB_color)
     plot2 <- plot2 + geom_point(data=Fdf_proj_FSD, aes(x=Year, y=Fmult), color=FSD_color)

     # add Yield projections
     Ydf_proj <- data.frame(Year = (nbaseyears + 1):ntotyears,
                            Yield = Yield[(nbaseyears + 1):ntotyears])
     Ydf_proj_FSD <- data.frame(Year = (nbaseyears + 1):ntotyears,
                                Yield = Yield_FSD[(nbaseyears + 1):ntotyears])
     plot3 <- plot3 + geom_point(data=Ydf_proj, aes(x=Year, y=Yield), color=PB_color)
     plot3 <- plot3 + geom_point(data=Ydf_proj_FSD, aes(x=Year, y=Yield), color=FSD_color)

     # add Survey projections
     plot4 <- ggplot(filter(avgS, Year <= nbaseyears), aes(x=Year, y=avg)) +
       geom_point() +
       geom_point(data=filter(avgS, Year > nbaseyears), aes(x=Year, y=avg), color=PB_color) +
       labs(y = "Average Survey") +
       ylim(0, NA) +
       theme_bw()
     avgS_FSD <- filter(surveys, Year > nbaseyears) %>%
       mutate(avg = (0.5 * (survey1 + survey2)))
     plot4 <- plot4 + geom_point(data=avgS_FSD, aes(x=Year, y=avg), color=FSD_color)
     
     # get the final PlanBsmooth fit
     plot5 <- PBres$tsplot
     
     # plot multipliers
     PBmultdf <- data.frame(Year = seq(1,ntotyears),
                            Multiplier = PBmult)
     FSDmultdf <- data.frame(Year = seq(1,ntotyears),
                             Multiplier = FSDmult)
     plot6 <- ggplot(PBmultdf, aes(x=Year, y=Multiplier)) +
       geom_line(color=PB_color, na.rm=TRUE) +
       geom_point(color=PB_color, na.rm=TRUE) +
       geom_hline(yintercept=1, color="black", linetype="dashed") +
       ylim(0, NA) +
       theme_bw()
     plot6 <- plot6 + geom_line(data=FSDmultdf, aes(x=Year, y=Multiplier), color=FSD_color, na.rm=TRUE)
     plot6 <- plot6 + geom_point(data=FSDmultdf, aes(x=Year, y=Multiplier), color=FSD_color, na.rm=TRUE)

     # add reference points to the plots
     ref.pt.color <- "dark green"
     ref.pt.type <- "solid"
     plot1 <- plot1 +
       geom_abline(slope=1/(ref.pts$spr[1]), intercept=0, color=ref.pt.color, linetype=ref.pt.type)
     
     plot2 <- plot2 +
       geom_hline(yintercept=ref.pts$Fval[1], color=ref.pt.color, linetype=ref.pt.type)
     
     plot3 <- plot3 +
       geom_hline(yintercept=ref.pts$Yield[1], color=ref.pt.color, linetype=ref.pt.type)
     
     plot4 <- plot4 +
       geom_hline(yintercept=ref.pts$SSB[1], color=ref.pt.color, linetype=ref.pt.type)
     
     # make the final combined plot
     multiplot(plot1, plot2, plot5, plot3, plot4, plot6, cols=2)

  }, height = 800, width = 600)
}

# Run the application 
shinyApp(ui = ui, server = server)

