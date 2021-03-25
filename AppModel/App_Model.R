#################################################################################
#################################################################################
#                                                                               #
#  Shiny App for the model of the titled manuscript "Lockdown, relaxation, and  #
#  ACME period in COVID-19: A study of disease dynamics on Hermosillo, Sonora,  #
#  Mexico." This app was used to manually fit the model solutions to the        #
#  initially reported data.         
#  https://doi.org/10.1371/journal.pone.0242957
#                                                                               #
#################################################################################
#################################################################################

#########################
#        PACKAGES       #
#########################
library(shiny)
library(deSolve)
library(extraDistr)

rm(list=ls())

###################################
#         COVID-19 MODEL          #
###################################

FunctionModel<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)),{
    if(Time<StartTimePropP)
    {
      valpha_ps<-0
      valpha_pa<-0
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0
    }
    if(Time>=StartTimePropP & Time<TimeGetPropP)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega20<-0
    }
    if(Time>=TimeGetPropP & Time<ReturnStartTime)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0  
    }
    if(Time>=ReturnStartTime & Time<TimeGetPropR)
    {
      vomega10<-0
    }
    if(Time>TimeGetPropR)
    {
      vomega10<-0
      vomega20<-0
    }
    N  <- S + E + Ia + Is + R + P + Pr
    dS <- -valpha_s*S*Is/N - valpha_a*S*Ia/N - vomega10 * S
    dE <- valpha_s*S*Is/N + valpha_a*S*Ia/N + valpha_ps*P*Is/N + valpha_pa*P*Ia/N + valpha_pls*Pr*Is/N + valpha_pla*Pr*Ia/N- vdelta*E
    dIa<- (1-vtheta)*vdelta*E - veta_a*Ia
    dIs<- vtheta*vdelta*E - vgamma_s*Is
    dH <- vbeta*vgamma_s*Is - vmu*H + vtau*vpsi*Q
    dD <- vpi*vmu*H
    dQ <- (1-vbeta)*vgamma_s*Is - vpsi*Q
    dR <- veta_a*Ia + (1-vpi)*vmu*H + (1-vtau)*vpsi*Q
    dP <- vomega10 * S + - valpha_ps*P*Is/N - valpha_pa*P*Ia/N - vomega20*P
    dPr <- vomega20*P - valpha_pls*Pr*Is/N - valpha_pla*Pr*Ia/N
    dIc <- vtheta*vdelta*E 
    dHc <- vbeta*vgamma_s*Is 
    dEc <- vdelta*E 
    dQc <- (1-vbeta)*vgamma_s*Is
    list(c(dS,dE,dIa,dIs,dH,dD,dQ,dR,dP,dPr,dIc,dHc,dEc,dQc))
  })
}

##############################################
#         FUNCTION TO SOLVE THE MODEL        #
##############################################

odesystem <- function(FinalDate,par){
  
  # DATA
  data<-read.csv("DataHillo200719.csv")[1:82,]
  data_obs<-cbind(data$DEATHS,data$INFEC_SYMP,data$HOSPITALIZED,data$AMBULATORY)
  
  # TIME GRID
  InitialTime<-data[1,2]
  InitialDate<-as.Date(data[1,1])
  FinalDate<-as.Date(FinalDate)
  diffDates<-as.numeric(FinalDate-InitialDate)
  FinalTime<-InitialTime+diffDates
  dt<-0.1
  TimeGrid<-seq(InitialTime,FinalTime,dt)
  
  # INITIAL CONDITIONS
  ValS<-930668 
  ValE<-0          
  ValIa<-0         
  ValIs<-1         
  ValH<-0
  ValD<-0
  ValQ<-0
  ValR<-0
  ValP<-0
  ValPr<-0
  ValIc<-ValIs
  ValHc<-0
  ValEc<-ValIs 
  ValQc<-ValQ  
  N<-ValS+ValE+ValIa+ValIs+ValH+ValD+ValQ+ValR+ValP+ValPr
  yini<-c(S=ValS,E=ValE,Ia=ValIa,Is=ValIs,H=ValH,D=ValD,Q=ValQ,R=ValR,P=ValP,
          Pr=ValPr,Ic=ValIc,Hc=ValHc,Ec=ValEc,Qc=ValQc) 
  
  # CONTROL PARAMETERS 
  ValStartTimePropP<-5  
  ValPeriodGetPropP<-30     
  ValReturnStartTime<-50 
  ValPeriodGetPropPr<-15 
  ValTimeGetPropP<-ValStartTimePropP+ValPeriodGetPropP
  ValTimeGetPropPr<-ValReturnStartTime+ValPeriodGetPropPr
  valomega10<- -(1/ValPeriodGetPropP)*log(1-par[16])
  valomega20<- -(1/ValPeriodGetPropPr)*log(1-par[17])
  
  # SYSTEM SOLUTION
  vecobsFinal<-seq(InitialTime,FinalTime,1)
  vecpar<-c(valpha_s=par[1],valpha_a=par[2],valpha_ps=par[3],valpha_pa=par[4],
            valpha_pls=par[5],valpha_pla=par[6],vdelta=par[7],vtheta=par[8],
            veta_a=par[9],vgamma_s=par[10],vbeta=par[11],vmu=par[12],vpi=par[13],
            vpsi=par[14],vtau=par[15],vomega10=valomega10,vomega20=valomega20,
            StartTimePropP=ValStartTimePropP,
            TimeGetPropP=ValTimeGetPropP,
            ReturnStartTime=ValReturnStartTime,
            TimeGetPropR=ValTimeGetPropPr)
  sol<-ode(yini,TimeGrid,FunctionModel,vecpar)
  Acum<-sol[,c(7,12,13,15)][sol[,1] %in% vecobsFinal,]
  NewCases_temp<-Acum[2:(length(vecobsFinal)),]-Acum[1:(length(vecobsFinal)-1),] 
  NewCases<-rbind(Acum[1,],NewCases_temp)
  NewCasesDc<-NewCases[,1]
  NewCasesIc<-NewCases[,2]
  NewCasesHc<-NewCases[,3]
  NewCasesQc<-NewCases[,4]

  # GRAPHICS   
  color=brewer.pal(n = 5, name = "Set1")
  par(mfrow= c(2,3))
  plot(vecobsFinal,NewCasesIc,type="l",col=color[1],lwd=7,xaxt="n",
       main="Daily New Cases",ylab="",xlab="",cex.main=2.5,cex.lab=2,cex.axis=1.7,
       col.axis="black",col.lab="blue",col.main=1)
  points(data$DAY,data$INFEC_SYMP,type="b",col=1,pch=15)
  axis(1, at = c(0,40,81,120,160,200),cex.axis=1.7,col.axis="black",
       labels=c("20/03/11","20/04/20","20/05/31","20/07/09","20/08/18","20/09/27"))
  grid()
  
  plot(vecobsFinal,NewCasesHc,type="l",col=color[2],lwd=7,xaxt="n",
       main="Daily New Hospitalizations",ylab="",xlab="",cex.main=2.5,cex.lab=2,cex.axis=1.7,
       col.axis="black",col.lab="blue",col.main=1)
  points(data$DAY,data$HOSPITALIZED,type="b",col=1,pch=15)
  axis(1, at = c(0,40,81,120,160,200),cex.axis=1.7,col.axis="black",
       labels=c("20/03/11","20/04/20","20/05/31","20/07/09","20/08/18","20/09/27"))
  grid()
  
  plot(vecobsFinal,NewCasesDc,type="l",col=color[3],lwd=7,xaxt="n",
       main="Daily New Deaths",ylab="",xlab="",cex.main=2.5,cex.lab=2,cex.axis=1.7,
       col.axis="black",col.lab="blue",col.main=1)
  points(data$DAY,data$DEATHS,type="b",col=1,pch=15)
  axis(1, at = c(0,40,81,120,160,200),cex.axis=1.7,col.axis="black",
       labels=c("20/03/11","20/04/20","20/05/31","20/07/09","20/08/18","20/09/27"))
  grid()

  plot(vecobsFinal,NewCasesQc,type="l",col=color[4],lwd=7,xaxt="n",
       main="Daily New Ambulatories",ylab="",xlab="",cex.main=2.5,cex.lab=2,cex.axis=1.7,
       col.axis="black",col.lab="blue",col.main=1)
  points(data$DAY,data$AMBULATORY,type="b",col=1,pch=15)
  axis(1, at = c(0,40,81,120,160,200),cex.axis=1.7,col.axis="black",
       labels=c("20/03/11","20/04/20","20/05/31","20/07/09","20/08/18","20/09/27"))
  grid()
  
  plot(sol[,1],sol[,14],type="l",col=color[5],lwd=7,xaxt="n",
       main="Cumulative Cases",ylab="",xlab="",cex.main=2.5,cex.lab=2,cex.axis=1.7,
       col.axis="black",col.lab="blue",col.main=1)
  axis(1, at = c(0,40,81,120,160,200),cex.axis=1.7,col.axis="black",
       labels=c("20/03/11","20/04/20","20/05/31","20/07/09","20/08/18","20/09/27"))
  grid()
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "topleft",inset = 0,
         legend = c("Data"),pch=15,
         col="black", cex=3, horiz = TRUE)
  return()
}

#####################################
#         DEFINE UI FOR APP         #
#####################################

ui <- fluidPage(
  theme = shinytheme("journal"),
  
  # Application title
  titlePanel("COVID-19 Model"),
  
  # Sidebar with a slider input for parameters
  fluidRow(
    column(2,
           # Incubation rate
           sliderInput("delta", step=0.01,
                       HTML("&delta;"),  
                       min = 0.07,
                       max = 1/2,
                       value = 1/5),
           br(),
           # Proportion of symptomatic individuals
           sliderInput("theta", step=0.01,
                       HTML("&theta;"),
                       min = 0,
                       max = 0.8,
                       value = 0.2),
           br(),
           # Transmission contact rate for susceptible class of symptomatic 
           # individuals
           sliderInput("alpha_s", step=0.001,
                       HTML("&alpha;<sub>s"),
                       min = 0.5,
                       max = 1.7,
                       value = 0.657)
    ),
    column(2,
           # Recovery rate for asymptomatic individuals
           sliderInput("eta_a", step=0.01,
                       HTML("&eta;<sub>a"),  
                       min = 0.06,
                       max = 0.13,
                       value = 1/10.5),
           br(),
           # Output rate from the symptomatic class by register
           sliderInput("gamma_s", step=0.01,
                       HTML("&gamma;<sub>s"), 
                       min = 0.12,
                       max = 1.25,
                       value = 1/3),
           br(),
           # Transmission contact rate for susceptible class of asymptomatic 
           # individuals
           sliderInput("alpha_a", step=0.001,
                       HTML("&alpha;<sub>a"),
                       min = 0.5,
                       max = 1.7,
                       value = 1.198)
    ),
    column(2,
           # Proportion of hospitalized individuals
           sliderInput("beta", step=0.01,
                       HTML("&beta;"),  
                       min = 0.1,
                       max = 0.3,
                       value = 0.126),
           br(),
           # Output rate from the hospitalized class by recovery/death
           sliderInput("mu", step=0.01,
                       HTML("&mu;"), 
                       min = 0.04,
                       max = 0.22,
                       value = 0.06),
           br(),
           # Transmission contact rate for susceptible protected class of 
           # symptomatic individuals
           sliderInput("alpha_ps", step=0.001,
                       HTML("&alpha;<sub>ps"),
                       min = 0,
                       max = 0.5,
                       value = 0.02)
    ),
    column(2,
           # Ï and (1âÏ) represent the proportion of hospitalized individuals that
           #recover or die, respectively
           sliderInput("pi", step=0.01,
                       HTML("&pi;"),
                       min = 0,
                       max = 0.4,
                       value = 0.3),
           br(),
           # Output rate from the quarantined class by hospitalization/recovery
           sliderInput("psi", step=0.01,
                       HTML("&psi;"), 
                       min = 0.06,
                       max = 0.13,
                       value = 0.5),
           br(),
           # Transmission contact rate for susceptible protected class of 
           # asymptomatic individuals
           sliderInput("alpha_pa", step=0.001,
                       HTML("&alpha;<sub>pa"),
                       min = 0,
                       max = 0.5,
                       value = 0.02)     
    ),
    column(2,
           # Ï and (1âÏ) are the ambulatory individuals proportions who are 
           # hospitalized and recovered, respectively
           sliderInput("tau", step=0.01,
                       HTML("&tau;"),  
                       min = 0,
                       max = 0.5,
                       value = 0.2),
           br(),
           # Final Time to graph
           sliderInput("TFinal", step=5,
                       HTML("Final Date"), 
                       min = as.Date("2020-05-31","%Y-%m-%d"),
                       max = as.Date("2020-09-27","%Y-%m-%d"),
                       value=as.Date("2020-09-27","%Y-%m-%d"),
                       timeFormat="%Y-%m-%d"),
           br(),
           # Transmission contact rate for susceptible protected released class 
           # of symptomatic individuals
           sliderInput("alpha_pls", step=0.001,
                       HTML("&alpha;<sub>pls"), 
                       min = 0,
                       max = 0.5,
                       value = 0.02)
    ),
    column(2,
           # Proportion of Protected 
           sliderInput("PropP", step=0.01,
                       "Proportion of Protected", 
                       min = 0.7,
                       max = 0.9,
                       value = 0.8),
           br(),
           # Proportion of Protected released
           sliderInput("PropRegreso", step=0.01,
                       "Proportion of Protected released",  
                       min = 0.1,
                       max = 0.35,
                       value = 0.2),
           br(),
           # Transmission contact rate for susceptible protected released class 
           # of asymptomatic individuals
           sliderInput("alpha_pla", step=0.001,
                       HTML("&alpha;<sub>pla"), 
                       min = 0,
                       max = 0.5,
                       value = 0.02)      
    ),
    # Show a plot of the generated distribution
    mainPanel(plotOutput("Is", height = "500px"), width = 25)
  )
)

#################################
#         DEFINE SERVER         #
#################################

server <- function(input, output) 
{
  
  output$Is <- renderPlot({
    FinalDate = as.Date(input$TFinal)
    par = c(input$alpha_s,input$alpha_a,input$alpha_ps,input$alpha_pa,
            input$alpha_pls,input$alpha_pla,input$delta,input$theta,
            input$eta_a,input$gamma_s,input$beta,input$mu,input$pi,
            input$psi,input$tau,input$PropP,input$PropRegreso)
    res = odesystem(FinalDate,par)
  })
}

#####################################
#        RUN THE APPLICATION        #
#####################################

shinyApp(ui = ui, server = server)



