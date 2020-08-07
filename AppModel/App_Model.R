#################
#   LIBRERIAS   #
#################
library(shiny)
library(deSolve)
library(extraDistr)
#################
#   FUNCIONES   #
#################

rm(list=ls())

# Graphical parameters
par(cex.axis = 3, cex.axis = 3, cex.sub = 3, cex.lab = 3, cex.main = 3, lwd=3)

FuncionModelo<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)),{
    # TiempoLograPropP<-TiempoInicioPropP+DuracionLograrPropP
    # TiempoLograrPropRegreso<-TiempoInicioRegreso+DuracionLograrRegreso
    # vomega10<- -(1/DuracionLograrPropP)*log(PropP)
    # vomega20<- -(1/DuracionLograrRegreso)*log(PropRegreso)
    if(Time<TiempoInicioPropP)
    {
      valpha_ps<-0
      valpha_pa<-0
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0
    }
    if(Time>=TiempoInicioPropP & Time<TiempoLograrPropP)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega20<-0
    }
    if(Time>=TiempoLograrPropP & Time<TiempoInicioRegreso)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0  
    }
    if(Time>=TiempoInicioRegreso & Time<TiempoLograrPropRegreso)
    {
      vomega10<-0
    }
    if(Time>TiempoLograrPropRegreso)
    {
      vomega10<-0
      vomega20<-0
    }
    N  <- S + E + Ia + Is + R + P + Pl
    dS <- -valpha_s*S*Is/N - valpha_a*S*Ia/N - vomega10 * S
    dE <- valpha_s*S*Is/N + valpha_a*S*Ia/N + valpha_ps*P*Is/N + valpha_pa*P*Ia/N + valpha_pls*Pl*Is/N + valpha_pla*Pl*Ia/N- vdelta*E
    dIa<- (1-vtheta)*vdelta*E - veta_a*Ia
    dIs<- vtheta*vdelta*E - vgamma_s*Is
    dH <- vbeta*vgamma_s*Is - vmu*H + vtau*vpsi*Q
    dD <- vpi*vmu*H
    dQ <- (1-vbeta)*vgamma_s*Is - vpsi*Q
    dR <- veta_a*Ia + (1-vpi)*vmu*H + (1-vtau)*vpsi*Q
    dP <- vomega10 * S + - valpha_ps*P*Is/N - valpha_pa*P*Ia/N - vomega20*P
    dPl <- vomega20*P - valpha_pls*Pl*Is/N - valpha_pla*Pl*Ia/N
    dIc <- vtheta*vdelta*E
    dHc <- vbeta*vgamma_s*Is 
    dEc <- vdelta*E
    dQc <- (1-vbeta)*vgamma_s*Is
    list(c(dS,dE,dIa,dIs,dH,dD,dQ,dR,dP,dPl,dIc,dHc,dEc,dQc))
  })
}

odesystem <- function(par){
  data<-read.csv("DatosHillo200719.csv")[1:82,]
  TInicial<-0
  #TFinal<-par[22]
  dt<-0.1
  RejillaTiempo<-seq(TInicial,par[22],dt)
  data_obs<-cbind(data$MUERTES,data$INFEC_SINTO,data$HOSPITALIZADOS,data$AMBULATORIOS)
  
  ValS<-930668 
  ValE<-0          
  ValIa<-0         
  ValIs<-1         
  ValH<-0
  ValD<-0
  ValQ<-0
  ValR<-0
  ValP<-0
  ValPl<-0
  ValIc<-ValIs
  ValHc<-0
  ValEc<-ValIs #Acumulado de Infectados Asintomaticos y Sintomaticos
  ValQc<-ValQ  #Acumulado de Ambulatorios
  N<-ValS+ValE+ValIa+ValIs+ValH+ValD+ValQ+ValR+ValP+ValPl
  yini<-c(S=ValS,E=ValE,Ia=ValIa,Is=ValIs,H=ValH,D=ValD,Q=ValQ,R=ValR,P=ValP,
          Pl=ValPl,Ic=ValIc,Hc=ValHc,Ec=ValEc,Qc=ValQc) 
  
  NuevosCasosDc<-c()
  NuevosCasosIc<-c()
  NuevosCasosHc<-c()
  NuevosCasosQc<-c()
  
  vecobsFinal<-seq(TInicial,par[22],1)
  valomega10<- -(1/par[19])*log(par[16])
  valomega20<- -(1/par[21])*log(par[17])
  ValTiempoLograrPropP<-par[18]+par[19]
  ValTiempoLograrPropRegreso<-par[20]+par[21]
  
  vecpar<-c(valpha_s=par[1],valpha_a=par[2],valpha_ps=par[3],valpha_pa=par[4],
            valpha_pls=par[5],valpha_pla=par[6],vdelta=par[7],vtheta=par[8],
            veta_a=par[9],vgamma_s=par[10],vbeta=par[11],vmu=par[12],vpi=par[13],
            vpsi=par[14],vtau=par[15],vomega10=valomega10,vomega20=valomega20,
            TiempoInicioPropP=par[18],TiempoLograrPropP=ValTiempoLograrPropP,
            TiempoInicioRegreso=par[20],TiempoLograrPropRegreso=ValTiempoLograrPropRegreso)
  sol<-ode(yini,RejillaTiempo,FuncionModelo,vecpar)
  Acumulado<-sol[,c(7,12,13,15)][sol[,1] %in% vecobsFinal,]
  NuevosCasos_temp<-Acumulado[2:(length(vecobsFinal)),]-Acumulado[1:(length(vecobsFinal)-1),] 
  NuevosCasos<-rbind(Acumulado[1,],NuevosCasos_temp)
  NuevosCasosDc<-NuevosCasos[,1]
  NuevosCasosIc<-NuevosCasos[,2]
  NuevosCasosHc<-NuevosCasos[,3]
  NuevosCasosQc<-NuevosCasos[,4]
  SumaError<-(NuevosCasos[1:82,]-data_obs)^2
  Error<-c(sum(SumaError[,1]),sum(SumaError[,2]),sum(SumaError[,3]),sum(SumaError[,4]))
  
  par(mfrow= c(2,3))
  plot(vecobsFinal,NuevosCasosIc,type="l",col="red",lwd=2,main="Nuevos Infectados Sintomáticos",cex=1.2)
  points(data$DIA,data$INFEC_SINTO,type="b",col=1,pch=15)
  plot(vecobsFinal,NuevosCasosHc,type="l",col="red",lwd=2,main="Nuevas Hospitalizados",cex=1.2)
  points(data$DIA,data$HOSPITALIZADOS,type="b",col=1,pch=15)
  plot(vecobsFinal,NuevosCasosDc,type="l",col="red",lwd=2,main="Nuevas Muertes",cex=1.2)
  points(data$DIA,data$MUERTES,type="b",col=1,pch=15)
  plot(vecobsFinal,NuevosCasosQc,type="l",col="red",lwd=2,main="Ambulatorios",cex=1.2)
  points(data$DIA,data$AMBULATORIOS,type="b",col=1,pch=15)
  
  plot(sol[,1],sol[,14],type="l",col="red",lwd=2,main="Acumulado Infectados",cex=1.2)
  abline(h=0.216*N, col="blue")

  return()
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Modelo para Covid-19"),
  
  # Sidebar with a slider input for parameters
  fluidRow(
    column(2,
           sliderInput("delta", step=0.001,
                       "delta",
                       min = 0,
                       max = 0.5,
                       value = 0.2),
           br(),
           sliderInput("theta", step=0.001,
                       "theta",
                       min = 0.17,
                       max = 0.25,
                       value = 0.20)      
    ),
    column(2,
           sliderInput("eta_a",
                       "eta_a",  step=0.001,
                       min = 0.001,
                       max = 0.2,
                       value = 1/10.5),
           br(),
           sliderInput("gamma_s",
                       "gamma_s", step=0.001,
                       min = 0.1,
                       max = 2,
                       value = 1/3)
    ),
    column(2,
           sliderInput("beta",
                       "beta",  step=0.001,
                       min = 0.01,
                       max = 0.4,
                       value = 0.126),
           br(),
           sliderInput("mu",
                       "mu", step=0.001,
                       min = 0.05,
                       max = 0.1,
                       value = 0.06)
    ),
    column(2,
           sliderInput("pi",
                       "pi",  step=0.001,
                       min = 0.2,
                       max = 0.4,
                       value = 0.3),
           br(),
           sliderInput("psi",
                       "psi", step=0.01,
                       min = 0.06,
                       max = 0.1,
                       value = 0.5)
    ),
    column(2,
           sliderInput("tau",
                       "tau",  step=0.001,
                       min = 0.1,
                       max = 0.3,
                       value = 0.2),
           br(),
           sliderInput("PropP",
                       "PropP", step=0.001,
                       min = 0.1,
                       max = 0.3,
                       value = 0.2)
    ),
    column(2,
           sliderInput("PropRegreso",
                       "PropRegreso",  step=0.01,
                       min = 0.65,
                       max = 0.9,
                       value = 0.8),
           br(),
           sliderInput("TiempoInicioPropP",
                       "TiempoInicioPropP", step=1,
                       min = 1,
                       max = 10,
                       value = 5)
    ),
    column(2,
           sliderInput("DuracionLograrPropP",
                       "DurLograrPropP",  step=1,
                       min = 30,
                       max = 40,
                       value = 30),
           br(),
           sliderInput("TiempoInicioRegreso",
                       "TiempoInicioRegreso", step=1,
                       min = 45,
                       max = 55,
                       value = 50)
    ),
    column(2,
           sliderInput("DuracionLograrRegreso",
                       "DurLograrRegreso",  step=1,
                       min = 10,
                       max = 30,
                       value = 15),
           br(),
           sliderInput("TFinal",
                       "TFinal", step=1,
                       min = 1,
                       max = 200,
                       value = 200)
    ),
    column(2,
           sliderInput("alpha_s", step=0.001,
                       "alpha_s",
                       min = 0,
                       max = 2,
                       value = 0.657),
           br(),
           sliderInput("alpha_a", step=0.001,
                       "alpha_a",
                       min = 0,
                       max = 2,
                       value = 1.116)
    ),
    column(2,           
           sliderInput("alpha_ps", step=0.001,
                       "alpha_ps",
                       min = 0,
                       max = 0.5,
                       value = 0.025),
           br(),
           sliderInput("alpha_pa", step=0.001,
                       "alpha_pa",
                       min = 0,
                       max = 0.5,
                       value = 0.021)      
    ),
    column(2,
           sliderInput("alpha_pls", step=0.001,
                       "alpha_pls",
                       min = 0,
                       max = 0.5,
                       value = 0.105),
           br(),
           sliderInput("alpha_pla", step=0.001,
                       "alpha_pla",
                       min = 0,
                       max = 0.5,
                       value = 0.044)      
    ),
    # Show a plot of the generated distribution
    mainPanel(plotOutput("Is", height = "400px"), width = 15)
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) 
{
  
  output$Is <- renderPlot({
    par = c(input$alpha_s,input$alpha_a,input$alpha_ps,input$alpha_pa,
            input$alpha_pls,input$alpha_pla,input$delta,input$theta,
            input$eta_a,input$gamma_s,input$beta,input$mu,input$pi,
            input$psi,input$tau,input$PropP,input$PropRegreso,
            input$TiempoInicioPropP,input$DuracionLograrPropP,     
            input$TiempoInicioRegreso,input$DuracionLograrRegreso,
            input$TFinal)
    res = odesystem(par)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)