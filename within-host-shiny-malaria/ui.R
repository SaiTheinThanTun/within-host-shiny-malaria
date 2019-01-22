library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Within Host Drug Action"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("initn", "Initial total parasites: ", min = 10, max=1000, value=10),
      sliderInput("mu", "Mean of initial age distribution: ", min = 1, max=48, value = 1),
      sliderInput("sig", "SD of initial age distribution: ", min = 1, max=48, value = 1),
      sliderInput("pmf", "Parasite Multiplication factor: ", min = 1, max=30, value=1),
      sliderInput("h", "h: ", min = 1, max = 10, value =1, step =1),
      sliderInput("initconc", "Initial drug A concentration: ", min = 25, max=150, value = 25, step = 5 ),
      sliderInput("halflife", "Drug A halflife: ", min = 5, max=240, value=5, step = 1),
      sliderInput("killrate", "Kill rate, drug A: ", min = .001, max = 1, value = .05, step = .05),
      sliderInput("ce50", "CE50, drug A: ", min = 10, max = 100, value = 10, step = 5),
      
      sliderInput("initconc_2", "Initial drug B concentration: ", min = 25, max=150, value = 50, step = 5 ),
      sliderInput("halflife_2", "Drug B halflife: ", min = 5, max=240, value=20, step = 1),
      sliderInput("killrate_2", "Kill rate, drug B: ", min = .001, max = 1, value = .001, step = .05),
      sliderInput("ce50_2", "CE50, drug B: ", min = 10, max = 100, value = 30, step = 5)
    ),
    mainPanel(
      plotOutput("combinedPlot"),
      plotOutput("drugeffPlot")
    )
  )
))


##initn=slider(10,1000),
##mu=slider(1,48),
##sig=slider(1,48),
##pmf=slider(1,30),
# k0=slider(0,1,step = 0.01),
# a=slider(0.0001,0.01,step = 0.001),
# tpar=slider(1000,100000,step = 10),
# delay=slider(0,4,step=1),
##initconc=slider(25,150,step=5),
##halflife=slider(5,240,step=1),
##killrate=slider(0.001,1,step=0.05),
##ce50=slider(10,100,step=5),
##h=slider(1,10,step=1),

