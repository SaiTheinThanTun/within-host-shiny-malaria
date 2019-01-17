library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Within Host Drug Action"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("initn", "Initial total parasites: ", min = 10, max=1000, value=10),
      sliderInput("mu", "Mean of initial age distribution: ", min = 1, max=48, value = 1),
      sliderInput("sig", "SD of initial age distribution: ", min = 1, max=48, value = 1),
      sliderInput("pmf", "Parasite Multiplication factor: ", min = 1, max=30, value=1),
      sliderInput("initconc", "Initial drug concentration: ", min = 25, max=150, value = 25, step = 5 ),
      sliderInput("halflife", "Drug halflife: ", min = 5, max=240, value=5, step = 1),
      sliderInput("killrate", "Kill rate: ", min = .001, max = 1, value = .001, step = .05),
      sliderInput("ce50", "CE50: ", min = 10, max = 100, value = 10, step = 5),
      sliderInput("h", "h: ", min = 1, max = 10, value =1, step =1)
    ),
    mainPanel(
      plotOutput("NJWimPlot"),
      plotOutput("drugfPlot"),
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

