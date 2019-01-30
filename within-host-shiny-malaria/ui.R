library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Within Host Drug Action"),
  
  sidebarLayout(
    sidebarPanel(
      #sliderInput("initn", "Initial total parasites: ", min = 1000, max=10^10, value=100),
      sliderInput("initn", "Initial total parasites (10^x): ", min = 3, max=10, value=8.93, step=.01),
      sliderInput("mu", "Mean of initial age distribution: ", min = 1, max=48, value = 28),
      sliderInput("sig", "SD of initial age distribution: ", min = 1, max=48, value = 7),
      sliderInput("pmf", "Parasite Multiplication factor: ", min = 8, max=10, value=1),
      sliderInput("h", "h: ", min = 1, max = 10, value =4, step =1),
      sliderInput("initconc", "Initial drug A concentration: ", min = 25, max=150, value = 72, step = 5 ),
      sliderInput("halflife", "Drug A halflife: ", min = 5, max=240, value=54, step = 1),
      sliderInput("killrate", "Kill rate, drug A: ", min = .05, max = .4, value = .2, step = .05),
      sliderInput("sen", "Sensitivity, drug A: ", min = .00, max = 1.00, value = 1, step = .01),
      sliderInput("ce50", "CE50, drug A: ", min = 10, max = 100, value = 15, step = 5),
      
      sliderInput("initconc_2", "Initial drug B concentration: ", min = 25, max=150, value = 50, step = 5 ),
      sliderInput("halflife_2", "Drug B halflife: ", min = 5, max=400, value=324, step = 1),
      sliderInput("killrate_2", "Kill rate, drug B: ", min = .05, max = .4, value = .1, step = .05),
      sliderInput("sen_2", "Sensitivity, drug B: ", min = .00, max = 1.00, value = 1, step = .01),
      sliderInput("ce50_2", "CE50, drug B: ", min = 10, max = 100, value = 30, step = 5)
    ),
    mainPanel(
      plotOutput("combinedPlot"),
      plotOutput("drugeffPlot")
    )
  )
))


