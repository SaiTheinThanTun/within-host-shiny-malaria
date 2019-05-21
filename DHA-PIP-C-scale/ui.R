library(shiny)

shinyUI(fluidPage(
  
  titlePanel("Within Host Drug Action: DHA & piperaquine only"),
  
  # sidebarLayout(
  #   sidebarPanel(
  fluidRow(
      #sliderInput("initn", "Initial total parasites: ", min = 1000, max=10^10, value=100),
    column(3,
      sliderInput("initn", "Log10 of initial total parasites: ", min = 3, max=14, value=8.93, step=.01),
      sliderInput("mu", "Mean of initial age distribution: ", min = 1, max=48, value = 28),
      sliderInput("sig", "SD of initial age distribution: ", min = 1, max=48, value = 7),
      sliderInput("pmf", "Parasite Multiplication factor: ", min = 8, max=10, value=8)
    ),
      
    column(4,
           strong("DHA pharmacodynamics"),
       # sliderInput("initconc", "Initial drug A concentration: ", min = 25, max=150, value = 72, step = 5 ),
       # sliderInput("halflife", "Drug A halflife: ", min = 5, max=240, value=16, step = 1), #value=54, step = 1),
      sliderInput("killrate", "Kill rate: ", min = 0, max = .4, value = .2, step = .01),
      sliderInput("sen", "Sensitivity: ", min = .00, max = 1.00, value = 1, step = .01),
      #sliderInput("ce50", "CE50: ", min = 10, max = 100, value = 15, step = 5),
      sliderInput("ce50", "Log10 of CE50: ", min = .5, max = 3, value = .4, step = .1),
      sliderInput("h", "h: ", min = 1, max = 10, value =4, step =1)
      ),
      column(4, strong("Piperaquine pharmacodynamics"),
      #drug2
      # sliderInput("initconc_2", "Initial drug B concentration: ", min = 25, max=150, value = 50, step = 5 ),
      # sliderInput("halflife_2", "Drug B halflife: ", min = 5, max=400, value=324, step = 1),
      sliderInput("killrate_2", "Kill rate: ", min = 0, max = .4, value = .1, step = .01),
      sliderInput("sen_2", "Sensitivity: ", min = .00, max = 1.00, value = 1, step = .01),
      #sliderInput("ce50_2", "CE50: ", min = 10, max = 100, value = 30, step = 5),
      sliderInput("ce50_2", "Log10 of CE50: ", min = .5, max = 3, value = .6, step = .1),
      sliderInput("h_2", "h: ", min = 1, max = 10, value =4, step =1)
      )
    # ,
    #   column(3,
    #   #drug3
    #   sliderInput("initconc_3", "Initial drug C concentration: ", min = 25, max=150, value = 60, step = 5 ),
    #   sliderInput("halflife_3", "Drug C halflife: ", min = 5, max=400, value=24, step = 1),
    #   sliderInput("killrate_3", "Kill rate, drug C: ", min = 0, max = .4, value = .15, step = .01),
    #   sliderInput("sen_3", "Sensitivity, drug C: ", min = .00, max = 1.00, value = .8, step = .01),
    #   sliderInput("ce50_3", "CE50, drug C: ", min = 10, max = 100, value = 25, step = 5)
    #   )
     ),
  fluidRow(p("Sensitivity to drug, toggle:"),
           actionButton("aSen", "DHA:100%, pip:50%"),
           actionButton("bSen", "DHA:50%, pip:100%")),
  # fluidRow(
  #   # mainPanel(
  #     column(4,plotOutput("paraPlot")),
  #     column(4,plotOutput("combinedPlot")),
  #     column(4,plotOutput("drugeffPlot"))
  #   ),
  fluidRow(
    column(4,plotOutput("paraPlot_DHApip")),
    column(4,plotOutput("DHA_PIP_Plot")),
    column(4,plotOutput("drugeffPlot_DHApip"))
  ),
  fluidRow(
    column(4,plotOutput(("killRate48")))
  )
  )
)
#)


