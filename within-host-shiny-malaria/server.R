library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
  output$NJWimPlot <- renderPlot({
    plot(NJWIm(initn,48,mu,sig,pmf,k0,a,tpar,delay,2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,245),ylim=c(0,8),type=show),
    initn=slider(10,1000),
    mu=slider(1,48),
    sig=slider(1,48),
    pmf=slider(1,30),
    # k0=slider(0,1,step = 0.01),
    # a=slider(0.0001,0.01,step = 0.001),
    # tpar=slider(1000,100000,step = 10),
    # delay=slider(0,4,step=1),
    initconc=slider(25,150,step=5),
    halflife=slider(5,240,step=1),
    killrate=slider(0.001,1,step=0.05),
    ce50=slider(10,100,step=5),
    h=slider(1,10,step=1),
    show=picker("line"="l", "points"="p")
  })
  
  output$drugfPlot <- renderPlot({
    plot(drugf(2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,245),ylim=c(0,80),type=show),
    initconc=slider(25,150,step=5),
    halflife=slider(5,240,step=1),
    killrate=slider(0.001,1,step=0.05),
    ce50=slider(10,100,step=5),
    h=slider(1,10,step=1),
    show=picker("line"="l", "points"="p")
  })
  
  output$drugeffPlot <- renderPlot({
    plot(drugeff(2400,initconc,0.693,halflife,killrate,ce50,h),xlim=c(0,125),ylim=c(0,100),type=show),
    initconc=slider(25,150,step=5),
    halflife=slider(5,240,step=1),
    killrate=slider(1,100,step=0.5),
    ce50=slider(10,100,step=5),
    h=slider(1,10,step=1),
    show=picker("line"="l", "points"="p")
  })
  
})
