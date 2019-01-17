library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$NJWimPlot <- renderPlot({
    plot(NJWIm(input$initn,48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
               input$killrate,input$ce50,input$h),xlim=c(0,245),ylim=c(0,8),type=show)
  })
  
  output$drugfPlot <- renderPlot({
    plot(drugf(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,245),ylim=c(0,80),
         type=show)
  })
  
  output$drugeffPlot <- renderPlot({
    plot(drugeff(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,125),ylim=c(0,100),
         type=show)
  })
  
})
