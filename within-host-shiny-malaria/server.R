library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

shinyServer(function(input, output) {
  
  output$combinedPlot <- renderPlot({
    par(mar=c(5, 4, 4, 6) + 0.1)
    plot(NJWIm(input$initn,48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
               input$killrate,input$ce50,input$h), axes=FALSE, xlab="", ylab="",xlim=c(0,245),ylim=c(0,8), type = 'l')
    axis(2, col="black", col.axis="black", las=1)
    mtext("Parasite density (log10)", side=2, col="black", line=3)
    #box()
    par(new=TRUE)
    plot(drugf(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h), axes=FALSE,xlab="", ylab="",xlim=c(0,245),ylim=c(0,80), type = 'l', col="blue")
    mtext("Drug concentration (log10)",side=4,col="blue",line=3) 
    axis(4, col="blue", col.axis="blue", las=1)
  })
  
  output$drugeffPlot <- renderPlot({
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l')
    plot(drugeff(2400,input$initconc,0.693,input$halflife,1000*input$killrate,input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l', ylab="Drug effect (1000*log10)")
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,245),ylim=c(0,10), type = 'l')
  })
  
})
