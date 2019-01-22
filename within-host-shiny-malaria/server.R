library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

shinyServer(function(input, output) {
  
  output$combinedPlot <- renderPlot({
    # parasiteDensity <- reactive(NJWIm(input$initn,48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
    #                                   input$killrate,input$ce50,input$h))
    parasiteDensity_2D <- reactive(NJWIm_2(input$initn,48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                      input$killrate,input$ce50,input$h, input$initconc_2,0.693,input$halflife_2,
                                      input$killrate_2,input$ce50_2))
    
    drugConc_A <- reactive(drugf(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h))
    drugConc_B <- reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,input$killrate_2,input$ce50_2,input$h))
    
    par(mar=c(5, 4, 4, 6) + 0.1)
    # plot(parasiteDensity(), axes=FALSE, xlab="", ylab="",xlim=c(0,245),ylim=c(0,8), type = 'l')
    plot(parasiteDensity_2D(), axes=FALSE, xlab="", ylab="",xlim=c(0,245),ylim=c(0,8), type = 'l')
    axis(2, col="black", col.axis="black", las=1)
    mtext("Parasite density (log10)", side=2, col="black", line=3)
    #box()
    par(new=TRUE)
    plot(drugConc_A(), axes=FALSE,xlab="", ylab="",xlim=c(0,245),ylim=c(0,max(c(80,max(drugConc_A()[,2]),max(drugConc_B()[,2])))), type = 'l', col="blue")
    mtext("Drug concentration (log10)",side=4,col="blue",line=3) 
    axis(4, col="blue", col.axis="blue", las=1)
    
    par(new=TRUE)
    plot(drugConc_B(), axes=FALSE,xlab="", ylab="",xlim=c(0,245),ylim=c(0,max(c(80,max(drugConc_A()[,2]),max(drugConc_B()[,2])))), type = 'l', col="red")
  })
  
  output$drugeffPlot <- renderPlot({
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l')
    plot(drugeff(2400,input$initconc,0.693,input$halflife,1000*input$killrate,input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l', ylab="Drug effect (1000*log10)", col="blue")
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,input$killrate,input$ce50,input$h),xlim=c(0,245),ylim=c(0,10), type = 'l')
    par(new=TRUE)
    plot(drugeff(2400,input$initconc_2,0.693,input$halflife_2,1000*input$killrate_2,input$ce50_2,input$h), axes=FALSE, xlab="", ylab="", xlim=c(0,125),ylim=c(0,100), type = 'l',  col="red")
    #ylab="Drug B effect (1000*log10)",
  })
  

})
