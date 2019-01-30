library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

shinyServer(function(input, output, session) {
  
  initn_R <- reactive(10^input$initn) #initial parasites in 10-power scale
  killrate_R <- reactive(input$killrate*input$sen) #killrate adjusted by sensitivity
  killrate_2_R <- reactive(input$killrate_2*input$sen_2)
  
  #toggle
  observeEvent(input$aSen, {
    updateSliderInput(session, "sen", value = 1)
    updateSliderInput(session, "sen_2", value = .5)
  })
  
  observeEvent(input$bSen, {
    updateSliderInput(session, "sen", value = .5)
    updateSliderInput(session, "sen_2", value = 1)
  })
  
  output$paraPlot <- renderPlot({
    parasiteDensity_2D <- reactive(NJWIm_2(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                           killrate_R(),input$ce50,input$h, input$initconc_2,0.693,input$halflife_2,
                                           killrate_2_R(),input$ce50_2))
    
    plot(parasiteDensity_2D(), xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,245),ylim=c(0,10), type = 'l')
    text(x = 120, y=7, paste("Sum of observable parasites: ", sum(round(parasiteDensity_2D()[,2]))))
    #axis(2, col="black", col.axis="black", las=1)
    #mtext("Parasite density (log10)", side=2, col="black", line=3)
  })
  
  output$combinedPlot <- renderPlot({
    # parasiteDensity <- reactive(NJWIm(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
    #                                   input$killrate,input$ce50,input$h))
    parasiteDensity_2D <- reactive(NJWIm_2(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                      killrate_R(),input$ce50,input$h, input$initconc_2,0.693,input$halflife_2,
                                      killrate_2_R(),input$ce50_2))
    
    drugConc_A <- reactive(drugf(2400,input$initconc,0.693,input$halflife,killrate_R(),input$ce50,input$h))
    drugConc_B <- reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),input$ce50_2,input$h))
    
    par(mar=c(5, 4, 4, 6) + 0.1)
    # plot(parasiteDensity(), axes=FALSE, xlab="", ylab="",xlim=c(0,245),ylim=c(0,8), type = 'l')
    plot(parasiteDensity_2D(), axes=FALSE, xlab="", ylab="",xlim=c(0,245),ylim=c(0,10), type = 'l')
    text(x = 120, y=7, paste("Sum of observable parasites: ", sum(round(parasiteDensity_2D()[,2]))))
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
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,killrate_R(),input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l')
    plot(drugeff(2400,input$initconc,0.693,input$halflife,1000*killrate_R(),input$ce50,input$h),xlim=c(0,125),ylim=c(0,100), type = 'l', ylab="Drug effect (1000*log10)", col="blue")
    #plot(drugeff(2400,input$initconc,0.693,input$halflife,killrate_R(),input$ce50,input$h),xlim=c(0,245),ylim=c(0,10), type = 'l')
    par(new=TRUE)
    plot(drugeff(2400,input$initconc_2,0.693,input$halflife_2,1000*killrate_2_R(),input$ce50_2,input$h), axes=FALSE, xlab="", ylab="", xlim=c(0,125),ylim=c(0,100), type = 'l',  col="red")
    #ylab="Drug B effect (1000*log10)",
  })
  

})
