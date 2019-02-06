library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

shinyServer(function(input, output, session) {
  
  initn_R <- reactive(10^input$initn) #initial parasites in 10-power scale
  killrate_R <- reactive(input$killrate*input$sen) #killrate adjusted by sensitivity
  killrate_2_R <- reactive(input$killrate_2*input$sen_2)
  #logical vector of where MIC is
  whereIsMIC_D1_R <- reactive(NJWIm(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                  killrate_R(),input$ce50,input$h)$MIC)
  whereIsMIC_D2_R <- reactive(NJWIm(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc_2,0.693,input$halflife_2,
                                  killrate_2_R(),input$ce50_2, input$h)$MIC)
  
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
    
    plot(x=parasiteDensity_2D()[,1],y=parasiteDensity_2D()[,2], xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,245),ylim=c(0,10), type = 'l')
    text(x = 120, y=7, paste("Log-Sum of observable parasites: ", sum(round(parasiteDensity_2D()[,2]))))
    abline(v=which(whereIsMIC_D1_R()==TRUE)[1], col="blue")
    abline(v=which(whereIsMIC_D2_R()==TRUE)[1], col="red")
    legend(150, 9.5, legend = c("MIC of drug A", "MIC of drug B"), col = c("blue", "red"), lty = 1)
  })
  
  output$combinedPlot <- renderPlot({
    drugConc_A <- reactive(drugf(2400,input$initconc,0.693,input$halflife,killrate_R(),input$ce50,input$h))
    drugConc_B <- reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),input$ce50_2,input$h))
    
    #par(mar=c(5, 4, 4, 6) + 0.1)
    plot(drugConc_A(),xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,245),ylim=c(0,max(c(80,max(drugConc_A()[,2]),max(drugConc_B()[,2])))), type = 'l', col="blue")
    par(new=TRUE)
    plot(drugConc_B(), axes=FALSE,xlab="", ylab="",xlim=c(0,245),ylim=c(0,max(c(80,max(drugConc_A()[,2]),max(drugConc_B()[,2])))), type = 'l', col="red")
    legend(150, 60, legend = c("A", "B"), col = c("blue", "red"), lty = 1)
  })
  
  output$drugeffPlot <- renderPlot({
    drugEff_A <- reactive(drugeff(2400,input$initconc,0.693,input$halflife,1000*killrate_R(),input$ce50,input$h))
    drugEff_B <- reactive(drugeff(2400,input$initconc_2,0.693,input$halflife_2,1000*killrate_2_R(),input$ce50_2,input$h))
    plot(drugEff_A(), xlim=c(0,245),ylim=c(0,max(c(100,max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l', ylab="Drug effect (1000*log10)", col="blue")
    par(new=TRUE)
    plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,245),ylim=c(0,max(c(100,max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l',  col="red")
  })
  

})
