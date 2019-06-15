library(shiny)

#utility functions go here
source("UtilityFunctions.R", local = TRUE)

shinyServer(function(input, output, session) {
  
  initn_R <- reactive(10^input$initn) #initial parasites in 10-power scale
  ce50_R <- reactive(10^input$ce50)
  ce50_2_R <- reactive(10^input$ce50_2)
  killrate_R <- reactive(input$killrate*input$sen) #killrate adjusted by sensitivity
  killrate_2_R <- reactive(input$killrate_2*input$sen_2)
  killrate_3_R <- reactive(input$killrate_3*input$sen_3)
  
  #TODO
  #MIC for 3rd drug hasn't be done as changes will still come to MIC #TODO
  #logical vector of where MIC is
  whereIsMIC_D1_R <- reactive(NJWIm(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                  killrate_R(),ce50_R(),input$h)$MIC)
  whereIsMIC_D2_R <- reactive(NJWIm(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc_2,0.693,input$halflife_2,
                                  killrate_2_R(),ce50_2_R(), input$h)$MIC)
  
  #toggle
  observeEvent(input$aSen, {
    updateSliderInput(session, "sen", value = 1)
    updateSliderInput(session, "sen_2", value = .5)
  })
  
  observeEvent(input$bSen, {
    updateSliderInput(session, "sen", value = .5)
    updateSliderInput(session, "sen_2", value = 1)
  })
  #end of #TODO
  output$paraPlot <- renderPlot({
    
    parasiteDensity_3D <- reactive(NJWIm_3(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                           killrate_R(),ce50_R(),input$h, input$initconc_2,0.693,input$halflife_2, 
                                           killrate_2_R(),ce50_2_R(), input$h_2, input$initconc_3,0.693,input$halflife_3,
                                           killrate_3_R(),input$ce50_3, input$h_3))
    
    plot(x=parasiteDensity_3D()[,1],y=parasiteDensity_3D()[,2], xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,1200),ylim=c(0,10), type = 'l')
    text(x = 520, y=7, paste("Log-Sum of observable parasites: ", round(log10(sum(parasiteDensity_3D()[,3]))))) #paste("Log-Sum of observable parasites: ", sum(round(parasiteDensity_3D()[,2]))))
    #abline(v=which(whereIsMIC_D1_R()==TRUE)[1], col="blue")
    #abline(v=which(whereIsMIC_D2_R()==TRUE)[1], col="red")
    
    #TODO
    # abline(v=TrueMIC(whereIsMIC_D1_R()), col="blue")
    # abline(v=TrueMIC(whereIsMIC_D2_R()), col="red")
    # legend(600, 9.5, legend = c("MIC of drug A", "MIC of drug B", "MIC is not correct yet!"), col = c("blue", "red", "black"), lty = 1)
  })
  
  output$paraPlot_obsolete <- renderPlot({
    
    parasiteDensity_2D <- reactive(NJWIm_2(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,input$initconc,0.693,input$halflife,
                                           killrate_R(),ce50_R(),input$h, input$initconc_2,0.693,input$halflife_2,
                                           killrate_2_R(),ce50_2_R(),input$h_2))
    
    plot(x=parasiteDensity_2D()[,1],y=parasiteDensity_2D()[,2], xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,1200),ylim=c(0,10), type = 'l')
    text(x = 520, y=7, paste("Log-Sum of observable parasites: ", round(log10(sum(parasiteDensity_2D()[,3]))))) #paste("Log-Sum of observable parasites: ", sum(round(parasiteDensity_2D()[,2]))))
    #abline(v=which(whereIsMIC_D1_R()==TRUE)[1], col="blue")
    #abline(v=which(whereIsMIC_D2_R()==TRUE)[1], col="red")
    abline(v=TrueMIC(whereIsMIC_D1_R()), col="blue")
    abline(v=TrueMIC(whereIsMIC_D2_R()), col="red")
    legend(600, 9.5, legend = c("MIC of drug A", "MIC of drug B"), col = c("blue", "red"), lty = 1)
  })
  
  output$combinedPlot <- renderPlot({
    drugConc_A <- reactive(drugf(2400,input$initconc,0.693,input$halflife,killrate_R(),ce50_R(),input$h))
    drugConc_B <- reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),ce50_2_R(),input$h_2))
    drugConc_C <- reactive(drugf(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h_3))
    
    #par(mar=c(5, 4, 4, 6) + 0.1)
    plot(drugConc_A(),xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A()[,2]),max(drugConc_B()[,2]),max(drugConc_C()[,2])))), type = 'l', col="blue")
    par(new=TRUE)
    plot(drugConc_B(), axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A()[,2]),max(drugConc_B()[,2]),max(drugConc_C()[,2])))), type = 'l', col="red")
    par(new=TRUE)
    plot(drugConc_C(), axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A()[,2]),max(drugConc_B()[,2]),max(drugConc_C()[,2])))), type = 'l', col="black")
    legend(150, 60, legend = c("A", "B", "C"), col = c("blue", "red", "black"), lty = 1)
  })
  
  output$drugeffPlot <- renderPlot({
    drugEff_A <- reactive(drugeff(2400,input$initconc,0.693,input$halflife,killrate_R(),ce50_R(),input$h))
    drugEff_B <- reactive(drugeff(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),ce50_2_R(),input$h_2))
    drugEff_C <- reactive(drugeff(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h_3))
    #plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(100,max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l', ylab="Drug effect", col="blue")
    plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l', xlab="Time (hours)", ylab="Drug effect", col="blue")
    par(new=TRUE)
    #plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(100,max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l',  col="red")
    plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="red")
    par(new=TRUE)
    plot(drugEff_C(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="black")
  })
  
  #new plots for DHApip calibrated model
  output$paraPlot_DHApip <- renderPlot({
    
    parasiteDensity_DHApip <- reactive(NJWIm_DHApip_C(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,
                                           killrate_R(),ce50_R(),input$h, 
                                           killrate_2_R(),ce50_2_R(), input$h_2,
                                           input$initconc_3,0.693,input$halflife_3,
                                           killrate_3_R(),input$ce50_3,input$h_3))
    
    plot(x=parasiteDensity_DHApip()[,1],y=parasiteDensity_DHApip()[,2], xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,1200),ylim=c(0,10), type = 'l')
    text(x = 520, y=7, paste("Log-Sum of observable parasites: ", round(log10(sum(parasiteDensity_DHApip()[,3]))))) #paste("Log-Sum of observable parasites: ", sum(round(parasiteDensity_DHApip()[,2]))))
    
    #TODO
    # abline(v=TrueMIC(whereIsMIC_D1_R()), col="blue")
    # abline(v=TrueMIC(whereIsMIC_D2_R()), col="red")
    # legend(600, 9.5, legend = c("MIC of drug A", "MIC of drug B", "MIC is not correct yet!"), col = c("blue", "red", "black"), lty = 1)
  })
  
  output$DHA_PIP_Plot <- renderPlot({
    drugConc_A <- DHAconcentration #reactive(drugf(2400,input$initconc,0.693,input$halflife,killrate_R(),ce50_R(),input$h))
    drugConc_A[,2] <- log10(drugConc_A[,2]) 
    drugConc_A[,2][is.infinite(drugConc_A[,2])] <- 1
    
    drugConc_B <- PIPconcentration #reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),ce50_2_R(),input$h))
    drugConc_B[,2] <- log10(drugConc_B[,2])
    drugConc_B[,2][is.infinite(drugConc_B[,2])] <- 1
    drugConc_C <- reactive(drugf(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h_3))
    
    #par(mar=c(5, 4, 4, 6) + 0.1)
    plot(drugConc_A,xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="blue")
    #plot(drugConc_A,xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2])))), type = 'l', col="blue")
    par(new=TRUE)
    plot(drugConc_B, axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="red")
    #plot(drugConc_B, axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2])))), type = 'l', col="red")
    par(new=TRUE)
    plot(drugConc_C(), axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="black")
    legend(150, 1, legend = c("DHA", "Piperaquine", "Drug C"), col = c("blue", "red", "black"), lty = 1)
    #legend(150, 1, legend = c("DHA", "Piperaquine"), col = c("blue", "red"), lty = 1)
  })
  
  output$drugeffPlot_DHApip <- renderPlot({
    drugEff_A <- reactive(SpecificDrugEffect(2400,"DHA",killrate_R(),ce50_R(),input$h))
    drugEff_B <- reactive(SpecificDrugEffect(2400,"pip",killrate_2_R(),ce50_2_R(),input$h_2))
    drugEff_C <- reactive(drugeff(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h_3))
    #plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l', ylab="Drug effect", col="blue")
    plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l', xlab="Time (hours)", ylab="Drug effect", col="blue")
    par(new=TRUE)
    #plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l',  col="red")
    plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="red")
    par(new=TRUE)
    plot(drugEff_C(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="black")
    #legend(150, .05, legend = c("DHA", "Piperaquine"), col = c("blue", "red"), lty = 1)
     legend(150, .05, legend = c("DHA", "Piperaquine", "Drug C"), col = c("blue", "red", "black"), lty = 1)
  })
  

})
