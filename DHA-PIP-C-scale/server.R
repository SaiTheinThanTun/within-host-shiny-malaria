library(shiny)
library(ggplot2)

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
    
    parasiteDensity_DHApip <- reactive(NJWIm_DHApip(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,
                                           killrate_R(),ce50_R(),input$h, 
                                           killrate_2_R(),ce50_2_R(), input$h_2))
    #output
    #1. time, 2. log10, 3. normal, 4. gam, 5. infect
    mycol <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue50")
    
    plot(x=parasiteDensity_DHApip()[,1],y=parasiteDensity_DHApip()[,2], xlab="Time (hours)", ylab="Parasite density (log10)",xlim=c(0,1200),ylim=c(0,max(c(8,parasiteDensity_DHApip()[,2]))), type = 'l')
    lines(x=parasiteDensity_DHApip()[,1],y=parasiteDensity_DHApip()[,4], lty="dotdash")
    auxillaryParasiteDensity <- parasiteDensity_DHApip()[,3]
    auxillaryParasiteDensity[auxillaryParasiteDensity==0] <- 1
    text(x = 520, y=7, paste("Log-Sum of observable parasites: ", round(sum(log10(auxillaryParasiteDensity)),3)))
    legend(650, 4.5, legend = c("Total countable parasites", "Gametocytes"), lty = c(1,4))
    par(new=TRUE)
    plot(x=parasiteDensity_DHApip()[,1],y=parasiteDensity_DHApip()[,5], type='l',col='blue', axes=FALSE, xlab="", ylab = "", ylim=c(0,1.2))
    #barplot(parasiteDensity_DHApip()[,5], axes=FALSE, col="red", border="red")
    polygon(x=parasiteDensity_DHApip()[,1],y=c(parasiteDensity_DHApip()[,5][-nrow(parasiteDensity_DHApip())],0), col=mycol)
    axis(4, at=c(0,1), col="blue",col.axis="blue",las=1)
    mtext("Probability of infectiousness",side=4,col="blue")
    
    #text(x = 520, y=7, paste("Log-Sum of observable parasites: ", round(log10(sum(parasiteDensity_DHApip()[,3])),3))) #paste("Log-Sum of observable parasites: ", sum(round(parasiteDensity_DHApip()[,2]))))
    
    #TODO
    # abline(v=TrueMIC(whereIsMIC_D1_R()), col="blue")
    # abline(v=TrueMIC(whereIsMIC_D2_R()), col="red")
    # legend(600, 9.5, legend = c("MIC of drug A", "MIC of drug B", "MIC is not correct yet!"), col = c("blue", "red", "black"), lty = 1)
  })
  
  output$DHA_PIP_Plot <- renderPlot({
    drugConc_A <- DHAconcentration #reactive(drugf(2400,input$initconc,0.693,input$halflife,killrate_R(),ce50_R(),input$h))
    drugConc_A[,2] <- log10(drugConc_A[,2]) 
    drugConc_A[,2][is.infinite(drugConc_A[,2])] <- 0
    
    drugConc_B <- PIPconcentration #reactive(drugf(2400,input$initconc_2,0.693,input$halflife_2,killrate_2_R(),ce50_2_R(),input$h))
    drugConc_B[,2] <- log10(drugConc_B[,2])
    drugConc_B[,2][is.infinite(drugConc_B[,2])] <- 0
    #drugConc_C <- reactive(drugf(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h))
    
    #par(mar=c(5, 4, 4, 6) + 0.1)
    #plot(drugConc_A,xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="blue")
    plot(drugConc_A,xlab="Time (hours)", ylab="Drug concentration (log10)",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2])))), type = 'l', col="blue")
    par(new=TRUE)
    #plot(drugConc_B, axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="red")
    plot(drugConc_B, axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2])))), type = 'l', col="red")
    par(new=TRUE)
    #plot(drugConc_C(), axes=FALSE,xlab="", ylab="",xlim=c(0,400),ylim=c(0,max(c(max(drugConc_A[,2]),max(drugConc_B[,2]),max(drugConc_C()[,2])))), type = 'l', col="black")
    #legend(150, 1, legend = c("DHA", "Piperaquine", "Drug C"), col = c("blue", "red", "black"), lty = 1)
    legend(150, 1, legend = c("DHA", "Piperaquine"), col = c("blue", "red"), lty = 1)
  })
  
  output$drugeffPlot_DHApip <- renderPlot({
    drugEff_A <- reactive(SpecificDrugEffect(2400,"DHA",killrate_R(),ce50_R(),input$h))
    drugEff_B <- reactive(SpecificDrugEffect(2400,"pip",killrate_2_R(),ce50_2_R(),input$h_2))
    #drugEff_C <- reactive(drugeff(2400,input$initconc_3,0.693,input$halflife_3,killrate_3_R(),input$ce50_3,input$h))
    plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l', ylab="Drug effect", col="blue")
    #plot(drugEff_A(), xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l', xlab="Time (hours)", ylab="Drug effect", col="blue")
    par(new=TRUE)
    plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2])))), type = 'l',  col="red")
    #plot(drugEff_B(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="red")
    #par(new=TRUE)
    #plot(drugEff_C(), axes=FALSE, xlab="", ylab="", xlim=c(0,400),ylim=c(0,max(c(max(drugEff_A()[,2]),max(drugEff_B()[,2]),max(drugEff_C()[,2])))), type = 'l',  col="black")
    legend(150, .05, legend = c("DHA", "Piperaquine"), col = c("blue", "red"), lty = 1)
    # legend(150, .05, legend = c("DHA", "Piperaquine", "Drug C"), col = c("blue", "red", "black"), lty = 1)
  })
  
  output$killRate48 <- renderPlot({
    paraDen <- reactive(NJWIm_DHApip(initn_R(),48,input$mu,input$sig,input$pmf,k0,a,tpar,delay,2400,
                                                    killrate_R(),ce50_R(),input$h,
                                                    killrate_2_R(),ce50_2_R(), input$h_2))
    #output
    #1. time, 2. log10, 3. normal, 4. gam, 5. infect
    every48 <- seq(from=1,to=2400, by=48)
    paraDen48 <- reactive(paraDen()[every48,3])
    killRate48 <- NA
    for(i in 1:(length(paraDen48())-1)){
      killRate48[i] <- paraDen48()[i]*8-paraDen48()[i+1]
    }

    killRate48[which(killRate48==0)] <- 1 #to prevent infinity values from log10

    barplot(log10(killRate48)[1:25], main = "Kill rate \n10*(# parasites at a timepoint)-(# parasites after 48 hours)", xlab="Time @48 hours interval", ylab = "Log10(parasite difference between 48 hr)", ylim=c(0,max(log10(killRate48))))
  })
})
