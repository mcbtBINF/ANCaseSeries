## Using Model of taxa ~ Energy Intake*Patient
## Outputs p-values, corrected p-values and graphs of log normalized abundance
## By: Matthew C. B. Tsilimigras, Revised Jan 2017

rm(list=ls())
setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal")

library("Kendall")
library("vegan")
library("lmtest")
library("pscl")
library("nlme")
library("gtools")

taxaLevels <- c("phylum","class","order","family","genus")

for(t in taxaLevels )
{
    pdf( paste(t, "_FINAL_EnergyIntakePatient_ANOVA_NoLow_plots.pdf", sep = ""))
    inFileName <- paste(t, "LogNormalwithMetadata_Edit.txt", sep="")
    myT <-read.table(inFileName, header=TRUE, sep="\t")
    numCols <- ncol(myT)
    myColClasses <- c(rep("character", 2), "numeric", "character", rep("numeric", numCols-4))
    myT <-read.table(inFileName, header=TRUE, sep="\t", colClasses=myColClasses)
    myT <- myT[ !is.na(myT[2]), ]

    names <-vector()

    EnergyIntakePatientpVal <- list()
    OLDEnergyIntakePatientpVal <- list()

    ## Removes samples of low sequencing depth.
    myT <- myT[-which(myT$Sample.ID %in% list(37, 45, 52, 58, 7, 9),arr.ind=TRUE),]

    ## Patients to colors
    colors <- vector()
    patient <- vector()
    cIndex <- 1

    for ( j in 1: nrow(myT))
    {
        if( substr(myT[j,]$Sample.ID,1,1) == "B")
        {
            colors[cIndex] <- "Blue"
            patient[cIndex] <- "B"
        } else if (substr(myT[j,]$Sample.ID,1,1) == "C" )
        {
            colors[cIndex] <- "Red"
            patient[cIndex] <- "C"
        } else
        {
            colors[cIndex] <- "Black"
            patient[cIndex] <- "A"
        }
        cIndex = cIndex + 1
    }

    index <-1

    myT <- myT[mixedorder(myT[,1]),]

    for( i in 2:(ncol(myT) - 14) )
    {
        ## Thresholding to require 25% of Samples that detects taxa to save statistical power
        if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
        {
            Day<-myT$Day
            ImputedBMI<-myT$Imputed.BMI
            BMI <- myT$BMI
            EnergyIntake<-myT$Energy.Intake..kcal.day.
            taxaType <- as.numeric(myT[,i])

            EnergyIntakePatient<-lm(taxaType ~  EnergyIntake*patient, x = TRUE)

            OLDEnergyIntakePatientpVal[index] <- list(summary(EnergyIntakePatient)$coefficients[,4][-1])
            EnergyIntakePatientpVal[index] <- list(anova(EnergyIntakePatient)$"Pr(>F)"[1:3])

            names[index] = names(myT)[i]
            index = index + 1
        }
    }

    ## Building the data.frames to eventually print out the p-values
    EnergyIntakePatientPV.df <- data.frame(EnergyIntakePatientpVal)
    EnergyIntakePatientPV.df <- t(EnergyIntakePatientPV.df)
    modeldf <- as.data.frame(matrix(unlist(OLDEnergyIntakePatientpVal), nrow=length(OLDEnergyIntakePatientpVal), byrow = TRUE))
    dFrameEnergyIntakePatient <- data.frame(names, EnergyIntakePatientPV.df, modeldf)
    colnames(dFrameEnergyIntakePatient) <- c("names", "ANOVA->Energy Intake", "ANOVA->patient", "ANOVA->Energy Intake:patient", "Day", "patientB", "patientC", "Energy Intake:patientB", "Energy Intake:patientC")

    ## Multiple hypothesis correction
    for (m in 2:dim(dFrameEnergyIntakePatient)[2])
    {
        dFrameEnergyIntakePatient[,dim(dFrameEnergyIntakePatient)[2] + 1] <- p.adjust(dFrameEnergyIntakePatient[,m], method = "BH")
        colnames(dFrameEnergyIntakePatient)[ncol(dFrameEnergyIntakePatient)]<-paste0("adj",colnames(dFrameEnergyIntakePatient)[m])
    }

    ## Repeat modeling so as to use corrected p-values for the graph display
    index <-1

    myT <- myT[mixedorder(myT[,1]),]

    for( i in 2:(ncol(myT) - 14) )
    {
        if( sum( myT[,i] >0 , na.rm=TRUE) > nrow(myT) /4 )
        {
            Day<-myT$Day
            ImputedBMI<-myT$Imputed.BMI
            BMI <- myT$BMI
            EnergyIntake<-myT$Energy.Intake..kcal.day.
            taxaType <- as.numeric(myT[,i])

            EnergyIntakePatient<-lm(taxaType ~ EnergyIntake*patient, x = TRUE)

            OLDEnergyIntakePatientpVal[index] <- list(summary(EnergyIntakePatient)$coefficients[,4][-1])
            EnergyIntakePatientpVal[index] <- list(anova(EnergyIntakePatient)$"Pr(>F)"[1:3])

            names[index] = names(myT)[i]

            ## Graphs for each of the models here...
            ## These are corrected p-values now
            graphMain = paste(names(myT)[i], "\n",
                              "pEnergyIntake=", format(dFrameEnergyIntakePatient[index, 13], digits=3), "\n",
                              "pPatientB=", format(dFrameEnergyIntakePatient[index, 14], digits=3),
                              "pPatientC=", format(dFrameEnergyIntakePatient[index, 15], digits=3), "\n",
                              "pEnergyIntake:PatientB=", format(dFrameEnergyIntakePatient[index,16], digits=3),
                              "pEnergyIntake:PatientC=", format(dFrameEnergyIntakePatient[index, 17], digits=3))
            par(mar = c(5, 4, 6, 2))

            plot(EnergyIntake, taxaType, col=colors, main=graphMain)
            legend("bottomright",
                   c("Patient A", "Patient B", "Patient C"),
                   pch = c(16, 16, 16),
                   col=c("BLACK", "BLUE", "RED"))

            abline(a = EnergyIntakePatient$coef[1], b = EnergyIntakePatient$coef[2])
            abline(a = EnergyIntakePatient$coef[1] + EnergyIntakePatient$coef[3], b = EnergyIntakePatient$coef[5] + EnergyIntakePatient$coef[2], col="BLUE")
            abline(a = EnergyIntakePatient$coef[1] + EnergyIntakePatient$coef[4], b = EnergyIntakePatient$coef[6] + EnergyIntakePatient$coef[2], col="RED")

            index = index + 1
        }
    }

    ## Finally, writing out the p-values and BH adjusted p-values
    write.table(dFrameEnergyIntakePatient, file = paste("FINAL_pValuesLongPatient_EnergyIntake_ANOVA_NoLow_", t, ".txt", sep=""), row.names=FALSE, sep="\t")

    dev.off()
}
