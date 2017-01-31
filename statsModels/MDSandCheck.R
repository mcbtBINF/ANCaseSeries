## Code for PCoA of AN microbiome samples passing cutoff threshold.
## By: Matthew C. B. Tsilimigras, Revised Jan 2017

rm(list=ls())
library("vegan")
library("calibrate")

setwd("/Users/mbrown67/Documents/Fodor/Datasets/CarrollData/Carroll_Longitudinal/")

taxaLevels <- c( "phylum", "class", "order", "family", "genus" )

countList <- list()
index<-1
for(taxa in taxaLevels )
{
    ## Get taxa counts
    rawFileName <- paste0("pivoted_", taxa, "asColumns.txt")
    myTCheck <- read.table(rawFileName, header=TRUE, sep="\t")

    myTCheck <- myTCheck[1:147,]
    rownames(myTCheck) <- myTCheck$sample
    myTCheck <- myTCheck[,-1]
    rownames(myTCheck) <- unlist(lapply(strsplit(rownames(myTCheck),split="r1_"),"[[",2))

    ## Drop bad samples based on sharp drops in reads relative to flanking samples
    ## r1_7, r1_9, r1_37, r1_45, r1_52, r1_58
    myTCheck<-myTCheck[!(rownames(myTCheck) %in% c("37", "45", "52", "58", "7", "9"
                                                   )),]
    mySums <- rowSums(myTCheck)

    ## Read in log normalized taxa abundance file
    inFileName <- paste(taxa, "LogNormalwithMetadata_Edit.txt", sep="")
    myT <-read.table(inFileName,header=TRUE,sep="\t")
    numCols <- ncol(myT)
    myColClasses <- c("character", rep("numeric", numCols-1))
    myT <-read.table(inFileName,header=TRUE,sep="\t", colClasses=myColClasses)
    myT <- myT[ !is.na(myT[2]), ]
    rownames(myT)<-myT[,1]
    myT<-myT[,-1]
    myT<-myT[!(rownames(myT) %in% c("37", "45", "52", "58", "7", "9"
                                    )),]

    colors <- rep("Black", dim(myT)[1])
    colors[grep("B", rownames(myT))] <- c("Blue")
    colors[grep("C", rownames(myT))] <- c("Red")

    ## capscale includes ONLY relative taxa abundances
    myMDS <- capscale(myT[,1:(ncol(myT)-14)]~1,distance="bray")

    pdf( paste(taxa, "CURRDROP_topMDS_R1_test.pdf",sep=""))
    for (xrun in 1:4) {
        for (yrun in 2:4) {
            if(xrun == yrun){
                break
            }
            plot(myMDS$CA$u[,xrun], myMDS$CA$u[,yrun],xlab=paste("MDS",xrun,sep=""), ylab=paste("MDS",yrun,sep=""), main=paste("PCoA at level:", taxa,sep=""),cex=1.0,
                 pch=NA,
                 col=colors
                 )
            text(myMDS$CA$u[,xrun],myMDS$CA$u[,yrun],labels=paste0(format(Shannon, digits=2),"_",mySums,"_",rownames(myTCheck)), cex=0.7, offset=0, col=colors)
            legend("topright",
                   c("Patient A", "Patient B", "Patient C"),
                   pch=c(16, 16, 16),
                   col=c("black", "blue", "red"))
        }
    }

    dev.off()
    write.table(myMDS$CA$u, sep="\t", file=paste("R1_test_pcoa_CURRDROP_", taxa, ".txt",sep=""))
    write.table(myMDS$CA$eig,file=paste("R1_test_eigenValues_CURRDROP_", taxa, ".txt", sep=""), sep="\t")

}
