#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019-2021 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2021.05.12.1530"

########################################################################

library(ggplot2)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'
library(sqldf)
library(grid)
library(gtable)
library(ggpubr)

########################################################################

# GETTING SCRIPT NAME
args = commandArgs(TRUE)
this_script = sub(".*=", "", commandArgs()[4])
script_name = file_path_sans_ext(basename(this_script))

########################################################################

#GLOBAL VARIABLES
#start_year = 2000
start_year = 2014

########################################################################

#SMALL HELPER FUNCTION
equalityCheck <- function(lenA, lenB, tolerance){
  if(missing(tolerance)){
    ifelse((lenA == lenB), TRUE, FALSE)
  }else{
    if ((tolerance <= 0) |  (tolerance > 1)){
      stop("Argument tolerance expects value between 0 and 1")
    }else{
      ifelse((lenA == lenB) | ((lenA < lenB) & ((lenA + lenA * tolerance) >= lenB)) | ((lenA > lenB) & ((lenB + lenB * tolerance) >= lenA)), TRUE, FALSE)
    }
  }
}

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")

## Load EXTENDED IR Stats Table (.csv-format)
ExtendedIRTableFnFn = tk_choose.files(caption = "Select the EXTENDED IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(ExtendedIRTableFnFn), "_", sep='')
out_fn = dirname(ExtendedIRTableFnFn)
ExtendedIRTableFnData = read.csv(ExtendedIRTableFnFn, sep = "\t")

# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, ExtendedIRTableFnData, by="ACCESSION")
# Alternative code via sqldf:
#library(sqldf)
#combinedDF = sqldf("SELECT AvailTableData.* FROM AvailTableData, ExtendedIRTableFnData WHERE AvailTableData.ACCESSION == ExtendedIRTableFnData.ACCESSION")

########################################################################

## Ensure that numeric values are recognized as numeric (i.e., transform from factor to numeric via character)
## See: https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
combinedDF = transform(combinedDF,
                       IRa_REPORTED_START = as.integer(as.character(IRa_REPORTED_START)),
                       IRb_REPORTED_START = as.integer(as.character(IRb_REPORTED_START)),
                       IRa_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),

                       IRb_REPORTED_END = as.integer(as.character(IRb_REPORTED_END)),
                       IRa_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH)),
                       IRb_REPORTED_LENGTH = as.integer(as.character(IRb_REPORTED_LENGTH)),

                       IRa_BLASTINFERRED_START = as.integer(as.character(IRa_BLASTINFERRED_START)),
                       IRa_BLASTINFERRED_END = as.integer(as.character(IRa_BLASTINFERRED_END)),
                       IRa_BLASTINFERRED_LENGTH = as.integer(as.character(IRa_BLASTINFERRED_LENGTH)),

                       IRb_BLASTINFERRED_START = as.integer(as.character(IRb_BLASTINFERRED_START)),
                       IRb_BLASTINFERRED_END = as.integer(as.character(IRb_BLASTINFERRED_END)),
                       IRb_BLASTINFERRED_LENGTH = as.integer(as.character(IRb_BLASTINFERRED_LENGTH))
)

# Convert column CREATE_DATE to dates
combinedDF$CREATE_DATE <- as.Date(combinedDF$CREATE_DATE)

########################################################################

## Extracting information for REPORTED IRs

## Selecting only genome records with validly reported IRs
## NOTE: To be validly reported, both IRs of a genome record must be specified and of equal length; this does NOT mean that the IRs are correctly identified!
validlyReportedIRs = combinedDF$IRa_REPORTED == "yes" &
                     combinedDF$IRb_REPORTED == "yes" &
                     equalityCheck(  combinedDF$IRa_REPORTED_LENGTH,
                                     combinedDF$IRb_REPORTED_LENGTH)

IRreportedDF = combinedDF[validlyReportedIRs,]
IRreported <- data.frame(IRreportedDF$ACCESSION)

## NOTES: In genome records with validly reported IRs that are correctly identified, self-BLASTING confirmed the presence and the same length of the IRs as specified in the annotations.
IRreported$CORR_IDENTI <- IRreportedDF$IRa_BLASTINFERRED == "yes" &              # Condition for CorrectlyReported
                          IRreportedDF$IRb_BLASTINFERRED == "yes" &              # Condition for CorrectlyReported
                          equalityCheck(  IRreportedDF$IRa_REPORTED_LENGTH,      # Condition for CorrectlyReported
                                          IRreportedDF$IRa_BLASTINFERRED_LENGTH) & # Condition for CorrectlyReported
                          equalityCheck(  IRreportedDF$IRb_REPORTED_LENGTH,      # Condition for CorrectlyReported
                                          IRreportedDF$IRb_BLASTINFERRED_LENGTH) # Condition for CorrectlyReported

IRreported$RELEASE_YEAR <- format(as.Date(IRreportedDF$CREATE_DATE, format="%d/%m/%Y"), "%Y")
colnames(IRreported) <- c("ACCESSION", "CORR_IDENTI", "RELEASE_YEAR")
IRreported$count <- 1


## Extracting information for UNreported IRs

## Selecting only genome records without reported IRs
unreportedIRs = combinedDF$IRa_REPORTED == "no" &
                combinedDF$IRb_REPORTED == "no"

IRunreportDF = combinedDF[unreportedIRs,]
IRunreport <- data.frame(IRunreportDF$ACCESSION)


## NOTES: In genome records with unreported IRs that do have IRs, self-BLASTING confirmed the presence and the same length of the inferred IRs.
IRunreport$IR_EXISTS <- IRunreportDF$IRa_BLASTINFERRED == "yes" &              # Condition for IR_EXISTS
                        IRunreportDF$IRb_BLASTINFERRED == "yes" &              # Condition for IR_EXISTS
                        equalityCheck(  IRunreportDF$IRa_BLASTINFERRED_LENGTH, # Condition for IR_EXISTS
                                        IRunreportDF$IRb_BLASTINFERRED_LENGTH) # Condition for IR_EXISTS

IRunreport$RELEASE_YEAR <- format(as.Date(IRunreportDF$CREATE_DATE, format="%d/%m/%Y"), "%Y")
colnames(IRunreport) <- c("ACCESSION", "IR_EXISTS", "RELEASE_YEAR")
IRunreport$count <- 1

########################################################################

# Obtain total numbers for REPORTED IRs
plotData_IRreported <- aggregate(IRreported$count,
                                 by=list(IRreported$RELEASE_YEAR,
                                         IRreported$CORR_IDENTI),
                                 FUN=sum)
colnames(plotData_IRreported) <- c("RELEASE_YEAR", "CORR_IDENTI", "TOTAL")

# Obtain total numbers for UNreported IRs
plotData_IRunreport <- aggregate(IRunreport$count,
                                 by=list(IRunreport$RELEASE_YEAR,
                                         IRunreport$IR_EXISTS),
                                 FUN=sum)
colnames(plotData_IRunreport) <- c("RELEASE_YEAR", "IR_EXISTS", "TOTAL")

########################################################################

# Obtain percentages for REPORTED IRs
reported_per_year <- aggregate(IRreported$count,
                               by=list(as.numeric(IRreported$RELEASE_YEAR)),
                               FUN=sum)
colnames(reported_per_year) <- c("RELEASE_YEAR", "TOTAL")
plotData_IRreported$PERCENTAGE <- as.numeric(
    sqldf("SELECT (plotData_IRreported.TOTAL/reported_per_year.TOTAL) FROM plotData_IRreported, reported_per_year 
          WHERE plotData_IRreported.RELEASE_YEAR = reported_per_year.RELEASE_YEAR")[,1]
)
# Round the precentages to three comma positions
plotData_IRreported = plotData_IRreported %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

# Obtain percentages for UNreported IRs
UNreported_per_year <- aggregate(IRunreport$count,
                               by=list(as.numeric(IRunreport$RELEASE_YEAR)),
                               FUN=sum)
colnames(UNreported_per_year) <- c("RELEASE_YEAR", "TOTAL")
plotData_IRunreport$PERCENTAGE <- as.numeric(
    sqldf("SELECT (plotData_IRunreport.TOTAL/UNreported_per_year.TOTAL) FROM plotData_IRunreport, UNreported_per_year 
          WHERE plotData_IRunreport.RELEASE_YEAR = UNreported_per_year.RELEASE_YEAR")[,1]
)
# Round the precentages to three comma positions
plotData_IRunreport = plotData_IRunreport %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

########################################################################

## PLOTTING FOR REPORTED IRs
plot_IRreported = ggplot(data=plotData_IRreported, aes(x=RELEASE_YEAR, y=PERCENTAGE), width=1) +
  geom_col(aes(fill=CORR_IDENTI), alpha=0.5) +
  geom_text(data=plotData_IRreported[which(plotData_IRreported$CORR_IDENTI==FALSE),], 
            aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
  geom_text(data=plotData_IRreported[which(plotData_IRreported$CORR_IDENTI==TRUE),], 
            aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3) + 
  ggtitle("(a) IRs reported") + 
  ylab("Percentage of Records\n") +
  scale_x_discrete(limits=factor(seq(start_year, 2020, 1)), 
                   labels=seq(start_year, 2020, 1)) +
  #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
  scale_fill_manual(values=c("grey0", "grey50"),
                    name="Reported IRs correct",
                    labels=c("No", "Yes")) +
  theme_minimal() + 
  theme(plot.title = element_text(size=14, hjust=0.5, face="bold"),
        plot.subtitle = element_text(size=14, face="italic"),
        axis.title.x = element_blank(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=14, face="bold"),
        #plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
        legend.key.width=unit(1,"cm"),
        legend.position = "bottom")

## PLOTTING FOR UNreported IRs
plot_IRunreport = ggplot(data=plotData_IRunreport, aes(x=RELEASE_YEAR, y=PERCENTAGE), width=1) +
  geom_col(aes(fill=IR_EXISTS), alpha=0.5) +
  geom_text(data=plotData_IRunreport[which(plotData_IRunreport$IR_EXISTS==FALSE),], 
            aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
  geom_text(data=plotData_IRunreport[which(plotData_IRunreport$IR_EXISTS==TRUE),], 
            aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3) + 
  ggtitle("(b) IRs not reported") +
  #ylab("Percentage of Records\n") +
  scale_x_discrete(limits=factor(seq(start_year, 2020, 1)), 
                   labels=seq(start_year, 2020, 1)) +
  #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
  scale_fill_manual(values=c("grey0", "grey50"),
                    name="IRs detected via self-BLASTing",
                    labels=c("No", "Yes")) +
  theme_minimal() +
  theme(plot.title = element_text(size=14, hjust=0.5, face="bold"),
        plot.subtitle = element_text(size=14, face="italic"),
        axis.text=element_text(size=13),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
        legend.key.width=unit(1,"cm"),
        legend.position = "bottom")


comb_plot <- ggarrange(plot_IRreported, plot_IRunreport, 
                       nrow=1, ncol=2, 
                       common.legend=FALSE, #legend="bottom", 
                       widths = c(5,5))
ggsave(file="./FalsePositivesFalseNegatives.pdf", plot=comb_plot, 
       width=11.5, height=5)

########################################################################

assign(script_name, comb_plot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
