#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019-2021 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2021.05.10.2000"

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
start_year = 2013

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

## Load IR Stats Table (.csv-format)
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(IRTableFn), "_", sep='')
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")

# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")
# Alternative code via sqldf:
#library(sqldf)
#combinedDF = sqldf("SELECT AvailTableData.* FROM AvailTableData, IRTableData WHERE AvailTableData.ACCESSION == IRTableData.ACCESSION")

########################################################################

## Ensure that numeric values are recognized as numeric (i.e., transform from factor to numeric via character)
## See: https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
combinedDF = transform(combinedDF,
                       IRa_REPORTED_START = as.integer(as.character(IRa_REPORTED_START)),
                       IRb_REPORTED_START = as.integer(as.character(IRb_REPORTED_START)),
                       IRa_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),
                       IRb_REPORTED_END = as.integer(as.character(IRb_REPORTED_END)),
                       IRa_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH)),
                       IRb_REPORTED_LENGTH = as.integer(as.character(IRb_REPORTED_LENGTH))
)

# Convert column CREATE_DATE to dates
combinedDF$CREATE_DATE <- as.Date(combinedDF$CREATE_DATE)

########################################################################


## Extracting information on RELEASE YEAR
DataOnReleaseYear <- data.frame(combinedDF$ACCESSION)
DataOnReleaseYear$IS_CONGRUENT <- combinedDF$IRa_REPORTED == "yes" &
  combinedDF$IRb_REPORTED == "yes" &
  equalityCheck(combinedDF$IRa_REPORTED_LENGTH, combinedDF$IRb_REPORTED_LENGTH)
DataOnReleaseYear$RELEASE_YEAR <- format(as.Date(combinedDF$CREATE_DATE, format="%d/%m/%Y"), "%Y")
colnames(DataOnReleaseYear) <- c("ACCESSION", "IS_CONGRUENT", "RELEASE_YEAR")
DataOnReleaseYear$count <- 1

## Extracting information on PUBLICATION STATUS
DataOnPublStatus <- data.frame(combinedDF$ACCESSION)
DataOnPublStatus$IS_PUBLISHED <- !(combinedDF$REFERENCE == "Unpublished")
DataOnPublStatus$IS_CONGRUENT <- combinedDF$IRa_REPORTED == "yes" &
  combinedDF$IRb_REPORTED == "yes" &
  equalityCheck(combinedDF$IRa_REPORTED_LENGTH, combinedDF$IRb_REPORTED_LENGTH)
DataOnPublStatus$count <- 1

## Extracting information on SEQUENCE VERSION
DataOnSeqVersion = data.frame(combinedDF$ACCESSION)
DataOnSeqVersion$IS_CONGRUENT <- combinedDF$IRa_REPORTED == "yes" &
  combinedDF$IRb_REPORTED == "yes" &
  equalityCheck(combinedDF$IRa_REPORTED_LENGTH, combinedDF$IRb_REPORTED_LENGTH)
DataOnSeqVersion$SEQVERSION <- as.numeric(combinedDF$VERSION)
colnames(DataOnSeqVersion) <- c("ACCESSION","IS_CONGRUENT", "SEQVERSION")
DataOnSeqVersion$count <- 1

# Obtain total numbers for year
plotData_year <- aggregate(DataOnReleaseYear$count,
                      by=list(DataOnReleaseYear$RELEASE_YEAR,
                              DataOnReleaseYear$IS_CONGRUENT),
                      FUN=sum)
colnames(plotData_year) <- c("RELEASE_YEAR", "IS_CONGRUENT", "TOTAL")

# Obtain total numbers for publication
plotData_publ <- aggregate(DataOnPublStatus$count,
                           by=list(DataOnPublStatus$IS_PUBLISHED,
                                   DataOnPublStatus$IS_CONGRUENT),
                           FUN=sum)
colnames(plotData_publ) <- c("IS_PUBLISHED", "IS_CONGRUENT", "TOTAL")

# Obtain total numbers for seq.number
plotData_seqn <- aggregate(DataOnSeqVersion$count,
                      by=list(DataOnSeqVersion$SEQVERSION,
                              DataOnSeqVersion$IS_CONGRUENT),
                      FUN=sum)
colnames(plotData_seqn) <- c("SEQVERSION", "IS_CONGRUENT", "TOTAL")

# Obtain percentages for year
acc_per_year <- aggregate(DataOnReleaseYear$count,
                          by=list(as.numeric(DataOnReleaseYear$RELEASE_YEAR)),
                          FUN=sum)
colnames(acc_per_year) <- c("RELEASE_YEAR", "TOTAL")
plotData_year$PERCENTAGE <- as.numeric(sqldf("SELECT (plotData_year.TOTAL/acc_per_year.TOTAL) FROM plotData_year,
                                        acc_per_year WHERE plotData_year.RELEASE_YEAR = acc_per_year.RELEASE_YEAR")[,1])
# Round the precentages to three comma positions
plotData_year = plotData_year %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

# Obtain percentages for publication
plotData_publ <- plotData_publ[order(plotData_publ$IS_PUBLISHED),]
plotData_publ$PERCENTAGE <- c(filter(plotData_publ, IS_PUBLISHED == FALSE)$TOTAL /
                           sum(filter(plotData_publ, IS_PUBLISHED == FALSE)$TOTAL),
                         filter(plotData_publ, IS_PUBLISHED == TRUE)$TOTAL /
                           sum(filter(plotData_publ, IS_PUBLISHED == TRUE)$TOTAL))
# Round the precentages to three comma positions
plotData_publ = plotData_publ %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

# Obtain percentages for seq.number
acc_per_version <- aggregate(DataOnSeqVersion$count,
                             by=list(DataOnSeqVersion$SEQVERSION),
                             FUN=sum)
colnames(acc_per_version) <- c("SEQVERSION", "TOTAL")
plotData_seqn$PERCENTAGE <- as.numeric(sqldf("SELECT (plotData_seqn.TOTAL/acc_per_version.TOTAL) FROM plotData_seqn,
                                        acc_per_version WHERE plotData_seqn.SEQVERSION = acc_per_version.SEQVERSION")[,1])
# Round the precentages to three comma positions
plotData_seqn = plotData_seqn %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

# Only data up to and including sequence version 3 is given
plotData_seqn = plotData_seqn[which(plotData_seqn$SEQVERSION<=3),]
plotData_seqn$SEQVERSION <- as.factor(plotData_seqn$SEQVERSION)

########################################################################

## PLOTTING OPERATIONS ##
year_plot = ggplot(data=plotData_year, aes(x=RELEASE_YEAR, y=PERCENTAGE), width=1) +
  geom_col(aes(fill=IS_CONGRUENT), alpha=0.5) +
  geom_text(data=plotData_year[which(plotData_year$IS_CONGRUENT==FALSE),], aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
  geom_text(data=plotData_year[which(plotData_year$IS_CONGRUENT==TRUE),], aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3) + 
  ggtitle("(a)") + 
  ylab("Percentage of Records\n") +
  scale_x_discrete(limits=factor(seq(start_year, 2020, 1)), labels=seq(start_year, 2020, 1)) +
  #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
  scale_fill_manual(values=c("grey0", "grey50"),
                    name="Complete and correct IR annotations present",
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

publication_plot = ggplot(data=plotData_publ, aes(x=IS_PUBLISHED, y=PERCENTAGE), width=1) +
  geom_col(aes(fill=IS_CONGRUENT), alpha=0.5) +
  geom_text(data=plotData_publ[which(plotData_publ$IS_CONGRUENT==FALSE),], aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
  geom_text(data=plotData_publ[which(plotData_publ$IS_CONGRUENT==TRUE),], aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3) + 
  ggtitle("(b)") +
  scale_x_discrete(labels=c("Unpubl.", "Publ.")) +
  #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
  scale_fill_manual(values=c("grey0", "grey50"),
                    name="Complete and correct IR annotations present",
                    labels=c("No", "Yes")) +
  theme_minimal() + 
  theme(plot.title = element_text(size=14, hjust=0.5, face="bold"),
        axis.text=element_text(size=13),
        axis.title=element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=14, face="italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
       legend.key.width=unit(1,"cm"),
       legend.position = "bottom")


seqn_plot = ggplot(data=plotData_seqn, aes(x=SEQVERSION, y=PERCENTAGE), width=1) +
  geom_col(aes(fill=IS_CONGRUENT), alpha=0.5) +
  geom_text(data=plotData_seqn[which(plotData_seqn$IS_CONGRUENT==FALSE),], aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
  geom_text(data=plotData_seqn[which(plotData_seqn$IS_CONGRUENT==TRUE),], aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3) + 
  ggtitle("(c)") +
  scale_x_discrete(labels=c("V1", "V2", "V3")) +
  scale_fill_manual(values=c("grey0", "grey50"),
                    name="Complete and correct IR annotations present",
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

comb_plot <- ggarrange(year_plot, publication_plot, seqn_plot, 
                       nrow=1, ncol=3, common.legend=TRUE, legend="bottom", 
                       widths = c(8,2,3))
ggsave(file="./CompositeIRequalityComparison.pdf", plot=comb_plot, 
       width=11, height=5)

########################################################################

assign(script_name, comb_plot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
