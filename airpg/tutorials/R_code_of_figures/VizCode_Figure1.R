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
library(forcats)

########################################################################

# GETTING SCRIPT NAME
args = commandArgs(TRUE)
this_script = sub(".*=", "", commandArgs()[4])
script_name = file_path_sans_ext(basename(this_script))

########################################################################

#GLOBAL VARIABLES
#start_year = 2000
start_year = 2010

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")

## Load IR Stats Table (.csv-format)
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(IRTableFn), "_", sep='')
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")

## Combine the data tables such that only accessions which exist in BOTH tables are maintained
# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")

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

########################################################################

## Seperate rows by criterion
posMatch = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="yes"),]
negMatch = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="no"),]

## Convert the column "CREATE_DATE" into a frequency table
# Convert to date
posDatesData = as.Date(posMatch$CREATE_DATE, format="%Y-%m-%d")
negDatesData = as.Date(negMatch$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
posTab = table(cut(posDatesData, 'month'))
negTab = table(cut(negDatesData, 'month'))
# Format as data.frame
posPlotData = data.frame(DATE=as.Date(names(posTab)), FREQ_RECORDS=as.vector(posTab), CRITERION="positive")
negPlotData = data.frame(DATE=as.Date(names(negTab)), FREQ_RECORDS=as.vector(negTab), CRITERION="negative")
## Concatenate plot data parts
plotData = rbind(posPlotData, negPlotData)
## Resort plot data
plotData = plotData[order(plotData$DATE),]
## Add a column that displays the growing cumulative sum for each criterion
plotData = plotData %>% group_by(CRITERION) %>% mutate(CUMFREQ=cumsum(FREQ_RECORDS))

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=CUMFREQ, fill=forcats::fct_rev(CRITERION)), width=1) +  #forcats::fct_rev inverts the order
    geom_bar(stat="identity", position="stack", alpha=0.5) # stacked barcharts
    #geom_bar(stat="identity", position="dodge", alpha=0.5) # side-by-side barcharts

myPlot = base_plot + 
    xlab("Year") + 
    ylab("Cumulative Number of Records\n") + 
    #ggtitle("Cumulative number of complete plastid\ngenomes on NCBI GenBank per year,\nseparated by presence of IR annotation",
    #    subtitle="Note: Only data after 2009 is displayed."
    #) + 
    scale_x_date(
        limits=c(as.Date(paste(start_year, "-01-01", sep='')), as.Date("2020-12-31")),
        date_breaks="1 year",
        minor_breaks=NULL,
        expand=expansion(0),
        date_labels="%Y"
    ) + 
    scale_y_continuous(breaks=seq(0, 9000, 1000), minor_breaks=seq(500, 9500, 1000)) +
    #scale_colour_grey(aesthetics = "fill") + 
    #scale_fill_brewer(palette="Dark2", name="Criterion positive/negative") + 
    scale_fill_manual(values=c("grey50", "grey0"), name="IR annotations present", labels=c("Yes", "No")) + 
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=20, face="italic"),
          axis.text=element_text(size=16),
          axis.text.x = element_text(size=16, angle=45, hjust=0, vjust=0.5),
          axis.title=element_text(size=18, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
          legend.key.width=unit(1,"cm"),
          legend.position = "bottom",
          legend.text = element_text(size=14),
          legend.title = element_text(size=14))


ggsave(file = "./PlastomeAccumulationWithWithoutIRAnnos.pdf", plot=myPlot)
################################

#    myPlot_transformed = ggplot(data=inData, aes(x=factor(year), y=accessions, fill=voucher)) +
#        geom_bar(stat="identity", position="dodge", alpha=0.5) + ##side-by-side barcharts
#        #geom_bar(stat="identity", position="identity", alpha=0.25) + ##overlaid barcharts
#        #geom_bar(stat="identity", alpha=0.5) + ##stacked barcharts
#        theme_bw() +
#        scale_fill_brewer(palette="Dark2", name="Specimen voucher") + 
#        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#        xlab("\nYear") + 
#        ylab("Number of accessions\n") +
#        labs(#tag="Log-transformed y-axis",
#             title="Accession numbers of new submissions of \ngenomic DNA to ENA per submission year",
#             subtitle="Log-transformed y-axis") #+
#        #ggtitle("Accession numbers of new submissions of genomic DNA to ENA per submission year") +
#        #ggsubtitle("Only submissions to the plant database have been counted.\nNote: Y-axis has been square-root transformed for better visibility")
#    myPlot_transformed = myPlot_transformed + annotation_logticks(base=10, sides="l")

########################################################################

assign(script_name, myPlot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
