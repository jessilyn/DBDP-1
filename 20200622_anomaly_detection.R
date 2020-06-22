############################################
#           Lymes Analysis                 #
############################################

library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)


is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

setwd("/Users/jessilyn/Documents/iPOP/Lyme")
factors <- c("skin_temp", "hr", "gsr", "accel_magnitude", "steps")
axesTitles <- c("Skin Temperature (F)", "Heart Rate (BPM)", "Galvanic Skin Response", "Acceleration", "Steps")
plotColors <- c("Sky Blue", "Red", "Grey40", "Green", "Purple")
plots <- list() 

# B1 data (multiple watches) from 5/23/2014 - 12/2/2014;  Peak data (one watch) from 5/20/2015 - 11/12/2015
# Dir1merged_head.csv:  mpsnyder_X.csv	                      5/23/14	-  8/25/14	96 files (daily)				
# Dir2merged_head.csv:  mpsnyder@stanford.edu_X.csv	          8/31/14 -	12/2/14	  38 files (daily)
# Dir3merged_head.csv:  mpsnyder_2015-01-02_data.csv	        1/2/15  -  6/7/15	  156
# MPSynder_Basis_2011-01-01T08-00-00_2015-11-12T01-48-00.csv	5/20/15	-  11/12/15	1 <-- will replace this with more current file soon				

#rectangle <- data.frame(xmin=as.Date("2014-05-24", format = "%Y-%m-%d"), xmax=as.Date("2015-06-06", format = "%Y-%m-%d"))
#, ymin=-Inf, ymax=Inf)
skinTempOutliers2 <- data.frame()


# * generate separate average and SD for each device
#updated latest data from Mike's watch with Dir5_merged
Mikesfile <- c("Dir1merged_head.csv", "Dir2merged_head.csv", "Dir3merged_head.csv", "Dir5merged_head.csv")

for (i in 1:3) {
  #i = 4 or 1:length(factors)
  print(factors[i])
  combo1<-data.frame()
  combo2<-data.frame()
  combo3<-data.frame()
  combo4<-data.frame()
  for (m in 1:length(Mikesfile)){
    print(Mikesfile[m])
    #m = 4
    rm(list=setdiff(ls(), c("combined", "factors", "combo1", "combo2", "combo3", "combo4", "axesTitles", "i", "m", "Mikesfile", "plotColors", "par", "plots", factors, "rectangle", "skinTempOutlierDays2")))
    combined <- read.csv(Mikesfile[m], header = TRUE, na.strings = "NA", stringsAsFactors=F) # may need to refine this to column types etc.
    combined$time <- as.Date(combined$time, "%Y-%m-%d")
    voi <- combined[[ factors[i] ]]
    
    a <- mean(voi, na.rm =TRUE) # average of voi
    s <- sd(voi, na.rm =TRUE) # stdev of voi
    write.table(data.frame(factors[i],a,s), paste(m, "summary.txt", sep=""), append=TRUE, sep = "\t", col.names = F, row.names = F)
    lowlimit <- a - 2*s  # values more than 2 stdevs below mean
    uplimit <-a + 2*s    # values more than 2 stdevs above mean
    midlowlimit <- a - 1*s # values more than 1 stdev below mean
    midupperlimit <- a + 1*s # values more than 1 stdev above mean
    breaks = c(0, lowlimit, midlowlimit, a, midupperlimit, uplimit, 1000)
    yticks = breaks[2:7]
    lab = c("lower", "lower-mid", "mid", "upper-mid", "upper", "max")
    
    # add if statement here to correct for negative values that are not physically possible
    if(lowlimit < 0) {
      lowlimit = 0
      if (midlowlimit > 0){
        breaks = c(0, midlowlimit, a, midupperlimit, uplimit, 1000)}
      yticks = breaks[2:6]
      if(midlowlimit <= 0){
        midlowlimit = 0
        breaks = c(0, a, midupperlimit, uplimit, 1000)
        yticks = breaks[2:5]
        lab = c("mid", "upper-mid", "upper", "max")}
    }
    
    mydf <- data.frame(Date = as.Date(c(combined$time), format = "%Y-%m-%d"), Value = c(voi), stringsAsFactors = FALSE) # define new data frame for plot
    foo <- expand.grid(Date = as.Date(c(unique(combined$time)), format = "%Y-%m-%d"),
                       group = lab) #expands out each date to create a data frame with 5 lines per each date (each corresponding to a different group e.g. lower) so that we can populate it later with the sums
    mutated = mutate(group_by(mydf, Date), group = cut(Value, breaks = breaks, labels = lab))
    c <- count(mutated, Date, group)  # count the number of times that particular date appears in each group
    
    ## ** should probably add in threshold here for min. total # reads in a day to be considered a day we should use.
    ndf <- left_join(foo, c, by = c("Date" = "Date", "group" = "group")) 
    names(ndf)[names(ndf)=="n"] <- "Total" # renames the values from above as Total
    out <- mutate(ndf, Total = replace(Total, which(Total %in% NA), 0)) #replaces NA values from the counting with 0 and creates df out
    
    proportion<-tapply(out$Total, out$Date, FUN=function(x){x/sum(x)}); # get proportions to make relative size/color of each dot
    propUnlist <- as.vector(unlist(proportion))
    sorted <- out[order(as.Date(out$Date, format = "%Y-%m-%d")),]
    combo <- cbind(sorted, propUnlist)
    combo<-combo[order(combo$propUnlist, decreasing = TRUE), ]
    combo[is.na(combo)] <- 0 #changed this from is.nan(combo)
    newname <-factors[i]
    assign(newname, combo)
    assign(paste("combo", m, sep=""), combo)
  }

    ####################
    #   OUTLIER DAYS   #
    ####################
  #if (factors[i] %in% "skin_temp"){
  if (factors[i] %in% "hr"){
    #rm(skinTempOutlierDays)
    skinTempOutlierDays <- structure(rep(NA_real_, 1 ), class="Date") # Make vector of "outlier" days for highlighting plots (based on proportion of points  (random threshold = n ) for that day falling 2 or more sds away from the mean): choose a threhsold here for what makes a day "strange." can make this informed based on the data (i.e. what are average values and SDs of proportions of daily values 1 or 2 sd out?)
    #print(factors[i])
    for(j in 1:length(combo$Date)){
      if (combo$propUnlist[j] > 0.04 & combo$group[j] != "lower" & combo$group[j] != "mid" & combo$group[j] != "lower-mid" & combo$group[j] != "upper-mid" & combo$group[j] != "upper"){
        print(combo[j,]) 
        skinTempOutlierDays <- c(skinTempOutlierDays, as.Date(combo$Date[j]))
      }
    }
    skinTempOutlierDays <- na.omit(skinTempOutlierDays)
    skinTempOutlierDays2 <- as.data.frame(skinTempOutlierDays)
    #skinTempOutlierDays <- as.Date(c("2015-07-04"))
  }

  
  
  ########
  # PLOT #
  ########
  rm(combo)
  combo <- rbind(combo1, combo2, combo3, combo4)
  #combo <- rbind(combo1, combo2, combo3)
  fil = paste(factors[i],".pdf", sep="")
  
  # will want to offset y axis so it captures the bin values (i.e. draw bins on outside of ticks)
  p1 <- ggplot(combo, aes(x = Date, y = as.integer(group)+0.5, size = propUnlist)) +
    geom_point(aes(color = combo$propUnlist), alpha = 0.8, show.legend = FALSE) +
    scale_x_date(breaks = date_breaks("2 weeks"), minor_breaks = date_breaks("1 day"), labels = date_format("%m-%d-%Y")) +  
    scale_y_continuous(name = axesTitles[i], breaks=c(2:(length(yticks)+1)), labels=round(yticks, digits=1), expand=c(0.1,0.1), minor_breaks = NULL) +
    theme(plot.title = element_text(face="bold", vjust=1),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(face="bold",angle=60, hjust = 1, size=10),
          axis.title.y = element_text(face="bold", size=10, vjust=1.8),
          axis.text.y = element_text(face="bold",size=10),
          panel.grid.major = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "grey"),
          panel.background = element_rect(fill = "white"),
          aspect.ratio=1/7) + scale_colour_continuous(high = plotColors[i], guide = FALSE) + 
          geom_rect(data=skinTempOutlierDays2, inherit.aes = FALSE, show.legend = FALSE, aes( xmin=skinTempOutlierDays - 0.5, xmax=skinTempOutlierDays + 0.5, ymin=-Inf, ymax=Inf), fill="yellow", linetype = 0, alpha=0.3)  #HIGHLIGHT OUTLIER DAYS IN PLOT # to pick a specific date to highlight: xmin=as.Date("2014-07-31", format = "%Y-%m-%d"), xmax=as.Date("2014-08-01", format = "%Y-%m-%d")
    
  #p2 = p1 + geom_rect(data=rectangle, inherit.aes = FALSE, aes(x = NULL,y = NULL, xmin=start, xmax=end, ymin=-Inf, ymax=+Inf), fill='pink', alpha=0.2)
  #geom_rect(data=skinTempOutlierDays, show.legend = FALSE, aes( xmin=as.Date(skinTempOutlierDays - 0.5), xmax=as.Date(skinTempOutlierDays + 0.5), ymin=-Inf, ymax=Inf), fill="yellow", linetype = 0, alpha=0.3, inherit.aes = FALSE)  #HIGHLIGHT OUTLIER DAYS IN PLOT # to pick a specific date to highlight: xmin=as.Date("2014-07-31", format = "%Y-%m-%d"), xmax=as.Date("2014-08-01", format = "%Y-%m-%d")

  ggsave(filename=fil, plot = p1, path = "/Users/jessilyn/Desktop", width = 15)
  graphics.off()

}





######################################################
# plot specific measurement type from specific dates #
######################################################
combined <- read.csv("Dir5merged_head.csv", header = TRUE, na.strings = "NA", stringsAsFactors=F) # may need to refine this to column types etc.


combined <- read.csv("Dir5merged_head.csv", header = TRUE, na.strings = "NA", stringsAsFactors=F) # may need to refine this to column types etc.
# combined$time <- as.Date(combined$time, "%Y-%m-%d")
# for plotting values over the course of 1 day, you need to preserve clock time
combined$time <- as.POSIXlt(combined$time, "%m/%d/%y %H:%M", tz="GMT")

# define 3 windows: a wide window (June-July 2015), a narrow window, and a single day (we expect Wednesday (Wed is 07-01-2015 or 07-08-2015) to be the greatest outlier)
t <-which(combined$time >= as.Date("2015-05-30") & combined$time <= as.Date("2015-07-30"))
wideWindowToPlot <- combined[t, ]
p <-which(combined$time >= as.Date("2015-06-30") & combined$time <= as.Date("2015-07-05"))
narrowWindowToPlot <- combined[p, ]


factors <- c("skin_temp", "hr", "gsr", "accel_magnitude", "steps")
axesTitles <- c("Skin Temperature (F)", "Heart Rate (BPM)", "Galvanic Skin Response", "Acceleration", "Steps")
plotColors <- c("Sky Blue", "Red", "Grey40", "Green", "Purple")
i=2
par(mar = rep(2, 4))
par(mfrow=c(1,5))
for (i in 3:6){
  boxplot(combined[i], main=names(combined[i]))
}



## Zoom in on specific days here to find hours where HR or ST is higher than the annual average (calculated during non activity)
# get average and SD for activity. define "bouts" of activity. 
#change time to be time of day "clock time", not sure if this is possible without ggplot
# will use this method below for the "daily average" plots to track cycles over the day later.
#j <-which(combined$time == as.POSIXct("2015-07-01"))
# make list fo dates
before <-c("2015-06-29", "2015-06-30", "2015-07-01", "2015-07-02", "2015-07-03")
lyme <- c("2015-07-04", "2015-07-05", "2015-07-06", "2015-07-07", "2015-07-08")
after <- c("2015-07-09", "2015-07-10", "2015-07-11", "2015-07-12", "2015-07-13")
dates <- c(before, lyme, after)
anno_dates <- c(rep("before", 5), rep("lyme", 5), rep("after", 5))
plots<-list()

i=2
  # for all dates
  for (k in 1:length(dates)){ 
    par(mfrow=c(5,3))
  #k=1
  #i=4
  lymeday <- dates[k]
  print(lymeday)
  #plotname <- paste(lymeday, factors[i], sep="_")
  j <-which(as.Date(combined$time) == lymeday)
  singleDay <- combined[j,]
  # check if steps is above certain threshold (10 steps) for the other sensor values
  correctedSingleDay<-singleDay[singleDay[,6]<10,]
  
  plots[[k]]<- ggplot(correctedSingleDay, aes(x=time,y=hr)) + 
  geom_point(size=0.5, color="darkred") +
  scale_x_datetime(breaks = date_breaks("3 hour"), labels=date_format("%H:%M"), name=NULL) + 
  scale_y_continuous(limits=c(50,140), name=NULL) +
  # fix y limits
  #scale_y_continuous(name = axesTitles[i]) +
  ggtitle(lymeday) +
  theme(plot.title = element_text(face="bold", vjust=1),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold",angle=90, hjust = 1, size=10),
        axis.title.y = element_text(face="bold", size=10, vjust=1.8),
        axis.text.y = element_text(face="bold",size=10))
#         panel.grid.major = element_line(colour = "grey"),
#         panel.grid.minor = element_line(colour = "grey"),
#         panel.background = element_rect(fill = "white"))
  }
 do.call(grid.arrange, c(plots[as.vector(gdata::interleave(1:5, 6:10, 11:15))], 
                                 list(nrow = 5))) 

 plots<-list()
 
 i=1
 # for all dates
 for (k in 1:length(dates)){ 
   par(mfrow=c(5,3))
   lymeday <- dates[k]
   print(lymeday)
   #plotname <- paste(lymeday, factors[i], sep="_")
   j <-which(as.Date(combined$time) == lymeday)
   singleDay <- combined[j,]
   # check if steps is above certain threshold (10 steps) for the other sensor values
   correctedSingleDay<-singleDay[singleDay[,6]<10,]
   
   plots[[k]]<- ggplot(correctedSingleDay, aes(x=time,y=skin_temp)) + 
     geom_point(size=0.5, color="darkblue") +
     scale_x_datetime(breaks = date_breaks("3 hour"), labels=date_format("%H:%M"), name=NULL) + 
     scale_y_continuous(limits=c(75,100), name=NULL) +
     # fix y limits
     #scale_y_continuous(name = axesTitles[i]) +
     ggtitle(lymeday) +
     theme(plot.title = element_text(face="bold", vjust=1),
           axis.title.x = element_blank(),
           axis.text.x  = element_text(face="bold",angle=90, hjust = 1, size=10),
           axis.title.y = element_text(face="bold", size=10, vjust=1.8),
           axis.text.y = element_text(face="bold",size=10))
 }
 do.call(grid.arrange, c(plots[as.vector(gdata::interleave(1:5, 6:10, 11:15))], 
                         list(nrow = 5))) 
  
  
  
# make these into box plots
boxplot(wideWindowToPlot[[ factors[i] ]] ~ as.Date(wideWindowToPlot$time), ylab=axesTitles[i], cex.names=0.5, las=3, outline=FALSE, col=plotColors[i])
#"skyblue" 
# ylim = c(80,100),
boxplot(combined[[ factors[i] ]] ~ combined$time, ylab=axesTitles[i], cex.names=0.5, outline=FALSE, col=plotColors[i])

boxplot(combined[[ factors[i] ]] ~ combined$time, ylab=axesTitles[i], xaxt="n", xlab="", outline=FALSE, col=plotColors[i])
axis(1, labels = FALSE)
text(x=combined$time, srt = 45, xpd = TRUE)





ggplot(data, aes(x=factor(Purpose), y=Rate)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))

ggplot(data = wideWindowToPlot, aes(x= as.Date(wideWindowToPlot$time)), y= wideWindowToPlot[[ factors[i] ]]) + geom_boxplot()

ggplot(wideWindowToPlot, aes(x = as.Date(wideWindowToPlot$time), y = wideWindowToPlot[[ factors[i] ]])) +
  geom_point(alpha = 0.8) +
  scale_x_date(breaks = "1 week", minor_breaks = "1 day", labels = date_format("%m-%d-%Y")) +  
  scale_y_continuous(name = axesTitles[i]) +
  theme(plot.title = element_text(face="bold", vjust=1),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold",angle=60, hjust = 1, size=10),
        axis.title.y = element_text(face="bold", size=10, vjust=1.8),
        axis.text.y = element_text(face="bold",size=10),
        panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "grey"),
        panel.background = element_rect(fill = "white"))

# make these into box plots
ggplot(narrowWindowToPlot, aes(x = as.Date(narrowWindowToPlot$time), y = narrowWindowToPlot[[ factors[i] ]])) +
  geom_point(alpha = 0.8) +
  scale_x_date(breaks = "1 day", labels = date_format("%m-%d-%Y")) +  
  scale_y_continuous(name = axesTitles[i]) +
  theme(plot.title = element_text(face="bold", vjust=1),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(face="bold",angle=60, hjust = 1, size=10),
        axis.title.y = element_text(face="bold", size=10, vjust=1.8),
        axis.text.y = element_text(face="bold",size=10),
        panel.grid.major = element_line(colour = "grey"),
        panel.grid.minor = element_line(colour = "grey"),
        panel.background = element_rect(fill = "white"))
plot(narrowWindowToPlot[[ factors[i] ]], type = "l", pch = 20, xlab= "dates", ylab = axesTitles[i]) 




# cyclical behavior of value over 1 day, 1 week, or 1 month
# values of peaks and troughs; velocity of change of peaks and troughs; etc. (other interesting patterns?) use previous binning method to account for noisy data. 
# if cyclical behavior is very slight, then we could shift the bins to something like 0.1*SD, 0.2*SD, etc.
