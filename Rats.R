# Libraries
library(sf)
library(tigris)
library(ggplot2)
library(patchwork)
library(lmerTest)
library(MASS)
library(lattice)
library(predictmeans) #residplot
library(plyr)

#####################################################
# Cleaning data
#####################################################

# Load the data
rats <- read.table("Data/Rats_Box.txt", header = TRUE, sep = ",")

# Reshaping to long format
rats_long <- reshape(
  rats,
  varying = list(c("y0","y1","y2","y3","y4")),
  v.names = "y",
  timevar = "week",
  times = 0:4,
  direction = "long"
)

rats_long <- rats_long[order(rats_long$Rat, rats_long$week), ]
row.names(rats_long) <- NULL

# Now storing the correct format
rats <- rats_long

# make treatment and cage factors 
rats$Trt <- as.factor(rats$Trt)
rats$Rat <- as.factor(rats$Rat)


# Investigate data
str(rats) 

# make two versions of the time variable 
# - one quantitative and one qualitative
#rats$monthQ <- rats$month
#rats$month <- factor(rats$month)
#-------

rats$weekQ <- as.numeric(rats$week)
rats$weekF <- as.factor(rats$week)

p1<-ggplot(rats, aes(x=rats$weekQ, y=y, group=rats$Rat, colour=Trt)) + 
  geom_line()
p1

mns <- ddply(rats, ~ Trt + weekF + weekQ, summarize, 
             y = mean(y))
p2<-ggplot(mns, aes(x=weekQ, y=y, group=Trt, colour=Trt)) + 
  geom_point() + geom_line()
p2


p1 / p2











