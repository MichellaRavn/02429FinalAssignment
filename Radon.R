# Libraries
library(sf)
library(tigris)
library(ggplot2)
library(patchwork)
library(lmerTest)
library(MASS)
library(lattice)
library(predictmeans) #residplot

#####################################################
# Cleaning data
#####################################################

# Load the data
df<-read.table("Data/Radon_MN.csv",header = T)

# Cleaning the data
df$county.name <- toupper(trimws(df$county.name))
df$county.name[df$county.name == "ST LOUIS"] <- "ST. LOUIS"
df$state2 <- NULL
names(df)[names(df) == "u.full"] <- "u"
df$x <- as.factor(df$x)

# Renaming basement, 0 for no basement, 1 for basement
df$b <- ifelse(df$basement == "Y", 1, ifelse(df$basement == "N", 0, NaN))

# If the measurement is in the basement, then b is true/1
df$b[df$x == 0] <- 1

# Consider assuming the ground floor measurements should assume no basement
#df$b[df$x == 1 & is.nan(df$b)] <- 0

# Drop NaNs for subanalysis
dfsub <- df[!is.nan(df$b), ]


# Creating a new grouping of basement
dfsub$bg<- ifelse(
  dfsub$b == 1 & dfsub$x == 0, "MB",
  ifelse(dfsub$b == 1 & dfsub$x == 1, "BNM",
         "NB"))

dfsub$bg<- as.factor(dfsub$bg)
df$b <- as.factor(df$b)

# Overview of county (un)balance
table(df$county)

#####################################################
# Exploratory plots
#####################################################

# Geographic placement of the counties
SW<-c("ROCK", "NOBLES","JACKSON","COTTONWOOD","MURRAY","PIPESTONE","LINCOLN","LYON","REDWOOD","YELLOW MEDICINE")
W<-c("LAC QUI PARLE","CHIPPEWA","SWIFT","BIG STONE","STEVENS","WILKIN","TRAVERSE","GRANT","CLAY","OTTER TAIL","BECKER", "POPE","DOUGLAS") #13
NW<-c("NORMAN","MAHNOMEN","POLK","RED LAKE","PENNINGTON","MARSHALL","KITTSON","ROSEAU","CLEARWATER") # 9
N<-c("LAKE OF THE WOODS","KOOCHICHING","BELTRAMI","ITASCA") # 4
NE<-c("ST. LOUIS","LAKE","COOK") # 3
E<-c("CARLTON","PINE","CHISAGO","ISANTI","AITKIN","KANABEC","MILLE LACS","CHISAGO","BENTON","SHERBURNE","ANOKA","WASHINGTON","RAMSEY","HENNEPIN","WRIGHT","CARVER")
SE<-c("DAKOTA","GOODHUE","WABASHA","WINONA","HOUSTON","FILLMORE","OLMSTED","MOWER","FREEBORN","DODGE","STEELE","RICE","WASECA","SCOTT","LE SUEUR")
S<-c("MARTIN","FARIBAULT","WATONWAN","BLUE EARTH","BROWN","NICOLLET","SIBLEY","RENVILLE")
M<-c("KANDIYOHI","MCLEOD","MEEKER","STEARNS","TODD","MORRISON","CROW WING","WADENA","CASS","HUBBARD")

regions <- list(
  SW = SW, W = W, NW = NW,
  N = N, NE = NE, E = E,
  SE = SE, S = S, M = M)

# Assigning region to each county
lookup <- rep(names(regions), sapply(regions, length))
names(lookup) <- unlist(regions)
df$region <- lookup[df$county.name]

# Check for unassigned counties:
#sum(is.na(df$region))
# Unassigned counties:
#unique(df$county_clean[is.na(df$region)])


# Storing color values
# Blue scale (north)
col_NW <- "#A6CEE3"   # light blue
col_N  <- "#1F78B4"   # medium blue
col_NE <- "#0B3C68"   # dark blue

# Green scale (middle)
col_W <- "#B2DF8A"    # light green
col_M <- "#33A02C"    # medium green
col_E <- "#145214"    # dark green

# Orange scale (south)
col_SW <- "#FDBF6F"   # light orange
col_S  <- "#FF7F00"   # medium orange
col_SE <- "#B15928"   # dark orange

# Combine into named vector
region_colors <- c(
  NW = col_NW,
  N  = col_N,
  NE = col_NE,
  W  = col_W,
  M  = col_M,
  E  = col_E,
  SW = col_SW,
  S  = col_S,
  SE = col_SE
)
region.order <- c("NW", "N", "NE", "W", "M", "E", "SW", "S", "SE")


# Provided code:
# Counties
(J <- length(unique(df$county.name)))

# Mean of response
(ybarbar = mean(df$y))

(sample.size <- as.vector(table(df$county)))
(sample.size.jittered <- sample.size*exp(runif (J, -.1, .1)))
(cty.mns = tapply(df$y,df$county,mean))
(cty.vars = tapply(df$y,df$county,var))
(cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size))
(cty.sds.sep = sqrt(tapply(df$y,df$county,var)/sample.size))


# Creating a map plot
options(tigris_use_cache = TRUE)


mn_map <- counties(state = "MN", year = 2020)
mn_map <- st_as_sf(mn_map)
mn_map$NAME <- toupper(trimws(mn_map$NAME))

names(mn_map)[names(mn_map) == "NAME"] <- "county.name"

county_region_map <- unique(df[, c("county.name", "region")])
county_number_map <- unique(df[, c("county.name", "county")])

mn_map_joined <- merge(
  mn_map,
  merge(county_region_map, county_number_map, by="county.name"),
  by = "county.name",
  all.x = TRUE
)

mn_map_joined$region <- factor(mn_map_joined$region,
                               levels = region.order)



centroids <- st_centroid(mn_map_joined)
highlight <- centroids[centroids$county %in% c(36, 37), ]

ggplot() +
  geom_sf(data = mn_map_joined, aes(fill = region), color = "white", size = 0.2) +
  
  geom_sf(
    data = highlight,
    color = "black",
    #fill = "red",
    size = 5,
    shape = 21,
    stroke = 1.2
  ) +
  
  # Uncomment for numbers added
  #geom_sf_text(
  #  data = centroids,
  #  aes(label = county),
  #  size = 3,
  #  color = "black"
  #) +
  
  scale_fill_manual(values = region_colors,
                    breaks = region.order ) +
  theme_void()



# Plotting means
# County names and regions aligned to tapply indexing for plotting
county_names <- tapply(df$county.name, df$county, function(x) x[1])
county_regions <- tapply(df$region, df$county, function(x) x[1])

plot_df <- data.frame(
  county.name = county_names,
  region = county_regions,
  mean_y = as.numeric(cty.mns),
  stringsAsFactors = FALSE
)

plot_df$region <- factor(plot_df$region, levels = region.order)
plot_df <- plot_df[order(plot_df$region, plot_df$mean_y), ]

plot_df$county.order <- factor(plot_df$county.name,
                               levels = plot_df$county.name)


ggplot(plot_df, aes(x = county.order, y = mean_y, color = region)) +
  geom_point(size = 3) +
  scale_color_manual(values = region_colors) +
  labs(
    title = "Observed mean radon per county",
    x = "County",
    y = "Mean log radon"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )




# The given plot about samplesize with a few alterations
par(mfrow=c(1,1))

# First an empty plot
plot(sample.size.jittered, cty.mns,
     type = "n",                       
     xlab="sample size in county j",
     ylab="Mean log radon in county j",
     log="x", cex.lab=.9, cex.axis=1, mgp=c(1.5,.5,0),
     ylim=c(0,3.2), yaxt="n", xaxt="n")

axis(1, c(1,3,10,30,100), cex.axis=.9, mgp=c(1.5,.5,0))
axis(2, seq(0,3), cex.axis=.9, mgp=c(1.5,.5,0))

# Draw error bars 
for (j in 1:J){
  lines(rep(sample.size.jittered[j],2),
        cty.mns[j] + c(-1,1)*cty.sds[j], lwd=.5)
}

# Adding points after to have error bars behind
points(sample.size.jittered, cty.mns,
       pch = 20, cex = 1.3,
       col = region_colors[ df$region[ match(names(cty.mns), df$county) ] ])

# Horizontal mean line
abline(h = mean(cty.mns), lwd = .5)

title("Observed mean response values per county", cex.main=.9, line=1)

# The max and min
cty.mns[cty.mns==max(cty.mns)]
cty.mns[cty.mns==min(cty.mns)]

county_names[cty.mns==max(cty.mns)] # Lac Qui Parle
county_names[cty.mns==min(cty.mns)] # Lake

# highlight max/min if needed
points(sample.size.jittered[36], cty.mns[36], cex=2)
points(sample.size.jittered[37], cty.mns[37], cex=2)


## Uranium plots
floor_colors <- c(
  "1" = "gold2",   # ground floor
  "0" = "grey37"      # basement
)

# Plot 1: Region
p1 <- ggplot(df, aes(x = u, y = y, color = region)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = region_colors, breaks=region.order) +
  labs(
    title = "Radon vs uranium - by region",
    x = "log uranium (u)",
    y = "log radon (y)"
  ) +
  theme_bw()

# Plot 2: Floor
p2 <- ggplot(df, aes(x = u, y = y, color = factor(x))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = floor_colors, name = "Floor (x)") +
  labs(
    title = "Radon vs uranium - by floor",
    x = "log uranium (u)",
    y = "log radon (y)"
  ) +
  theme_bw()

# Combine
p1 | p2

table(df$x)


## Investigate floor status
p_strip <- ggplot(df, aes(x = x, y = y, color = x)) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  stat_summary(fun = mean, geom = "point",
               size = 3, shape = 21, fill = "black", color = "black") +
  scale_color_manual(values = floor_colors) +
  labs(
    title = "Log radon - by floor",
    x = "Floor (x)",
    y = "Log radon (y)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )

p_box <- ggplot(df, aes(x = x, y = y, fill = x)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = floor_colors) +
  labs(
    title = "Log radon – boxplot by floor",
    x = "Floor (x)",
    y = "Log radon (y)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )

p_box / p_strip

#####################################################
# Modelling
#####################################################

# Testing basement for personal interest
m0 <- lmer(y ~ u:x + u + x + ( 1 | county), data = dfsub) 
mb <- lmer(y ~ u:b + u + b + ( 1 | county), data = dfsub)
mbg <- lmer(y ~ u:bg+ u + bg+( 1 | county), data = dfsub)
m1 <- lmer(y ~ u:x + u + I(u^2) + x + ( 1 | county), data = dfsub)

ranova(m0)
ranova(mb)
ranova(mbg)
ranova(m1)

m0 <- update(m0,REML=F)
mb <- update(mb,REML=F)
mbg <- update(mbg,REML=F)


AIC(m0, mb, mbg, m1)
BIC(m0, mb, mbg, m1)

anova(mbg,mb) # Prefer mbg
anova(mbg,m0) # Prefer m0 (basement no effect)
anova(m0,m3) # Prefer m0

# Moving on with m0, full data
m0 <- lmer(y ~ u*x + u + x + ( 1 | county), data = df, REML=TRUE) 
ranova(m0)
m0 <- update(m0,REML=F)
drop1(m0) # interaction significant 

m0 <- update(m0,REML=T)
summary(m0)

residplot(m0) 
residplot(mb) 
residplot(mbg) 
residplot(m2) 


m0<-update(m0,REML=T)
summary(m0) 


# Cooks distance
cd <- cooks.distance(m0)
plot(cd, type = "h", main = "Cook's distance", xlab = "Observation", ylab = "Distance")
df[cd>0.15,c("y","county","county.name","region","u","b")]







# Final model results
mfinal <- m0


residplot(mfinal)

#ranef(mfinal)

confint(mfinal)

