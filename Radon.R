# Libraries
library(sf)
library(tigris)
library(ggplot2)

# Load the data
df<-read.table("Data/Radon_MN.csv",header = T)

# Cleaning the data
df$county.name <- toupper(trimws(df$county.name))
df$county.name[df$county.name == "ST LOUIS"] <- "ST. LOUIS"
df$state2 <- NULL

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

mn_map_joined <- merge(
  mn_map,
  county_region_map,
  by = "county.name",
  all.x = TRUE
)

ggplot(mn_map_joined) +
  geom_sf(aes(fill = region), color = "white", size = 0.2) +
  scale_fill_manual(values = region_colors) +
  theme_void()


# Plotting means
# County names and regions aligned to tapply indexing
county_names <- tapply(df$county.name, df$county, function(x) x[1])
county_regions <- tapply(df$region, df$county, function(x) x[1])

plot_df <- data.frame(
  county.name = county_names,
  region = county_regions,
  mean_y = as.numeric(cty.mns),
  stringsAsFactors = FALSE
)

region.order <- c("NW", "N", "NE", "E", "M", "W", "SW", "S", "SE")

plot_df$region <- factor(plot_df$region, levels = region.order)
plot_df <- plot_df[order(plot_df$region, plot_df$mean_y), ]

plot_df$county.order <- factor(plot_df$county.name,
                               levels = plot_df$county.name)


ggplot(plot_df, aes(x = county.order, y = mean_y, color = region)) +
  geom_point(size = 3) +
  scale_color_manual(values = region_colors) +
  labs(
    title = "Mean radon per county",
    x = "County",
    y = "Mean radon"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )









