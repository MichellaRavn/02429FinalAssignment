# Libraries
library(sf)
library(tigris)
library(ggplot2)
library(patchwork)
library(lmerTest)
library(MASS)
library(lattice)
library(performance)
library(predictmeans) #residplot
library(emmeans)
library(scales)


## Cleaning data

# Load the data
df<-read.table("Data/Radon_MN.csv",header = T)

# Cleaning the data
df$county.name <- toupper(trimws(df$county.name))
df$county.name[df$county.name == "ST LOUIS"] <- "ST. LOUIS"
df$state2 <- NULL
names(df)[names(df) == "u.full"] <- "u"
df$x <- as.factor(df$x)
df$county <- as.factor(df$county)

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

# Overview of county (un)balance (sample sizes)
table(df$county)


## Exploratory plots


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
#dev.off()   


mn_map <- counties(state = "MN", year = 2020)
mn_map <- st_as_sf(mn_map)
mn_map$NAME <- toupper(trimws(mn_map$NAME))

names(mn_map)[names(mn_map) == "NAME"] <- "county.name"

county_region_map <- unique(df[, c("county.name", "region")])
county_number_map <- unique(df[, c("county.name", "county")])

mn_map_joined <- merge(mn_map,
  merge(county_region_map, county_number_map, by="county.name"),
  by = "county.name",
  all.x = TRUE)

mn_map_joined$region <- factor(mn_map_joined$region,
                               levels = region.order)
centroids <- st_centroid(mn_map_joined)
highlight <- centroids[centroids$county %in% c(36, 37), ]

ggplot() +
  geom_sf(data = mn_map_joined, aes(fill = region), color = "white", size = 0.2) +
  
  # Uncomment for min/max circles
  #geom_sf(
  #  data = highlight,
  #  color = "black",
  #  size = 5,
  #  shape = 21,
  #  stroke = 1.2
  #) +
  
  # Uncomment for numbers added
  geom_sf_text(
    data = centroids,
    aes(label = county),
    size = 3,
    color = "black" ) +
  
  scale_fill_manual(values = region_colors,
                    breaks = region.order ) +
  theme_void()


# Plotting means
# County names and regions aligned to tapply indexing for plotting
county_names <- unique(df$county.name)#tapply(df$county.name, df$county, function(x) x[1])
county_regions <- tapply(df$region, df$county, function(x) x[1])

plot_df <- data.frame(
  county = county_names,
  region = factor(county_regions, levels = region.order),
  mean_y = as.numeric(cty.mns),
  sd_y   = as.numeric(cty.sds.sep))

plot_df$lower <- plot_df$mean_y - plot_df$sd_y
plot_df$upper <- plot_df$mean_y + plot_df$sd_y
plot_df <- plot_df[order(plot_df$region, plot_df$mean_y), ]
plot_df$county.order <- factor(plot_df$county, levels = plot_df$county)


ggplot(plot_df, aes(x = county.order, y = mean_y, color = region)) +
  
  # Error bars first (behind points)
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.0,       # no horizontal bar
                color = "black",   # or same color as points if preferred
                alpha = 0.6,       # slightly transparent
                linewidth = 0.5) +
  
  # Points on top
  geom_point(size = 3) +
  geom_hline(yintercept = mean(cty.mns), linetype = "dashed", linewidth = 0.5) +
  
  scale_color_manual(values = region_colors) +
  labs(
    title = "Observed mean radon per county",
    x = "County",
    y = "Mean log radon") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
        cty.mns[j] + c(-1,1)*cty.sds[j], lwd=.5)}

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

#### Modelling

# Testing basement factor
m0 <- lmer(y ~ u:x + u + x + ( 1 | county), data = dfsub) 
mb <- lmer(y ~ u:b + u + b + ( 1 | county), data = dfsub)
mbg <- lmer(y ~ u:bg+ u + bg + ( 1 | county), data = dfsub)


m1 <- lmer(y ~ u:x + u + I(u^2) + x + ( 1 | county), data = dfsub)

mrs<-lmer(y ~ u:x + u + x + ( 1 + x | county), data = df) 
ranova(mrs)

ranova(m0)
ranova(mb)
ranova(mbg)
ranova(m1)
ranova(mrs)

m0 <- update(m0,REML=F)
mb <- update(mb,REML=F)
mbg <- update(mbg,REML=F)


AIC(m0, mb, mbg)
BIC(m0, mb, mbg)

anova(mbg,mb) # Prefer mbg
anova(mbg,m0) # Prefer m0 (basement no effect)
anova(m0,m1) # Prefer m0


residplot(m0) 
residplot(mb) 
residplot(mbg) 
residplot(m1) 

# Moving on with m0, full data
m0 <- lmer(y ~ u*x + u + x + ( 1 | county), data = df, REML=TRUE) 
ranova(m0)
m0 <- update(m0,REML=F)
drop1(m0) # interaction significant

AIC(m0)
BIC(m0)
r2_nakagawa(m0)
var_county   <- as.numeric(VarCorr(m0)$county)          # 9.323
var_resid <- sigma(m0)^2                         # 22.875
icc <- var_county / (var_county + var_resid)
icc


m0<-update(m0,REML=T)
summary(m0) 

## predictions 
df$pred <- predict(m0)

# Cooks distance
cd <- cooks.distance(m0)
plot(cd, type = "h", main = "Cook's distance", xlab = "Observation", ylab = "Distance")
df[cd>0.2,c("y","pred","county","county.name","region","u","x")]

# Random effect
dotplot(ranef(m0, condVar = TRUE),
        strip = FALSE,
        ylab = "",
        scales = list(cex = 0.5))   # smaller labels

## Model parameters
# Raw
est  <- fixef(m0)                   
ci   <- confint(m0, parm = names(est))  # 2.5% and 97.5%

# Back-transform (exp)
est_bt <- exp(est)
ci_bt  <- exp(ci)

#--- Combine into table ---
param_table <- data.frame(
  Parameter = names(est),
  Estimate_raw = round(est, 2),
  CI_low_raw = round(ci[,1], 2),
  CI_high_raw = round(ci[,2], 2),
  Estimate_bt = round(est_bt, 2),
  CI_low_bt = round(ci_bt[,1],2),
  CI_high_bt = round(ci_bt[,2], 2)
)

param_table




 ## 0.0399 ~ 4% of variance described by random effect of county
############## Results plots 

plot_df <- data.frame(
  county = county_names,
  region = factor(county_regions, levels = region.order),
  mean_y = as.numeric(cty.mns),
  sd_y   = as.numeric(cty.sds.sep)
)

# Mixed predictions
df$pred_y <- predict(m0)
# Fixed predictions
df$pred_y_fix <- predict(m0,re.form=NA)
# Random effect
df$pred_y_rand <- df$pred_y - df$pred_y_fix

(cty.mns.pred = tapply(df$pred_y,df$county,mean))
plot_df$mean_pred_y <- as.numeric(cty.mns.pred)

(cty.new.mns.pred = tapply(df$pred_y_fix,df$county,mean))
plot_df$mean_pred_y_fix <- as.numeric(cty.new.mns.pred)

plot_df <- plot_df[order(plot_df$region, plot_df$mean_y), ]
plot_df$county.order <- factor(plot_df$county, levels = plot_df$county)

ggplot(plot_df, aes(x = county.order)) +
  
  # Error bar for observed
  #geom_errorbar(
  #  aes(ymin = lower, ymax = upper, color = region),
  #  width = 0, alpha = 0.6, linewidth = 0.5
  #) +
  
  geom_hline(yintercept = mean(cty.mns), linetype = "dashed", linewidth = 0.5) +
  
  # Observed points (region colors)
  geom_point(aes(y = mean_y, color = region), size = 3) +
  
  # Predicted mixed (RED) – color mapped!
  geom_point(
    aes(y = mean_pred_y, color = "Predicted mixed"),
    size = 2
  ) +
  
  geom_errorbar(
    aes(ymin = pred_lower, ymax = pred_upper, color = "Predicted mixed"),
    width = 0, alpha = 0.4
  ) +
  
  
  # Predicted fixed (BLACK) – color mapped!
  #geom_point(
  #  aes(y = mean_pred_y_fix, color = "Predicted fixed"),
  #  size = 1.5
  #) +
  
  # One single merged legend
  scale_color_manual(
    name = "Legend",
    values = c(
      region_colors,
      "Predicted mixed"  = "red",
      "Predicted fixed"  = "black"
    )
  ) +
  
  labs(
    title = "Observed and predicted mean radon per county",
    x = "County",
    y = "Mean log radon"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )




### Emmeans

# 25%, 50%, 75% quantiles
u_vals <-quantile(df$u, probs = c(0.25, 0.50, 0.75))

# emmeans table
em <- emmeans(m0, "x", by = "u", at = list(u = u_vals),adjust="tukey")
em

# Back-transform to radon scale
em_bt <- transform(
  as.data.frame(em),
  em_bt      = exp(emmean),
  lower_bt   = exp(lower.CL),
  upper_bt   = exp(upper.CL)
)

em_bt


emmip(m0,x ~ u, at = list(u = u_vals) , CIs = TRUE) + theme(legend.position="top")+
  #ggtitle("Emmeans interaction")+
  ylab("Log radon") + 
  labs(color="Floor") + 
  scale_color_manual(values=floor_colors)

emtrends(m0, pairwise ~x, var="u",infer=TRUE,adjust="tukey")
# low change in radon c for higher val of u, 
# significant slopes, strong evidence for x=0
# significant difference


slopes <- emtrends(m0, ~ x, var = "u")  # dy/du
slopes_df <- as.data.frame(slopes)

ggplot(slopes_df, aes(x = x, y = u.trend)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(
    title = "Estimated slope",
    x = "Floor (x)",
    y = "u slope (dy/du)"
  ) +
  theme_bw()



### New county prediction example
new_data <- data.frame(
  u      = -1.2,
  x      = factor(0, levels = levels(df$x)),
  county = NA     # not in data -> random effect = 0 -> using re.form = NA
)

# Population‐average prediction (new county)
pred_log  <- predict(m0, newdata = new_data, re.form = NA) # E(d(j))=0
pred_rad  <- exp(pred_log)
pred_rad






### Estimated intercept
# Random effects with conditional variance
re  <- ranef(m0, condVar = TRUE) # condVar to get posterior variance

# random intercept deviations (b0_j)
ri  <- re$county
ri$county <- rownames(ri)
u_cty <- tapply(df$u, df$county, mean)   # same value of u in each county, mean = same level
ri$u <- u_cty[ri$county]


# intercept
beta0 <- fixef(m0)["(Intercept)"]
alpha<- fixef(m0)["u"]

# county-specific regression intercepts = mu + alpha * u + b0_j
ri$intercept <- beta0 + ri[,"(Intercept)"] + alpha * ri$u

# extract posterior variances for b0_j, store as standard error
post_var <- attr(ranef(m0, condVar = TRUE)$county, "postVar")[1,1,]
ri$SE <- sqrt(post_var)
ri$lower <- ri$intercept - qnorm(1-(0.05/2)) * ri$SE
ri$upper <- ri$intercept + qnorm(1-(0.05/2)) * ri$SE

# adding region for coloring 
ri$region <- df$region[ match(ri$county, df$county) ]
ri$region <- factor(ri$region, levels = region.order)


ggplot(ri, aes(x = u, y = intercept, color = region)) +
  
  # error bars first (behind points)
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.6) +
  
  # colored points
  geom_point(size = 2) +
  
  geom_abline(
    intercept = beta0, slope = alpha,
    color     = "black", alpha=0.4, linewidth = 0.2) +
  
  scale_color_manual(values = region_colors, name = "Region") +
  
  labs(
    title = "Estimated intercepts for counties, x=0",
    x = "County uranium (u)",
    y = "Estimated intercept (log radon)"
  ) +
  theme_bw()


## At ground level
# fixed basement effect
beta_x1 <- fixef(m0)["x1"]   # name may differ: check your model
alpha_x1 <- fixef(m0)["u:x1"] 

# county-specific mean at x = 1
ri$intercept_x1 <- (beta0 + beta_x1) + ri$"(Intercept)" + (alpha+alpha_x1)*ri$u #+ ri$b1

# CI for x=1 as well
ri$lower_x1 <- ri$intercept_x1 - qnorm(.975)*ri$SE
ri$upper_x1 <- ri$intercept_x1 + qnorm(.975)*ri$SE


ggplot(ri, aes(x = u, y = intercept_x1, color = region)) +
  
  # error bars first (behind points)
  geom_errorbar(aes(ymin = lower_x1, ymax = upper_x1), width = 0, alpha = 0.6) +
  
  # colored points
  geom_point(size = 2) +
  
  geom_abline(
    intercept = beta0+beta_x1, slope = alpha+alpha_x1,
    color     = "black", alpha=0.4, linewidth = 0.2) +
  
  scale_color_manual(values = region_colors, name = "Region") +
  
  labs(
    title = "Estimated intercepts for counties, x=1",
    x = "County uranium (u)",
    y = "Estimated intercept (log radon)"
  ) +
  theme_bw()





######### Prediction for a new house 
# Range of uranium levels

u_seq <- seq(min(df$u), max(df$u), length.out = 100)

new_house <- data.frame(
  u = u_seq,
  x = factor(0, levels = levels(df$x)),     
  county = "newCounty"  #
)
new_house_ground <- data.frame(
  u = u_seq,
  x = factor(1, levels = levels(df$x)),
  county = "newCounty"
)


# Make predictions E(d_j)=0
pred0 <- predict(m0, newdata = new_house, re.form = NA)
pred1 <- predict(m0, newdata = new_house_ground, re.form = NA)

new_house$pred  <- pred0
new_house_ground$pred <- pred1


ggplot(new_house, aes(x = u, y = pred)) +
  geom_line() +
  labs(x = "Uranium level", y = "Predicted y", 
       title = "Predicted exposure for new county")

boot_fun0 <- function(fit) {
  predict(fit, newdata = new_house, re.form = NA)
}

set.seed(2201)
b0 <- bootMer(m0, boot_fun0, nsim = 500)

CI0 <- apply(b0$t, 2, quantile, c(0.025, 0.975))

new_house$lwr <- CI0[1, ]
new_house$upr <- CI0[2, ]

boot_fun1 <- function(fit) {
  predict(fit, newdata = new_house_ground, re.form = NA)
}

set.seed(2201)
b1 <- bootMer(m0, boot_fun1, nsim = 500)

CI1 <- apply(b1$t, 2, quantile, c(0.025, 0.975))

new_house_ground$lwr <- CI1[1, ]
new_house_ground$upr <- CI1[2, ]



ggplot(new_house, aes(x = u, y = pred)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Uranium level", y = "Predicted y (CI)",
       title = "Population-level prediction with bootstrap CIs")

## Back transformed 
new_house$pred_bt <- exp(new_house$pred)
new_house$lwr_bt  <- exp(new_house$lwr)
new_house$upr_bt  <- exp(new_house$upr)

new_house_ground$pred_bt <- exp(new_house_ground$pred)
new_house_ground$lwr_bt  <- exp(new_house_ground$lwr)
new_house_ground$upr_bt  <- exp(new_house_ground$upr)

# Back-transform uranium predictor
new_house$u_orig        <- exp(new_house$u)
new_house_ground$u_orig <- exp(new_house_ground$u)

ggplot() +
  
  # --- Basement (x = 0) CI ribbon ---
  geom_ribbon(
    data = new_house,
    aes(x = u_orig, ymin = lwr_bt, ymax = upr_bt, fill = "0"),
    alpha = 0.25
  ) +
  
  # --- Basement (x = 0) line ---
  geom_line(
    data = new_house,
    aes(x = u_orig, y = pred_bt, color = "0"),
    size = 1
  ) +
  
  # --- Ground floor (x = 1) CI ribbon ---
  geom_ribbon(
    data = new_house_ground,
    aes(x = u_orig, ymin = lwr_bt, ymax = upr_bt, fill = "1"),
    alpha = 0.25
  ) +
  
  # --- Ground floor (x = 1) line ---
  geom_line(
    data = new_house_ground,
    aes(x = u_orig, y = pred_bt, color = "1"),
    size = 1
  ) +
  
  scale_color_manual(
    name = " ",
    values = floor_colors,
    labels = c("0" = "Basement", "1" = "Ground floor")
  ) +
  
  scale_fill_manual(
    name = "",
    values = floor_colors,
    labels = c("0" = "Basement CI", "1" = "Ground floor CI")
  ) +
  
  labs(
    x = "Uranium level (ppm)",
    y = "Predicted radon level (pCi/L)",
    title = "Predicted radon level for a new county"
  ) +
  
  theme_bw() +
  theme(
    plot.title  = element_text(size = 20, face = "bold"),
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    legend.position = "top"
  )













####### Region as fixed effect 

m1 <- lmer(y ~ u*x + u + x + region + ( 1 | county), data = df, REML=TRUE) 
summary(m1)
ranova(m1) # As expected county insignificant because it is explained in region

m2 <- lm(y ~ u*x + region , data = df) 

drop1(m2,test="F")

m3 <- lm(y ~ u + x + region , data = df) 

drop1(m3,test="F")


anova(m0,m3)

residplot(m0)
plot(m3,which=1:4)

