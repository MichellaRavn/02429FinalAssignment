# Libraries
library(ggplot2)
library(patchwork)
library(lmerTest)
library(MASS)
library(lattice)
library(predictmeans) 
library(plyr)
library(dplyr)
library(gridExtra)
library(emmeans)
library(lattice)
library(performance)


### Cleaning data
# Load the data
rats <- read.table("Data/Rats_Box.txt", header = TRUE, sep = ",")
treat_names <- c("Control", "Thyroxin", "Thioracil")

# Finding weight parameter, cumulative weight gain
rats <- within(rats, {
  y1 <- y0 + y1        
  y2 <- y1 + y2       
  y3 <- y2 + y3          
  y4 <- y3 + y4        
})

# Reshaping to long format
rats_long <- reshape(
  rats,
  varying = list(c("y1","y2","y3","y4")),
  v.names = "y",
  timevar = "week",
  times = 1:4,
  direction = "long"
)

rats_long <- rats_long[order(rats_long$Rat, rats_long$week), ]
row.names(rats_long) <- NULL

# Now storing the correct format
rats <- rats_long

# make treatment and cage factors 
rats$Trt <- as.factor(rats$Trt)
rats$Rat <- as.factor(rats$Rat)
rats$y <- as.numeric(rats$y)

# Data structure
str(rats) 

# make two versions of the time variable 
# - one quantitative and one qualitative
rats$weekQ <- as.numeric(rats$week)
rats$weekF <- as.factor(rats$week)

# Baseline weight for plotting
baseline <- rats |>
  distinct(Rat, Trt, y0) |>
  mutate(
    weekQ = 0,
    weekF = factor(0),
    y     = y0
  )

# for plotting
rats_plot <- bind_rows(
  rats[, c("Rat", "Trt", "weekQ", "weekF", "y")],
  baseline[, c("Rat", "Trt", "weekQ", "weekF", "y")]
)

# Individual profiles
p1 <- ggplot(rats_plot, aes(x = weekQ, y = y, group = Rat, colour = Trt)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 0.7) +
  labs(
    title = "Weight over time",
    x = "Week",
    y = "Weight (gram)"
  ) +
  theme_bw()

# Mean profiles 
mns <- rats_plot |>
  group_by(Trt, weekQ) |>
  summarise(y = mean(y), .groups = "drop")

p2 <- ggplot(mns, aes(x = weekQ, y = y, group = Trt, colour = Trt)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean weight per treatment",
    x = "Week",
    y = "Mean weight (gram)"
  ) +
  theme_bw()

p1 | p2


###### Log profiles
#baseline <- rats |>
#  distinct(Rat, Trt, y0) |>
#  mutate(
#    weekQ = 0,
#    weekF = factor(0),
#    logy     = log(y0)
#  )

# for plotting
#rats_plot <- bind_rows(
#  rats[, c("Rat", "Trt", "weekQ", "weekF", "logy")],
#  baseline[, c("Rat", "Trt", "weekQ", "weekF", "logy")]
#)

# Individual profiles
#p1 <- ggplot(rats_plot, aes(x = weekQ, y = logy, group = Rat, colour = Trt)) +
#  geom_line(alpha = 0.8) +
#  labs(
#    title = "Log weight over time",
#    x = "Week",
#    y = "Log weight"
#  ) +
#  theme_bw()
#p1

# Mean profiles 
#mns <- rats_plot |>
#  group_by(Trt, weekQ) |>
#  summarise(y = mean(logy), .groups = "drop")

#p2 <- ggplot(mns, aes(x = weekQ, y = y, group = Trt, colour = Trt)) +
#  geom_line(size = 1) +
#  geom_point(size = 2) +
#  labs(
#    title = "Mean log weight per treatment",
#    x = "Week",
#    y = "Mean log weight"
#  ) +
#  theme_bw()

#p1 / p2


# I. Fit a random effects model using function lmer() from lme4/lme4Test
#    Covariance with structure corresponding to Compound symmetry.

m1 <- lmer(y ~ weekQ + Trt + weekQ:Trt + y0 + (1 | Rat), data = rats)
m2 <- lmer(y ~ weekF + Trt + weekF:Trt + y0 + (1 | Rat), data = rats)
m3 <- lmer(y ~ weekQ + Trt + weekQ:Trt + y0 + (1 + weekQ | Rat), data = rats)
m4 <- lmer(y ~ I(weekQ^2) + weekQ + Trt + weekQ:Trt + y0 + (1 + weekQ | Rat), data = rats)
m5 <- lmer(y ~ I(weekQ^3) + I(weekQ^2) + weekQ + Trt + weekQ:Trt + y0 + (1 + weekQ | Rat), data = rats)

# Testing random effect
ranova(m1)
ranova(m2)
ranova(m3) 
ranova(m4) 
ranova(m5) 

# Residual plots
residplot(m1) # looks like adding curvature is needed
residplot(m2)
residplot(m3) # funnel shaped
residplot(m4) # almost identical to m3
residplot(m5) # slightly more flat res vs fit

# Comparing nested models
anova(m1,m2) # prefers m1, time as a numeric
anova(m1,m3) # prefers m3, random slope
anova(m3,m4) # quadratic is insignificant
anova(m3,m5) # cubic is significant

# Continuing with m3
m3 <- update(m3,REML=F)
drop1(m3) # Interaction significant
drop1(m5) # Interactions significant

# AIC, BIC - prefer m3, simple model with random slope
AIC(m1,m2,m3,m4)
BIC(m1,m2,m3,m4)

# Linear models for boxcox
m3_lm <- lm(y ~ weekQ + Trt + weekQ:Trt + y0, data = rats)
m5_lm <- lm(y ~ I(weekQ^3) + I(weekQ^2) + weekQ + Trt + weekQ:Trt + y0, data = rats)

# Checking for transformation (no indication, but making sure)
aux <- boxcox(m5_lm, lambda = seq(-1, 2, by = 0.05))
(lambda<-aux$x[which.max(aux$y)]) # 0 is in the CI, 1 is not = log transform

# Log transforming
rats$logy <- log(rats$y)
rats$logy0 <- log(rats$y0)

# Fitting new models
m6 <- lmer(logy ~ weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ | Rat), data = rats)
# Residuals are curved, consider squaring week
m7 <- lmer(logy ~ I(weekQ^2) + weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ | Rat), data = rats)
m8 <- lmer(logy ~ I(weekQ^3)+I(weekQ^2) + weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ + weekQ^2 + weekQ^3 | Rat), data = rats)

# Residual plots
residplot(m6)
residplot(m7)
residplot(m8)

anova(m4,m5) # prefer m5
anova(m5,m6) # prefer m6

AIC(m1,m2,m3,m4,m5,m6,m7,m8)
BIC(m1,m2,m3,m4,m5,m6,m7,m8) # m8 has lowest, but already m6 and m7 are much better

# Variance explained (fixed and mixed)
r2_nakagawa(m1)
r2_nakagawa(m2)
r2_nakagawa(m3)
r2_nakagawa(m4)
r2_nakagawa(m5)
r2_nakagawa(m6)
r2_nakagawa(m7)
r2_nakagawa(m8)

# Acceptable results from m3 in residplot, thus recommending for final model
mfinal <- m3


#### Checking variogram structures
## Correlation structures
cor_structures <- list(
  COMP   = corCompSymm(form = ~ weekQ | Rat),
  GAUS   = corGaus(form = ~  weekQ | Rat),
  EXP    = corExp(form = ~ weekQ | Rat),
  AR1    = corAR1(form = ~  weekQ | Rat), # Serial, order 1
  CAR1   = corCAR1(form = ~ weekQ | Rat), # Continous AR1
  LIN    = corLin(form = ~ weekQ | Rat),
  RATIO  = corRatio(form = ~ weekQ | Rat),
  SPHER  = corSpher(form = ~ weekQ | Rat)
)

# Fit one model for each structure
M <- list()

for (nm in names(cor_structures)) {
  M[[nm]] <- try(lme(
    logy ~ weekQ + Trt + weekQ:Trt + y0 ,
    random = ~ 1  | Rat, # random slope does not converge
    correlation = cor_structures[[nm]],
    data = rats,
    control = lmeControl(msMaxIter = 5000, niterEM = 50))
  )
}
# Check names
names(M)

V <- list()
for (nm in names(M)) {
  
  var_obj <- Variogram(M[[nm]], form = ~ weekQ | Rat, data = rats)
  
  V[[nm]] <- plot(
    var_obj,
    main = paste("Variogram -", nm),
    pch = 16,
    cex = 1,
    ylim = c(0, 2),
    smooth = FALSE,
    xlab = "Distance (week)"
  )
}

# Display all in a grid
grid.arrange(
  grobs = V,
  ncol = 4
)

# With only three datapoints it is hard to tell

# This breaks, does not converge:
#M8<-lme(logy ~ weekQ + Trt + weekQ:Trt + y0 , random = ~ 1 + weekQ  | Rat, 
#        correlation = corExp(form = ~ weekQ | Rat), data = rats,
#        control = lmeControl(msMaxIter = 500000, niterEM = 50))




# --- Compare structures ------------------------------------------------------
anova(M[[1]],M[[2]],M[[3]],M[[6]],M[[7]])
# Not really better AIC or BIC
# Even though AR1 or exp might be theoretically most correct because of time dist dependence


m9<-M[[2]]
      
anova(m3,m9)

residplot(M[[1]])
residplot(M[[2]])
residplot(M[[3]])
residplot(M[[6]])
residplot(M[[7]])


# comparing model fit
#m3=m5=m8?
par(mfrow=c(1,3))
plot(predict(m3, level=1),predict(m5, level=1), pch=20,
     main = "m3 random slope vs m5 poly3")
abline(0,1, col="red")

plot(predict(m3, level=1),exp(predict(m7, level=1)), pch=20,
     main = "m3 random slope vs m7 log poly2")
abline(0,1, col="red")

plot(predict(m3, level=1),exp(predict(m8, level=1)), pch=20,
     main = "m3 random slope vs m8 log poly3")
abline(0,1, col="red")
# m3 seems just as fine as the complicated m8


# Comparing fitted vs observed
par(mfrow=c(2,4), mar=c(3,3,2,1), oma=c(0,0,3,0))  

plot(predict(m1, level=1), rats$y, pch=20,
     main="m1 ")
abline(0,1, col="red")

plot(predict(m2, level=1), rats$y, pch=20,
     main="m2")
abline(0,1, col="red")

plot(predict(m3, level=1), rats$y, pch=20,
     main="m3 ")
abline(0,1, col="red")

plot(predict(m4, level=1), rats$y, pch=20,
     main="m4")
abline(0,1, col="red")

plot(predict(m5, level=1), rats$y, pch=20,
     main="m5")
abline(0,1, col="red")

plot(exp(predict(m6, level=1)), rats$y, pch=20,
     main="m6")
abline(0,1, col="red")

plot(exp(predict(m7, level=1)), rats$y, pch=20,
     main="m7")
abline(0,1, col="red")

plot(exp(predict(m8, level=1)), rats$y, pch=20,
     main="m8")
abline(0,1, col="red")

# Add outer title
mtext("Fitted vs Observed for all 8 models", outer=TRUE, cex=1.5)




### BLUPS 
dotplot(ranef(m3, condVar = TRUE),
        strip = FALSE,
        ylab = "",
        scales = list(cex = 0.8))   # smaller labels



### Post hoc
round(fixef(m3),digits=2)
ci_m3 <- confint(m3, level = 0.95)
round(ci_m3, 2)

# Predictions
# Mixed predictions
rats$m3 <- predict(m3)
# Fixed predictions
rats$m3_fix <- predict(m3,re.form=NA)

# Mixed predictions
rats$m7 <- exp(predict(m7))
# Fixed predictions
rats$m7_fix <- exp(predict(m7,re.form=NA))

newdat <- expand.grid(
  weekQ = 0:4,
  Trt = levels(rats$Trt),
  y0 = mean(rats$y0),
  logy0 = mean(rats$logy0)
)

newdat$pred_m3 <- predict(m3, newdata = newdat, re.form = NA)
newdat$pred_m7 <- exp(predict(m7, newdata = newdat, re.form = NA))

p2 +
  geom_line(
    data = newdat,
    aes(x = weekQ, y = pred_m3, group = Trt, linetype = Trt),
    linewidth = 0.8, alpha = 0.8
  )



### Emmeans
wvals <- c(1,2.5,4)

# emmeans table
em <- emmeans(m3, "Trt", by = "weekQ", at = list(weekQ=wvals))
em

# For log
em7 <- emmeans(m7, "Trt", by = "weekQ", at = list(weekQ=wvals))
em7

# Back-transform to radon scale
em_bt <- transform(
  as.data.frame(em7),
  em_bt      = round(exp(emmean),2),
  lower_bt   = round(exp(lower.CL),2),
  upper_bt   = round(exp(upper.CL),2)
)

em_bt

emmip(m3,Trt ~ weekQ, at = list(weekQ = wvals) , CIs = TRUE) + theme(legend.position="top")+
  #ggtitle("Emmeans interaction")+
  ylab("Weight") + 
  labs(color="Trt") 

tr <- emtrends(m3, pairwise ~Trt, var="weekQ",infer=TRUE, adjust="Tukey")  
# plot 
library(multcomp)
L <- contrast(tr, "pairwise")@linfct
rownames(L) <- c("1 - 2", "1 - 3", "2 - 3")
mult_slope <- glht(m3, linfct = L)
summary(mult_slope)
confint(mult_slope)
plot(mult_slope, col = 2:7)


slopes <- emtrends(m3, ~ Trt, var = "weekQ")  
slopes_df <- as.data.frame(slopes)
ggplot(slopes_df, aes(x = Trt, y = weekQ.trend)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  labs(
    title = "Estimated slope",
    x = "Trt",
    y = "Weight slope (dy/dweek)"
  ) +
  theme_bw()





# Cooks distance
cd <- cooks.distance(m3)
plot(cd, type = "h", main = "Cook's distance", xlab = "Observation", ylab = "Distance")
rats[cd>0.6,]


# Time dependent ICC
vc <- VarCorr(m3)

# Variances
sigma0_sq <- vc$Rat["(Intercept)", "(Intercept)"]    # random intercept variance
sigma1_sq <- vc$Rat["weekQ", "weekQ"]                # random slope variance

# Covariance between intercept and slope
sigma01 <- vc$Rat["(Intercept)", "weekQ"]

# Residual variance
sigma_e_sq <- attr(vc, "sc")^2

ICC <- function(t) {
  # between-rat variance at time t
  var_between <- sigma0_sq + 2*t*sigma01 + t^2*sigma1_sq
  
  # ICC(t)
  var_between / (var_between + sigma_e_sq)
}

ICC(1.0)
ICC(2.5)
ICC(4.0)







### Individual slopes
rats_plot2 <- rats_plot %>%
  arrange(Rat, weekQ) %>%        # ensure time is sorted within rats
  group_by(Rat) %>%
  mutate(y0_rat = first(y)) %>%  # baseline = first observed weight
  ungroup()

# Plot function with color = baseline weight
plot_trt <- function(trt_value) {
  ggplot(
    rats_plot2 %>% filter(Trt == trt_value),
    aes(x = weekQ, y = y, group = Rat, colour = y0_rat)
  ) +
    geom_line(alpha = 0.9, linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_colour_viridis_c(
      option = "C",
      begin = 0.1,   # light
      end   = 0.9    # dark
    ) +
    labs(
      title = paste("Treatment", trt_value),
      x = "Week",
      y = "Weight (gram)",
      colour = "Initial\nWeight"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )
}

p1 <- plot_trt(1)
p2 <- plot_trt(2)
p3 <- plot_trt(3)

grid.arrange(p1, p2, p3, ncol = 3)











