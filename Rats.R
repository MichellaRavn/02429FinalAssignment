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
library(dplyr)

#####################################################
# Cleaning data
#####################################################

# Load the data
rats <- read.table("Data/Rats_Box.txt", header = TRUE, sep = ",")

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

# Data structure
str(rats) 

# make two versions of the time variable 
# - one quantitative and one qualitative

rats$weekQ <- as.numeric(rats$week)
rats$weekF <- as.factor(rats$week)

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
  geom_line(alpha = 0.5) +
  labs(
    title = "Weight per rat over time",
    x = "Week",
    y = "Weight"
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
    x = "Week (0 = baseline y0)",
    y = "Mean weight"
  ) +
  theme_bw()

p1 / p2


# I. Fit a random effects model using function lmer() from lme4/lme4Test
#    Covariance with structure corresponding to Compound symmetry.

m1 <- lmer(y ~ weekQ + Trt + weekQ:Trt + y0 + (1 | Rat), data = rats)
m2 <- lmer(y ~ weekF + Trt + weekF:Trt + y0 + (1 | Rat), data = rats)
m3 <- lmer(y ~ weekQ + Trt + weekQ:Trt + y0 + (1 + weekQ | Rat), data = rats)

ranova(m1)
ranova(m2)
ranova(m3) # random slope significant

anova(m1,m2) # prefers m1, time as a numeric
anova(m1,m3) # prefers m3, random slope

# Checking for transformation (no indication, but making sure)
m3_lm <- lm(y ~ weekQ + Trt + weekQ:Trt + y0, data = rats)
aux <- boxcox(m3_lm, lambda = seq(-1, 2, by = 0.05))
(lambda<-aux$x[which.max(aux$y)]) # 0 is in the CI, 1 is not = log transform

# Log transforming

rats$logy <- log(rats$y)
rats$logy0 <- log(rats$y0)

m4 <- lmer(logy ~ weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ | Rat), data = rats)

# AIC, BIC - prefer m4
AIC(m1,m2,m3,m4)
BIC(m1,m2,m3,m4)

# Residuals are curved, consider squaring week
m5 <- lmer(logy ~ I(weekQ^2) + weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ | Rat), data = rats)
m6 <- lmer(logy ~ I(weekQ^3)+I(weekQ^2) + weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ + weekQ^2 + weekQ^3 | Rat), data = rats)

residplot(m4)
residplot(m5)
residplot(m6)

anova(m4,m5) # prefer m5
anova(m5,m6) # prefer m6, but suggest m5 for simplicity

AIC(m4,m5,m6)
BIC(m4,m5,m6) # m6 has lowest, but not much more than m5

m7 <- lmer(logy ~ I(weekQ^2) + weekQ + Trt + weekQ:Trt + logy0 + (1 + weekQ + weekQ^2 | Rat), data = rats)


m_final <- m5
# Variance
# Extract variance components
VC     <- VarCorr(m_final)
Sigma  <- as.matrix(VC$Rat)        # 2×2 matrix: intercept & slope
sigres <- VC$Residual[1]           # residual variance

# ICC function for random-intercept + random-slope model
icc_at <- function(x) {
  rand_var <- Sigma[1,1] + 2*Sigma[1,2]*x + Sigma[2,2]*x^2
  rand_var / (rand_var + sigres)
}

# ICC at weeks 1–4
ICC <- sapply(1:4, icc_at)
ICC


residplot(m1)
residplot(m2)
residplot(m3)
residplot(m4)

## Checking correlation structures

# Define correlation structures
cor_structures <- list(
  COMP   = corCompSymm(form = ~ weekF | Rats),      # OK: constant correlation
  GAUS   = corGaus(form = ~ weekF | Rats),     # OK: continuous spatial
  EXP    = corExp(form = ~ weekF | Rats),      # OK: continuous spatial
  AR1    = corAR1(form = ~ weekF | Rats),      # OK: discrete AR(1), order only
  CAR1   = corCAR1(form = ~ weekF | Rats),     # OK: continuous-time AR(1)
  LIN    = corLin(form = ~ weekF | Rats),       # OK: continuous spatial
  # ARMA   = corARMA(form = ~ weekF | Rats, p=1,q=1),  # FAIL: requires equally spaced time
  #SYMM   = corSymm(form = ~ weekF | Rats)          # FAIL: covariate must be 1,2,...,T
  RATIO  = corRatio(form = ~ weekF | Rats),          # FAIL: zero distances not allowed
  SPHER  = corSpher(form = ~ weekF | Rats)           # FAIL: zero distances not allowed
)



# Fit one model for each structure
M <- list()

for (nm in names(cor_structures)) {
  M[[nm]] <- lme(
    logph ~ logmin * treatm,
    random = ~1 | pigno,
    correlation = cor_structures[[nm]],
    data = dfph,
    control = lmeControl(msMaxIter = 200, niterEM = 50)
  )
}

# Check names
names(M)















baseline <- rats |>
  distinct(Rat, Trt, y0) |>
  mutate(
    weekQ = 0,
    weekF = factor(0),
    logy     = log(y0)
  )

# for plotting
rats_plot <- bind_rows(
  rats[, c("Rat", "Trt", "weekQ", "weekF", "logy")],
  baseline[, c("Rat", "Trt", "weekQ", "weekF", "logy")]
)


# Individual profiles
p1 <- ggplot(rats_plot, aes(x = weekQ, y = logy, group = Rat, colour = Trt)) +
  geom_line(alpha = 0.5) +
  labs(
    title = "Log weight per rat over time",
    x = "Week",
    y = "Log weight"
  ) +
  theme_bw()
p1

# Mean profiles 
mns <- rats_plot |>
  group_by(Trt, weekQ) |>
  summarise(y = mean(logy), .groups = "drop")

p2 <- ggplot(mns, aes(x = weekQ, y = y, group = Trt, colour = Trt)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean log weight per treatment",
    x = "Week",
    y = "Mean log weight"
  ) +
  theme_bw()

p1 / p2

