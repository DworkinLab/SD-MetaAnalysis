############ NT Bio4C12 SSD and SShD Meta-analysis #######################
# Libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)
library(tidyr)
library(HardyWeinberg)
library(emmeans)
library(bayesplot)
library(reshape2)
library(RColorBrewer)


# Setting working directory
setwd("~/OneDrive - McMaster University/SSD and SShD Meta-analysis/Data")

# Reading in data-set
metadata <- read.csv("metadata_4c12.csv")

# Changing the phyla for Bonduriansky
metadata$Phyla <- with(metadata, ifelse(Study == "Bonduriansky", "Arthropoda", Phyla))

# Strata Type for Houle et al.
metadata$Strata_Type <- with(metadata, ifelse(is.na(Strata_Type), "None", Strata_Type))

# Make male-biased positive
metadata$SSD <-metadata$SSD*(-1)

# Add new column with info about female or male bias
metadata$Biased <- with(metadata, ifelse(SSD*(-1) < 0, "Male", "Female"))

# Convert to factor
metadata[,c(3, 4, 5, 6, 7, 8, 9, 13, 12, 25)] <- lapply(metadata[,c(3, 4, 5, 6, 7, 8, 9, 13, 12, 25)], factor)

# Taking the abolute value of SSD
metadata$SSD <- abs(metadata$SSD)

# Sample sizes
with(metadata, table(Phyla, Biased))

# Remove one Chordata female-based species for now (till more data-sets are added)
metadata <- metadata %>%
  filter((Phyla == "Chordata" & Biased == "Male")| Phyla == "Arthropoda")

# Create interaction term for phyla and biased
metadata$phyla_biased <- with(metadata, interaction(Phyla, Biased))

##################### Analyzing the bootstrap samples ######################
# Calculate sample sizes (back from Fisher's standard errors)
metadata$Sample.Sizes <- with(metadata, 1/(Fisher_SE)^2 - 3)

# Bootstrap standard errors for SSD

png(filename = "SSDBoot.png", width = 20, height = 15, units = "cm", res = 300)

ggplot(metadata, aes(x = Sample.Sizes, y = SSD.SE)) + geom_point() + labs(x = "Sample Sizes", y = "SSD Bootstrapped Standard Error") + theme_bw() + xlim(c(0, 250))

dev.off()

png(filename = "SShDBoot.png", width = 20, height = 15, units = "cm", res = 300)

ggplot(metadata, aes(x = Sample.Sizes, y = SShD.SE)) + geom_point() + labs(x = "Sample Sizes", y = "SShD Bootstrapped Standard Error") + theme_bw() + xlim(c(0, 250))

dev.off()

ggplot(metadata, aes(x = Sample.Sizes, y = Mag_Diff_SE)) + geom_point() + labs(x = "Sample Sizes", y = "SShD Bootstrapped Standard Error") + theme_bw() + xlim(c(0, 1000))


#############################################################################
# Some exploratory plots

# SSD vs Phyla
ggplot(metadata, aes(x = Biased, y = SSD)) + geom_boxplot(aes(colour= Phyla)) + labs(x = "Biased", y = "SSD") + theme_bw()

# SSD in male and female biased organisms (Sexually Selected as Strata)
ggplot(metadata, aes(x = Biased, y = SSD)) + geom_boxplot(aes( colour = Sexually_Selected, fill = Phyla)) + labs(x = "Bias", y = "SSD") + theme_bw()

# SSD in different phyla by bias
ggplot(metadata, aes(x = Phyla, y = SSD)) + geom_boxplot(aes(colour= Biased)) + labs(x = "Phyla", y = "SSD") + theme_bw()

ggplot(metadata, aes(x = Sexually_Selected, y = SSD)) + geom_boxplot(aes(colour= interaction(Phyla, Biased))) + labs(x = "Strata Type", y = "SSD") + theme_bw()


# SSD vs SShD
ggplot(metadata, aes(x = Biased, y = SShD)) + geom_boxplot(aes(colour = Phyla)) + labs(x = "Bias", y = "SShD") + theme_bw()

ggplot(metadata, aes(x = Biased, y = SShD)) + geom_point(aes(colour= Sexually_selected)) + labs(x = "Biased", y = "SShD") + theme_bw()


# SShD vs Strata-Type
ggplot(metadata, aes(x = Sexually_Selected, y = SShD)) + geom_boxplot(aes(colour= interaction(Phyla, Biased))) + labs(x = "Strata Type", y = "SShD") + theme_bw()

ggplot(metadata, aes(x = SSD, y = SShD)) + geom_boxplot(aes(colour= interaction(Phyla, Biased))) + labs(x = "Strata Type", y = "SShD") + theme_bw()

# SShD sexually selected
ggplot(metadata, aes(x = Sexually_Selected, y = SShD)) + geom_boxplot(aes(colour= Phyla)) + labs(x = "Sexually Selected", y = "SShD") + theme_bw()

# Vector correlations by strata
# Phyla

png(filename = "CorrOB.png", width = 20, height = 15, units = "cm", res = 300)

ggplot(metadata, aes(x = Corr)) + geom_density(aes(fill = Phyla), alpha = 0.5) + 
  labs(title = "Observed Distribution of Vector Correlations")

dev.off()

png(filename = "CorrOB.png", width = 20, height = 15, units = "cm", res = 300)

ggplot(metadata, aes(x = Corr)) + geom_histogram(aes(fill = Phyla),  alpha = 0.5, 
                                                 binwidth = 0.1)


dev.off()

ggplot(metadata, aes(x = SShD)) + geom_density() 

# Phyla biased
ggplot(metadata, aes(x = Corr)) + geom_density(aes(fill = phyla_biased)) + labs(x = "Strata Type", y = "Corr") 

ggplot(metadata, aes(x = Corr)) + geom_density() 

ggplot(metadata, aes(x = FisherZ)) + geom_density(aes(fill = Strata_Type)) + labs(x = "Strata Type", y = "Corr") 

ggplot(metadata, aes(x = FisherZ)) + geom_density() 

ggplot(metadata, aes(x = Sexually_Selected, y = Mag_Diff)) + geom_boxplot(aes(colour= Phyla)) + labs(x = "Strata Type", y = "Corr") + theme_bw()

ggplot(metadata, aes(x = Biased, y = Corr)) + geom_boxplot(aes(colour= Biased)) + labs(x = "Biased", y = "Corr") + theme_bw()

ggplot(metadata, aes(x = Corr)) + geom_density()

# Log SSD vs bootstrapped standard error
with(metadata, plot(SSD.SE, SSD))

# Extracting the approx sample sizes from Fisher SE
metadata$Sample.Sizes <- (1/(metadata$Fisher_SE)^2) + 3

# Sample sizes vs SSD.SE
with(metadata, plot(SSD.SE, Sample.Sizes))

# Sample sizes vs SShD.SE
with(metadata, plot(SShD.SE, Sample.Sizes))

ggplot(metadata, aes(x = Corr)) + geom_density()

metadata$SSD <- abs(metadata$SSD)


# Some basic models


##################### Model for SSD #################################
# Setting priors
priors <- c(prior(normal(0, 2), class = Intercept),
            prior(normal(0, 5), class = "b"))

SSD_fit <- brm(data = metadata, family = student(),
                formula = SSD|trunc(lb = 0)| se(SSD.SE, sigma = TRUE) ~ Phyla + Sexually_Selected + (1|Study),  seed = 801,  iter = 20000, warmup = 1000, inits = "random",  chains =  1,
               prior = priors)

# Posterior predictive check
png(filename = "pp_checkSSD.png", width = 20, height = 15, units = "cm", res = 300)

pp_check(SSD_fit, nsamples = 100) + xlim(c(-0.3, 0.3))

dev.off()

#creating conditional effects object
c_eff <- conditional_effects(SSD_fit)


# Plot of phyla and Sexually_Selected effects
dat_1 <- conditional_effects(SSD_fit)[[1]]




png(filename = "SShDss.png", width = 20, height = 15, units = "cm", res = 300)
p1 <-
  conditional_effects(SSD_fit,
                      effects = "Phyla", points = TRUE)

plot(p1,
     plot = F)[[1]] +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(y = "SShD", colour = "Sexually Selected", x = "Phyla and Direction of Bias") + theme_bw() 

dev.off()

conditions <- data.frame(Biased = c("Female", "Male"))

# Grouped posterior distributions
y <- metadata$SSD

yrep1 <- posterior_predict(SSD_fit, nsamples = 1000)

png(filename = "SSD_Biased.png", width = 20, height = 15, units = "cm", res = 300)

color_scheme_set("purple")
ppc_stat_grouped(metadata$SSD, yrep1, stat = "median", group = metadata$Biased)
dev.off()


conditional_effects(SSD_fit, effects = "Phyla:Sexually_Selected")


################### Model for SShD  #################################
# Setting priors
priors <- c(prior(normal(0, 10), class = Intercept),
            prior(normal(0,5), class = "b"))


# Student's t
# Recovers kurtosis beautifully!!
SShD_fit2 <- brm(data = metadata, family = student(),
                 formula = SShD|trunc(lb = 0) ~ SSD + (1|Study),seed = 800,  iter = 15000, warmup = 1000, inits = "random",  chains =  1)



summary(SShD_fit2)

# Posterior predictive check
png(filename = "pp_checkSShD.png", width = 20, height = 15, units = "cm", res = 300)

pp_check(SShD_fit2, nsamples = 10) + xlim(c(0,0.1))

dev.off()

# Grouped posterior predictive plots
y <- metadata$SShD

yrep1 <- posterior_predict(SShD_fit2, nsamples = 1000)

png(filename = "SShD_Phyla1.png", width = 20, height = 15, units = "cm", res = 300)

color_scheme_set("blue")
ppc_stat_grouped(metadata$SShD, yrep1, stat = 'median', group = metadata$phyla_biased) +
  legend_none()

dev.off()

ggsave(filename = "plots/ppc_skew1.png", width = 4.5, height = 3.75)

# Pareto k values
SShD_loo <- loo(SShD_fit1, save_psis = TRUE, cores= 2)

# No bad pareto k values
plot(SShD_loo)

conditions <- data.frame(phyla_biased = c("Arthropoda.Female","Arthropoda.Male","Chordata.Male"))

conditions <- data.frame(Phyla= c("Chordata", "Arthropoda"))

# Figures
png(filename = "SShDss.png", width = 20, height = 15, units = "cm", res = 300)
p1 <-
  conditional_effects(SShD_fit2,
                      effects = "phyla_biased:Sexually_Selected")

plot(p1,
     plot = F)[[1]] +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(y = "Shape Dimorphism", x = "Size Dimorphism")

dev.off()


png(filename = "SShDfit.png", width = 20, height = 15, units = "cm", res = 300)
p1 <-
  conditional_effects(SShD_fit2,
                      effects = "SSD:phyla_biased", spaghetti = TRUE, nsamples = 200)

plot(p1, points = TRUE, point_args = c(alpha = 1/2, size = 1),
     plot = F)[[1]] +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(y = "SShD", colour = "Phyla and Direction of Bias")

dev.off()

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

png(filename = "SSDvsSShDStudy.png", width = 20, height = 15, units = "cm", res = 300)

metadata %>%
  group_by(phyla_biased)  %>%
  add_fitted_draws(SShD_fit2, n = 100) %>%
  ggplot(aes(x = SSD, y = SShD, color = ordered(phyla_biased))) +
  geom_line(aes(y = .value, group = paste(phyla_biased, .draw)), alpha = .1) +
  geom_point(data = metadata, aes(alpha = 0.0000000001, shape = ".")) +
  scale_colour_manual(values = mycolors) + theme_bw()

dev.off()

################## Model for correlation ############################
corr_fit <- brm(data = metadata, family = skew_normal(),
                bf(FisherZ|se(Fisher_SE, sigma = TRUE) ~  phyla_biased + Sexually_Selected + (1|Study)), seed =800, iter = 20000, warmup = 1000, inits = "random", chains = 1)


corr_fit <- brm(data = metadata, family = Beta(),
                Corr ~  phyla_biased + Sexually_Selected + (1|Study), seed =800, iter = 20000, warmup = 1000, inits = "random", chains = 1)

plot(corr_fit)

summary(corr_fit)

png(filen)

y <- metadata$FisherZ

yrep <- posterior_predict(corr_fit, draws = 1000)

yrep <- ifisherz(yrep)

# Model fit checks
png(filename = "pp_checkCorrfit.png", width = 20, height = 15, units = "cm", res = 300)

pp_check(corr_fit, nsamples = 100)

dev.off()

png(filename = "DistCorrfit3.png", width = 20, height = 15, units = "cm", res = 300)

color_scheme_set("viridis")
ppc_stat(metadata$Corr, yrep, stat = "median")

dev.off()

pp <- posterior_predict(SShD_fit2)

pp <- transpose(data.frame(pp)

hist(pp$X3)


ppc_intervals(metadata$Corr, yrep)

# Pareto k values
corr_loo <- loo(corr_fit, save_psis = TRUE, cores= 2)

# No bad pareto k values
plot(corr_loo)

data <- ppc_data(y, yrep)

# Some conditional effects plots
posterior <- as.matrix(corr_fit)

mcmc_areas(posterior,
           pars = c("b_phyla_biasedArthropoda.Male", "b_phyla_biasedChordata.Male"),
           prob = 0.8) 

plot(conditional_effects(corr_fit, effects = "Sexually_Selected:Strata_Type"))

draws <- metadata %>%
  group_b(Study) %>%
  add_fitted_draws(SShD_fit1, dpar = TRUE) %>%
  ggplot(aes(y = condition)) +
  stat_dist_slab(aes(dist = "norm", arg1 = mu, arg2 = sigma),
                 slab_color = "gray65", alpha = 1/10, fill = NA
  ) +
  geom_point(aes(x = response), data = ABC, shape = 21, fill = "#9ECAE1", size = 2)

save.image("thesis.RData")

# Contrasts
tt <- corr_fit %>%
  emmeans( ~ Strata_Type) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  median_qi()

png(filename = "ContSShDfit1.png", width = 20, height = 15, units = "cm", res = 300)

SShD_fit2 %>%
  emmeans( ~ Sexually_Selected) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x = .value, y = contrast)) +
  stat_halfeye(colour = "black") + theme_bw() +
  labs(x = "Values", y = "Posterior Distribution of Contrasts")

dev.off()

# Attempt to plot the posterior distributions
pp <- melt(yrep)

pp$index <- rep(seq(1, 546, 1), 19000)

pp$Phyla <- rep(metadata$Phyla, 19000)

# Stratified samples
ss <- pp %>%
  group_by(Phyla) %>%
  sample_n(1000)

ss <- data.frame(ss)

ggplot(ss, aes(x = value)) + geom_histogram(aes(colour = Phyla, alpha = 0.1))
