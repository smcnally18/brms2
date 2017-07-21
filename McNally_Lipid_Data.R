library(brms)
library(dplyr)
library(ggplot2)
library(readr)
library(car)
library(broom)
library(lsmeans)
library(visreg)
library(MASS)
library(forcats)
library(dplyr)
library(gplots)
library(sjPlot)
library(arm)
library(visreg)
library(DHARMa)
library(maps)
library(mapdata)

#load data
lipids_brms <- read.csv("UPDATED_Normalized_Energetic_Lipid_Content.csv") %>%
  na.omit(lipids_brms)

#what is my data
hist(lipids_brms$Lipid)

#log normal distrubtion, lipids are bound at 0

#set my priors
prior <- c(set_prior("normal(0,10)", class = "b"), 
           set_prior("student_t(10,0,1)", class = "sigma"),
           set_prior("student_t(10,0,1)", class = "sd"))

######Full Modelw/intercept Area and Colony Effect
lipid_brms_1 <- brm(Lipid  ~ 1 + Date + Bleaching.Level + Upwelling + 
                      (1 + Date + Bleaching.Level + Upwelling | Area) +
                      (1 + Date + Bleaching.Level + Upwelling | Colony), data=lipids_brms, family = lognormal, iter=4000, chains=4) 

#model check
pp_check(lipid_brms)

#plot 
plot(marginal_effects(lipid_brms, probs = c(0.05, 0.95)))

#####Model with only Colony Effect
lipid_brms_2 <- brm(Lipid  ~ 1 + Date + Bleaching.Level + Upwelling +
                    (1 | Area) +
                    (1 + Date + Bleaching.Level + Upwelling | Colony), data=lipids_brms, family = lognormal, prior = prior, sample_prior = TRUE, iter=4000, chains=4) 

pp_check(lipid_brms_2)
plot(marginal_effects(lipid_brms_2, probs = c(0.05, 0.95)))


#####Model with only Area Effect
lipid_brms_3 <- brm(Lipid  ~ 1 + Date + Bleaching.Level + Upwelling +
                      (1 + Date + Bleaching.Level + Upwelling | Area) +
                      (1 | Colony), data=lipids_brms, family = lognormal, prior = prior, sample_prior = TRUE, iter=4000, chains=4) 

pp_check(lipid_brms_3)
plot(marginal_effects(lipid_brms_3, probs = c(0.05, 0.95), ask =FALSE))

####3Model with neither
lipid_brms_4 <- brm(Lipid  ~ 1 + Date + Bleaching.Level + Upwelling + 
                      (1 | Area) +
                      (1 | Colony), data=lipids_brms, family = lognormal, prior = prior, sample_prior = TRUE, iter=4000, chains=4)

pp_check(lipid_brms_4)
plot(marginal_effects(lipid_brms_4, probs = c(0.05, 0.95)))

#leave on out (LOO) cross validation to compare model fits
LOO <- LOO(lipid_brms, lipid_brms_2, lipid_brms_3, lipid_brms_4)

#LOO doesnt really make sense :/ why are the values negative? Do I chose the model with the more negative values or the least negative value? All LOO info online deals with postive values. 

summary(lipid_brms_3)

#Hypopthesis Testing https://www.rdocumentation.org/packages/brms/versions/1.4.0/topics/hypothesis.brmsfit

hypothesis(lipid_brms_3, "UpwellingYes > Intercept", class = "b", group = "Upwelling")
hypothesis(lipid_brms_3, "Bleaching.Level > Intercept", class = "b", group = "Bleaching")

#this doesnt make much sense what is "inf" why is there no value? https://groups.google.com/forum/#!topic/brms-users/4rRYECc6b0c

#marginal effects by Area

conditions <- data.frame(Area = unique(lipids_brms$Area))
rownames(conditions) <- unique(lipids_brms$Area)
plot(marginal_effects(lipid_brms_3, conditions = conditions, 
                      re_formula = NULL, method = "predict"), 
     ncol = 2, points = TRUE)



############# brms Examples 
help("brm") 
vignette("brms_distreg")
help("update.brmsfit")
help("brmsformula")
vignette("brms_nonlinear")
vignette("brms_overview")
