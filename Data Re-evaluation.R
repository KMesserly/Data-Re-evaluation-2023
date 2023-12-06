#install packages
install.packages("lme4")
install.packages("Matrix")
install.packages("dplyr")
install.packages("optimx")
install.packages("performance")

#read in packages 
library(lme4)
library(Matrix)
library(dplyr)
library(optimx)
library(performance)
library(ggplot2)

# Read in data 
telomere <- read.csv("Data for telomere analyses.csv")

#Data Cleaning 
# Remove rows with missing data in the 'tel' column - matches what is published
# with 127 observations remaining for 12 variables
telomere <- telomere[complete.cases(telomere$tel), ]

# Remove rows with missing data in the 'sex' column - matches what is published
# with 116 observations remaining for 12 variables
telomere$sex[telomere$sex == ""] <- NA
telomere <- telomere[complete.cases(telomere$sex), ]

#Summary - output matches what is published in the paper
summary(telomere$tel)

# Set sum-to-0 contrasts for the categorical variable
telomere$TREAT <- as.factor(telomere$TREAT)
telomere$TREAT <- relevel(telomere$TREAT, ref = "CTR")

telomere$sex <- as.factor(telomere$sex)
telomere$sex <- relevel(telomere$sex, ref = "F")

telomere$ageclass <- as.factor(telomere$ageclass)
telomere$ageclass <- relevel(telomere$ageclass, ref = "A")

# Standardize the continuous variable
# datawizard::standardize is used because it is compatable with the prformance
#check used later while the function scale is not
telomere$brood_std <- datawizard::standardize(telomere$BS)
telomere$hatchorder_std <- datawizard::standardize(telomere$hatchorder)
telomere$K_std <- datawizard::standardize(telomere$K)
telomere$mass_std <- datawizard::standardize(telomere$Mass)

#Initial GLM model as in supplemental- Best GVIF Values of the 2 attempted models
fit_tel_init <- glmer(tel ~ TREAT * ageclass +brood_std + hatchorder_std + 
                        sex + K_std + mass_std + (1 |band), data = telomere, 
                      family = Gamma(link = "log"))
#Model fails to converge
car::vif(fit_tel_init)

#Initial GLM model as stated in published paper
#fit_tel_init <- glmer(tel ~ TREAT * ageclass + BS + hatchorder + sex + K_std + 
#mass_std + (1 |band), data = telomere, family = Gamma(link = "log"))

#car::vif(fit_tel_init)

#Solution for collinearity is to create a variable called Mass index
# Calculate the mean and standard deviation of mass for each age class
age_class_stats <- telomere %>%
  group_by(ageclass) %>%
  summarise(mean_mass = mean(Mass), sd_mass = sd(Mass))

# Merge the calculated statistics back into the original data frame
telomere <- merge(telomere, age_class_stats, by = "ageclass", all.x = TRUE)

# Calculate the mass index using the formula
telomere$mass_index <- with(telomere, (Mass - mean_mass) / sd_mass)

#Refit the model using mass_index
fit_tel_init2 <- glmer(tel ~ TREAT * ageclass +brood_std + hatchorder_std + 
                         sex + K_std + mass_index + (1 |band),
                      data = telomere, family = Gamma(link = "log"),
                      control = glmerControl(optimizer ='optimx',
                                             optCtrl=list(method='L-BFGS-B')))
summary(fit_tel_init2)

#Check for colinearity issues - GVIF is below 1o for all variables -Good
car::vif(fit_tel_init2)

#opens viewing pane for graphic outside R for windows computers
windows()

#Check the performance of the data graphically
performance::check_model(fit_tel_init2)
