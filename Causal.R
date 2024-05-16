## Load packages 

library(tableone) 
library(table1) 
library(ipw) 
library(sandwich) 
library(survey) 
library(ggplot2) 
library(tidyverse)

# Part 1: IPW for point treatment

## import the data

rhc <- read.csv("rhc.csv")

## step 1: select variables for propensity score model

### explore the data
treatment <- as.numeric(rhc$swang1=='RHC') 
died <- as.numeric(rhc$death=='Yes')
age <- rhc$age
female <- as.numeric(rhc$sex=='Female') 
ARF <- as.numeric(rhc$cat1=='ARF')
CHF <- as.numeric(rhc$cat1=='CHF')
cirr <- as.numeric(rhc$cat1=='Cirrhosis')
colcan <- as.numeric(rhc$cat1=='Colon Cancer') 
coma <- as.numeric(rhc$cat1=='Coma')
COPD <- as.numeric(rhc$cat1=='COPD')
lungcan <- as.numeric(rhc$cat1=='Lung Cancer') 
MOSF <- as.numeric(rhc$cat1=='MOSF w/Malignancy') 
sepsis <- as.numeric(rhc$cat1=='MOSF w/Sepsis') 
meanbp1 <- rhc$meanbp1
apache <- rhc$aps1

mydata <- as.data.frame(cbind(treatment, died, age, female, ARF, CHF, cirr, colcan, coma, COPD, lungcan, MOSF, sepsis, meanbp1, apache))

## step 2: calculate the propensity score 
psmodel <- glm(treatment ~ age + female + meanbp1 + apache + ARF + CHF + cirr + colcan + coma + lungcan + MOSF + sepsis,
               family = binomial(link = "logit"))
summary(psmodel)$coefficients

### calculate the ps
ps <- predict(psmodel, type = "response") 
propensity <- cbind(mydata, ps)

### plot the ps
ggplot(propensity, aes(x=ps, fill=as.factor(treatment), group=as.factor(treatment))) + 
  geom_histogram(aes(y=-1*..density..) , data = ~ subset(., treatment =="0")) + 
  geom_histogram(aes(y=..density..), data = ~ subset(., treatment =="1")) +
  theme_classic() + scale_fill_discrete(name = "Treatment group", labels=c("Control", "RHC")) + 
  labs(x="Propensity Score", y="") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

## step 3: calculate weights

### create weights
weight <- ifelse(treatment==1, 1/(ps), 1/(1-ps))
propensity <- propensity %>% mutate(weight=ifelse(treatment==1, 1/ps, 1/(1-ps)))

### apply weights to the data
weighted_data <- svydesign(ids = ~1, data= mydata, weights = ~weight)

## step 4: check assumptions of IPW

### conditional exchangeability
weighted_table <- svyCreateTableOne(vars = c("age","female", "meanbp1", "apache", "ARF","CHF", 
                                             "cirr", "colcan", "coma", "COPD", "lungcan", "MOSF",
                                             "sepsis"), strata = "treatment", 
                                    data= weighted_data, test= FALSE, smd= TRUE)
as.data.frame(print(weighted_table, smd = TRUE)) # check SMD of all confounders

#### compare SMD before and after weighting
raw_data <- svydesign(ids = ~1, data= mydata)
table_raw <- svyCreateTableOne(vars = c("age","female", "meanbp1", "apache", "ARF","CHF",
                                        "cirr", "colcan", "coma", "COPD", "lungcan", "MOSF",
                                        "sepsis"), strata = "treatment", 
                               data= raw_data, test= FALSE, smd= TRUE)

smd <- as.data.frame(cbind(SMD_After= ExtractSmd(weighted_table), SMD_Before = ExtractSmd(table_raw))) 
colnames(smd) <- c("SMD_after_weighting", "SMD_before_weighting")
smd$Factor <- rownames(smd)
smd_plot <- ggplot(smd, aes(y = Factor, x =SMD_after_weighting)) +
  geom_point(color="indianred") + xlim(0,0.6) +
  geom_point(aes(y = Factor, x =SMD_before_weighting), color="royalblue") + 
  theme_bw() + labs(x="Standardized mean difference \n Red = after weighting, Blue = before weighting") + 
  geom_vline(xintercept = 0.1, linetype = "dashed")
smd_plot

### positivity--there are both exposed and unexposed individuals at each level of every confounder
table_raw # check baseline characteristics in the original data

#### check weight of individual

## Bin width defaults to 1/30 of the range of the data. Pick better value with ## ‘binwidth‘.
ggplot(propensity, aes(weight)) + geom_dotplot() + theme_bw()

### consistency--about clear definition of treatment/exposure

## step 5: Fir marginal structural model
glm_model <- glm(died ~ treatment, weights=weight, family = binomial(link = "log"))
beta_ipw <- coef(glm_model) 
beta_ipw
### use asymptotic sandwich variance
SE <- sqrt(diag(vcovHC(glm_model, type = "HC0")))

### get estimate and 95% CI
as.data.frame(cbind(causal_RR = exp(beta_ipw[2]),
                    lower_CI = exp(beta_ipw[2] - 1.96*SE[2]),
                    upper_CI = exp(beta_ipw[2] + 1.96*SE[2])))
### causal risk difference
glm_model2 <- glm(died ~ treatment, weights=weight, family = binomial(link = "identity")) 
beta_ipw2 <- coef(glm_model2)
SE2 <- sqrt(diag(vcovHC(glm_model2, type = "HC0")))

as.data.frame(cbind(causal_RD = beta_ipw2[2],
                    lower_CI = beta_ipw2[2] - 1.96*SE2[2],
                    upper_CI = beta_ipw2[2] + 1.96*SE2[2]))

## Alternative: use IPW package

### calculate weights
weight_model <- ipwpoint(exposure = treatment, family = "binomial", link = "logit", 
                         denominator = ~ age + female + meanbp1 + apache + ARF + CHF + cirr +
                           colcan + coma + lungcan + MOSF + sepsis, data=mydata)
summary(weight_model$ipw.weights)

### plot weights
ipwplot(weights = weight_model$ipw.weights, logscale = FALSE, xlim=c(0,18), xlab="Weight")

## Fit MSM using svyglm from survey package (get sandwich estimators directly)
mydata$wt <- weight_model$ipw.weights
### msm for risk difference
msm <- svyglm(died ~ treatment, design = svydesign(~ 1, weight = ~wt, data = mydata))
coef(msm)
confint(msm)

## IPW using truncated weights at 1st and 99th percentiles
weight_model_trunc <- ipwpoint(exposure = treatment, family = "binomial", link = "logit", 
                               denominator = ~ age + female + meanbp1 + apache + ARF + CHF + cirr +
                                 colcan + coma + lungcan + MOSF + sepsis, data=mydata, trunc = 0.01)
mydata$wt2 <- weight_model_trunc$ipw.weights
msm2 <- svyglm(died ~ treatment, design = svydesign(~ 1, weight = ~wt2, data = mydata)) 
coef(msm2)
confint(msm2)

# Part 2: IPW for time-dependent confounders and informative censoring 

## use the haartdat dataset from ipw
data("haartdat")

## step 1: calculate weights
temp <- ipwtm(exposure = haartind, family = "survival", numerator = ~ sex + age, 
              denominator = ~cd4.sqrt + sex + age, id = patient, tstart = tstart, timevar = fuptime, 
              type = "first", data = haartdat)
summary(temp$ipw.weights)

### compare with stabilized weights
temp.unstab <- ipwtm(exposure = haartind, family = "survival",
                     denominator = ~cd4.sqrt + sex + age, id = patient, tstart = tstart, timevar = fuptime, 
                     type = "first", data = haartdat)
summary(temp.unstab$ipw.weights)

## step 2: investigating weights
ipwplot(weights = temp$ipw.weights, timevar = haartdat$fuptime, binwidth = 100, ylim = c(-1.5, 1.5), 
        main="Stabilized weights", xaxt="n", yaxt="n")
axis(side = 1, at = c(0, 5, 10, 15, 20, 25, 30, 35),
     labels = as.character(c(0, 5, 10, 15, 20, 25, 30, 35)*100))
axis(side = 2, at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
     labels = as.character(c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)))

## step 3: Fit the MSM using robust variance estimator
mort_model <- coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient), 
                    data = haartdat, weights = temp$ipw.weights)
summary(mort_model)

## Estimate the inverse probability of censoring weights (IPCW)
temp2 <- ipwtm(exposure = dropout, family = "survival", numerator = ~ sex + age, 
               denominator = ~cd4.sqrt + sex + age, id = patient, tstart = tstart, timevar = fuptime, 
               type = "first", data = haartdat)
summary(temp2$ipw.weights)

## Fit MSM for IPCW
mort_model2 <- coxph(Surv(tstart, fuptime, event) ~ haartind + cluster(patient), 
                     data = haartdat, weights = temp$ipw.weights*temp2$ipw.weights)
summary(mort_model2)

## comparing with time-varying cox model without IPW
mort_model3 <- coxph(Surv(tstart, fuptime, event) ~ haartind +
                       cluster(patient) + cd4.sqrt + sex + age, data = haartdat)
summary(mort_model3)
