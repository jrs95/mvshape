###################################################################
## Regression in R                                               ##
##                                                               ##
## James Staley                                                  ##
## University of Bristol                                         ##
## Email: js16174@bristol.ac.uk                                  ##
###################################################################

## Working directory
setwd("E:/UoB/Presentations/R 18-04-16")

## Read data
data<-read.csv("data.csv")
data[1:5,]

## Summarise data
summary(data$bmi)
summary(data[,c("age", "bmi")]) 

## Simple linear regression
model <- lm(bmi~age, data=data)
summary(model)

## Multiple linear regression
model <- lm(bmi~age+sex, data=data)
summary(model)
names(summary(model))
coef <- summary(model)$coefficients

## Generalised linear models: linear regression
model <- glm(bmi~age, data=data, family=gaussian)
summary(model)

## Generalised linear models: logistic regression
model <- glm(disease~age, data=data, family=binomial)
summary(model)

## Survival analysis
library(survival)
km <- survfit(Surv(time, disease) ~ 1, data=data)
plot(km, ylim=c(0.85,1))
model <- coxph(Surv(time, disease) ~ age, data=data)
summary(model)

## Meta-analysis
model1 <- lm(bmi~age, data=data[data$cohort==0,])
model2 <- lm(bmi~age, data=data[data$cohort==1,])
beta1 <- summary(model1)$coefficients[2,1]
beta2 <- summary(model2)$coefficients[2,1]
beta <- c(beta1, beta2)
se1 <- summary(model1)$coefficients[2,2]
se2 <- summary(model2)$coefficients[2,2]
se <- c(se1, se2)
library(metafor)
model <- rma(yi=beta, sei=se, method="DL")
summary(model)
forest(model)




