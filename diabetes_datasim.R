# diabetes data simulation project

load(file = "data2.rdata")
summary(data2)


logistic <- function(t) 1 / (1 + exp(-t))

# simulate factor/categorical variable

x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
x <- sample( c("A","B","C","D"), 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
prop.table(table(x))

x <- rep( c("A","B","C","D"), 10000*c(0.1,0.2,0.65,0.05) )
x <- sample(x, 10000) 

sample(1:4,10000,rep=TRUE,prob=c(.1,.2,.3,.4))


mycat <- cut(runif(10000), c(0, 0.1, 0.3, 0.6, 1))
prop.table(table(mycat))
mycat <- cut(runif(10000), c(0, 0.1, 0.3, 0.6, 1), labels=FALSE)

# simulate random multinomial categorical correlated data
library(GenOrd)
set.seed(1)
# Sets the marginals.
# The values are cumulative so for the first variable the first marginal will be .1, the second is .2, the third is .3, and the fourth is .4
marginal <- list(c(0.1,0.3,0.6),c(0.4,0.7,0.9))
# Checks the lower and upper bounds of the correlation coefficients.
corrcheck(marginal)
# Sets the correlation coefficients
R <- matrix(c(1,-0.6,-0.6,1),2,2) # Correlation matrix
n <- 1000
##Selects and ordinal sample with given correlation R and given marginals.
m <- ordsample(n, marginal, R)
##compare it with the pre-defined R
cor(m)
prop.table(table(m[,1],m[,2]))
margin.table(table(m[,1],m[,2]),1)

chisq.test(m)

gbar < - tapply(m[,1], list(m[,1], m[,2]), length)

par(mfrow=c(1,1))
barplot(gbar, beside=T, col=cm.colors(4), main="Example Bar Chart of Counts by Group",xlab="Group",ylab="Frequency")



library(manipulate)
N <- 1e3
logistic <- function(t) 1 / (1 + exp(-t))

simLOS <- function(N) {
  AGE  <- sample(1:94, N, replace=T)
  COPD <- runif(N)< .12*logistic((AGE-50)/10)
  ASTH <- runif(N)< .095
  HRTD <- runif(N)< .19*logistic((AGE-50)/10)
  DIAB <- runif(N)< .08*logistic((AGE-50)/10)
  CANC <- runif(N)< .19*logistic((AGE-40)/10)
  LOWR <- !(COPD+ASTH+HRTD+DIAB+CANC)
  FLUSHOT <- runif(N)< .34 + .5*(COPD+ASTH+HRTD+DIAB+CANC)
  pflu  <- ifelse(FLUSHOT,.05,.20)
  FLU  <- runif(N)<pflu
  LOS <- 4 + COPD*8 + ASTH*6 +HRTD*6 + DIAB*4 + CANC*5 +FLU*2 +round(rnorm(N,sd=2,mean=0),digits=0)
  LOS[LOS<1]<- 1
  data.frame(AGE,COPD,ASTH,HRTD,DIAB,CANC,LOWR,FLUSHOT,FLU,LOS)    
}
datafr<- (simLOS(N))


head(datafr)



# simulate binomial variable

a <- 1
amx <- ifelse(runif(10000)<.3,-a,a)
sum(amx==1)

a <- 10
amx <- sample(c(-a,a),12,replace=TRUE)
sum(amx==10)

a <- 6
ifelse(rbinom(4,1,.5),-a,a)



library(bindata)
## Construct a binary correlation matrix
rho <- 0.7905694
m <- matrix(c(1,rho,rho,1), ncol=2)   
## Simulate 10000 x-y pairs, and check that they have the specified
## correlation structure
x <- rmvbin(1e5, margprob = c(0.5, 0.5), bincorr = m) 
cor(x)


# simulate Random Multivariate Correlated Data (Continuous Variables)

R = matrix(cbind(1,.80,.2,  .80,1,.7,  .2,.7,1),nrow=3)
U = t(chol(R))
nvars = dim(U)[1]
numobs = 100000
set.seed(1)
random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
orig.raw = as.data.frame(t(random.normal))
names(raw) = c("response","predictor1","predictor2")
cor(raw)
plot(head(raw, 100))
plot(head(orig.raw,100))


library(polycor) 
sim1 <- function(thresh=0.5, r=0.3) { 
  x <- rmvnorm(1000,c(0,0),matrix(c(1,r,r,1), nr=2)) 
  x[x>thresh] <- 2 
  x[x<2] <- 1 
  polychor(x[,1], x[,2]) 
} 

tr <- double(100) 
for(i in 1:100) tr[i] <- sim1() 
summary(tr) 

?manipulate

# CORRELATION OF CATEGORICALS
# example data
set.seed(1)
DF <- data.frame(x=sample(c("Y","N"),100,T),y=sample(c("Y","N"),100,T))

# how to get correlation
DF[] <- lapply(DF,as.integer)
cor(DF)
#            x          y
# x  1.0000000 -0.0369479
# y -0.0369479  1.0000000

# visualize it
library(corrplot)
corrplot(cor(DF))


?aov
?ICC 

t1 = c(164, 172, 168, 177, 156, 195)
t2 = c(178, 191, 197, 182, 185, 177)
t3 = c(175, 193, 178, 171, 163, 176)
t4 = c(155, 166, 149, 164, 170, 168)

val = c(t1, t2, t3, t4)
fac = gl(n=4, k=6, labels=c('type1', 'type2', 'type3', 'type4'))

aov1 = aov(val ~ fac)
summary(aov1)



df <- transform(df, score = 0.01 * (aerobicCap - mean(aerobicCap)) + 0.04 * 
                  (energyIntake - mean(energyIntake))+ 481 * (bmi - mean(bmi)) - 2.1)
df$prob <- logistic(df$score)

# The outcome is then masked into a range of outcome for high and low scores
df$scoreMath <- round(ifelse(df$prob > runif(N), rnorm(N,mean = 94, sd= 4),
                             rnorm(N,mean = 45, sd= 3)),2)

## Helper function to calculate cooefficients
library(manipulate)
manipulate({
  df1 <- transform(df, score = 10^a * (aerobicCap - mean(aerobicCap)) + 10^b * 
                     (energyIntake - mean(energyIntake)) + 10^c * (bmi - mean(bmi)) - 2.1)
  df1$prob <- logistic(df$score)
  hist(df1$prob, breaks=50)
}, a=slider(-9, 9, step=0.1, initial = 0), b=slider(-9, 9, step=0.1, initial = 0), 
c=slider(-9, 9, step=0.1, initial = 0))


# correlated continuous variables
R = matrix(cbind(1,.80,.2,  .80,1,.7,  .2,.7,1),nrow=3)
U = t(chol(R))
nvars = dim(U)[1]
numobs = 100000
set.seed(1)
random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
orig.raw = as.data.frame(t(random.normal))
names(raw) = c("response","predictor1","predictor2")
cor(raw)
plot(head(raw, 100))
plot(head(orig.raw,100))



R = matrix(cbind(1,.32,.19,.47,-.01,-0.01,.07,.22,  
                 .32,1,.06,.27,-.01,-.002,.04,.15,  
                 .19,.06,1,.39,-.02,-.04,-.07,.07,
                 .47,.27,.39,1,.05,.01,.06,.26,
                 -0.01,-0.01,-0.02,0.05,1,0.09,0.11,0.09,
                 -0.01,-0.002,-0.04,0.01,0.09,1,0.27,0.06,
                 0.07,0.034,-0.07,0.06,0.11,0.27,1,0.1,
                 0.22,0.15,0.07,0.26,0.09,0.06,0.1,1
                 ),
           nrow=8)
U = t(chol(R))
nvars = dim(U)[1]
numobs = 10000
set.seed(1)
random.normal = matrix(rpois(nvars*numobs,1.5), nrow=nvars, ncol=numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
orig.raw = as.data.frame(t(random.normal))
names(raw) = c('time_in_hospital','num_lab_procedures','num_procedures','num_medications',
               'number_outpatient','number_emergency','number_inpatient','number_diagnoses')
cor(raw)
plot(head(raw, 100))
plot(head(orig.raw,100))

num_procedures <- rpois(10000,0.8)

set.seed(1)
num_medications <- floor(raw$num_medications*8.9*logistic((age-3.5)/.8))
hist(num_medications)
plot(age,num_medications)
hist(data2$num_medications)
hist(raw$num_medications)
sd(data2$num_medications)
mean(data2$num_medications)
g <- ggplot(data2, aes(x=age, y=num_medications))
g + geom_boxplot()

num_medications 
number_outpatient 
number_emergency  
number_inpatient  
number_diagnoses

########################################################

# start datasim

N <- 1e6
set.seed(123)

generate_dataset <- function(N) {
  age <- sample(c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"),
                N, replace=T, prob=c(0.001, 0.007, 0.016, 0.037, 0.095, 0.17, 0.221, 0.257, 0.169, 0.027))
  gender <- sample(c('Female','Male'), N, replace = T, prob = c(0.538,0.462))
  race <- sample(c('AfricanAmerican','Asian','Caucasian','Hispanic','Other'), N,
                 replace = T,prob = c(0.189,0.006,0.748,0.02,0.037))
  #time_in_hospital <- floor(age * rpois(10000,1.8))/6
  num_lab_procedures <- rnorm(10000,43.1,19.67)
  R = matrix(cbind(1,.32,.19,.47,-.01,-0.01,.07,.22,  
                   .32,1,.06,.27,-.01,-.002,.04,.15,  
                   .19,.06,1,.39,-.02,-.04,-.07,.07,
                   .47,.27,.39,1,.05,.01,.06,.26,
                   -0.01,-0.01,-0.02,0.05,1,0.09,0.11,0.09,
                   -0.01,-0.002,-0.04,0.01,0.09,1,0.27,0.06,
                   0.07,0.034,-0.07,0.06,0.11,0.27,1,0.1,
                   0.22,0.15,0.07,0.26,0.09,0.06,0.1,1
  ),
  nrow=8)
  U = t(chol(R))
  nvars = dim(U)[1]
  numobs = N
  set.seed(1)
  random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
  X = U %*% random.normal
  newX = t(X)
  raw = as.data.frame(newX)
  names(raw) = c('time_in_hospital','num_lab_procedures','num_procedures','num_medications',
                 'number_outpatient','number_emergency','number_inpatient','number_diagnoses')
  time_in_hospital <- raw$time_in_hospital*1.8*logistic((age-4)/.7)
  num_lab_procedures <- raw$num_lab_procedures*17.5*logistic((age-0)/.2)
  num_procedures  <- ifelse(raw$num_procedures>0, floor(raw$num_procedures*0.9*logistic((age-4)/.8)),raw$num_procedures+1)
  max_glu_serum <- sample(c('>200','>300','None','Norm'), N, replace = T, prob = c(0.015,0.012,0.947,0.026))
  A1Cresult <- sample(c('>7','>8','None','Norm'), N, replace = T, prob = c(0.037,0.081,0.833,0.049))
  insulin <- sample(c('Down','No','Steady','Up'), N, replace = T, prob = c(0.12,0.466,0.303,.111))
  change <- sample(c('Ch','No'), N, replace = T, prob = c(0.462,0.538))
  diabetesMed <- sample(c('No','Yes'), N, replace = T, prob = c(0.23,0.77))
  diag1 <- sample(c('circulatory','respiratory','digestive','diabetes','injury',
                    'musculoskeletal','genitourinary','neoplasms','other'), N, replace = T, 
                  prob = c(.299,.142,.093,.086,.068,.049,.05,.138,.074))
  diag2 <- sample(c('circulatory','respiratory','digestive','diabetes','injury',
                    'musculoskeletal','genitourinary','neoplasms','other'), N, replace = T, 
                  prob = c(.313,.107,.041,.126,.024,.017,.082,.185,.105))
  admission_source <- sample(c('clinic_referral','emergency','other'), N, replace = T, prob = c(0.101,0.565,.334))
  discharged_to <- sample(c('home','transferred','left_AMA'), N, replace = T, prob = c(0.241,0.753,.006))
  payer_code <- sample(c('Insured','Self_pay'), N, replace = T, prob = c(0.951,0.049))
  readmitted <- sample(c('<30','>30','NO'), N, replace = T, prob = c(0.538,0.462))
  data.frame(age,gender,race,time_in_hospital,num_lab_procedures,num_procedures,
             max_glu_serum,A1Cresult,insulin,change,diabetesMed,diag1,diag2,
             admission_source,discharged_to,payer_code,readmitted)
}

time_in_hospital 
num_lab_procedures
num_procedures 
num_medications 
number_outpatient 
number_emergency  
number_inpatient  
number_diagnoses



prop.table(table(data2$diag2))

sd(data2$num_procedures)
mean(data2$num_procedures)
plot(data2$discharged_to) # transferred to another facility 70%
hist(data2$num_procedures)
hist(raw$num_procedures)
hist(num_procedures)

set.seed(123)
num_proceduress <- rnorm(10000,43.1,19.67)
hist(num_procedures)
num_lab_procedures <- floor(age * rpois(10000,34))/4.5
hist(num_lab_procedures)
num_lab_procedures <- floor(age * runif(10000,1,132))/12
hist(num_lab_procedures)
set.seed(123)
num_procedures <- rpois(10000,1.4)*1.2*logistic((age-5)/.8)
hist(num_procedures)
plot(age,num_procedures)

df <- data.frame(age,time_in_hospital)

summary(data2)

mean(age)
age <- sample(1:10, 10000, replace=T, prob=c(0.001, 0.007, 0.016, 0.037, 0.095, 0.17, 0.221, 0.257, 0.169, 0.027))
space <- runif(1000, 1, 10)
plot(data2$age)


library(ggplot2)

g <- ggplot(data2, aes(x=age, y=num_procedures))
g + geom_boxplot()

aov1 <- aov(data2$time_in_hospital ~ data2$age)
summary(aov1)

fit <- lm(time_in_hospital~age, data = df)
summary(fit)
plot(df$age,df$time_in_hospital)
