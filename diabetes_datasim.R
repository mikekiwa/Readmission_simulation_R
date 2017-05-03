# diabetes data simulation project

load(file = "data2.rdata")
summary(data2)


logistic <- function(t) 1 / (1 + exp(-t))


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
table(marginal)
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



########################################################

# start datasim
