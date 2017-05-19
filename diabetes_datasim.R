# diabetes data simulation project

load(file = "data2.rdata")
summary(data2)


logistic <- function(t) 1 / (1 + exp(-t))


########################################################

# start datasim

N <- 1e5

generate_dataset <- function(N) {
  set.seed(123)
  age <- sample(c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"),
                N, replace=T, prob=c(0.001, 0.007, 0.016, 0.037, 0.095, 0.17, 0.221, 0.257, 0.169, 0.027))
  age_score <- sample(1:10, N, replace=T, prob=c(0.001, 0.007, 0.016, 0.037, 0.095, 0.17, 0.221, 0.257, 0.169, 0.027))
  age <- ifelse(age_score==1,"0-10",(ifelse(age_score==2,"11-20",(ifelse(age_score==3,"21-30",
                    (ifelse(age_score==4,"31-40",(ifelse(age_score==5,"41-50",
                    (ifelse(age_score==6,"51-60",(ifelse(age_score==7,"61-70",
                    (ifelse(age_score==8,"71-80",(ifelse(age_score==9,"81-90","91-100")))))))))))))))))
  gender <- sample(c('Female','Male'), N, replace = T, prob = c(0.538,0.462))
  #gender_score <- sample(1:2, N, replace = T, prob = c(0.538,0.462))
  #gender <- ifelse(gender_score==1,'Female','Male')
  race <- sample(c('AfricanAmerican','Asian','Caucasian','Hispanic','Other'), N,
                 replace = T,prob = c(0.189,0.006,0.748,0.02,0.037))
  #### generate correlation matrix of 8 numeric variables
  # https://www.r-bloggers.com/simulating-random-multivariate-correlated-data-continuous-variables/
  R = matrix(cbind(1,.32,.19,.47,-.01,-0.01,.07,.22,  
                   .32,1,.06,.27,-.01,-.002,.04,.15,  
                   .19,.06,1,.39,-.02,-.04,-.07,.07,
                   .47,.27,.39,1,.05,.01,.06,.26,
                   -0.01,-0.01,-0.02,0.05,1,0.09,0.11,0.09,
                   -0.01,-0.002,-0.04,0.01,0.09,1,0.27,0.06,
                   0.07,0.034,-0.07,0.06,0.11,0.27,1,0.1,
                   0.22,0.15,0.07,0.26,0.09,0.06,0.1,1
  ), nrow=8)
  U = t(chol(R))
  nvars = dim(U)[1]
  numobs = N
  random.normal = matrix(rpois(nvars*numobs,1.2), nrow=nvars, ncol=numobs);
  X = U %*% random.normal
  newX = t(X)
  raw = as.data.frame(newX)
  names(raw) = c('time_in_hospital','num_lab_procedures','num_procedures','num_medications',
                 'number_outpatient','number_emergency','number_inpatient','number_diagnoses')
  time_in_hospital <- ceiling(raw$time_in_hospital*1.8*logistic((age_score-3.5)/.7))
  num_lab_procedures <- ceiling(raw$num_medications*17.5*logistic((age_score-0)/.2))
  num_procedures  <- ceiling(ifelse(raw$num_procedures>0, floor(raw$num_procedures*0.9*logistic((age_score-4)/.8)),raw$num_procedures+1))
  num_medications <- ceiling(raw$num_medications*9*logistic((age_score-2.5)/.8))
  number_outpatient <- ceiling(rpois(N,0.2)*14*logistic((age_score-3)/1.2))
  number_emergency <- ceiling(rpois(N,0.2)*24*logistic((age_score-5)/2))
  number_inpatient <- ceiling(rpois(N,0.18)*8*logistic((age_score-1)/4))
  number_diagnoses <- ceiling(ifelse(raw$number_diagnoses<10,raw$number_diagnoses*2,floor(rpois(N,0.12)*8*logistic((age_score-1)/4))))
  max_glu_serum <- sample(c('>200','>300','None','Norm'), N, replace = T, prob = c(0.015,0.012,0.947,0.026))
  A1Cresult <- sample(c('>7','>8','None','Norm'), N, replace = T, prob = c(0.037,0.081,0.833,0.049))
  insulin <- sample(c('Down','No','Steady','Up'), N, replace = T, prob = c(0.12,0.466,0.303,.111))
  insulin_score <- sample(1:4,N,replace=T,prob = c(0.12,0.466,0.303,.111) )
  insulin <- ifelse(insulin_score==1,'Down',(ifelse(insulin_score==2,'No',
                       (ifelse(insulin_score==3,'Steady','Up')))))
  change <- sample(c('Ch','No'), N, replace = T, prob = c(0.462,0.538))
  change_score <- sample(1:2,N,replace=T,prob = c(0.462,0.538))
  change <- ifelse(change_score==1,'Ch','No')
  diabetesMed <- sample(c('No','Yes'), N, replace = T, prob = c(0.23,0.77))
  diabetesMed_score <- sample(1:2,N,replace=T,prob = c(0.23,0.77))
  diabetesMed <- ifelse(diabetesMed_score==1,'No','Yes')
  diag1 <- sample(c('circulatory','respiratory','digestive','diabetes','injury',
                    'musculoskeletal','genitourinary','neoplasms','other'), N, replace = T, 
                  prob = c(.299,.142,.093,.086,.068,.049,.05,.138,.074))
  #diag1_score <- sample(1:9,N,replace=T,prob = c(.299,.142,.093,.086,.068,.049,.05,.138,.074))
  admission_source <- sample(c('clinic_referral','emergency','other'), N, replace = T, prob = c(0.101,0.565,.334))
  discharged_to <- sample(c('home','transferred','left_AMA'), N, replace = T, prob = c(0.241,0.753,.006))
  payer_code <- sample(c('Insured','Self_pay'), N, replace = T, prob = c(0.951,0.049))
  score <- 0.67*age_score + 0.5*time_in_hospital + 0.5 * num_lab_procedures -
    0.49 *num_procedures + 0.5 *num_medications + 0.52 *number_outpatient +
    0.56 *number_emergency + 0.6 *number_inpatient + 0.5 *number_diagnoses - 
    1.43 * insulin_score + 0.9 * change_score + 0.6 * diabetesMed_score
  class <- runif(length(score))< .64*logistic((score-30)/2)
  readmitted <- ifelse(class==T,"YES","NO")
  data.frame(gender,race,age,time_in_hospital,num_lab_procedures,num_procedures,
             num_medications,number_outpatient,number_emergency,number_inpatient,number_diagnoses,
             max_glu_serum,A1Cresult,insulin,change,diabetesMed,diag1,
             admission_source,discharged_to,payer_code,readmitted)
}


data <- generate_dataset(N)

plot(data$readmitted, main = "readmissions") # readmission: >50% no readmission


fit_sim <- glm(readmitted~.,data=train, family=binomial)
summary(fit_sim)
logisticPseudoR2s(fit_sim)
fit_main <- glm(readmitted~age+time_in_hospital+num_lab_procedures+num_medications+
                number_outpatient+number_emergency+number_inpatient+number_diagnoses+
                  insulin+diag1,
                data=train, family=binomial)
summary(fit_main)
logisticPseudoR2s(fit_main)

set.seed(123)
inTrain <- createDataPartition(y = data$readmitted, p = .66,list = FALSE)
train <- data[ inTrain,]
test <- data[-inTrain,]
nrow(train) # 67167
nrow(test) # 3459
plot(train$readmitted)

nnet_model <- nnet(formula = readmitted~age+time_in_hospital+num_lab_procedures+num_medications+
                     number_outpatient+number_emergency+number_inpatient+number_diagnoses+
                     insulin+diag1, 
                   data=train, size = 10, maxit = 100)
test$pred <- predict(nnet_model, test, type = "class")
prop.table(table(test$readmitted, test$pred),1)
confusionMatrix(test$pred, test$readmitted)

# Predict labels on test
ypred <- predict(nnet_model,test,type = "class")
# Compute at the prediction scores
test$ypredscore = predict(nnet_model,test,type="raw")
# Check that the predicted labels are the signs of the scores
table(test$ypredscore > 0,test$pred)
# compute ROC curve, precision-recall etc...
pred <- prediction(test$ypredscore,test$readmitted)
# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)
# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec")
plot(perf)
# Plot accuracy as function of threshold
perf <- performance(pred, measure = "acc")
plot(perf)



