
setwd("/Users/Lsrex15/Desktop/Project Things") 
#### import data/data split ####
data = read.csv("10000sample.csv")
data = data[,5:315] #getting rid of qualitative data

#Visualization of the Target variable
dim(data)
table(data$target_kazutsugi) #about the same number of response for each 
hist(data$target_kazutsugi,main = "Plot of Target Variable: Target_kazutsugi", xlab = "Values")

#next we decided that to split data in training and test set we will use proportion 80/20
#we chose set of 1000 observations for the simplicity of calculations
set.seed(123)
N=sample(nrow(data),nrow(data)*0.8)
train_set=data[N,]#N -rows
test_set=data[-N,]# Not N rows

#### Least Squares Regression Model ####

#We ran CV for least squares regression for educational purposes
#install.packages("DAAG") 
library(DAAG) #We use the cv.lm function of the DAAG library to run the k-fold validation (reference added in final report)
linear_model=cv.lm(data=train_set,form.lm=formula(target_kazutsugi~.),m=5)
attributes(linear_model)
#as we se in the result of output MSE for CV is 0.275
# Also in the output we can see that only a few variables are significant.

# Now, we run a least squares regression on the simple train split

leastsquares_train=lm(target_kazutsugi~.,data=train_set)
leastsquares_train_predict=predict(leastsquares_train)
leastsquares_predict=predict(leastsquares_train,newdata=test_set)
summary(leastsquares_train)
plot(leastsquares_train)
MSE_Train_lsr=mean((train_set$target_kazutsugi-leastsquares_train_predict)^2)
MSE_Test_lsr=mean((test_set$target_kazutsugi-leastsquares_predict)^2)
anovatable = anova(leastsquares_train)
bar <- subset(anovatable, `Pr(>F)`<0.025) #Stores variables that are most statistically significant with very low threshhold.

#### Forward Stepwise Selection ####
library(MASS)

#Suppress all the step outputs 
sink("/dev/null")
forward_train = stepAIC(leastsquares_train,direction ='forward')
sink()
summary(forward_train)
#F statistic - " 1.05 on 310" means it chose all 310 variables
#No subset selection was performed by the Forward Stepwise Regression
#Therefore, the results for this model will be the same as the least squares regression

#### Backward Stepwise Regression ####
library(MASS)

#Suppress all the step outputs - This computation was particularly expensive - ran for over 5 minutes 
#sink("/dev/null")
backward_train = stepAIC(leastsquares_train,direction ='backward',steps = 10)
#sink()
summary(backward_train)
pred_back = predict(backward_train, newdata = train_set)
MSE_train_back = mean((pred_back - data$target_kazutsugi)^2)
backward_train_predict= predict(backward_train,newdata =test_set)
MSE_Test_Backward = mean((backward_train_predict - data$target_kazutsugi)^2) #output is 0.0985822




#### Lasso Regression ####
#install.packages("lars") 
library(lars)
lasso.tr=lars(x=as.matrix(train_set[,1:310]),y=train_set$target_kazutsugi,type="lasso") # note this function takes x and y separately, unlike other functions that take formula y~x
#plot(lasso.tr, plottype="coefficients") #beware the computation expense
coef(lasso.tr)[109,coef(lasso.tr)[109,]!=0] # Print of the first ten nonzero predictors

plot(lasso.tr)

# select tuning parameter s by CV, the lower the better but one standard deviation rule applies
cvlasso.tr=cv.lars(x=as.matrix(train_set[,1:310]),y=train_set$target_kazutsugi,type="lasso")
bestfrac = cvlasso.tr$index[which.min(cvlasso.tr$cv)]
bestfrac


##### Now, for the CV-Lasso
x=as.matrix(train_set[,1:310])
y=as.matrix(train_set$target_kazutsugi)
K=5
n=nrow(x)
sam=sample(1:n,n) #create a random list of subject indices
S=list(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
for(i in 1:K)
{
  ind = sam[round((i-1)*n/K + 1):round(i*n/K)] # partition the random list of indices into K equal parts and pick the ith part to be validation
  Xin = x[-ind,]; Yin = y[-ind]
  Xout = x[ind,]; Yout = y[ind]
  for(s in S)
  {
    lasso.in=lars(x=as.matrix(Xin),y=Yin,type="lasso")
    pred.out=predict(lasso.in,newx=Xout,type="fit",mode="fraction",s=s)$fit
    err.te=mean((Yout-pred.out)^2)
    print(i)
    print(s)
    print(err.te)
  }
}

las.coef = predict(lasso.tr, type="coefficients", mode="fraction", s=0)
las.coef
#well actually it turns out that the CV tells us that the best s is 0
#However in this case all the coefficients are shrinked to 0
#In Our case it actually turns out that The less parameter s is the smaller is MSE.
#So after 0 the best parameter is 0.1 with mean CV MSE = 0,14548894.
# Appying one standart deviation rule we do obtain new MSE 0,158621133 - wich is new trashhold.
#starting from the smallest s=0.1 we ca see it is our trashhold, so optimal S=0.1

las.coef = predict(lasso.tr, type="coefficients", mode="fraction", s=0.1)
las.coef
sum(las.coef$coefficients!=0)
#102 variables sellected
pred.te=predict(lasso.tr,newx=test_set[,1:310],type="fit",mode="fraction",s=0.1)$fit
MSE_Test_Lasso=mean((test_set$target_kazutsugi-pred.te)^2)


#### Principal Component Regression ####

#Implementing Bartlett's Sphericity Test
#install.packages("REdaS")
library(REdaS)
B = bart_spher(data[,1:310])
B$p.value #yield 0 which ensures heteroskedacity in data

#implementing PCR
library(pls)
pcr.tr = pcr(target_kazutsugi~.,data=train_set, validation="CV")
summary(pcr.tr)
#well actually i did run it several times and it seems like intercept is the smallest
#after the intercept 6 components are the smallest
validationplot(pcr.tr)
pred.tr = predict(pcr.tr,data=train_set , ncomp=6)
MSE_TR_PCR=mean((train_set$target_kazutsugi-pred.tr)^2)
pred.test=predict(pcr.tr,newdata=test_set, ncomp=6)
MSE_TEST_PCR=mean((test_set$target_kazutsugi-pred.test)^2)

#60 comps for 80 percent variance explained
pred.tr.60 = predict(pcr.tr,data=train_set , ncomp=60)
MSE_TR_PCR.60=mean((train_set$target_kazutsugi-pred.tr.60)^2)
pred.test.60=predict(pcr.tr,newdata=test_set, ncomp=60)
MSE_TEST_PCR.60=mean((test_set$target_kazutsugi-pred.test)^2)


#### Ridge Regression ####
#install.packages("MASS") 
library(MASS)
ridreg.tr=lm.ridge(target_kazutsugi~., data=train_set,lambda=seq(0,5000,len=10))
ridreg.tr$coef
plot(ridreg.tr)
select(ridreg.tr)
# i run ridreg.tr=lm.ridge(target_kazutsugi~., data=train_set,lambda=seq(0,5000,len=10)) for different values
# of len for 50, 1000, and 5000 it tends to choose the highest lamda possible
#prediction of ridge for lambda=500
ridreg.tr=lm.ridge(target_kazutsugi~., data=train_set,lambda=5000)
X.tr=as.matrix(cbind(rep(1,nrow(train_set)),train_set[,1:310]))#original X with vector of 1's for intercept on training data
pred.tr=X.tr%*%coef(ridreg.tr)
MSE_Train_Ridge=mean((train_set$target_kazutsugi-pred.tr)^2)
X.te=as.matrix(cbind(rep(1,nrow(test_set)),test_set[,1:310])) #original X with vector of 1's for intercept on test data
pred.te=X.te%*%coef(ridreg.tr)
MSE_Test_Ridge=mean((test_set$target_kazutsugi-pred.te)^2) #test error


#### Logistic Regression ####
install.packages("nnet")
library(nnet)

g = train_set$target_kazutsugi*4 + 1

# Create Indicator Matrix Y #
Y_data = matrix(0,nrow=dim(train_set)[1],ncol =5)
for (i in 1:5){
  Y_data[g == i , i ]=1
}

newset = cbind(train_set[,1:310],Y_data)

multinomModel <- multinom(Y_data~ ., data=newset,maxit=1000,MaxNWts=1585) # multinom Model
summary (multinomModel) # model summary

g_t = test_set$target_kazutsugi*4 + 1

# Create Indicator Matrix Y #
Y_data_t = matrix(0,nrow=dim(test_set)[1],ncol =5)
for (i in 1:5){
  Y_data_t[g_t == i , i ]=1
}

newset_t = cbind(test_set[,1:310],Y_data_t)

predicted_scores <- predict (multinomModel, newset_t, "probs") # predict on new data
predicted_scores

predicted_class <- predict (multinomModel, newset_t)
predicted_class

table(predicted_class, g_t)

mean(as.character(predicted_class) != as.character(g_t))

oldset = read.csv("1000sample.csv")
oldset = oldset[,5:315]
g_0 = oldset$target_kazutsugi*4 + 1

# Create Indicator Matrix Y #
Y_data_0 = matrix(0,nrow=dim(oldset)[1],ncol =5)
for (i in 1:5){
  Y_data_0[g_0 == i , i ]=1
}

oldset_t = cbind(oldset[,1:310],Y_data_0)

predicted_class <- predict (multinomModel, oldset_t)
predicted_class

table(predicted_class, g_0)

mean(as.character(predicted_class) != as.character(g_0))

xc = train_set[,1:310]
sv = svd(xc)
U = sv$u # normalized scores
U[,1] # PC1 normalized scores
V = sv$v # pc directions
pc1 = V[,1] #PC1 direction
pc2 = V[,2]
plot = (pc1)
D = sv$d # d_j's
U%*%diag(D) # unnormalized score
Z=xc%*%V # unnormalized score, same as U*D
Z[,1] #PC1 unnormalized scores
D^2/nrow(xc)

ord = order(abs(V[,5]),decreasing=TRUE)
x = as.matrix(xc[,ord[1:10]])
colnames(x)
