lambda = 30
lambdaX = seq(lambda,lambda/10,by = -1)
featuresY = sample(1:100, length(lambdaX))
accuracyZ = sample(1:100, length(lambdaX))

dfx = data.frame(lambdaX, featuresY, accuracyZ)

with(dfx, symbols(x=lambdaX, y=featuresY, circles=accuracyZ, inches=1/3,
                  ann=F, bg="steelblue2", fg=NULL))




lambdaX = seq(lambda,lambda/10,by = -1)
featuresX <- c()
accuracyZ <- c()

for (lambda in lambdaX){
  
  # calculate the linear regression relation between y and x1
  model <- glmnet(x,y, family = "gaussian", alpha =1 , lambda = lambda , standardize = FALSE)
  
  # print the linear regression relation summary
  # print(summary(model))
  
  coef.fit <- coef(model, s = lambda)[2:nrow(x)]
  
  features.in <- which(abs(coef.fit) > 0 )
  
  featuresX = append(featuresX , length(features.in))
  
  x = x[, features.in]
  
  relation <- lm(y~x)
  S = summary(relation)
  accuracyZ = append(accuracyZ , S$r.squared)
}

dfx = data.frame(lambdaX, featuresY, accuracyZ)

with(dfx, symbols(x=lambdaX, y=featuresY, circles=accuracyZ, inches=1/3,
                  ann=F, bg="steelblue2", fg=NULL))
