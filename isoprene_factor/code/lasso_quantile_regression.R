install.packages(c("quantreg", "glmnet","caret"))
library(quantreg)
library(glmnet)
library(caret)

# Set the seed for reproducibility
set.seed(123)

df <- read.csv("../data/cleaned_regress_for_R.csv")
correlation <- cor(df$J, df$F)
print(correlation)

# Create a partition
partition <- createDataPartition(df$F, p = 0.1, list = FALSE)

# Split the data using the partition
train_data <- df[partition, ]
test_data <- df[-partition, ]

tau_value <- 0.95  # Adjust this for the desired quantile

qr_model <- rq(F ~ gpp + jmax + lue + J + temp + swdown + vcmax, data=train_data, tau=tau_value)
summary(qr_model)  # View model summary














df <- read.csv("../data/cleaned_regress_for_R.csv")

# Create a partition
partition <- createDataPartition(df$F, p = 0.15, list = FALSE)

# Split the data using the partition
train_data <- df[partition, ]
test_data <- df[-partition, ]

#+ biome + sm_stress_binary + sm_stress_continuous + vcmax

# Separate independent and dependent variables
X <- model.matrix(F ~ gpp + jmax + lue + J + temp + swdown + vcmax , data=train_data)[,-1]
y <- train_data$F

# Quantile regression with Lasso
tau_value <- 0.95  # Adjust this for the desired quantile
qr_lasso <- rq.fit.lasso(X, y, tau=tau_value)

# Get coefficients
print(coef(qr_lasso))




# Predict on test data
X_test <- model.matrix(F ~ gpp + jmax + lue + J + temp + swdown + vcmax, data=test_data)[,-1]
y_test <- test_data$F
y_pred <- predict(qr_lasso, newx = X_test)

# Compute MAE (Mean Absolute Error)
mae <- mean(abs(y_test - y_pred))
print(paste("Mean Absolute Error (MAE):", round(mae, 3)))

# Plot predicted vs actual values for visual inspection
plot(y_test, y_pred, main="Actual vs Predicted", xlab="Actual values", ylab="Predicted values")
abline(a=0, b=1, col="red")  # Adds a y=x line for reference


