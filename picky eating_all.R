###Orca prey slection study
library(randomForest)
library(caret)
library(MASS)
library(ggordiplots)
library(klaR)

###randomForest
# Load the CSV file into a data frame
data <- read.csv("reduced_targeted_Merged_Data.csv")

# Convert Targeted to numeric (0/1)
data$Targeted <- ifelse(data$Targeted == "Y", 1, 0)

# Select relevant columns, missing values retained
data_clean <- (data[, c("Lipid", "Targeted", "MUFA", "PUFA","SFA", "Calories.kj.g", "Protein")])

# Convert MUFA and PUFA to numeric if they are factors or characters
data_clean$MUFA <- as.numeric(as.character(data_clean$MUFA))
data_clean$PUFA <- as.numeric(as.character(data_clean$PUFA))

# Ensure Targeted is numeric
data_clean$Targeted <- as.factor(data_clean$Targeted)

table(data_clean$Targeted)

####random forest analysis

#impute missing data using rf algorithm
data.imputed <- rfImpute(Targeted ~ ., data_clean, ntree=1000)

#create training and testing data sets
set.seed(222)
ind <- sample(2, nrow(data_clean), replace = TRUE, prob = c(0.7, 0.3))
train <- data.imputed[ind==1,]
test <- data.imputed[ind==2,]

#run randomforest 
set.seed(333)

#rf <- randomForest(Targeted~.,data=train, proximity=TRUE,na.action = na.roughfix) #alterantive way of dealing with missing values - gives very similar result

rf <- randomForest(Targeted~.,data=train, proximity=TRUE) 

#show the basic results
print(rf)

#get the confusion matrix (training data)
p1 <- predict(rf, train)
confusionMatrix(p1, train$Targeted)

#get the confusion matrix testing data
p2 <- predict(rf, test)
confusionMatrix(p2, test$Targeted)

#plot the resuklts after x trees
plot(rf)

#MDS plot of the daat
MDSplot(rf, train$Targeted)

#calculate the importance values
imp<-importance(rf)

#plot the importance values
varImpPlot(rf)

# Make predictions on the test set
predictions <- predict(rf, test)

# Evaluate accuracy
accuracy <- sum(predictions == test$Targeted) / nrow(test)
cat("Accuracy:", accuracy, "\n")

##GLMM

##Lipid:
# Load necessary libraries
library(patchwork)
library(dplyr)
library(MuMIn)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(gridExtra) 


# Fixed effects: Fixed effects are constant across individuals or groups. 
# Random effects: Random effects account for variability and differences between different entities or subjects within a larger group.
### every new variable replaced and rerun

# Load the CSV file into a data frame
data <- read.csv("reduced_targeted_Merged_Data.csv")

# Ensure the relevant columns are treated as factors
data$Species2 <- as.factor(data$Species2)
data$Targeted <- as.factor(data$Targeted)

# Check and convert Lipid column if necessary
if (!is.numeric(data$Lipid)) {
  data$Lipid <- as.numeric(as.character(data$Lipid))
}

# Plot histograms for the response variable in original and transformed forms
# Original data
hist(data$Lipid, main = "Histogram of Lipid (Original)", xlab = "Calories kj/g", breaks = 20)

# Log-transformed data
hist(log(data$Lipid), main = "Histogram of Lipid (Log Transformed)", xlab = "Log Calories kj/g", breaks = 20)

# Square root-transformed data
hist(sqrt(data$Lipid), main = "Histogram of Calories (Square Root Transformed)", xlab = "Square Root Calories kj/g", breaks = 20)

##### clean data
# Remove rows with missing values in relevant columns
data_clean <- na.omit(data[, c( "Lipid", "Targeted", "Species2")])

##### Run GLM model #####
# Fit GLMM model
glmm_model <- lmer(sqrt(Lipid) ~ Targeted + (1| Species2), data = data_clean)

# Diagnostic plots

# Q-Q plot of residuals
qqnorm(residuals(glmm_model))

residuals_df <- data.frame(Fitted = runif(100), Residuals = rnorm(100))

# Plot residuals vs fitted
residual_plot <- ggplot(residuals_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Plot") +
  theme_bw()

# Display the plot
print(residual_plot)

# Summary of the GLM model
summary(glmm_model)

# Set options to fail on NA values
options(na.action = "na.fail")

# Perform dredging
glmm_model_dredge <- dredge(glmm_model)

# Print the dredge results
glmm_model_dredge # Ranks models based on AICc values. Lower AIC = better model fit!

### Calculate R-squared values ###

# Get the first model from the dredged models
first_model <- get.models(glmm_model_dredge, 1)[[1]]

# Calculate the overall R2 for the full model (first model)
overall_r2_first_model <- r.squaredGLMM(first_model)

# Print the overall R-squared
print(overall_r2_first_model)

#### predicted values
glmm_bestmod <- get.models(glmm_model_dredge, 1)[[1]]

# Predictions for interaction of Targetted and Species
glmm_predicted_targeted <- ggpredict(glmm_model, terms = c("Targeted"))

# Plot mean predicted values with error bars for SE
ggplot(glmm_predicted_targeted, aes(x = x, y = predicted)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(
    x = "Targeted",
    y = "Predicted log(Lipid)",
    title = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  guides(color = FALSE)  # Remove legend for color mapping if not needed


##Calories
# Load necessary libraries
library(patchwork)
library(dplyr)
library(MuMIn)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(gridExtra) 


# Fixed effects: Fixed effects are constant across individuals or groups. 
# Random effects: Random effects account for variability and differences between different entities or subjects within a larger group.
### every new variable replaced and rerun

# Load the CSV file into a data frame
data <- read.csv("reduced_targeted_Merged_Data.csv")

# Ensure the relevant columns are treated as factors
data$Species2 <- as.factor(data$Species2)
data$Targeted <- as.factor(data$Targeted)

# Check and convert Calories.kj.g column if necessary
if (!is.numeric(data$Calories.kj.g)) {
  data$Calories.kj.g <- as.numeric(as.character(data$Calories.kj.g))
}

# Plot histograms for the response variable in original and transformed forms
# Original data
hist(data$Calories.kj.g, main = "Histogram of Calories.kj.g (Original)", xlab = "Calories kj/g", breaks = 20)

# Log-transformed data
hist(log(data$Calories.kj.g), main = "Histogram of Calories.kj.g (Log Transformed)", xlab = "Log Calories kj/g", breaks = 20)

# Square root-transformed data
hist(sqrt(data$Calories.kj.g), main = "Histogram of Calories (Square Root Transformed)", xlab = "Square Root Calories kj/g", breaks = 20)

##### clean data
# Remove rows with missing values in relevant columns
data_clean <- na.omit(data[, c( "Calories.kj.g", "Targeted", "Species2")])

##### Run GLM model #####
# Fit GLMM model
glmm_model <- lmer(log(Calories.kj.g) ~ Targeted + (1| Species2), data = data_clean)

# Diagnostic plots

# Q-Q plot of residuals
qqnorm(residuals(glmm_model))

residuals_df <- data.frame(Fitted = runif(100), Residuals = rnorm(100))

# Plot residuals vs fitted
residual_plot <- ggplot(residuals_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Plot") +
  theme_bw()

# Display the plot
print(residual_plot)

# Summary of the GLM model
summary(glmm_model)

# Set options to fail on NA values
options(na.action = "na.fail")

# Perform dredging
glmm_model_dredge <- dredge(glmm_model)

# Print the dredge results
glmm_model_dredge # Ranks models based on AICc values. Lower AIC = better model fit!

### Calculate R-squared values ###

# Get the first model from the dredged models
first_model <- get.models(glmm_model_dredge, 1)[[1]]

# Calculate the overall R2 for the full model (first model)
overall_r2_first_model <- r.squaredGLMM(first_model)

# Print the overall R-squared
print(overall_r2_first_model)

#### predicted values
glmm_bestmod <- get.models(glmm_model_dredge, 1)[[1]]

# Predictions for interaction of Targetted and Species
glmm_predicted_targeted <- ggpredict(glmm_model, terms = c("Targeted"))

# Plot mean predicted values with error bars for SE
ggplot(glmm_predicted_targeted, aes(x = x, y = predicted)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(
    x = "Targeted",
    y = "Predicted log(Calories.kj.g)",
    title = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  guides(color = FALSE)  # Remove legend for color mapping if not needed


##Protein
# Load necessary libraries
library(patchwork)
library(dplyr)
library(MuMIn)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(gridExtra) 


# Fixed effects: Fixed effects are constant across individuals or groups. 
# Random effects: Random effects account for variability and differences between different entities or subjects within a larger group.
### every new variable replaced and rerun

# Load the CSV file into a data frame
data <- read.csv("reduced_targeted_Merged_Data.csv")

# Ensure the relevant columns are treated as factors
data$Species2 <- as.factor(data$Species2)
data$Targeted <- as.factor(data$Targeted)

# Check and convert Protein column if necessary
if (!is.numeric(data$Protein)) {
  data$Protein <- as.numeric(as.character(data$Protein))
}

# Plot histograms for the response variable in original and transformed forms
# Original data
hist(data$Protein, main = "Histogram of Protein (Original)", xlab = "Calories kj/g", breaks = 20)

# Log-transformed data
hist(log(data$Protein), main = "Histogram of Protein (Log Transformed)", xlab = "Log Calories kj/g", breaks = 20)

# Square root-transformed data
hist(sqrt(data$Protein), main = "Histogram of Calories (Square Root Transformed)", xlab = "Square Root Calories kj/g", breaks = 20)

##### clean data
# Remove rows with missing values in relevant columns
data_clean <- na.omit(data[, c( "Protein", "Targeted", "Species2")])

##### Run GLM model #####
# Fit GLMM model
glmm_model <- lmer(log(Protein) ~ Targeted + (1| Species2), data = data_clean)

# Diagnostic plots

# Q-Q plot of residuals
qqnorm(residuals(glmm_model))

residuals_df <- data.frame(Fitted = runif(100), Residuals = rnorm(100))

# Plot residuals vs fitted
residual_plot <- ggplot(residuals_df, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Plot") +
  theme_bw()

# Display the plot
print(residual_plot)

# Summary of the GLM model
summary(glmm_model)

# Set options to fail on NA values
options(na.action = "na.fail")

# Perform dredging
glmm_model_dredge <- dredge(glmm_model)

# Print the dredge results
glmm_model_dredge # Ranks models based on AICc values. Lower AIC = better model fit!

### Calculate R-squared values ###

# Get the first model from the dredged models
first_model <- get.models(glmm_model_dredge, 2)[[1]]

# Calculate the overall R2 for the full model (first model)
overall_r2_first_model <- r.squaredGLMM(first_model)

# Print the overall R-squared
print(overall_r2_first_model)

#### predicted values
glmm_bestmod <- get.models(glmm_model_dredge, 1)[[1]]

# Predictions for interaction of Targetted and Species
glmm_predicted_targeted <- ggpredict(glmm_model, terms = c("Targeted"))

# Plot mean predicted values with error bars for SE
ggplot(glmm_predicted_targeted, aes(x = x, y = predicted)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = predicted - std.error, ymax = predicted + std.error), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(
    x = "Targeted",
    y = "Predicted log(Protein)",
    title = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  guides(color = FALSE)  # Remove legend for color mapping if not needed


