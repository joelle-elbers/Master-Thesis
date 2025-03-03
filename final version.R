

##Set working directory and packages, open dataset
setwd("~/M thesis")
library(haven)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(lavaan)
library(lmtest)
library(car)
library(MASS)
library(ggseg3d)
library(ggseg)
library(psych)

plot(dk)


setwd("~/M thesis")

source("Process.R")


MCC_data_JolleElbers_encrypted <- read_sav("MCC_data_JolleElbers_encrypted.sav")
View(MCC_data_JolleElbers_encrypted)



# Set a seed for reproducibility
set.seed(2785552)  



##Brain plots for figures
##Create OFC plot
ofc_plot <- ggseg(atlas = dk, 
                  mapping = aes(fill = region)) +
  scale_fill_manual(values = c(
    "medial orbitofrontal" = "#1f78b4", 
    "lateral orbitofrontal" = "#a6cee3"  
  ), 
  na.value = "#fce4ec") +               
  theme_void() +                        
  theme(legend.position = "right", 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  ggtitle("Orbitofrontal Cortex")

# Create Amygdala plot
amygdala_plot <- ggseg(atlas = aseg, 
                       mapping = aes(fill = label)) +
  scale_fill_manual(values = c(
    "Left-Amygdala" = "#F4A261",  
    "Right-Amygdala" = "#E63946" 
  ), 
  na.value = "#fce4ec") +                
  theme_void() +                        
  theme(legend.position = "right", 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  ggtitle("Amygdala")

# Combine the two plots side by side
combined_plot <- ofc_plot + amygdala_plot + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)






#####Data exploration 
#copying the data so I can clean it without changing the original file
data_a <- MCC_data_JolleElbers_encrypted
View(data_a)

summary(data_a)

# Get all variables and their labels
value_labels <- sapply(data_a, function(x) attr(x, "labels"))
value_labels



####Combining variables to make the variables for primary analyses

##Filtering wave 3 that I intend to use to see if we remove the brain missings how many cases I have left
brain_variables <- c(
  "Left.Amygdala_3", "Left.Amygdala_1", 
  "Right.Amygdala_3", "Right.Amygdala_1", 
  "EstimatedTotalIntraCranialVol_3", "EstimatedTotalIntraCranialVol_1", 
  "lh_lateralorbitofrontal_area_3", "lh_lateralorbitofrontal_area_1", 
  "lh_medialorbitofrontal_area_3", "lh_medialorbitofrontal_area_1", 
  "rh_lateralorbitofrontal_area_3", "rh_lateralorbitofrontal_area_1", 
  "rh_medialorbitofrontal_area_3", "rh_medialorbitofrontal_area_1", 
  "lh_lateralorbitofrontal_thickness_3", "lh_lateralorbitofrontal_thickness_1", 
  "lh_medialorbitofrontal_thickness_3", "lh_medialorbitofrontal_thickness_1", 
  "rh_lateralorbitofrontal_thickness_3", "rh_lateralorbitofrontal_thickness_1", 
  "rh_medialorbitofrontal_thickness_3", "rh_medialorbitofrontal_thickness_1", 
  "lh_medialorbitofrontal_volume_1", "lh_medialorbitofrontal_volume_3", 
  "lh_lateralorbitofrontal_volume_1", "lh_lateralorbitofrontal_volume_3", 
  "rh_medialorbitofrontal_volume_1", "rh_medialorbitofrontal_volume_3", 
  "rh_lateralorbitofrontal_volume_1", "rh_lateralorbitofrontal_volume_3", 
  "lh_OFC_volume_1", "lh_OFC_volume_3", 
  "rh_OFC_volume_1", "rh_OFC_volume_3"
)


# List of the anxiety measure variables
anxiety_variables <- c("mcc_w05_bsi_anxiety_q4_pp", 
                       "mcc_w05_bsi_anxiety_q9_pp", 
                       "mcc_w05_bsi_anxiety_q14_pp", 
                       "mcc_w05_bsi_anxiety_q18_pp", 
                       "mcc_w05_bsi_anxiety_q20_pp")


### Calculating Structural covariance 
data_a <- data_a %>%
  mutate(
    # Calculate volumes for medial and lateral orbitofrontal cortices for each hemisphere
    lh_medialorbitofrontal_volume_1 = lh_medialorbitofrontal_thickness_1 * lh_medialorbitofrontal_area_1,
    lh_medialorbitofrontal_volume_3 = lh_medialorbitofrontal_thickness_3 * lh_medialorbitofrontal_area_3,
    lh_lateralorbitofrontal_volume_1 = lh_lateralorbitofrontal_thickness_1 * lh_lateralorbitofrontal_area_1,
    lh_lateralorbitofrontal_volume_3 = lh_lateralorbitofrontal_thickness_3 * lh_lateralorbitofrontal_area_3,
    
    rh_medialorbitofrontal_volume_1 = rh_medialorbitofrontal_thickness_1 * rh_medialorbitofrontal_area_1,
    rh_medialorbitofrontal_volume_3 = rh_medialorbitofrontal_thickness_3 * rh_medialorbitofrontal_area_3,
    rh_lateralorbitofrontal_volume_1 = rh_lateralorbitofrontal_thickness_1 * rh_lateralorbitofrontal_area_1,
    rh_lateralorbitofrontal_volume_3 = rh_lateralorbitofrontal_thickness_3 * rh_lateralorbitofrontal_area_3,
    
    # Combine medial and lateral volumes to create the overall OFC volume for each hemisphere
    lh_OFC_volume_1 = lh_medialorbitofrontal_volume_1 + lh_lateralorbitofrontal_volume_1,
    lh_OFC_volume_3 = lh_medialorbitofrontal_volume_3 + lh_lateralorbitofrontal_volume_3,
    rh_OFC_volume_1 = rh_medialorbitofrontal_volume_1 + rh_lateralorbitofrontal_volume_1,
    rh_OFC_volume_3 = rh_medialorbitofrontal_volume_3 + rh_lateralorbitofrontal_volume_3
  )

#Calculating total PDS score (with the 3 items I have)
data_a <- data_a %>%
  mutate(
    pds_total = mcc_w03_pds_q1_c + mcc_w03_pds_q2_c + mcc_w03_pds_q3_c,
  )




#Filtering data:
primary_analyses_set <- data_a %>%
  # Step 1: Exclude rows with missing DHEA hair data for wave 3
  filter(mcc_w03_hair_dhea_level_c != 999) %>%
  
  # Step 2: Group by family ID and randomly select one child per family
  group_by(familyid_rand) %>%
  sample_n(1) %>%
  ungroup() %>%
  
  # Step 3: Exclude rows with missing brain data
  filter(if_all(all_of(brain_variables), ~ !is.na(.))) %>%
  
  # Step 4: Exclude rows with 999s in anxiety measures
  filter(if_all(all_of(anxiety_variables), ~ . != 999))

# View the filtered dataset and its summary
summary(primary_analyses_set)
View(primary_analyses_set)


#Steps for flow chart:
# Step 1: Exclude rows with missing DHEA hair data for wave 3
step1_removed <- nrow(data_a) - nrow(filter(data_a, mcc_w03_hair_dhea_level_c != 999))
cat("Rows removed in Step 1 (missing DHEA):", step1_removed, "\n")

# Step 2: Group by family ID and randomly select one child per family
step2_removed <- nrow(filter(data_a, mcc_w03_hair_dhea_level_c != 999)) - nrow(
  data_a %>%
    filter(mcc_w03_hair_dhea_level_c != 999) %>%
    group_by(familyid_rand) %>%
    sample_n(1) %>%
    ungroup()
)
cat("Rows removed in Step 2 (family selection):", step2_removed, "\n")

# Step 3: Exclude rows with missing brain data
step3_removed <- nrow(data_a %>%
                        filter(mcc_w03_hair_dhea_level_c != 999) %>%
                        group_by(familyid_rand) %>%
                        sample_n(1) %>%
                        ungroup()) - nrow(
                          data_a %>%
                            filter(mcc_w03_hair_dhea_level_c != 999) %>%
                            group_by(familyid_rand) %>%
                            sample_n(1) %>%
                            ungroup() %>%
                            filter(if_all(all_of(brain_variables), ~ !is.na(.)))
                        )
cat("Rows removed in Step 3 (missing brain data):", step3_removed, "\n")

# Step 4: Exclude rows with 999s in anxiety measures
step4_removed <- nrow(data_a %>%
                        filter(mcc_w03_hair_dhea_level_c != 999) %>%
                        group_by(familyid_rand) %>%
                        sample_n(1) %>%
                        ungroup() %>%
                        filter(if_all(all_of(brain_variables), ~ !is.na(.)))) - nrow(
                          primary_analyses_set
                        )
cat("Rows removed in Step 4 (999 in anxiety):", step4_removed, "\n")



##creating SC variables:
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    # Calculate volume change scores for OFC (left and right hemispheres)
    lh_OFC_volume_change = lh_OFC_volume_3 - lh_OFC_volume_1,
    rh_OFC_volume_change = rh_OFC_volume_3 - rh_OFC_volume_1,
    
    # Calculate volume change scores for Amygdala (left and right hemispheres)
    Left.Amygdala_volume_change = Left.Amygdala_3 - Left.Amygdala_1,
    Right.Amygdala_volume_change = Right.Amygdala_3 - Right.Amygdala_1
  )


# Normalize the change scores for each brain region
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    # Normalized change scores for OFC (left and right hemispheres)
    lh_OFC_volume_normalized_change = lh_OFC_volume_change / lh_OFC_volume_1,
    rh_OFC_volume_normalized_change = rh_OFC_volume_change / rh_OFC_volume_1,
    
    # Normalized change scores for Amygdala (left and right hemispheres)
    Left.Amygdala_normalized_change = Left.Amygdala_volume_change / Left.Amygdala_1,
    Right.Amygdala_normalized_change = Right.Amygdala_volume_change / Right.Amygdala_1
  )


# Add direction and magnitude comparison for co-development scores
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    # Direction for OFC and Amygdala changes (increase or decrease)
    direction_lh_OFC = if_else(lh_OFC_volume_normalized_change > 0, "increase", "decrease"),
    direction_rh_OFC = if_else(rh_OFC_volume_normalized_change > 0, "increase", "decrease"),
    direction_Left_Amygdala = if_else(Left.Amygdala_normalized_change > 0, "increase", "decrease"),
    direction_Right_Amygdala = if_else(Right.Amygdala_normalized_change > 0, "increase", "decrease"),
    
    # Magnitude comparison (absolute change) between regions
    magnitude_left_OFC_left_Amygdala = abs(lh_OFC_volume_normalized_change) - abs(Left.Amygdala_normalized_change),
    magnitude_right_OFC_left_Amygdala = abs(rh_OFC_volume_normalized_change) - abs(Left.Amygdala_normalized_change),
    magnitude_right_Amygdala_left_OFC = abs(Right.Amygdala_normalized_change) - abs(lh_OFC_volume_normalized_change),
    magnitude_right_Amygdala_right_OFC = abs(Right.Amygdala_normalized_change) - abs(rh_OFC_volume_normalized_change),
    
    # Combined co-development measure (direction + magnitude comparison)
    left_Amygdala_left_OFC_co_development_combined = Left.Amygdala_normalized_change - lh_OFC_volume_normalized_change,
    left_Amygdala_right_OFC_co_development_combined = Left.Amygdala_normalized_change - rh_OFC_volume_normalized_change,
    right_Amygdala_left_OFC_co_development_combined = Right.Amygdala_normalized_change - lh_OFC_volume_normalized_change,
    right_Amygdala_right_OFC_co_development_combined = Right.Amygdala_normalized_change - rh_OFC_volume_normalized_change
  )

##Have to correct for the time in between measures
primary_analyses_set$time_interval <- primary_analyses_set$mcc_w03_age_c - primary_analyses_set$mcc_w01_age_c

# Time-correct the normalized volume changes
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    lh_OFC_volume_normalized_change_time_corrected = lh_OFC_volume_normalized_change / time_interval,
    rh_OFC_volume_normalized_change_time_corrected = rh_OFC_volume_normalized_change / time_interval,
    Left.Amygdala_normalized_change_time_corrected = Left.Amygdala_normalized_change / time_interval,
    Right.Amygdala_normalized_change_time_corrected = Right.Amygdala_normalized_change / time_interval
  )


# Time-correct the normalized volume changes
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    lh_OFC_volume_normalized_change_time_corrected = lh_OFC_volume_normalized_change / time_interval,
    rh_OFC_volume_normalized_change_time_corrected = rh_OFC_volume_normalized_change / time_interval,
    Left.Amygdala_normalized_change_time_corrected = Left.Amygdala_normalized_change / time_interval,
    Right.Amygdala_normalized_change_time_corrected = Right.Amygdala_normalized_change / time_interval
  )


# Recalculate co-development scores with time-corrected changes
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    left_Amygdala_left_OFC_co_dev_time_corrected = 
      Left.Amygdala_normalized_change_time_corrected - lh_OFC_volume_normalized_change_time_corrected,
    left_Amygdala_right_OFC_co_dev_time_corrected = 
      Left.Amygdala_normalized_change_time_corrected - rh_OFC_volume_normalized_change_time_corrected,
    right_Amygdala_left_OFC_co_dev_time_corrected = 
      Right.Amygdala_normalized_change_time_corrected - lh_OFC_volume_normalized_change_time_corrected,
    right_Amygdala_right_OFC_co_dev_time_corrected = 
      Right.Amygdala_normalized_change_time_corrected - rh_OFC_volume_normalized_change_time_corrected)
    



# Anxiety measure: Creating the average anxiety score variable by adding it up and dividing by total items
primary_analyses_set <- primary_analyses_set %>%
  mutate(
    anxiety_mean_w05 = (mcc_w05_bsi_anxiety_q4_pp + mcc_w05_bsi_anxiety_q9_pp + mcc_w05_bsi_anxiety_q14_pp + mcc_w05_bsi_anxiety_q18_pp + mcc_w05_bsi_anxiety_q20_pp) / 5
  )


## correcting hair_dhea for age and pds total
# Fit the model
model <- lm(mcc_w03_hair_dhea_level_c ~ mcc_w03_age_c + pds_total, data = primary_analyses_set)

primary_analyses_set$dhea_corrected <- resid(model)


# Calculate IQR, Q1, and Q3 for outliers for DHEA
iqr <- IQR(primary_analyses_set$dhea_corrected, na.rm = TRUE)
q1 <- quantile(primary_analyses_set$dhea_corrected, 0.25, na.rm = TRUE)
q3 <- quantile(primary_analyses_set$dhea_corrected, 0.75, na.rm = TRUE)

# Define the lower and upper bounds for outliers
lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

# Identify outliers
primary_analyses_set$outlier <- primary_analyses_set$mcc_w03_hair_dhea_level_c < lower_bound | 
  primary_analyses_set$dhea_corrected > upper_bound

# Count and view outliers
outlier_count <- sum(primary_analyses_set$outlier, na.rm = TRUE)
cat("Number of outliers:", outlier_count, "\n")


# View rows with outliers
outliers <- primary_analyses_set[primary_analyses_set$outlier, ]
print(outliers)

# Visualize outliers using a boxplot
library(ggplot2)
ggplot(primary_analyses_set, aes(y = dhea_corrected)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red") +
  labs(title = "Boxplot of Corrected DHEA Levels",
       y = "Corrected DHEA Levels") +
  theme_minimal()


# Create a new dataset excluding outliers
primary_analyses_set_no_DHEA_outliers <- primary_analyses_set %>%
  filter(!outlier)

# View the filtered dataset
View(primary_analyses_set_no_DHEA_outliers)

# Summary of the new dataset
summary(primary_analyses_set_no_DHEA_outliers)

# Count rows before and after removing outliers
cat("Number of cases before removing outliers:", nrow(primary_analyses_set), "\n")
cat("Number of cases after removing outliers:", nrow(primary_analyses_set_no_DHEA_outliers), "\n")




##3 outliers verwijderd van DHEA (deze passen ook minder binnen de reference values van de andere studies)


####Assumption checks + data transformations + data transformation checks
##DHEA

# Boxplot for mcc_w03_hair_dhea_level_c --> 3 outliers 
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = dhea_corrected)) +
  geom_boxplot() +
  labs(title = "Boxplot of mcc_w03_hair_dhea_level_c", y = "dhea_corrected") +
  theme_minimal() 

# normally distributed
hist(primary_analyses_set_no_DHEA_outliers$dhea_corrected, main = "Distribution of Variable", xlab = "Hair DHEA", breaks = 20)



##Anxiety
# Boxplot for anxiety_mean_w05 --> 3 outliers
library(ggplot2)
ggplot(primary_analyses_set_no_DHEA_outliers, aes(y = anxiety_mean_w05)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red") +
  labs(title = "Boxplot of Anxiety Mean (Wave 5)",
       y = "Anxiety Mean (Wave 5)") +
  theme_minimal()


##calculating the quartiles
lower_bound_anxiety <- quantile(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05, 0.25) - 1.5 * IQR(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05)
upper_bound_anxiety <- quantile(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05, 0.75) + 1.5 * IQR(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05)

primary_analyses_set_no_DHEA_outliers$outlier_anxiety <- primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05 < lower_bound_anxiety | 
  primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05 > upper_bound_anxiety


outlier_count_anxiety <- sum(primary_analyses_set_no_DHEA_outliers$outlier_anxiety, na.rm = TRUE)
cat("Number of outliers for anxiety_mean_w05:", outlier_count_anxiety, "\n")


outliers_anxiety <- primary_analyses_set_no_DHEA_outliers[primary_analyses_set_no_DHEA_outliers$outlier_anxiety, ]
print(outliers_anxiety)


#histograms: --> skewed to the right (very much so, but I dont want to lose the variance completely as 1.5 is not an insane value to have, maybe log-transform? but then I use interpretability? )
hist(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05, main = "Distribution of Variable", xlab = "anxiety mean", breaks = 20) 

##Let's log transform it anyway ---  #okay that did not help
primary_analyses_set_no_DHEA_outliers$log_anxiety_mean_w05 <- log(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05 + 1) 
hist(primary_analyses_set_no_DHEA_outliers$log_anxiety_mean_w05)

##Trying the square root transformation
primary_analyses_set_no_DHEA_outliers$sqrt_anxiety_mean_w05 <- sqrt(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05)
hist(primary_analyses_set_no_DHEA_outliers$sqrt_anxiety_mean_w05)


# Model with original anxiety_mean_w05
model_original <- lm(anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_left_OFC_co_development_combined + mcc_sex_c, 
                     data = primary_analyses_set_no_DHEA_outliers)

# Model with transformed anxiety_sqrt ##later edit: found out I cant do this because it is used for count data and this is not count data 
model_sqrt <- lm(sqrt_anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_left_OFC_co_development_combined + mcc_sex_c, 
                 data = primary_analyses_set_no_DHEA_outliers)

# Check diagnostics for the original model
residuals_original <- model_original$residuals
shapiro.test(residuals_original)  # Shapiro-Wilk test for normality
par(mfrow = c(1, 2))  # Plot two graphs side by side
hist(residuals_original, breaks = 20, main = "Original Residuals", xlab = "Residuals")
qqnorm(residuals_original)
qqline(residuals_original)

# Check diagnostics for the transformed model ##histogram is better
residuals_sqrt <- model_sqrt$residuals
shapiro.test(residuals_sqrt)  
hist(residuals_sqrt, breaks = 20, main = "Transformed Residuals", xlab = "Residuals")
qqnorm(residuals_sqrt)
qqline(residuals_sqrt)


# Plot residuals vs fitted values for original model
plot(model_original$fitted.values, residuals_original, 
     main = "Residuals vs Fitted (Original)", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs fitted values for transformed model ##ehh tbf idk if this improves it
plot(model_sqrt$fitted.values, residuals_sqrt, 
     main = "Residuals vs Fitted (Transformed)", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Compare model fit for the original model
summary(model_original)

# Compare model fit for the transformed model
summary(model_sqrt)

# Compare AIC/BIC for both models ##model fit is better for transformed
AIC(model_original, model_sqrt)
BIC(model_original, model_sqrt)

# Compare coefficients for both models ##they are not THAT different, especially the covariance variable changes a bit
coef(model_original)
coef(model_sqrt)

##overall,there is improvement in the data with the sqrt transformation when it comes to tions. I choose to sqrt transform anxiety_mean


###box-cox transformation to see which of the transformations makes the data better meet the assumptions

# Shift the data to make all values positive
shifted_anxiety <- primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05 + 1  # Adding 1 to shift values

# Fit a linear model to the shifted anxiety variable
model_boxcox <- lm(shifted_anxiety ~ dhea_corrected + left_Amygdala_left_OFC_co_development_combined + mcc_sex_c, 
                   data = primary_analyses_set_no_DHEA_outliers)

# Apply the Box-Cox transformation
boxcox_result <- boxcox(model_boxcox, lambda = seq(-2, 2, by = 0.1))

# Find the best lambda (the value that maximizes the log-likelihood)
best_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
cat("Best lambda:", best_lambda, "\n")

# Apply the Box-Cox transformation to the shifted anxiety_mean_w05 variable with the best lambda
primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05 <- (shifted_anxiety^best_lambda - 1) / best_lambda

# Check the distribution of the transformed data
hist(primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05, main = "Box-Cox Transformed Anxiety Mean", xlab = "Transformed Anxiety Mean")

hist(primary_analyses_set_no_DHEA_outliers$anxiety_mean_w05)


# Model with original anxiety_mean_w05
model_original <- lm(anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_left_OFC_co_dev_time_corrected + mcc_sex_c, 
                     data = primary_analyses_set_no_DHEA_outliers)

# Model with Box-Cox transformed anxiety_mean_w05
model_boxcox <- lm(boxcox_anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_left_OFC_co_dev_time_corrected + mcc_sex_c, 
                   data = primary_analyses_set_no_DHEA_outliers)

# Check diagnostics for the original model
residuals_original <- model_original$residuals
shapiro.test(residuals_original)
par(mfrow = c(1, 2))  
hist(residuals_original, breaks = 20, main = "Original Residuals", xlab = "Residuals")
qqnorm(residuals_original)
qqline(residuals_original)

# Check diagnostics for the Box-Cox transformed model
residuals_boxcox <- model_boxcox$residuals
shapiro.test(residuals_boxcox)  
hist(residuals_boxcox, breaks = 20, main = "Box-Cox Transformed Residuals", xlab = "Residuals")
qqnorm(residuals_boxcox)
qqline(residuals_boxcox)

# Plot residuals vs fitted values for original model
plot(model_original$fitted.values, residuals_original, 
     main = "Residuals vs Fitted (Original)", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Plot residuals vs fitted values for Box-Cox transformed model
plot(model_boxcox$fitted.values, residuals_boxcox, 
     main = "Residuals vs Fitted (Box-Cox Transformed)", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Compare model fit for the original model
summary(model_original)

# Compare model fit for the Box-Cox transformed model
summary(model_boxcox)

# Compare AIC/BIC for both models
AIC(model_original, model_boxcox)
BIC(model_original, model_boxcox)

# Compare coefficients for both models
coef(model_original)
coef(model_boxcox)

##boxcox improves the quality of the data, as does sqrt 


###structural covariance 

# Create a list of variables to plot
variables <- c("left_Amygdala_left_OFC_co_dev_time_corrected",
               "left_Amygdala_right_OFC_co_dev_time_corrected",
               "right_Amygdala_left_OFC_co_dev_time_corrected",
               "right_Amygdala_right_OFC_co_dev_time_corrected")

# Melt the data to make it long format for easier plotting
library(reshape2)
data_long <- melt(primary_analyses_set_no_DHEA_outliers, id.vars = NULL, measure.vars = variables)



# Create boxplot
ggplot(data_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Structural Covariance Measures", y = "Value", 
       title = "Boxplots of Structural Covariance Measures") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




###Checking the Process Macro assumpttions of the data: independence of observations, linearity  of relationships among variables, error values are homoscedastic, multicollinearity does not exist among multiple IVs, error values are normally distributed
##Building the multiple regression models
# Model 1: Mediator - left_Amygdala_left_OFC_co_development_combined
model_1 <- lm(boxcox_anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_left_OFC_co_dev_time_corrected + mcc_sex_c, 
              data = primary_analyses_set_no_DHEA_outliers)

# Model 2: Mediator - left_Amygdala_right_OFC_co_development_combined
model_2 <- lm(boxcox_anxiety_mean_w05 ~ dhea_corrected + left_Amygdala_right_OFC_co_dev_time_corrected + mcc_sex_c, 
              data = primary_analyses_set_no_DHEA_outliers)

# Model 3: Mediator - right_Amygdala_left_OFC_co_development_combined
model_3 <- lm(boxcox_anxiety_mean_w05 ~ dhea_corrected + right_Amygdala_left_OFC_co_dev_time_corrected + mcc_sex_c, 
              data = primary_analyses_set_no_DHEA_outliers)

# Model 4: Mediator - right_Amygdala_right_OFC_co_development_combined
model_4 <- lm(boxcox_anxiety_mean_w05 ~ dhea_corrected + right_Amygdala_right_OFC_co_dev_time_corrected + mcc_sex_c, 
              data = primary_analyses_set_no_DHEA_outliers)


# Residual diagnostics #pretty much normally distributed! ((they  weren't before I sqrt transformed anxiety scores))
residuals <- model_1$residuals
hist(residuals, breaks = 20, main = "Residuals", xlab = "Residuals")
shapiro.test(residuals)

residuals <- model_2$residuals
hist(residuals, breaks = 20, main = "Residuals", xlab = "Residuals")
shapiro.test(residuals)

residuals <- model_3$residuals
hist(residuals, breaks = 20, main = "Residuals", xlab = "Residuals")
shapiro.test(residuals)

residuals <- model_4$residuals
hist(residuals, breaks = 20, main = "Residuals", xlab = "Residuals")
shapiro.test(residuals)

##Q-Q plots look fine (deviations not too big, multiple regression is robust against non-severe violations of normality according to Hayes 2018)
# Q-Q plot for model_1
residuals_1 <- model_1$residuals
qqnorm(residuals_1, main = "Q-Q Plot for Model 1 Residuals")
qqline(residuals_1, col = "red")

# Q-Q plot for model_2
residuals_2 <- model_2$residuals
qqnorm(residuals_2, main = "Q-Q Plot for Model 2 Residuals")
qqline(residuals_2, col = "red")

# Q-Q plot for model_3
residuals_3 <- model_3$residuals
qqnorm(residuals_3, main = "Q-Q Plot for Model 3 Residuals")
qqline(residuals_3, col = "red")

# Q-Q plot for model_4
residuals_4 <- model_4$residuals
qqnorm(residuals_4, main = "Q-Q Plot for Model 4 Residuals")
qqline(residuals_4, col = "red")


##independence --> was completely okay in this data
# Durbin-Watson test for model_1 ##fine!
dwtest(model_1)

# Model 2 ##fine!
dwtest(model_2) 

# Model 3 ##fine!!
dwtest(model_3)

# Model 4 ##fine! 
dwtest(model_4)


##Checking linearity 
# Model 1: Residuals vs Fitted values
plot(model_1$fitted.values, model_1$residuals, 
     main = "Residuals vs Fitted Values (Model 1)", xlab = "Fitted Values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

# Model 2: Residuals vs Fitted values
plot(model_2$fitted.values, model_2$residuals, 
     main = "Residuals vs Fitted Values (Model 2)", xlab = "Fitted Values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

# Model 3: Residuals vs Fitted values
plot(model_3$fitted.values, model_3$residuals, 
     main = "Residuals vs Fitted Values (Model 3)", xlab = "Fitted Values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")

# Model 4: Residuals vs Fitted values
plot(model_4$fitted.values, model_4$residuals, 
     main = "Residuals vs Fitted Values (Model 4)", xlab = "Fitted Values", ylab = "Residuals", pch = 16)
abline(h = 0, col = "red")



##Blue line looks okay-ish (has to be a straight line), purple (residuals) looks off (we knew that)
# Model 1: Partial Residual Plot 
crPlots(model_1)

#Model 2: Partial residuals plot
crPlots(model_2)

#model 3: partial residuals plot
crPlots(model_3)

#model 4: partial residuals plot
crPlots(model_4)


# Model 1: Scatter plot of DHEA and Anxiety
plot(primary_analyses_set_no_DHEA_outliers$dhea_corrected, primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05,
     main = "Scatter plot of DHEA vs Anxiety (Model 1)", xlab = "DHEA", ylab = "Anxiety", pch = 16)

# Model 1: Scatter plot of left_Amygdala_left_OFC_co_development_combined and Anxiety
plot(primary_analyses_set_no_DHEA_outliers$left_Amygdala_left_OFC_co_dev_time_corrected, 
     primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05,
     main = "Scatter plot of left_Amygdala_left_OFC_co_development_combined vs Anxiety (Model 1)", 
     xlab = "left_Amygdala_left_OFC_co_dev_time_corrected", ylab = "Anxiety", pch = 16)

# Model 1: Scatter plot of Sex and Anxiety
plot(primary_analyses_set_no_DHEA_outliers$mcc_sex_c, primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05,
     main = "Scatter plot of Sex vs Anxiety (Model 1)", xlab = "Sex", ylab = "Anxiety", pch = 16)

###conclusion for linearity: I am pretty sure the data does not meet the tion for linearity....


##Multicollinearity ##No problem!!!
vif(model_1)
vif(model_2)
vif(model_3)
vif(model_4)


##homoscedasticity of error #also fine! no heteroscedasticity

library(lmtest)
bptest(model_1)
bptest(model_2)
bptest(model_3)
bptest(model_4)



# Create boxplot
ggplot(data_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Structural Covariance Measures", y = "Value", 
       title = "Boxplots of Structural Covariance Measures") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###Histograms and qq-plots for the structural covariance measures ##Look pretty okay
hist_left_Amygdala_left_OFC <- ggplot(primary_analyses_set, aes(x = left_Amygdala_left_OFC_co_dev_time_corrected)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of left_Amygdala_left_OFC_co__dev_time_corrected", x = "left_Amygdala_left_OFC_co__dev_time_corrected", y = "Frequency") +
  theme_minimal()
print(hist_left_Amygdala_left_OFC)

qq_left_Amygdala_left_OFC <- ggplot(primary_analyses_set, aes(sample = left_Amygdala_left_OFC_co_dev_time_corrected)) +
  geom_qq() +
  geom_qq_line(color = "blue") +
  labs(title = "Q-Q plot of left_Amygdala_left_OFC_co__dev_time_corrected") +
  theme_minimal()
print(qq_left_Amygdala_left_OFC)


hist_left_Amygdala_right_OFC <- ggplot(primary_analyses_set, aes(x = left_Amygdala_right_OFC_co_dev_time_corrected)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of left_Amygdala_right_OFC_co_dev_time_corrected", x = "left_Amygdala_right_OFC_co_dev_time_corrected", y = "Frequency") +
  theme_minimal()
print(hist_left_Amygdala_right_OFC)


qq_left_Amygdala_right_OFC <- ggplot(primary_analyses_set, aes(sample = left_Amygdala_right_OFC_co_dev_time_corrected)) +
  geom_qq() +
  geom_qq_line(color = "blue") +
  labs(title = "Q-Q plot of left_Amygdala_right_OFC_co_dev_time_corrected") +
  theme_minimal()
print(qq_left_Amygdala_right_OFC)


hist_right_Amygdala_left_OFC <- ggplot(primary_analyses_set, aes(x = right_Amygdala_left_OFC_co_dev_time_corrected)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of right_Amygdala_left_OFC_co_dev_time_corrected", x = "right_Amygdala_left_OFC_co_development_combined", y = "Frequency") +
  theme_minimal()
print(hist_right_Amygdala_left_OFC)

qq_right_Amygdala_left_OFC <- ggplot(primary_analyses_set, aes(sample = right_Amygdala_left_OFC_co_development_combined)) +
  geom_qq() +
  geom_qq_line(color = "blue") +
  labs(title = "Q-Q plot of right_Amygdala_left_OFC_co_dev_time_corrected") +
  theme_minimal()
print(qq_right_Amygdala_left_OFC)


hist_right_Amygdala_right_OFC <- ggplot(primary_analyses_set, aes(x = right_Amygdala_right_OFC_co_development_combined)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of right_Amygdala_right_OFC_co_dev_time_corrected", x = "right_Amygdala_right_OFC_co_dev_time_corrected", y = "Frequency") +
  theme_minimal()
print(hist_right_Amygdala_right_OFC)


qq_right_Amygdala_right_OFC <- ggplot(primary_analyses_set, aes(sample = right_Amygdala_right_OFC_co_dev_time_corrected)) +
  geom_qq() +
  geom_qq_line(color = "blue") +
  labs(title = "Q-Q plot of right_Amygdala_right_OFC_co_development_combined") +
  theme_minimal()
print(qq_right_Amygdala_right_OFC)



# Filter the dataset for pds_total < 12
  pds_filtered <- primary_analyses_set_no_DHEA_outliers[primary_analyses_set_no_DHEA_outliers$pds_total < 12, ]

# Boxplot for pds_total (only values < 12)
ggplot(pds_filtered, aes(x = "", y = pds_total)) +
  geom_boxplot() +
  labs(title = "Boxplot of pds_total (values < 12)", y = "pds_total") +
  theme_minimal()


# Boxplot for lh_OFC_volume_normalized_change
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = lh_OFC_volume_normalized_change_time_corrected)) +
  geom_boxplot() +
  labs(title = "Boxplot of lh_OFC_volume_normalized_change", y = "lh_OFC_volume_normalized_change") +
  theme_minimal()

# Boxplot for rh_OFC_volume_normalized_change
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = rh_OFC_volume_normalized_change_time_corrected)) +
  geom_boxplot() +
  labs(title = "Boxplot of rh_OFC_volume_normalized_change", y = "rh_OFC_volume_normalized_change") +
  theme_minimal()

# Boxplot for Left.Amygdala_normalized_change
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = Left.Amygdala_normalized_change_time_corrected)) +
  geom_boxplot() +
  labs(title = "Boxplot of Left.Amygdala_normalized_change", y = "Left.Amygdala_normalized_change") +
  theme_minimal()

# Boxplot for Right.Amygdala_normalized_change
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = Right.Amygdala_normalized_change_time_corrected)) +
  geom_boxplot() +
  labs(title = "Boxplot of Right.Amygdala_normalized_change", y = "Right.Amygdala_normalized_change") +
  theme_minimal()



# Boxplot for anxiety_total_w05
ggplot(primary_analyses_set_no_DHEA_outliers, aes(x = "", y = boxcox_anxiety_mean_w05)) +
  geom_boxplot() +
  labs(title = "Boxplot of anxiety_total_w05", y = "anxiety_total_w05") +
  theme_minimal()


###RELIABILITIES

##My data testosterone
# Correlate mcc_w03_hair_dhea_level_c with mcc_w03_hair_testo_testosterone_c
cor_dhea_testo <- cor(
  primary_analyses_set$mcc_w03_hair_dhea_level_c[primary_analyses_set$mcc_w03_hair_testo_testosterone_c != 999],
  primary_analyses_set$mcc_w03_hair_testo_testosterone_c[primary_analyses_set$mcc_w03_hair_testo_testosterone_c != 999],
  use = "complete.obs"
)

cor_dhea_testo


# Plot the correlation
ggplot(primary_analyses_set[primary_analyses_set$mcc_w03_hair_testo_testosterone_c != 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = mcc_w03_hair_testo_testosterone_c)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  
  labs(x = "DHEA Level", y = "Testosterone Level", 
       title = "Correlation between DHEA and Testosterone") +
  theme_minimal()



##Full data testosterone
# Correlate DHEA with Testosterone in data_a
cor_dhea_testo_a <- cor(
  data_a$mcc_w03_hair_dhea_level_c[data_a$mcc_w03_hair_testo_testosterone_c != 999],
  data_a$mcc_w03_hair_testo_testosterone_c[data_a$mcc_w03_hair_testo_testosterone_c != 999],
  use = "complete.obs"
)

# Print the correlation result
print(cor_dhea_testo_a)

# Plot the correlation for data_a
ggplot(data_a[!is.na(data_a$mcc_w03_hair_dhea_level_c) & 
                !is.na(data_a$mcc_w03_hair_testo_testosterone_c) & 
                data_a$mcc_w03_hair_dhea_level_c != 999 &
                data_a$mcc_w03_hair_testo_testosterone_c != 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = mcc_w03_hair_testo_testosterone_c)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +  
  labs(x = "DHEA Level", y = "Testosterone Level", 
       title = "Correlation between DHEA and Testosterone in data_a") +
  theme_minimal()



##My data dhea with previous DHEA


# Correlate mcc_w03_hair_dhea_level_c with mcc_w02_hair_dhea_level_c
cor_dhea_prev_dhea <- cor(
  primary_analyses_set$mcc_w03_hair_dhea_level_c[primary_analyses_set$mcc_w02_hair_dhea_level_c != 999],
  primary_analyses_set$mcc_w02_hair_dhea_level_c[primary_analyses_set$mcc_w02_hair_dhea_level_c != 999],
  use = "complete.obs"
)

cor_dhea_prev_dhea


# Plot the correlation between DHEA levels in wave 2 and wave 3
ggplot(primary_analyses_set[primary_analyses_set$mcc_w02_hair_dhea_level_c != 999 & 
                                               primary_analyses_set$mcc_w03_hair_dhea_level_c != 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = mcc_w02_hair_dhea_level_c)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  # Linear regression line
  labs(x = "DHEA Level (Wave 3)", y = "DHEA Level (Wave 2)", 
       title = "Correlation between DHEA Levels at Wave 2 and Wave 3") +
  theme_minimal()





##full data
# Correlate mcc_w03_hair_dhea_level_c with mcc_w02_hair_dhea_level_c for data_a
# Calculate correlation while excluding values 999 and values above 999 in both waves
cor_dhea_prev_dhea_data_a <- cor(
  data_a$mcc_w03_hair_dhea_level_c[data_a$mcc_w02_hair_dhea_level_c != 999 & data_a$mcc_w02_hair_dhea_level_c <= 999 & 
                                     data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999],
  data_a$mcc_w02_hair_dhea_level_c[data_a$mcc_w02_hair_dhea_level_c != 999 & data_a$mcc_w02_hair_dhea_level_c <= 999 & 
                                     data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999],
  use = "complete.obs"
)

# View the correlation result
cor_dhea_prev_dhea_data_a

# Plot the correlation between DHEA levels in wave 2 and wave 3 for data_a, excluding invalid values in wave 2
ggplot(data_a[data_a$mcc_w02_hair_dhea_level_c != 999 & data_a$mcc_w02_hair_dhea_level_c <= 999 & 
                data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = mcc_w02_hair_dhea_level_c)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  # Linear regression line
  labs(x = "DHEA Level (Wave 3)", y = "DHEA Level (Wave 2)", 
       title = "Correlation between DHEA Levels at Wave 2 and Wave 3 for data_a") +
  theme_minimal()





###PDS total
# Correlate mcc_w03_hair_dhea_level_c with pds_total for the full dataset
# Exclude values 999 and values above 999 in mcc_w03_hair_dhea_level_c
# Exclude values above 12 in pds_total
cor_dhea_pds_full <- cor(
  primary_analyses_set$mcc_w03_hair_dhea_level_c[primary_analyses_set$pds_total <= 12 & primary_analyses_set$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set$mcc_w03_hair_dhea_level_c <= 999],
  primary_analyses_set$pds_total[primary_analyses_set$pds_total <= 12 & primary_analyses_set$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set$mcc_w03_hair_dhea_level_c <= 999],
  use = "complete.obs"
)

# View the correlation result for the full dataset
cor_dhea_pds_full

# Plot the correlation between DHEA levels and PDS total for the full dataset, excluding invalid values in both
ggplot(primary_analyses_set[primary_analyses_set$pds_total <= 12 & primary_analyses_set$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set$mcc_w03_hair_dhea_level_c <= 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = pds_total)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  # Linear regression line
  labs(x = "DHEA Level (Wave 3)", y = "PDS Total", 
       title = "Correlation between DHEA Levels (Wave 3) and PDS Total for Full Dataset") +
  theme_minimal()



##Zonder de outliers op DHEA...
cor_dhea_pds_full <- cor(
  primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$pds_total <= 12 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c <= 999],
  primary_analyses_set_no_DHEA_outliers$pds_total[primary_analyses_set_no_DHEA_outliers$pds_total <= 12 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c <= 999],
  use = "complete.obs"
)

cor_dhea_pds_full


ggplot(primary_analyses_set_no_DHEA_outliers[primary_analyses_set_no_DHEA_outliers$pds_total <= 12 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c != 999 & primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c <= 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = pds_total)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  # Linear regression line
  labs(x = "DHEA Level (Wave 3)", y = "PDS Total", 
       title = "Correlation between DHEA Levels (Wave 3) and PDS Total for Full Dataset") +
  theme_minimal()



# Correlate mcc_w03_hair_dhea_level_c with pds_total for data_a
# Exclude values 999 and values above 999 in mcc_w03_hair_dhea_level_c
# Exclude values above 12 in pds_total
cor_dhea_pds_data_a <- cor(
  data_a$mcc_w03_hair_dhea_level_c[data_a$pds_total <= 12 & data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999],
  data_a$pds_total[data_a$pds_total <= 12 & data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999],
  use = "complete.obs"
)

# View the correlation result for data_a
cor_dhea_pds_data_a

# Plot the correlation between DHEA levels and PDS total for data_a, excluding invalid values in both
ggplot(data_a[data_a$pds_total <= 12 & data_a$mcc_w03_hair_dhea_level_c != 999 & data_a$mcc_w03_hair_dhea_level_c <= 999, ], 
       aes(x = mcc_w03_hair_dhea_level_c, y = pds_total)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", color = "blue") +  # Linear regression line
  labs(x = "DHEA Level (Wave 3)", y = "PDS Total", 
       title = "Correlation between DHEA Levels (Wave 3) and PDS Total for data_a") +
  theme_minimal()





##Making descriptive statistics

# Continuous variables
cat("Descriptive Statistics for Age:\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for DHEA (unit):\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c, na.rm = TRUE), "\n")

# Filter out values greater than 12 for PDS
pds_filtered <- primary_analyses_set_no_DHEA_outliers$pds_total[primary_analyses_set_no_DHEA_outliers$pds_total <= 12]

# Calculate descriptive statistics for the filtered PDS variable
cat("\nDescriptive Statistics for PDS (Puberty) (filtered, ≤12):\n")
summary(pds_filtered)
cat("\nMean:", mean(pds_filtered, na.rm = TRUE), 
    "SD:", sd(pds_filtered, na.rm = TRUE), "\n")



cat("\nDescriptive Statistics for Anxiety:\n")
summary(primary_analyses_set_no_DHEA_outliers$anxiety_total_w05)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$anxiety_total_w05, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$anxiety_total_w05, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for LH OFC Change:\n")
summary(primary_analyses_set_no_DHEA_outliers$lh_OFC_volume_normalized_change_time_corrected)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$lh_OFC_volume_normalized_change_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$lh_OFC_volume_normalized_change_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for RH OFC Change:\n")
summary(primary_analyses_set_no_DHEA_outliers$rh_OFC_volume_normalized_change_time_corrected)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$rh_OFC_volume_normalized_change_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$rh_OFC_volume_normalized_change_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for Left Amygdala Change:\n")
summary(primary_analyses_set_no_DHEA_outliers$Left.Amygdala_normalized_change_time_corrected)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$Left.Amygdala_normalized_change_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$Left.Amygdala_normalized_change_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for Right Amygdala Change:\n")
summary(primary_analyses_set_no_DHEA_outliers$Right.Amygdala_normalized_change_time_corrected)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$Right.Amygdala_normalized_change_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$Right.Amygdala_normalized_change_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for left_Amygdala_left_OFC_co_development_combined:\n")
summary(primary_analyses_set_no_DHEA_outliers$left_Amygdala_left_OFC_co_dev_time_corrected)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$left_Amygdala_left_OFC_co_dev_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$left_Amygdala_left_OFC_co_dev_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for left_Amygdala_right_OFC_co_development_combined:\n")
summary(primary_analyses_set_no_DHEA_outliers$left_Amygdala_right_OFC_co_dev_time_corrected, na.rm = TRUE)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$left_Amygdala_right_OFC_co_dev_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$left_Amygdala_right_OFC_co_dev_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for right_Amygdala_left_OFC_co_development_combined:\n")
summary(primary_analyses_set_no_DHEA_outliers$right_Amygdala_left_OFC_co_dev_time_corrected, na.rm = TRUE)
cat("\nMean:", mean(primary_analyses_set_no_DHEA_outliers$right_Amygdala_left_OFC_co_dev_time_corrected, na.rm = TRUE),
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$right_Amygdala_left_OFC_co_dev_time_corrected, na.rm = TRUE), "\n")

cat("\nDescriptive Statistics for right_Amygdala_right_OFC_co_development_combined:\n")
summary(primary_analyses_set$right_Amygdala_right_OFC_co_dev_time_corrected)
cat("\nMean:", mean(primary_analyses_set$right_Amygdala_right_OFC_co_dev_time_corrected, na.rm = TRUE), 
    "SD:", sd(primary_analyses_set$right_Amygdala_right_OFC_co_dev_time_corrected, na.rm = TRUE), "\n")

# Categorical variable
cat("\nFrequency Table for Sex:\n")
table(primary_analyses_set_no_DHEA_outliers$mcc_sex_c)


##Feedback: do hormone values separately for males and females:


# Boys (mcc_sex_c == 0)
cat("\n--- Descriptive Statistics for Boys (mcc_sex_c == 0) ---\n")

cat("\nAge:\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0])
cat("Mean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0], na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0], na.rm = TRUE), "\n")

cat("\nDHEA (unit):\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0])
cat("Mean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0], na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0], na.rm = TRUE), "\n")

cat("\nPDS (filtered, ≤12):\n")
pds_filtered_boys <- primary_analyses_set_no_DHEA_outliers$pds_total[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 0 & primary_analyses_set_no_DHEA_outliers$pds_total <= 12]
summary(pds_filtered_boys)
cat("Mean:", mean(pds_filtered_boys, na.rm = TRUE), 
    "SD:", sd(pds_filtered_boys, na.rm = TRUE), "\n")

# Girls (mcc_sex_c == 1)
cat("\n--- Descriptive Statistics for Girls (mcc_sex_c == 1) ---\n")

cat("\nAge:\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1])
cat("Mean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1], na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_age_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1], na.rm = TRUE), "\n")

cat("\nDHEA (unit):\n")
summary(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1])
cat("Mean:", mean(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1], na.rm = TRUE), 
    "SD:", sd(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1], na.rm = TRUE), "\n")

cat("\nPDS (filtered, ≤12):\n")
pds_filtered_girls <- primary_analyses_set_no_DHEA_outliers$pds_total[primary_analyses_set_no_DHEA_outliers$mcc_sex_c == 1 & primary_analyses_set_no_DHEA_outliers$pds_total <= 12]
summary(pds_filtered_girls)
cat("Mean:", mean(pds_filtered_girls, na.rm = TRUE), 
    "SD:", sd(pds_filtered_girls, na.rm = TRUE), "\n")




variables <- c("dhea_corrected", "mcc_w03_hair_dhea_level_c", 
               "left_Amygdala_left_OFC_co_dev_time_corrected", 
               "left_Amygdala_right_OFC_co_dev_time_corrected",
               "right_Amygdala_left_OFC_co_dev_time_corrected", 
               "right_Amygdala_right_OFC_co_dev_time_corrected", 
               "anxiety_total_w05", "boxcox_anxiety_mean_w05", "anxiety_mean_w05",
               "lh_OFC_volume_normalized_change_time_corrected", 
               "rh_OFC_volume_normalized_change_time_corrected", 
               "Left.Amygdala_normalized_change_time_corrected", 
               "Right.Amygdala_normalized_change_time_corrected")

# Check if all columns exist in the dataset
missing_vars <- setdiff(variables, colnames(primary_analyses_set_no_DHEA_outliers))
if (length(missing_vars) > 0) {
  print(paste("These variables are missing:", paste(missing_vars, collapse = ", ")))
}

# Create a subset of the data with the variables of interest
subset_data <- primary_analyses_set_no_DHEA_outliers[, variables, drop = FALSE]



#Correlations:

# Calculate the correlation matrix
cor_matrix <- cor(subset_data, use = "complete.obs")

# Melt the correlation matrix for ggplot2
cor_matrix_melted <- melt(cor_matrix)

# Create a clearer heatmap with diverging color scale and correlation values
ggplot(cor_matrix_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, 
                       name = "Correlation", limits = c(-1, 1)) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +  # Add correlation values
  theme_minimal() +
  labs(title = "Correlation Matrix", x = "Variables", y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1))


# Load necessary library
library(Hmisc)

# Calculate correlation matrix and p-values
rcorr_results <- rcorr(as.matrix(subset_data), type = "pearson")

# Extract correlation coefficients and p-values
cor_matrix <- rcorr_results$r
p_values <- rcorr_results$P

# Create a dataframe for visualization of significance
library(reshape2)
cor_matrix_melted <- melt(cor_matrix)
p_values_melted <- melt(p_values)

# Combine correlations and p-values
cor_data <- cbind(cor_matrix_melted, p_value = p_values_melted$value)
colnames(cor_data) <- c("Var1", "Var2", "Correlation", "P_Value")

# Add significance levels
cor_data$Significance <- cut(cor_data$P_Value,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", "ns"),
                             right = FALSE)

# Plot the correlation matrix with significance
library(ggplot2)

ggplot(cor_data, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, 
                       name = "Correlation", limits = c(-1, 1)) +
  geom_text(aes(label = paste0(round(Correlation, 2), "\n", Significance)), 
            color = "black", size = 3) +  # Add correlation values with significance
  theme_minimal() +
  labs(title = "Correlation Matrix with Significance", x = "Variables", y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1))





library(corrplot)

# Create the correlation matrix
cor_matrix <- cor(subset_data, use = "complete.obs")

# Visualize it with a circular plot
corrplot(cor_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.cex = 0.8, 
         addCoef.col = "black", number.cex = 0.7,
         title = "Correlation Matrix - Circular Plot")






##Running the moderated-mediation with process model 7 

#Model 1: mediator: left_amygdaka_left_OFC
process(data = primary_analyses_set_no_DHEA_outliers,
        y = "boxcox_anxiety_mean_w05",
        x = "dhea_corrected",
        m = "left_Amygdala_left_OFC_co_dev_time_corrected", 
        w = "mcc_sex_c",
        model = 7)

#Model 2: mediator: left_amygdala_right OFC
process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "left_Amygdala_right_OFC_co_dev_time_corrected", 
        w = "mcc_sex_c", 
        model = 7)


#Model 3: mediator: right_amygdala_left_OFC
process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "right_Amygdala_left_OFC_co_dev_time_corrected", 
        w = "mcc_sex_c", 
        model = 7)

#Model 4: mediator: right amygdala right ofc

process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "right_Amygdala_right_OFC_co_dev_time_corrected", 
        w = "mcc_sex_c", 
        model = 7)


##Making the mediations
# Model 1: Mediator: left_amygdala_left_OFC
process(data = primary_analyses_set_no_DHEA_outliers,
        y = "boxcox_anxiety_mean_w05",
        x = "dhea_corrected",
        m = "left_Amygdala_left_OFC_co_dev_time_corrected", 
        model = 4)  # Model 4 is for simple mediation (no moderator)

# Model 2: Mediator: left_amygdala_right_OFC
process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "left_Amygdala_right_OFC_co_dev_time_corrected", 
        model = 4)  # Model 4 is for simple mediation (no moderator)

# Model 3: Mediator: right_amygdala_left_OFC
process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "right_Amygdala_left_OFC_co_dev_time_corrected", 
        model = 4)  # Model 4 is for simple mediation (no moderator)

# Model 4: Mediator: right_amygdala_right_OFC
process(data = primary_analyses_set_no_DHEA_outliers, 
        y = "boxcox_anxiety_mean_w05", 
        x = "dhea_corrected", 
        m = "right_Amygdala_right_OFC_co_dev_time_corrected", 
        model = 4)  # Model 4 is for simple mediation (no moderator)




##PSychometrics 

# Extract the anxiety items from your dataset
anxiety_items <- primary_analyses_set_no_DHEA_outliers[, c(
  "mcc_w05_bsi_anxiety_q4_pp",
  "mcc_w05_bsi_anxiety_q9_pp",
  "mcc_w05_bsi_anxiety_q14_pp",
  "mcc_w05_bsi_anxiety_q18_pp",
  "mcc_w05_bsi_anxiety_q20_pp"
)]

# Compute Cronbach's alpha
alpha_result <- psych::alpha(anxiety_items)

# View the results
print(alpha_result)


# Select the columns for DHEA at waves 2 and 3
data_for_icc <- MCC_data_JolleElbers_encrypted[, c("mcc_w02_hair_dhea_level_c", "mcc_w03_hair_dhea_level_c")]

# Remove rows with missing data
data_for_icc <- na.omit(data_for_icc)

data_for_icc$mcc_w02_hair_dhea_level_c <- as.numeric(data_for_icc$mcc_w02_hair_dhea_level_c)
data_for_icc$mcc_w03_hair_dhea_level_c <- as.numeric(data_for_icc$mcc_w03_hair_dhea_level_c)

icc_result <- ICC(data_for_icc)
print(icc_result)

# Select the columns for DHEA at waves 2 and 3
data_for_icc <- MCC_data_JolleElbers_encrypted[, c("mcc_w02_hair_dhea_level_c", "mcc_w03_hair_dhea_level_c")]

# Remove rows with missing data or 999 values
data_for_icc <- na.omit(data_for_icc)
data_for_icc <- data_for_icc[data_for_icc$mcc_w02_hair_dhea_level_c != 999 & 
                               data_for_icc$mcc_w03_hair_dhea_level_c != 999, ]

# Ensure values are numeric
data_for_icc$mcc_w02_hair_dhea_level_c <- as.numeric(data_for_icc$mcc_w02_hair_dhea_level_c)
data_for_icc$mcc_w03_hair_dhea_level_c <- as.numeric(data_for_icc$mcc_w03_hair_dhea_level_c)

# Calculate ICC
icc_result <- ICC(data_for_icc)
print(icc_result)



# Load the required package
library(psych)

# Select the three items for the PDS total
pds_items <- data_a[, c("mcc_w03_pds_q1_c", "mcc_w03_pds_q2_c", "mcc_w03_pds_q3_c")]

# Calculate Cronbach's alpha for the three items
alpha_result <- alpha(pds_items)

# Print the result
print(alpha_result)



##Posthoc test mod


# Prepare the dataset as before (ensure no missing values)
data_for_analysis <- primary_analyses_set_no_DHEA_outliers[!is.na(primary_analyses_set_no_DHEA_outliers$mcc_w03_hair_dhea_level_c) &
                                                             !is.na(primary_analyses_set_no_DHEA_outliers$boxcox_anxiety_mean_w05) &
                                                             !is.na(primary_analyses_set_no_DHEA_outliers$mcc_sex_c), ]

# Fit the moderation model using PROCESS (Model 1 for moderation)
results <- process(data = data_for_analysis,
                   y = "boxcox_anxiety_mean_w05",  # Dependent variable
                   x = "mcc_w03_hair_dhea_level_c", # Independent variable
                   w = "mcc_sex_c",  # Moderator variable (sex)
                   model = 1,  # Model 1 is the simple moderation model
                   boot = 5000)  # Bootstrapping with 5000 iterations

# View the results
summary(results)





####Chwecking final descriptives for last version:
## Filter for pds_total under 12 for both boys and girls
primary_analyses_set_no_DHEA_outliers_filtered <- primary_analyses_set_no_DHEA_outliers[primary_analyses_set_no_DHEA_outliers$pds_total < 12,]

# Boys (mcc_sex_c == 0)
cat("\n--- Descriptive Statistics for Boys (mcc_sex_c == 0) ---\n")

# Descriptive for all variables (excluding pds_total > 12)
vars_boys <- c("mcc_w03_age_c", "mcc_w03_hair_dhea_level_c", "pds_total", "anxiety_mean_w05",
               "left_Amygdala_left_OFC_co_dev_time_corrected", "left_Amygdala_right_OFC_co_dev_time_corrected",
               "right_Amygdala_left_OFC_co_dev_time_corrected", "right_Amygdala_right_OFC_co_dev_time_corrected")

for (var in vars_boys) {
  cat("\n", var, ":\n")
  var_data <- primary_analyses_set_no_DHEA_outliers_filtered[[var]][primary_analyses_set_no_DHEA_outliers_filtered$mcc_sex_c == 0]
  summary(var_data)
  cat("Mean:", mean(var_data, na.rm = TRUE), 
      "SD:", sd(var_data, na.rm = TRUE), "\n")
}

# Girls (mcc_sex_c == 1)
cat("\n--- Descriptive Statistics for Girls (mcc_sex_c == 1) ---\n")

# Descriptive for all variables (excluding pds_total > 12)
vars_girls <- c("mcc_w03_age_c", "mcc_w03_hair_dhea_level_c", "pds_total", "anxiety_mean_w05",
                "left_Amygdala_left_OFC_co_dev_time_corrected", "left_Amygdala_right_OFC_co_dev_time_corrected",
                "right_Amygdala_left_OFC_co_dev_time_corrected", "right_Amygdala_right_OFC_co_dev_time_corrected")

for (var in vars_girls) {
  cat("\n", var, ":\n")
  var_data <- primary_analyses_set_no_DHEA_outliers_filtered[[var]][primary_analyses_set_no_DHEA_outliers_filtered$mcc_sex_c == 1]
  summary(var_data)
  cat("Mean:", mean(var_data, na.rm = TRUE), 
      "SD:", sd(var_data, na.rm = TRUE), "\n")
}

# T-tests for all variables
cat("\n--- T-tests for Boys vs Girls ---\n")

for (var in vars_boys) {
  cat("\nT-test for", var, ":\n")
  t_test_result <- t.test(primary_analyses_set_no_DHEA_outliers_filtered[[var]] ~ primary_analyses_set_no_DHEA_outliers_filtered$mcc_sex_c)
  print(t_test_result)
  cat("\n")
}




# Descriptive Statistics for the entire sample (both boys and girls together)
cat("\n--- Descriptive Statistics for the Entire Sample ---\n")

# Descriptive for all variables (excluding pds_total > 12)
vars_all <- c("mcc_w03_age_c", "mcc_w03_hair_dhea_level_c", "pds_total", "anxiety_mean_w05",
              "left_Amygdala_left_OFC_co_dev_time_corrected", "left_Amygdala_right_OFC_co_dev_time_corrected",
              "right_Amygdala_left_OFC_co_dev_time_corrected", "right_Amygdala_right_OFC_co_dev_time_corrected")

# Calculate and display the descriptive statistics for the entire sample
for (var in vars_all) {
  print(paste("\n", var, ":"))
  var_data <- primary_analyses_set_no_DHEA_outliers_filtered[[var]]
  print(summary(var_data))
  print(paste("Mean:", mean(var_data, na.rm = TRUE), 
              "SD:", sd(var_data, na.rm = TRUE)))
}


