# Load Required Libraries
library(dplyr)
library(forecast)
library(tidyr)
library(knitr)
library(stats)
library(ggplot2)
library(boot)
library(mice)
library(MASS)
require(JointAI)
library(brms)
library(loo)
library(sf)
library(spdep)
library(sp)
library(sf)
require(sf)
library(tmap)
library(CARBayes)
library(coda)
# Load the dataset
lung_cancer_data <- read.csv("Lung Cancer.csv")
smoking_data <- read.csv("Smoking.csv")
manuf_data <- read.csv("Manufacturing.csv")
constr_data <- read.csv("Construction.csv")
emi_data <- read.csv("CO2.csv")
earning_data <- read.csv("Earnings.csv")
pop_data <- read.csv("Population.csv")
# Try Imputation Method using backwardcast by ARIMA (Times Series)
smoking_data_t <- t(smoking_data[-1])
smoking_ts <- ts(smoking_data_t, start = 2012, end = 2019, frequency = 1)
backcasts <- matrix(NA, nrow = 4, ncol = ncol(smoking_ts))
rownames(backcasts) <- 2008:2011
for(i in 1:ncol(smoking_ts)) {
  # Reverse the time series data for backcasting
  ts_data_reversed <- rev(smoking_ts[, i])
  
  
  
  

  
  
  
#################################################  
  # Adjust the start time for the reversed series
  start_year <- 2019 # Assuming your data goes up to 2019
  ts_data_reversed <- ts(ts_data_reversed, start = start_year, frequency = 1)
  # Fit an ARIMA model to the reversed data
  fit <- auto.arima(ts_data_reversed)
  # Forecast (which is effectively backcasting) the missing years
  bc <- forecast(fit, h = 4)
  # Store the backcasted values in reverse order to align with original time
  backcasts[, i] <- rev(bc$mean)
}
# Convert backcasts to a data frame for easier viewing/manipulation
backcast_df <- as.data.frame(backcasts)
# Initialize an empty data frame to store backcasted values
backcast_df <- data.frame(
  council_area = rep(smoking_data$Council.Areas, each = 4),
  year = rep(2008:2011, times = ncol(smoking_ts)),
  backcasted_value = as.vector(backcasts)
)
colnames(backcast_df) <- c('council_area','year','backcasted_value')
names(backcast_df) <- c('Council.Areas','Year','Smoking_Rate')
backcast_df$Year <- paste0("X", backcast_df$Year)
# Standardise the Earning Data
earning_data[, -1] <- earning_data[, -1] / 1000
# Remove the last 15 rows
emi_data_cleaned <- emi_data %>%
  slice(1:(n() - 15)) %>%
  filter(Year >= 2008 & Year <= 2017)
emi_data_selected <- emi_data_cleaned %>%
  dplyr::select(Council.Areas, Year, Emissions_per_km2 )
# Pivot lung_cancer_data
Smoking_long <- pivot_longer(smoking_data,
                             cols = -Council.Areas,
                             names_to = "Year",
                             values_to = "Smoking_Rate")
# We only need the data from 2008-2017
Smoking_1 <- Smoking_long %>%
  filter(Year >= 'X2008' & Year <= 'X2017')
merged_df <- merge(Smoking_1, backcast_df, by = c("Council.Areas",
                                                  "Year",
                                                  'Smoking_Rate'),
                   all = TRUE)
 
  
   
#########################################  
# Pivot lung_cancer_data 
lung_long <- pivot_longer(lung_cancer_data,
                          cols = -Council.Areas,
                          names_to = "Year",
                          values_to = "Lung_cancer_count")
# Pivot manuf_data
manuf_long <- pivot_longer(manuf_data,
                           cols = -Council.Areas,
                           names_to = "Year",
                           values_to = "Manufacturing_Employment_count")
# Pivot constr_data
constr_data_long <- pivot_longer(constr_data,
                                 cols = -Council.Areas,
                                 names_to = "Year",
                                 values_to = "Construction_Employment_count")
# Pivot earning_data
earning_data_long <- pivot_longer(earning_data,
                                  cols = -Council.Areas,
                                  names_to = "Year",
                                  values_to = "Earnings")
# Pivot pop_data
pop <- pivot_longer(pop_data,
                    cols = -Council.Areas,
                    names_to = "Year",
                    values_to = "Pop_count")
pop_data <- pop %>%
  filter(Year >= 'X2008' & Year <= 'X2017')
# Merge the data frames by Council.Areas and Year
merged_data <- merge(constr_data_long, pop_data, by = c("Council.Areas", "Year"))
# Perform the rates calculation on the merged data frame
merged_data$Construction_Rate <- merged_data$Construction_Employment_count *
  1000 / merged_data$Pop_count *100
merged_data_1 <- merge(manuf_long, pop_data, by = c("Council.Areas", "Year"))
merged_data_1$Manuf_Rate <- merged_data_1$Manufacturing_Employment_count *
  1000 / merged_data$Pop_count *100
# Final Merged Data
merged_reg_data <- cbind(lung_long[,1:3], Construction = merged_data[,5],
                         Manufacturing = merged_data_1[,5],
                         Emissions = emi_data_selected[,3],
                         
                         Earnings = earning_data_long[,3],
                         Smoking_Rate = merged_df[,3],
                         Pop = pop_data[,3])
# Assuming your dataframe is named 'df' and the year column is named 'Year'
merged_reg_data$Smoking_Rate[merged_reg_data$Year %in% c("X2008",
                                                         "X2009",
                                                         "X2010",
                                                         "X2011")] <- NA
data <- merged_reg_data[, !colnames(merged_reg_data) %in% c('Council.Areas',
                                                            'Year')]
md.pattern(merged_reg_data)
# Perform MICE imputation
imputed_data <- mice(merged_reg_data, m=40, method='norm', maxit=10, seed=500)
plot(imputed_data)
stripplot(imputed_data, data = Smoking_Rate ~ .imp)
# Randomly choose the imputation data for the bayesian analysis
completed_data <- complete(imputed_data, action = 1)
# Fit Poisson and Negative Binomial Models on Imputed Datasets
Poi_model <- lapply(1:40, function(i) {
  completed_data <- complete(imputed_data, action = i)
  glm(Lung_cancer_count ~ Manufacturing + Construction + Emissions
      + Earnings + Smoking_Rate,
      data = completed_data, family = "poisson",
      offset = log(completed_data$Pop))
})
### Pool Model Results According to Rubin's Rules
pooled_results <- pool(Poi_model)
summary(pooled_results)
### Calculate and Print Mean AIC for Each Set of Models
aic_model <- sapply(Poi_model, AIC)
mean_aic_model <- mean(aic_model)
print(paste("Mean AIC for Poisson models with offset:", mean_aic_model))
### Sensitivity Analysis
# Sensitivity analysis by varying imputation parameters
imputed_data_varied <- mice(data, m=20, method='norm', maxit=5, seed=501)
# Comparing regression coefficients across varied imputation parameters
models_varied <- list()
for (i in 1:20) {
  completed_data_varied <- complete(imputed_data_varied, action = i)
  models_varied[[i]] <- glm(Lung_cancer_count ~ Manufacturing + Construction
                            + Emissions + Earnings + Smoking_Rate,
                            data = completed_data_varied, family = "poisson",
                            offset=log(Pop_count))
}
pooled_results_varied <- pool(models_varied)
summary(pooled_results_varied)                         
                         
                         
#####################################                         
### Bootstrap Analysis on Imputed Dataset
# Fitting Poisson model to bootstrap samples
fit_model <- function(data, indices) {
  boot_data <- data[indices, ]
  model <- glm(Lung_cancer_count ~ Manufacturing + Construction + Emissions
               + Earnings + Smoking_Rate,
               data = boot_data, family = "poisson", offset=log(Pop_count))
  coef(model)
}
num_coefficients <- 6
bootstrap_results <- matrix(NA, nrow = 40, ncol = num_coefficients)
# Iterates over each imputed dataset to perform bootstrap analysis.
for (i in 1:40) {
  # Retrieves the i-th imputed dataset.
  completed_data <- complete(imputed_data, action = i)
  # Performs bootstrap analysis on the imputed dataset.
  boot_res <- boot(data = completed_data, statistic = fit_model, R = 1000)
  # Calculates the mean of coefficients across bootstrap samples
  bootstrap_results[i, ] <- apply(boot_res$t, 2, mean)
}
# Computes the overall mean of the bootstrap results to estimate the variability
overall_coef_sd <- colMeans(bootstrap_results)
# Outputs the overall variability (standard deviation) of each coefficient.
print(overall_coef_sd)
# Calculates the 95% confidence intervals for the coefficients
ci_95 <- apply(bootstrap_results,
               2, function(x) quantile(x, probs = c(0.025, 0.975)))
# Converts the confidence intervals into a more readable data frame format.
ci_95_df <- as.data.frame(t(ci_95))
# Assigns names to the columns
colnames(ci_95_df) <- c("2.5%", "97.5%")
# Assigns names to the rows based on the predictors and intercept.
rownames(ci_95_df) <- c("Intercept", "Manufacturing", "Construction",
                        "Emissions", "Earnings", "Smoking_Rate")
# Displays the confidence intervals for each coefficient.
print(ci_95_df)
#########################################





# WAIC
# Convert the model to a Bayesian framework using brms
bayesian_model <- brm(
  formula = Lung_cancer_count ~ Manufacturing + Construction + Emissions
  + Earnings + Smoking_Rate + offset(log(Pop_count)),
  family = poisson(),
  data = completed_data,
  prior <- c(
    set_prior("normal(0, 10)", class = "b"),
    set_prior("normal(0, 10)", class = "Intercept")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE) # Ensure all parameters are saved
)
# Calculate WAIC with moment matching
waic_result <- loo(bayesian_model, save_psis = TRUE)
# Print WAIC result
print(waic_result)
WAIC_Poi <- waic_result$estimates["looic", "Estimate"]
### Spatial Analysis Section
# Ratio Data
shape <- read_sf(dsn = "/home/baoqizhang/Desktop/Projects/Shape")
shape <- st_read("/home/baoqizhang/Desktop/Projects/Shape")
# Unified council areas name for shape file and the lung cancer data
shape$local_auth <- gsub("Eilean Siar", "Na h-Eileanan Siar", shape$local_auth)
lung_cancer_data <- read.csv("Lung Cancer.csv")
pop_data <- read.csv("Population.csv")
calculate_ratio <- function(lung_cancer_data, pop_data) {
  # Extract the year columns from both data frames
  years <- colnames(lung_cancer_data)[-1]
  # Initialize an empty data frame to store the results
  result <- data.frame(council.areas = pop_data$Council.Areas)
  # Loop through each year and calculate the ratio
  for (year in years) {
    lung_cancer_col <- as.numeric(lung_cancer_data[, year])
    pop_col <- as.numeric(pop_data[, year])
    # Check for missing or non-numeric values
    if (any(is.na(lung_cancer_col)) || any(is.na(pop_col))) {
      warning(paste("Skipping year", year,
                    "due to missing or non-numeric values"))
      result[[year]] <- NA
      
      
#################################
    } else {
      # Calculate the ratio (lung cancer cases / population)
      ratio <- (lung_cancer_col / pop_col)*100
      result[[year]] <- ratio
    }
  }
  return(result)
}
ratio_data <- calculate_ratio(lung_cancer_data, pop_data)
# Calculate the average ratio across specified columns, ignoring NA values
ratio_data$mean_ratio = rowMeans(ratio_data[, 2:ncol(ratio_data)], na.rm = TRUE)
# Display the first few rows of the updated ratio dataset for verification
head(ratio_data)
lung_cancer_data <- read.csv("Lung Cancer.csv")
lung_cancer_data$mean_count = rowMeans(lung_cancer_data[, 2:ncol(lung_cancer_data)])
# Display the first few rows of the lung cancer data to ensure accuracy
head(lung_cancer_data)
# Merge ratio and lung cancer datasets based on council areas
comb <- merge(ratio_data, lung_cancer_data, by.x = "council.areas",
              by.y = "Council.Areas")
# Further merge the combined dataset with spatial boundary data
merged_data_1 <- merge(shape, comb, by.x = "local_auth", by.y = "council.areas")
# Converting shape file to the Spatial File
shape_sp <- as(merged_data_1, "Spatial")
# Identification and Analysis of Spatial Neighbors
# Conversion of the Spatial*DataFrame to an 'sf' object
shape_sf <- st_as_sf(shape_sp)
# Calculation of centroids for each spatial feature
centroids <- st_centroid(shape_sf)
# Extraction of x and y coordinates from the calculated centroids
centroids_coords <- st_coordinates(centroids)
# Inclusion of centroid coordinates back into the 'sf' object
# This enriches the spatial dataset with precise location data for each area
shape_sf$x <- centroids_coords[, 1]
shape_sf$y <- centroids_coords[, 2]
# Identification of spatial neighbors based on the Queen contiguity criterion
# This method determines adjacency based on shared borders or vertices
nb_q_sp <- poly2nb(shape_sf)
# Manual adjustments for specific spatial relationships
nb_q_sp[[20]] <- c(nb_q_sp[[20]], 4)
nb_q_sp[[4]] <- c(nb_q_sp[[4]], 20)
nb_q_sp[[27]] <- c(nb_q_sp[[27]], 1)
nb_q_sp[[1]] <- c(nb_q_sp[[1]], 27)
nb_q_sp[[27]] <- c(nb_q_sp[[27]], 23)
nb_q_sp[[23]] <- c(nb_q_sp[[23]], 27)
nb_q_sp[[1]] <- c(nb_q_sp[[1]], 23)
nb_q_sp[[23]] <- c(nb_q_sp[[23]], 1)
nb_q_sp[[16]] <- c(nb_q_sp[[16]], 23)
nb_q_sp[[23]] <- c(nb_q_sp[[23]], 16)
nb_q_sp[[16]] <- c(nb_q_sp[[16]], 20)
nb_q_sp[[20]] <- c(nb_q_sp[[20]], 16)
# Refinement of the neighbors list to exclude non-neighbors
nb_q_sp <- lapply(nb_q_sp, function(neighbors) neighbors[neighbors != 0])
# Summary of the neighbors list to validate the adjustments
summary(nb_q_sp)
# Construction of an adjacency matrix from the neighbors list
# The adjacency matrix represents the neighbor relationships between areas
n_neigh <- nrow(shape_sp)
adj_mat <- matrix(data = 0, nrow = n_neigh, ncol = n_neigh)
for(i in 1:n_neigh) {
  adj_mat[i, nb_q_sp[[i]]] <- 1
}
# Conversion of the neighbors list into a 'nb' class object
nb_q_sp <- lapply(nb_q_sp, as.integer)
class(nb_q_sp) <- "nb"
attr(nb_q_sp, "region.id") <- as.character(1:32)
# Creation of a 'listw' object with a zero policy
lw_q_B <- nb2listw(nb_q_sp, style="B", zero.policy=TRUE)
# Perform Shapiro-Wilk normality test
shapiro_result <- shapiro.test(shape_sf$mean_ratio)
# Print the results
print(shapiro_result)
# Execution of the Global Moran's I Test using the mean ratio data
moran_result <- moran.test(shape_sp$mean_ratio, lw_q_B)
# Presentation of the Moran's I Test results
print(moran_result)
# Moran's I scatterplot to illustrate the spatial autocorrelation
moran.plot(shape_sp$mean_ratio, lw_q_B, return_df = TRUE)
### Visualization of Neighborhood Structures
# Calculate the number of neighbors for each spatial unit
num_neighbors <- sapply(nb_q_sp, length)
# Generate a color palette to visually differentiate spatial units
num_breaks <- length(unique(num_neighbors))
colors <- colorRampPalette(c("white", "cornflowerblue"))(num_breaks)
# Map each unique neighbor count to a specific color index
color_mapping <- setNames(seq_len(num_breaks), sort(unique(num_neighbors)))
# Assign colors to each spatial unit according to its number of neighbors
area_colors <- colors[color_mapping[as.character(num_neighbors)]]
# Plot the spatial framework
plot(st_geometry(shape_sf), col = area_colors, border = "black")
# Create legend labels
legend_labels <- paste(names(color_mapping), "neighbors")
legend_colors <- colors
# Display a legend on the plot to interpret the color scheme
legend("topright", legend = legend_labels, fill = legend_colors, cex = 0.7)
# Illustrate the neighborhood relationships by connecting them
for (i in 1:length(nb_q_sp)) {
  # Iterate through each spatial unit to identify its neighbors
  neighbors <- nb_q_sp[[i]]
  for (j in neighbors) {
    # Draw lines connecting each unit to its neighbors, emphasizing the
    # neighborhood structure
    lines(rbind(centroids_coords[i,], centroids_coords[j,]), col="black")
  }
}
# Define a semi-transparent color for highlighting centroids
transparent_red <- rgb(1, 0, 0, alpha = 0.5)
# Centroids represent the geometric centers of spatial units,
# providing a focal point for visualizing neighborhood connections
points(centroids_coords[,1], centroids_coords[,2], pch = 21, col = "black",
       bg = transparent_red, cex = 1.5)
### Local Moran's I
# Calculate Local Moran's I
local_moran <- localmoran(shape_sf$mean_ratio, lw_q_B)
print(local_moran)
tmap_mode("plot")
shape_sp_1 <- st_as_sf(shape_sp)
# Add the Local Moran's I statistics to the sf object
shape_sp_1$lmI <- local_moran[, "Ii"] # local Moran's I
shape_sp_1$lmZ <- local_moran[, "Z.Ii"] # z-scores
shape_sp_1$lmp <- local_moran[, "Pr(z != E(Ii))"] # p-values
# Map 1: Lung Cancer Mean Ratio
p1 <- tm_shape(shape_sp_1) +
  tm_polygons(col = "mean_ratio", title = "Lung Cancer Rate",
              style = "quantile") +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 0.9)
p1
# Map 2: Local Moran's I
p2 <- tm_shape(shape_sp_1) +
  tm_polygons(col = "lmI", title = "Local Moran's I",
              style = "quantile") +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 0.9)
p2
# Map 3: Z-score
p3 <- tm_shape(shape_sp_1) +
  tm_polygons(col = "lmZ", title = "Z-score",
              # Use qnorm(0.975) for two-sided test
              breaks = c(-Inf, qnorm(0.95), Inf),
              # Red for significant positive autocorrelation
              palette = c("red", "white")) +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 0.9)
p3
# Map 4: p-value
p4 <- tm_shape(shape_sp_1) +
  tm_polygons(col = "lmp", title = "p-value",
              # Use 0.025 and 0.975 for two-sided test
              breaks = c(-Inf, 0.05, Inf),
              # Red for p-value less than 0.05
              palette = c("red", "white")) +
  tm_layout(legend.outside = TRUE,
            legend.text.size = 0.9)
p4
# LISA Cluster Plot
mp <- moran.plot(as.vector(scale(shape_sp_1$mean_ratio)), lw_q_B)
shape_sp_1$lmp <- local_moran[, 5]
shape_sp_1$quadrant <- NA
# high-high
shape_sp_1[(mp$x >= 0 & mp$wx >= 0) & (shape_sp_1$lmp <= 0.05), "quadrant"]<- 1
# low-low
shape_sp_1[(mp$x <= 0 & mp$wx <= 0) & (shape_sp_1$lmp <= 0.05), "quadrant"]<- 2
# high-low
shape_sp_1[(mp$x >= 0 & mp$wx <= 0) & (shape_sp_1$lmp <= 0.05), "quadrant"]<- 3
# low-high
shape_sp_1[(mp$x <= 0 & mp$wx >= 0) & (shape_sp_1$lmp <= 0.05), "quadrant"]<- 4
# non-significant
shape_sp_1[(shape_sp_1$lmp > 0.05), "quadrant"] <- 5
tm_shape(shape_sp_1) + tm_fill(col = "quadrant", title = "",
                               breaks = c(1, 2, 3, 4, 5, 6),
                               palette = c("red", "blue", "lightpink", "skyblue2",
                                           "white"),
                               labels = c("High-High", "Low-Low", "High-Low",
                                          "Low-High", "Non-significant")) +
  tm_legend(text.size = 1) + tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE, title = "Clusters") +
  tm_layout(legend.outside = TRUE)
# Define the necessary variables
weights <- nb2listw(nb_q_sp, style="B", zero.policy=TRUE)
neighbors_list <- as(weights, "listw")$neighbours
W_matrix <- nb2mat(nb_q_sp, style="B", zero.policy=TRUE)
# MCMC parameters
M.burnin <- 10000 # Number of burn-in iterations (discarded)
M <- 100000 # Number of iterations retained
# use 2016 model to predict 2017
# car model for 2016 since need to predict for 2017
combined_data_2016 <- filter(completed_data, Year == "X2016")
# BYM
set.seed(444) # For reproducability
MCMC_bym_all <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate
  + Manufacturing +
    Construction + Emissions + Earnings,
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_bym_all$summary.results)
# Additional results on model fit criteria
MCMC_bym_all$modelfit
# model check
beta_bym_all <- MCMC_bym_all$samples$beta
gelman_results_bym_all <- gelman.diag(beta_bym_all)
print(gelman_results_bym_all) #1 converges
plot(beta_bym_all)
autocorr.plot(beta_bym_all)
gelman.plot(beta_bym_all)
# Check posterior predict model
psi <- MCMC_bym_all$samples$psi
samples_beta_bym_all <- do.call(rbind, lapply(MCMC_bym_all$samples$beta,
                                              function(x) x))
samples_psi_bym_all <- do.call(rbind, lapply(MCMC_bym_all$samples$psi,
                                             function(x) x))
Nsamples = 1000
postpred_Y <- array(NA,dim=c(32,Nsamples))
postpred_mean <- array(NA,dim=c(32,Nsamples))
for (i_s in 1:Nsamples) {
  # randomly chose a set of samples after the combining samples
  intercept_s <- samples_beta_bym_all[s, 1]
  beta_smoking_s <- samples_beta_bym_all[s, 2]
  beta_construction_s <- samples_beta_bym_all[s, 3]
  beta_earnings_s <- samples_beta_bym_all[s, 4]
  beta_manufacturing_s <- samples_beta_bym_all[s, 5]
  beta_emissions_s <- samples_beta_bym_all[s,6]
  psi_means_s <- samples_psi_bym_all[s, ]
  # calculate the average for the prediction
  postpred_mean[, i_s] <- with(combined_data_2016,
                               exp(intercept_s +
                                     beta_smoking_s * Smoking_Rate +
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings +
                                     beta_emissions_s * Emissions + psi_means_s +
                                     offset(log(Pop_count))
                               ))
  # generating prediction counts
  postpred_Y[, i_s] <- rpois(32, lambda = postpred_mean[, i_s])
}
# Compare actual and prediction
hist(colMeans(postpred_Y))
abline(v = mean(combined_data_2016$Lung_cancer_count), col = "red")
hist(apply(postpred_Y, 2, sd))
abline(v = sd(combined_data_2016$Lung_cancer_count), col = "red")
# predict
car_2017 <- filter(completed_data, Year == "X2017")
# Extract posterior means
intercept_bym <- MCMC_bym_all$summary.results["(Intercept)", "Mean"]
beta_smoking_bym <- MCMC_bym_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_bym <- MCMC_bym_all$summary.results["Construction", "Mean"]
beta_earnings_bym <- MCMC_bym_all$summary.results["Earnings", "Mean"]
beta_manufacturing_bym <- MCMC_bym_all$summary.results["Manufacturing", "Mean"]
beta_emissions_bym <- MCMC_bym_all$summary.results["Emissions", "Mean"]
psi_means_bym = colMeans(samples_psi_bym_all)
# Predicting lung cancer cases for 2017
car_2017$predicted_lung_cancer_bym <- with(car_2017,
                                           exp(intercept_bym +
                                                 beta_smoking_bym * Smoking_Rate +
                                                 beta_construction_bym * Construction +
                                                 beta_manufacturing_bym * Manufacturing +
                                                 beta_earnings_bym * Earnings +
                                                 beta_emissions_bym * Emissions + psi_means_bym +
                                                 offset(log(Pop_count))
                                           )
)
# Add a column to compare the predicted cases to actual cases
car_2017$comparison_b <- with(car_2017, predicted_lung_cancer_bym - Lung_cancer_count)
mean(car_2017$comparison_bˆ2) # 489.4742 BYM(1 year)
# Create a scatter plot
ggplot(car_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_bym)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count",
       title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()
# LEROUX
set.seed(444)
MCMC_leroux_all <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate
  + Manufacturing + Construction + Emissions + Earnings,
  data = combined_data_2016,
  family = "poisson",
  W = W_matrix,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_leroux_all$summary.results)
# Additional results on model fit criteria
MCMC_leroux_all$modelfit
# model check
beta_le_all <- MCMC_leroux_all$samples$beta
gelman_results_le_all <- gelman.diag(beta_le_all)
print(gelman_results_le_all) #1 converge
plot(beta_le_all)
autocorr.plot(beta_le_all)
gelman.plot(beta_le_all)
# Check posterior predict model
phi <- MCMC_leroux_all$samples$phi
samples_beta_le_all <- do.call(rbind, lapply(MCMC_leroux_all$samples$beta,
                                             function(x) x))
samples_phi_le_all <- do.call(rbind, lapply(MCMC_leroux_all$samples$phi,
                                            function(x) x))
Nsamples = 1000
postpred_Y <- array(NA,dim=c(32,Nsamples))
postpred_mean <- array(NA,dim=c(32,Nsamples))
for (i_s in 1:Nsamples) {
  # Choose randomly a set of samples
  s = 10*i_s
  intercept_s_le <- samples_beta_le_all[s, 1]
  beta_smoking_s_le <- samples_beta_le_all[s, 2]
  beta_construction_s_le <- samples_beta_le_all[s, 3]
  beta_earnings_s_le <- samples_beta_le_all[s, 4]
  beta_manufacturing_s_le <- samples_beta_le_all[s, 5]
  beta_emissions_s_le <- samples_beta_le_all[s, 6]
  phi_means_s_le <- samples_phi_le_all[s, ]
  # Calculate the prediction average
  postpred_mean[, i_s] <- with(combined_data_2016,
                               exp(intercept_s_le +
                                     beta_smoking_s_le * Smoking_Rate +
                                     beta_construction_s_le * Construction +
                                     beta_manufacturing_s_le * Manufacturing +
                                     beta_earnings_s_le * Earnings +
                                     beta_emissions_s_le * Emissions +
                                     phi_means_s_le +
                                     offset(log(Pop_count))
                               ))
  # Generating prediction counts
  postpred_Y[, i_s] <- rpois(32, lambda = postpred_mean[, i_s])
}
# Compare actual and prediction
hist(colMeans(postpred_Y))
abline(v = mean(combined_data_2016$Lung_cancer_count), col = "red")
hist(apply(postpred_Y, 2, sd))
abline(v = sd(combined_data_2016$Lung_cancer_count), col = "red")
# predict for 2017
# Extract posterior means
intercept_le <- MCMC_leroux_all$summary.results["(Intercept)", "Mean"]
beta_smoking_le <- MCMC_leroux_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_le <- MCMC_leroux_all$summary.results["Construction", "Mean"]
beta_earnings_le <- MCMC_leroux_all$summary.results["Earnings", "Mean"]
beta_manufacturing_le <- MCMC_leroux_all$summary.results["Manufacturing", "Mean"]
beta_emissions_le <- MCMC_leroux_all$summary.results["Emissions", "Mean"]
phi_means_le <- colMeans(samples_phi_le_all)
# Predicting lung cancer cases for 2017
car_2017$predicted_lung_cancer_le <- with(car_2017,
                                          exp(intercept_le +
                                                beta_smoking_le * Smoking_Rate +
                                                beta_construction_le * Construction +
                                                beta_manufacturing_le * Manufacturing +
                                                beta_earnings_le * Earnings +
                                                beta_emissions_le * Emissions +
                                                phi_means_le +
                                                offset(log(Pop_count))))
# Add a column to compare the predicted cases to actual cases
car_2017$comparison_le <- with(car_2017, predicted_lung_cancer_le - Lung_cancer_count)
mean(car_2017$comparison_leˆ2) # 498.3557 Leroux(1 year)
# Create a scatter plot
ggplot(car_2017, aes(x = Lung_cancer_count, y = predicted_lung_cancer_le)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count",
       title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()
# 9 Year predict 2017
W_big_2 <- matrix(0, nrow = 288, ncol = 288)
# Replace diagonal with W_matrix
for (i in 0:8) {
  rows <- (1:32) + i * 32
  cols <- (1:32) + i * 32
  W_big_2[rows, cols] <- W_matrix
}
library(dplyr)
filtered_data <- completed_data %>%
  filter(!grepl("X2017", Year))
completed_data_car <-completed_data[order(completed_data$Year), ]
filtered_data_car <- completed_data_car[completed_data_car$Year != "X2017", ]
M.burnin <- 10000 # Number of burn-in iterations (discarded)
M <- 100000 # Number of iterations retained
# BYM
set.seed(444) # For reproducability
MCMC_big_bym_all <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) + Smoking_Rate
  + Manufacturing +
    Construction + Emissions + Earnings,
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_big_bym_all$summary.results)
# Additional results on model fit criteria
MCMC_big_bym_all$modelfit
# model check
beta_bym_all_big <- MCMC_big_bym_all$samples$beta
gelman_results_bym_all_big <- gelman.diag(beta_bym_all_big)
print(gelman_results_bym_all_big) # 1.01
plot(beta_bym_all_big)
autocorr.plot(beta_bym_all_big)
gelman.plot(beta_bym_all_big)
# Check posterior predict model
psi <- MCMC_big_bym_all$samples$psi
samples_beta_bym_all_big <- do.call(rbind, lapply(MCMC_big_bym_all$samples$beta,
                                                  function(x) x))
samples_psi_bym_all_big <- do.call(rbind, lapply(MCMC_big_bym_all$samples$psi,
                                                 function(x) x))
Nsamples = 1000
postpred_Y <- array(NA,dim=c(288,Nsamples))
postpred_mean <- array(NA,dim=c(288,Nsamples))
for (i_s in 1:Nsamples) {
  # Randomly select a set of parameter samples
  s = 10*i_s
  intercept_s <- samples_beta_bym_all_big[s, 1]
  beta_smoking_s <- samples_beta_bym_all_big[s, 2]
  beta_construction_s <- samples_beta_bym_all_big[s, 3]
  beta_earnings_s <- samples_beta_bym_all_big[s, 4]
  beta_manufacturing_s <- samples_beta_bym_all_big[s, 5]
  beta_emissions_s <- samples_beta_bym_all_big[s,6]
  psi_means_s <- samples_psi_bym_all_big[s, ]
  # Calculate the mean value of the prediction
  postpred_mean[, i_s] <- with(filtered_data_car,
                               exp(intercept_s +
                                     beta_smoking_s * Smoking_Rate +
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings +
                                     beta_emissions_s * Emissions + psi_means_s +
                                     offset(log(Pop_count))
                               ))
  # Generate predictive counts
  postpred_Y[, i_s] <- rpois(288, lambda = postpred_mean[, i_s])
}
# Comparison of predicted and actual observed data
hist(colMeans(postpred_Y))
abline(v = mean(filtered_data_car$Lung_cancer_count), col = "red")
hist(apply(postpred_Y, 2, sd))
abline(v = sd(filtered_data_car$Lung_cancer_count), col = "red")
# predict!
car_big_2017 <- filter(completed_data, Year == "X2017")
# Extract posterior means
intercept_big <- MCMC_big_bym_all$summary.results["(Intercept)", "Mean"]
beta_smoking_big <- MCMC_big_bym_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_big <- MCMC_big_bym_all$summary.results["Construction", "Mean"]
beta_earnings_big <- MCMC_big_bym_all$summary.results["Earnings", "Mean"]
beta_manufacturing_big <- MCMC_big_bym_all$summary.results["Manufacturing", "Mean"]
beta_emissions_big <- MCMC_big_bym_all$summary.results["Emissions", "Mean"]
psi_all_chains <- do.call(rbind, MCMC_big_bym_all$samples$psi)
psi_means <- colMeans(psi_all_chains)
psi_means_big_chain <- matrix(psi_means, nrow = 32, ncol = 9, byrow = FALSE)
psi_row_averages <- rowMeans(psi_means_big_chain)
# Predicting lung cancer cases for 2017
car_big_2017$predicted_lung_cancer_bym_big <- with(car_big_2017,
                                                   exp(intercept_big +
                                                         beta_smoking_big * Smoking_Rate +
                                                         beta_construction_big * Construction +
                                                         beta_manufacturing_big * Manufacturing +
                                                         beta_earnings_big * Earnings +
                                                         beta_emissions_big * Emissions +
                                                         psi_row_averages +
                                                         offset(log(Pop_count))))
# Add a column to compare the predicted cases to actual cases
car_big_2017$comparison_big_B <- with(car_big_2017,
                                      predicted_lung_cancer_bym_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_Bˆ2) # 671.8069 BYM multi_year
# Create a scatter plot
ggplot(car_big_2017, aes(x = Lung_cancer_count,
                         y = predicted_lung_cancer_bym_big)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Actual Lung Cancer Count", y = "Predicted Lung Cancer Count",
       title = "Actual vs Predicted Lung Cancer Cases") +
  theme_minimal()
# LEROUX
set.seed(444) # For reproducability
MCMC_big_le_all <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) +
    Smoking_Rate + Manufacturing +
    Construction + Emissions + Earnings,
  data = filtered_data_car,
  family = "poisson",
  W = W_big_2,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_big_le_all$summary.results)
# Additional results on model fit criteria
MCMC_big_le_all$modelfit
# Check posterior predict model
phi <- MCMC_big_le_all$samples$phi
samples_beta_le_all_big <- do.call(rbind,
                                   lapply(MCMC_big_le_all$samples$beta, function(x) x))
samples_phi_bym_all_big <- do.call(rbind,
                                   lapply(MCMC_big_le_all$samples$phi, function(x) x))
Nsamples = 1000
postpred_Y <- array(NA,dim=c(288,Nsamples))
postpred_mean <- array(NA,dim=c(288,Nsamples))
for (i_s in 1:Nsamples) {
  s = 10*i_s
  intercept_s <- samples_beta_le_all_big[s, 1]
  beta_smoking_s <- samples_beta_le_all_big[s, 2]
  beta_construction_s <- samples_beta_le_all_big[s, 3]
  beta_earnings_s <- samples_beta_le_all_big[s, 4]
  beta_manufacturing_s <- samples_beta_le_all_big[s, 5]
  beta_emissions_s <- samples_beta_le_all_big[s,6]
  phi_means_s <- samples_phi_bym_all_big[s, ]
  postpred_mean[, i_s] <- with(filtered_data_car,
                               exp(intercept_s +
                                     beta_smoking_s * Smoking_Rate +
                                     beta_construction_s * Construction +
                                     beta_manufacturing_s * Manufacturing +
                                     beta_earnings_s * Earnings +
                                     beta_emissions_s * Emissions + phi_means_s +
                                     offset(log(Pop_count))
                               ))
  postpred_Y[, i_s] <- rpois(288, lambda = postpred_mean[, i_s])
}
# Extract posterior means
intercept_big <- MCMC_big_le_all$summary.results["(Intercept)", "Mean"]
beta_smoking_big <- MCMC_big_le_all$summary.results["Smoking_Rate", "Mean"]
beta_construction_big <- MCMC_big_le_all$summary.results["Construction", "Mean"]
beta_earnings_big <- MCMC_big_le_all$summary.results["Earnings", "Mean"]
beta_manufacturing_big <- MCMC_big_le_all$summary.results["Manufacturing", "Mean"]
beta_emissions_big <- MCMC_big_le_all$summary.results["Emissions", "Mean"]
phi_all_chains <- do.call(rbind, MCMC_big_le_all$samples$phi)
phi_means <- colMeans(phi_all_chains)
phi_means_big_chain <- matrix(phi_means, nrow = 32, ncol = 9, byrow = FALSE)
phi_row_averages <- rowMeans(phi_means_big_chain)
# Predicting lung cancer cases for 2017
car_big_2017$predicted_lung_cancer_le_big <- with(car_big_2017,
                                                  exp(intercept_big +
                                                        beta_smoking_big * Smoking_Rate +
                                                        beta_construction_big * Construction +
                                                        beta_manufacturing_big * Manufacturing +
                                                        beta_earnings_big * Earnings +
                                                        beta_emissions_big * Emissions + phi_row_averages +
                                                        offset(log(Pop_count))
                                                  )
)
# Add a column to compare the predicted cases to actual cases
car_big_2017$comparison_big_L <- with(car_big_2017,
                                      predicted_lung_cancer_le_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_Lˆ2) # 666.2226 Leroux(multi-year)
poisson_model_4 <- glm(Lung_cancer_count ~ Manufacturing + Construction
                       + Emissions + Earnings + Smoking_Rate,
                       offset=offset(log(Pop_count)),
                       data = completed_data,
                       family = "poisson")
summary(poisson_model_4)s
# prepare the model for 2008-2016 since we need to predict for 2017
filtered_data <- completed_data %>%
  filter(completed_data$Year >= "X2008" & completed_data$Year <= "X2016")
# final model for glm in poisson
poisson_model_4 <- glm(Lung_cancer_count ~ Manufacturing + Construction
                       + Emissions + Earnings + Smoking_Rate,
                       offset=offset(log(Pop_count)),
                       data = filtered_data,
                       family = "poisson")
summary(poisson_model_4)
plot(poisson_model_4)
# predict for 2017
y_2017 <- completed_data %>% filter(Year == "X2017")
# Predicting lung cancer cases in 2017
predicted_counts <- predict(poisson_model_4, y_2017, type = "response")
# Add the predicted counts to your y_2017 dataset
y_2017$predicted_lung_cancer = predicted_counts
# Compute absolute value
y_2017$absolute_error <- abs(y_2017$Lung_cancer_count - y_2017$predicted_lung_cancer)
# View the results
View(y_2017)
# calculate the MSE
squared_errors <- (y_2017$Lung_cancer_count - y_2017$predicted_lung_cancer)ˆ2
mse <- mean(squared_errors)
print(mse) ## 1311.372
# prepare the model for 2008-2016 since we need to predict for 2017
filtered_data <- completed_data %>%
  filter(completed_data$Year >= "X2008" & completed_data$Year <= "X2016")
#make a 320x320 matrix filled with 0
W_big_3 <- matrix(0, nrow = 320, ncol = 320)
# Replace diagonal with W_matrix
for (i in 0:9) {
  rows <- (1:32) + i * 32
  cols <- (1:32) + i * 32
  W_big_3[rows, cols] <- W_matrix
}
# BYM plus all covariates
set.seed(444) # For reproducability
MCMC_bym_allyear <- S.CARbym(
  formula = Lung_cancer_count ~ offset(log(Pop_count))
  + Smoking_Rate + Manufacturing +
    Construction + Emissions + Earnings,
  data = completed_data_car,
  family = "poisson",
  W = W_big_3,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_bym_allyear$summary.results)
# Additional results on model fit criteria
MCMC_bym_allyear$modelfit
WAIC_BYM <- MCMC_bym_allyear$modelfit["WAIC"]
# LEROUX plus all covariates
set.seed(444) # For reproducability
MCMC_leroux_allyear <- S.CARleroux(
  formula = Lung_cancer_count ~ offset(log(Pop_count)) +
    Smoking_Rate + Manufacturing +
    Construction + Emissions + Earnings,
  data = completed_data_car,
  family = "poisson",
  W = W_big_3,
  burnin = M.burnin,
  n.sample = M.burnin + M, # Total iterations
  n.chains = 3,
  n.cores = 3,
  verbose = FALSE
)
# Broad summary of results
print(MCMC_leroux_allyear$summary.results)
# Additional results on model fit criteria
MCMC_leroux_allyear$modelfit
WAIC_LEROUX <- MCMC_leroux_allyear$modelfit["WAIC"]
# RMSE calculated
# bym 1
car_2017$comparison_b <- with(car_2017,
                              predicted_lung_cancer_bym - Lung_cancer_count)
mean(car_2017$comparison_bˆ2) # 489.4742
rmse1 <- mean((car_2017$comparison_b/car_2017$Lung_cancer_count)ˆ2)
sqrt(rmse1)
# le 1
car_2017$comparison_le <- with(car_2017,
                               predicted_lung_cancer_le - Lung_cancer_count)
mean(car_2017$comparison_leˆ2) # 498.3557
rmse2 <- mean((car_2017$comparison_le/car_2017$Lung_cancer_count)ˆ2)
sqrt(rmse2)
# bym multi
car_big_2017$comparison_big_B <- with(car_big_2017,
                                      predicted_lung_cancer_bym_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_Bˆ2) # 671.8069
rmse3 <- mean((car_big_2017$comparison_big_B/car_big_2017$Lung_cancer_count)ˆ2)
sqrt(rmse3)
# le multi
car_big_2017$comparison_big_L <- with(car_big_2017,
                                      predicted_lung_cancer_le_big - Lung_cancer_count)
mean(car_big_2017$comparison_big_Lˆ2) # 666.2226
rmse4 <- mean((car_big_2017$comparison_big_L/car_big_2017$Lung_cancer_count)ˆ2)
sqrt(rmse4)
# poisson
squared_errors <- (y_2017$Lung_cancer_count - y_2017$predicted_lung_cancer)ˆ2
mse <- mean(squared_errors)
print(mse)
y_2017$error <- with(y_2017, Lung_cancer_count - predicted_lung_cancer)
rmse5 <- mean((y_2017$error/y_2017$Lung_cancer_count)ˆ2)
sqrt(rmse5)




  
  