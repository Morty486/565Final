###############################################################
# We have four distributions
# 1. Standard Normal Distribution
# 2. Weibull scale = 1, shape = 1, mode at the boundary
# 3. Weibull scale = 1, shape = 1.2, positively skewed.
# 4. Weibull scale = 1, shape = 3.6, nearly symmetric
# 5. Weibull scale = 1, shape = 25, negatively skewed shapes
###############################################################

# Load necessary library
library(GenOrd)
library(ggplot2)

# Define the Weibull density function
weibull_pdf <- function(y, gamma, delta) {
  (delta / gamma) * (y / gamma)^(delta - 1) * exp(-(y / gamma)^delta)
}

# Define parameter sets
parameters <- list(
  "Mode at the boundary" = c(1, 1),
  "Positively skewed" = c(1, 1.2),
  "Nearly symmetric" = c(1, 3.6),
  "Negatively skewed" = c(1, 25)
)

# Generate data for plotting
y <- seq(0, 3, length.out = 1000)
data <- data.frame(
  y = rep(y, times = length(parameters)),
  Density = unlist(lapply(parameters, function(p) weibull_pdf(y, p[1], p[2]))),
  Type = rep(names(parameters), each = length(y))
)

# Plot the Weibull distributions
ggplot(data, aes(x = y, y = Density, color = Type)) +
  geom_line() +
  labs(
    title = "Weibull Distributions for Different Parameter Sets",
    x = "y",
    y = "Density",
    color = "Parameter Set"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



# 1. For each variable, compute the first four moments and store them
result_gen <- function(distribution, shape = NULL, scale = NULL) {
  n <- 100000
  if (distribution == "normal") {
    data <- rnorm(n, mean = 0, sd = 1)
  } else if (distribution == "weibull") {
    if (is.null(shape) | is.null(scale)) {
      stop("For Weibull distribution, both 'shape' and 'scale' must be specified.")
    }
    data <- rweibull(n, shape = shape, scale = scale)  # Weibull distribution
  } else {
    stop("Choose 'normal' or 'weibull'.")
  }

  first_moment <- mean(data)         # E(X)
  second_moment <- mean(data^2)      # E(X^2)
  third_moment <- mean(data^3)       # E(X^3)
  fourth_moment <- mean(data^4)      # E(X^4)


# 2. For each variable, compute all 99 percentiles and store them.
percentiles <-  quantile(data, probs = seq(0.01, 1, by = 0.01))

  list(
    first_moment = first_moment,
    second_moment = second_moment,
    third_moment = third_moment,
    fourth_moment = fourth_moment,
    percentiles = percentiles,
    data = data
  )
}

set.seed(123)
result_normal <- result_gen("normal")
result_weibull1<- result_gen("weibull", shape = 1, scale = 1)
result_weibull2<- result_gen("weibull", shape = 1.6, scale = 1)
result_weibull3<- result_gen("weibull", shape = 3.6, scale = 1)
result_weibull4<- result_gen("weibull", shape = 25, scale = 1)

result_normal_2 <- result_gen("normal")
result_weibull1_2<- result_gen("weibull", shape = 1, scale = 1)
result_weibull2_2<- result_gen("weibull", shape = 1.6, scale = 1)
result_weibull3_2<- result_gen("weibull", shape = 3.6, scale = 1)
result_weibull4_2<- result_gen("weibull", shape = 25, scale = 1)



# 3. For each pair of variables, compute Demirtas and Hedeker bounds
all_result <- list(result_normal,result_weibull1,result_weibull2,
                   result_weibull3,result_weibull4)

all_result_2 <- list(result_normal_2,result_weibull1_2,result_weibull2_2,
                   result_weibull3_2,result_weibull4_2)

## The code here is modified based on the library GenCorSeqSort
Pair_Bounds <- function(lst,lst2){
  dim = length(lst)
  low_bdd = matrix(0, dim, dim)
  up_bdd = matrix(0, dim, dim)
  for (i in 1:(dim - 1)) {
    for (j in (i + 1):dim) {
      x = lst[[i]]$data
      y = lst[[j]]$data
      low_bdd[i, j] = stats::cor(x[order(x)], rev(y[order(y)]))
      up_bdd[i, j] = stats::cor(x[order(x)], y[order(y)])
    }
  }

  for (i in 1:dim){
    x = lst[[i]]$data
    y = lst2[[i]]$data
    low_bdd[i, i] = stats::cor(x[order(x)], rev(y[order(y)]))
    up_bdd[i, i] = stats::cor(x[order(x)], y[order(y)])
  }

  low_bdd[lower.tri(low_bdd)] <- t(low_bdd)[lower.tri(low_bdd)]
  up_bdd[lower.tri(up_bdd)] <- t(up_bdd)[lower.tri(up_bdd)]


  list(low_bdd = low_bdd, up_bdd = up_bdd)
}

DH_Bounds <- Pair_Bounds(all_result,all_result_2)


# 4. Specify correlation values that are within these bounds with an increment of 0.01.
specified_corr <- list()

for (i in 1:5) {
  for (j in i:5) {
    # Extract bounds
    lower <- DH_Bounds$low_bdd[i, j]
    upper <- DH_Bounds$up_bdd[i, j]
   specified_corr[[paste0("Pair_", i, "_", j)]] <- seq(lower+0.01,
                                                       upper-0.01, by = 0.01)

  }
}

# 5. For each correlation within each continuous pair (δspec), generate bivariate ordinal data with 100 categories that are equally likely
probs <- seq(1/100, 99/100, by = 1/100)
marginal <- list(probs, probs)  # Two variables
corrcheck(marginal)

# # This is only one example, try to generalize it
#
# # Example correlation values
# delta_spec <- specified_corr$Pair_1_1[1] # Replace with your value
# delta_max <- DH_Bounds$up_bdd[1]   # Replace with your computed max bound
# delta_min <- DH_Bounds$low_bdd[1]  # Replace with your computed min bound
#
# # Adjust correlation
# if (delta_spec > 0) {
#   input_corr <- delta_spec / delta_max
# } else {
#   input_corr <- -delta_spec / delta_min
# }
#
# # Create correlation matrix
# Sigma <- matrix(c(1, input_corr, input_corr, 1), nrow = 2)
#
# n <- 100
# data <- ordsample(n, marginal, Sigma)

## Function to generate bivariate ordinal data based on specified correlation matrix
generate_bivord_data <- function(index) {
all_data <- list()
num <- index
for (i in 1:5){
  for (j in i:5){
    delta_spec <- specified_corr[[paste0("Pair_", i, "_", j)]][num] ## Select a random index from range of plausible values
    delta_max <- DH_Bounds$up_bdd[i,j]
    delta_min <- DH_Bounds$low_bdd[i,j]
    if (delta_spec > 0) {
      input_corr <- delta_spec / delta_max
    } else {
      input_corr <- -delta_spec / delta_min
    }
    Sigma <- matrix(c(1, input_corr, input_corr, 1), nrow = 2)
    data <- ordsample(100000, marginal, Sigma)
    all_data[[paste0("Pair_", i, "_", j)]] <- data
  }
}
  return(all_data)
}


# 6. Revert to continuous equivalent

# ordinal_to_continuous <- function(data, percentiles) {
#   apply(data, 2, function(column) sapply(column, function(ordinal_value) {
#     percentiles[ordinal_value]
#   }))
# }
#
# # Apply the mapping
# continuous_data_nn <- ordinal_to_continuous(data, result_normal$percentiles)
#
# continuous_data_all <- lapply(all_data, function(data) {
#   ordinal_to_continuous(data, result_normal$percentiles)
# })


percentiles_list <- list(
  result_normal$percentiles,
  result_weibull1$percentiles,
  result_weibull2$percentiles,
  result_weibull3$percentiles,
  result_weibull4$percentiles
)


# Updated function to handle column-specific percentiles
ordinal_to_continuous_multi <- function(data, percentiles_list) {
  # Iterate over each column and apply the corresponding percentiles
  result <- sapply(seq_along(percentiles_list), function(col_index) {
    sapply(data[, col_index], function(ordinal_value) {
      percentiles_list[[col_index]][ordinal_value]
    })
  })
  return(result)
}

# Function to convert back to original scale using percentiles
convert_to_percentiles <- function(ord_data) {
# Initialize a list to store continuous data
cts_data <- list()
input_data <- ord_data

# Loop through all pairs
for (i in 1:5) {
  for (j in i:5) {
    # Extract ordinal data for the current pair
    ordinal_data <- input_data[[paste0("Pair_", i, "_", j)]]

    # Create a list of percentiles for the current pair
    percentiles_list_pair <- list(
      percentiles_list[[i]],  # First column uses percentiles from index i
      percentiles_list[[j]]   # Second column uses percentiles from index j
    )

    # Convert ordinal data to continuous data
    continuous_data <- ordinal_to_continuous_multi(ordinal_data, percentiles_list_pair)

    # Store the continuous data in cts_data
    cts_data[[paste0("Pair_", i, "_", j)]] <- continuous_data
  }
}
  return(cts_data)
}

# Check the results
names(cts_data)  # List of pairs
head(cts_data[["Pair_1_2"]])  # Check one dataset


# 7&8. Design simulation study to evaluate performance of mixed data generation for different marginal distributions
# by comparing ratios and differences between specified and empirical correlations


# In the previous code, I have already finished one specified correlation
# I randomly choose the 10th value in the specified correlation value.

# Next, you need to do it for others

## Run simulation study with 3 different correlation matrices
set.seed(123)
bivord_data1 <- generate_bivord_data(10)
cts_data1 <- convert_to_percentiles(bivord_data1)

set.seed(123)
bivord_data2 <- generate_bivord_data(33)
cts_data2 <- convert_to_percentiles(bivord_data2)

set.seed(123)
bivord_data3 <- generate_bivord_data(50)
cts_data3 <- convert_to_percentiles(bivord_data3)

## Compute empirical correlations for each dataset
empirical_corr_list1 <- lapply(cts_data1, cor)
specified_corr_list1 <- sapply(names(specified_corr), function(pair_name) {
  specified_corr[[pair_name]][10]
})

empirical_corr_list2 <- lapply(cts_data2, cor)
specified_corr_list2 <- sapply(names(specified_corr), function(pair_name) {
  specified_corr[[pair_name]][33]
})

empirical_corr_list3 <- lapply(cts_data3, cor)
specified_corr_list3 <- sapply(names(specified_corr), function(pair_name) {
  specified_corr[[pair_name]][50]
})

## Function to calculate differences and ratios
# compute_diff_ratio <- function(empirical_corr, specified_corr) {
#   difference <- empirical_corr[1, 2] - specified_corr
#   ratio <- empirical_corr[1, 2] / specified_corr
#   return(list(difference = difference, ratio = ratio))
# }
#
#
#
# # Extract results
# differences <- sapply(results, function(x) x$difference)
# ratios <- sapply(results, function(x) x$ratio)

## Function to compute differences and ratios for all pairs
compute_diff_ratios_all <- function(empirical_corr_list, specified_corr_vector) {
  # Use mapply to process each pair of empirical correlation and specified correlation
  results <- mapply(function(empirical_corr, specified_corr) {
    difference <- empirical_corr[1, 2] - specified_corr
    ratio <- empirical_corr[1, 2] / specified_corr
    return(list(difference = difference, ratio = ratio))
  }, empirical_corr_list, specified_corr_vector, SIMPLIFY = FALSE)

  return(results)
}

## Apply the function to empirical correlations and specified correlations
results1 <- compute_diff_ratios_all(empirical_corr_list1, specified_corr_list1)
results2 <- compute_diff_ratios_all(empirical_corr_list2, specified_corr_list2)
results3 <- compute_diff_ratios_all(empirical_corr_list3, specified_corr_list3)

## Extract differences and ratios
differences1 <- sapply(results1, function(x) x$difference)
ratios1 <- sapply(results1, function(x) x$ratio)

differences2 <- sapply(results2, function(x) x$difference)
ratios2 <- sapply(results2, function(x) x$ratio)

differences3 <- sapply(results3, function(x) x$difference)
ratios3 <- sapply(results3, function(x) x$ratio)


## Check results
print(differences1)
print(ratios1)

print(differences2)
print(ratios2)

print(differences3)
print(ratios3)


# Generate plots
## Function to visualize simulation performance
plot_corr <- function(diffs, ratios) {
diff_ratio_df <- data.frame(Pair = factor(paste0("Pair", 1:15),
                                          levels = paste0("Pair", 1:15)),
                            Difference = diffs,
                            Ratio = ratios)

diff_ratio_long <- reshape2::melt(diff_ratio_df, id.vars = "Pair", variable.name = "Metric", value.name = "Value")

p <- ggplot(diff_ratio_long, aes(x = Pair, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparison of Differences and Ratios between specified and empirical correlations",
       x = "Pair", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
  scale_fill_manual(values = c("Difference" = "blue", "Ratio" = "green")) +
  scale_y_continuous(limits = c(-0.2, 1.5)) +
  scale_x_discrete(expand = expansion(mult = 0.2))

return(p)

}


## Display graphs
plot_corr(differences1, ratios1)
plot_corr(differences2, ratios2)
plot_corr(differences3, ratios3)



