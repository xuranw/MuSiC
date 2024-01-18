# Load necessary libraries
library(ggplot2)

# Using built-in dataset 'mtcars'
data(mtcars)

# Print the first few rows of the dataset
head(mtcars)

# Basic summary statistics
summary(mtcars)``

# Create a new column 'kpl' (kilometers per liter) for fuel efficiency
mtcars$kpl <- mtcars$mpg * 0.425144

# Simple plot: Miles per Gallon vs. Horsepower
ggplot(mtcars, aes(x=hp, y=mpg)) +
  geom_point() +
  ggtitle("Miles per Gallon vs. Horsepower") +
  xlab("Horsepower") +
  ylab("Miles per Gallon")

# Basic linear regression model
fit <- lm(mpg ~ wt + cyl + hp, data=mtcars)
summary(fit)
