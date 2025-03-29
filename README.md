This page hosts various functions for accessing rblimp features. The functions can be accessed as follows.

source("https://raw.githubusercontent.com/craigenders/rblimp-adds/main/rblimp-functions.R")

The chibar_test function performs a chi-bar test of random slopes in a multilevel model with approximate mixture weights. The function currently supports models with two random slopes. The function is called as follows.
   
chibar_test(my_model, raneff = c("var1"))

The plot_posteriors function creates density plots of posterior distributions from an rblimp model object. The function is called as follows, var "var1" is the name of a variable that appears in the model (e.g., an outcome). Specifying "All" prints all estimated parameters.
   
plot_posteriors(my_model, var = "var1")
