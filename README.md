This page hosts various functions for accessing rblimp features.

1. ChiBarWaldTest.R performs a chi-bar test of random slopes in a multilevel model. The function currently supports models with two random slopes. The function is called as follows.
   
chibar_test_slopes(my_model, testvars = c("var1"))

2. PlotPosteriors.R creates density plots of posterior distributions from an rblimp model object. The function is called as follows, var "var1" is the name of a variable that appears in the model (e.g., an outcome).
   
plot_posteriors(my_model, var = "var1")
