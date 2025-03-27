This page hosts various functions for accessing rblimp features.

1. ChiBarWaldTest.R performs a chi-bar test of random slopes in a multilevel model. The function currently supports models with two random slopes. The function is called as follows.
   
chibar_test_slopes(my_model, testvars = c("var1"))
