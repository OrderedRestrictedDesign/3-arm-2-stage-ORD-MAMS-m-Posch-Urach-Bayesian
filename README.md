# Code to reproduce the results on Section 4 of "A Bayesian multi-arm multi-stage clinical trial design incorporating information about treatment ordering"

# The file "functionsORDMAMSUPBAYESIAN_3arms2stages.R" contains 3 functions used for the MAMS(m) design:
1. bounds_3armJstagem that finds the critical bounds for a 3-arm or 4-arm J-stage MAMS design
2. boundaries_3armJstagem that calls the function at 1. and find the sample size to reject all for a 3-arm 2-stage or 3-arm 3-stage MAMS(m) design
3. simulAll that simulates a 3-arm J-stage MAMS(m) design

# 3 functions for the ORD design:
1. bounds_3armsORD_differentbounds that finds the critical bounds for a 3-arm 2-stage ORD design
2. boundariesSample_3armsORD_differentbounds that calls the function at 1. and find the sample size to reject all or at least one for a 3-arm J-stage ORD design
3. simul that simulates a 3-arm J-stage ORD design

# 1 functions for the Urach & Posch design:
1. boundpowerUrach3arm2stage that finds the critical bounds and the sample size for a 3-arm 2-stage design

# 2 functions for the Bayesian design:
1. sim that computes the FWER under the global and partial nulls and the power for a 3-arm 2-stage Bayesian design
2. integ_function that calls the function at 1. for each value of the grid

# The file "simulationsORDMAMSUPBAYESIAN_3arms2stages.R":
contains the call of the functions and perform simulations for different scenarios
