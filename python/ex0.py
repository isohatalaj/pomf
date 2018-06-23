
import pypomf as pp

# ******************************************************************************
# Example 0 -- Simple sampling of process distribution
# This example demonstrates how to evaluate the process distribution
# functions in simple cases.
#
# Ex1 shows how to better or more conveniently do repeated evaluations
# of the pdfs, and how to plot the functions.
#
print("Example 0 {}".format(50*'*'))


# Evaluate the (a_tilde, r) joint distribution

# First three arguments are the time, a_tilde, r point at which to
# evaluate the pdf, next two are the initial, time zero (a_tilde, r)
# values. Rest are the model parameters.
f = pp.mckde_eval_pdf(10.0,
                      0.15, 0.03,
                      0.05, 0.02,
                      a_bar = 0.075, gamma = 0.15,
                      sigma_c = 0.025, sigma_c_bar = 0.025,
                      kappa = -0.0, r_bar = 0.010,
                      rho = 0.070, rho_bar = 0.050)

print("f = ", f)

# If the pdf needs to be evaluated for different (a_tilde, r), but
# fixed time and parameters, then it is very much advisable to first
# generate the pd function itself, and then evaluate it. This avoids
# needing the do the costly numerical sampling of the distribution,
# and the subsequent kernel density estimation of the pdf.

# The following gets the pdf at time 10.0 for initial (a_tilde = 0.05,
# r = 0.02), and the stated parameters.
pdf = pp.mckde_pdf(10.0, 0.05, 0.02,
                   a_bar = 0.075, gamma = 0.15,
                   sigma_c = 0.025, sigma_c_bar = 0.025,
                   kappa = -0.0, r_bar = 0.010,
                   rho = 0.070, rho_bar = 0.050)

print("f = ", pdf((0.15, 0.03)))


# The unconditional a_tilde pdf is also available. It has an
# analytic formula, and thus can be always be evaluated relatively
# cheaply.
#
# First argument is the time at which to eval the pdf, second is the
# a_tilde point where to evaluate, and the third is the time zero
# initial value of a_tilde.
f_a_tilde = pp.a_tilde_pdf(1.0, 0.15, 0.05,
                           a_bar = 0.075, gamma = 0.15,
                           sigma_c = 0.025, sigma_c_bar = 0.025,
                           kappa = -0.0, r_bar = 0.010,
                           rho = 0.070, rho_bar = 0.050)

print("f_a_tilde = ", f_a_tilde)
