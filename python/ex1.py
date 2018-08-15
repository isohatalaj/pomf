
import pypomf as pp

# ******************************************************************************
# Example 1 -- Running basic simulations
#
print("Example 1 {}".format(50*'*'))

# ******************************************************************************
# Define simulation parameters, and construct simulation object
#

# Probability densities are computed through simulation of the
# stochastic differential equations. The algorithm uses the simple
# Euler method with `dt` as its time-step length. 
dt = 0.01

# This parametre selects the number of sample runs to perform, i.e.
# how many samples of the process pdf to generate.
n_samples = 1000

# Construct a simulation object
mcsimu = pp.MCSimu(dt=dt, n_samples=n_samples)


# ******************************************************************************
# Set model parameters, and generate a sampling of the pdfs
#

# Upon construction, the `mcsimu` object has been set to somewhat
# arbitrary initial values. The parameters can now be assigned new
# values..
mcsimu.pars = pp.Pars(a_bar = 0.15,
                      gamma = 0.1,
                      sigma_c = 0.30,
                      sigma_c_bar = 0.30,
                      kappa = 0.0,
                      r_bar = 0.010,
                      rho = 0.040,
                      rho_bar = 0.030,
                      sin_theta = 0.1,
                      sin_theta_bar = 0.1,
                      a_tilde_init = 0.15,
                      r_init = 0.03,
                      N0_over_y0 = 1.0,
                      Ntot0_over_ytot0 = 1.0)

# ..or alternatively/additionally directly adjusted:
mcsimu.pars.gamma = 0.15
mcsimu.pars.kappa = -0.5

# We can print the parameter values to screen as follows:
mcsimu.pars.print_values()

# Let us say we're interested in the outcome distribution at time
# t_end:
t_end = 5.0


# ******************************************************************************
# Generate samples of the process trajectories, and construct an
# estimate of the pdf.
#

# When initial values for the unobserved variable are not given, its
# values are drawn from its equilibrium distribution.
res = mcsimu.run_step(t_end, None)

# The `res` object now contains all the relevant information for us to
# construct the pdf of the observables (a_tilde, r) at time t_end.


# ******************************************************************************
# Display results
print("\nResults:")

# Obtain the (a_tilde, r) probability density function
pdf = res.get_obs_pdf()

# Test evaluating it at the initial observation
# print("a_tilde = {}, r = {}, f = {}".format(a_tilde_init,
#                                             r_init,
#                                             pdf((a_tilde_init, r_init))))

# Plot the pdf of observables
res.plot_obs_pdf()

# Plot the pdf of the hidden variables, i.e. eta, delta
res.plot_hidden_pdf()
