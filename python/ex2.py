
import pypomf as pp
import numpy as np

# ******************************************************************************
# Example 2 -- Max likelihood parameter estimation
# More of a test than example. The simulations offer a way of
# evaluating the loglikelihood of parameters, which in turn can be
# used in simple max likelihood parameter estimation.

print("Example 2 {}".format(50*'*'))

# ******************************************************************************
# Generate a single sample of the process trajectory.  This serves as
# our artificial data for the parameter estimation.

# Sequence of times, first element needs to be zero.
ts = [0] + list(np.arange(0.1, 10.0, 1.0))

# Initial observation values
a_tilde_init = 0.05
r_init = 0.025

# Generate a single sample of the process in `mcres`.
mcsimu = pp.MCSimu(n_samples=1, dt=0.01)

# Perturb some of the parameters slightly from defaults
mcsimu.pars.a_bar *= 1.1
mcsimu.pars.kappa *= 1.25
mcsimu.pars.rho *= 1.15
mcres = mcsimu.run_step(ts[1:], (a_tilde_init, r_init))

# Print used parameters for reference
print("SAMPLE DATA PARAMS:")
mcres.pars.print_values()

# Extract the observables from `mcres`, and including the initial
# condition to be the first element of the list (this is required)
obs = [[a_tilde_init, r_init]] + [[v[0][0], v[0][1]] for v in mcres.snapshots]

# Construct a likelihood model instance that's been customised for the
# process we're studying here.
#
# NOTE: It may be important here to increase n_samples, to make sure
# the method converges. The max likelihood algorithm simply seeks a
# maximum of a sum of log pdf evaluations.  But since the pdf is in
# fact here a random variable, it slightly fluctuates from one
# evaluation to the next. If these fluctations are larger than the
# convergence set of the algorithm, the root finder will never reach
# its stopping condition (it will instead forever keep chasing the
# moving maximum).  Increasing the parameter n_samples will make the
# pdf less random, and thus improves convergence. The computational
# complexity of pdf evaluations should scale roughly linearly with
# n_samples, so multiplying it by ten makes the process about ten
# times slower.
mlm = pp.MCLikelihoodModel(obs, ts, dt=0.01, n_samples=1000)

# We can now do the MLE fit. Nelder-Mead algorithm used here to max
# the loglikelihood function, xtol and ftol are tolerance parameters
# for position and value of the optimal point, respectively.  We skip
# the parameter cov estimation (skip_hessian=True), since the Hessian
# of the loglikelihood appear singular (I think this is probably
# because we have at least one redundant parameter in the model: the
# equations depend only on the difference rho - rho_bar, not the
# actual values).
#
# NOTE: This may take a very long time, minutes or hours depending on
# parameters!
mlmres = mlm.fit(maxiter=1000, maxfun=1000, method='nm', xtol=1e-3, ftol=1e-3,
                 skip_hessian=True)

# Print summary.
print("MLE ESTIMATED PARAMETERS:")
print(mlmres.summary())
