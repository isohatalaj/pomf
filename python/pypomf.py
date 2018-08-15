
import itertools
import ctypes
import ctypes.util

import numpy as np

import statsmodels.api as sm
from statsmodels.base.model import GenericLikelihoodModel

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

# ******************************************************************************
# INTERFACE FOR THE UNDERLYING C LIBRARY
# See below for main Python API.

try:
    libpomf = ctypes.CDLL("libpomf.dll")
except:
    libpomf = ctypes.CDLL("./libpomf.so")

_r_star = libpomf.pomf_r_starx
_r_star.argtypes = [ctypes.c_double, ctypes.c_double,
                    ctypes.c_void_p]
_r_star.restype = ctypes.c_double

_eta_star = libpomf.pomf_eta_starx
_eta_star.argtypes = [ctypes.c_double, ctypes.c_double,
                      ctypes.c_void_p]
_eta_star.restype = ctypes.c_double

_delta_equil_stdev = libpomf.pomf_delta_equil_stdevx
_delta_equil_stdev.argtypes = [ctypes.c_void_p]
_delta_equil_stdev.restype = ctypes.c_double

_sigma_a_tilde_star = libpomf.pomf_sigma_a_tilde_starx
_sigma_a_tilde_star.argtypes = [ctypes.c_void_p]
_sigma_a_tilde_star.restype = ctypes.c_double

_mcsimu_make = libpomf.pomf_mcsimux_make
_mcsimu_make.argtypes = [ctypes.c_int]
_mcsimu_make.restype = ctypes.c_void_p

_mcsimu_free = libpomf.pomf_mcsimux_free
_mcsimu_free.argtypes = [ctypes.c_void_p]
_mcsimu_free.restype = None

_mcsimu_get_samples = libpomf.pomf_mcsimux_get_samples
_mcsimu_get_samples.argtypes = [ctypes.c_void_p]
_mcsimu_get_samples.restype = ctypes.POINTER(ctypes.c_double)

_mcsimu_init = libpomf.pomf_mcsimux_init
_mcsimu_init.argtypes = [ctypes.c_void_p,   # self
                         ctypes.c_void_p]   # px
_mcsimu_init.restype = ctypes.c_int

_mcsimu_run = libpomf.pomf_mcsimux_run
_mcsimu_run.argtypes = [ctypes.c_void_p,   # self
                        ctypes.c_double,   # t_end
                        ctypes.c_double]   # dt
_mcsimu_run.restype = ctypes.c_int

_mcsimu_init_wobs = libpomf.pomf_mcsimux_init_wobs
_mcsimu_init_wobs.argtypes = [ctypes.c_void_p,  # self
                              ctypes.c_double,  # r
                              ctypes.c_double,  # y
                              ctypes.c_double]  # ytot


# ******************************************************************************
# PYTHON API BEGINS HERE

# This flag controls whether the output `y` variable is
# interpreted as the output of FIRMS or the TOTAL output,
# which also contains the output of banks.  Setting this
# `True` means that output `y` is interpreted as the total
# output.
use_total_y = True

xdim = ctypes.c_int.in_dll(libpomf, "pomf_xdim").value
xwdim = ctypes.c_int.in_dll(libpomf, "pomf_xwdim").value
xsim_ia = ctypes.c_int.in_dll(libpomf, "pomf_xsim_ia").value
xsim_idelta = ctypes.c_int.in_dll(libpomf, "pomf_xsim_idelta").value
xsim_ieta = ctypes.c_int.in_dll(libpomf, "pomf_xsim_ieta").value
xsim_intot = ctypes.c_int.in_dll(libpomf, "pomf_xsim_intot").value
xsim_iy = ctypes.c_int.in_dll(libpomf, "pomf_xsim_iy").value
xsim_iytot = ctypes.c_int.in_dll(libpomf, "pomf_xsim_iytot").value
xelems = ctypes.c_int.in_dll(libpomf, "pomf_xelems").value
xiatilde = ctypes.c_int.in_dll(libpomf, "pomf_xiatilde").value
xir = ctypes.c_int.in_dll(libpomf, "pomf_xir").value
xidelta = ctypes.c_int.in_dll(libpomf, "pomf_xidelta").value
xintot = ctypes.c_int.in_dll(libpomf, "pomf_xintot").value
xiy = ctypes.c_int.in_dll(libpomf, "pomf_xiy").value
xiytot = ctypes.c_int.in_dll(libpomf, "pomf_xiytot").value

if use_total_y:
    xiyuse = xiytot
else:
    xiyuse = xiy

def gaussian(x, m, v):
    return np.exp(-0.5*(x-m)**2/v) / np.sqrt(2.0*np.pi*v)

class Pars(ctypes.Structure):
    """Encapsulates static model paramaters and also interfaces with the
    underlying C code.

    """
    
    _fields_ = [
        ("a_bar", ctypes.c_double),
        ("gamma", ctypes.c_double),
        ("sigma_c", ctypes.c_double),
        ("sigma_c_bar", ctypes.c_double),
        ("sigma_a", ctypes.c_double),
        ("kappa", ctypes.c_double),
        ("r_bar", ctypes.c_double),
        ("rho", ctypes.c_double),
        ("rho_bar", ctypes.c_double),
        ("sin_theta", ctypes.c_double),
        ("sin_theta_bar", ctypes.c_double),
        ("a_tilde_init", ctypes.c_double),
        ("r_init", ctypes.c_double),
        ("N0_over_y0", ctypes.c_double),
        ("Ntot0_over_ytot0", ctypes.c_double),
    ]


    def __init__(self, *args, **kwargs):
        self.set_defaults()
        super().__init__(*args, **kwargs)

    def from_list(self, vals):
        """Set parameter values from a list."""

        for val, (key, _) in zip(vals, self._fields_):
            setattr(self, key, val)

    def to_list(self):
        """Return a list of parameter values"""
        return [getattr(self, key) for key, _ in self._fields_]
            
    def set_defaults(self):
        """Initialize the object to the 'default' values of the parameters."""

        self.a_bar = 0.15
        self.gamma = 0.10
        self.sigma_a = 0.03
        self.sigma_c = 0.30
        self.sigma_c_bar = 0.30
        self.kappa = -0.5
        self.r_bar = 0.010
        self.rho = 0.040
        self.rho_bar = 0.030
        self.sin_theta = 0.1
        self.sin_theta_bar = 0.1
        self.a_tilde_init = 0.15
        self.r_init = 0.05
        self.N0_over_y0 = 1.0         # ???
        self.Ntot0_over_ytot0 = 1.0   # ???


    def print_values(self):
        """Dump the current values to standard output"""
        
        for key, _ in self._fields_:
            if key[0] != '_':
                val = getattr(self, key)
                print("{:12s} = {}".format(key, val))

    def r_range(self, a_tilde):
        """Compute the range of the interest rate given `a_tilde`."""

        r0 = self.r_star(a_tilde, 0.0)
        r1 = self.r_star(a_tilde, 1.0)
        return (min(r0, r1), max(r0, r1))

    def a_tilde_pdf(self, t, a_tilde):
        """Evaluate the unconditional a_tilde pdf at `a_tilde` at time `t`
        given initial value `a_tilde_init` at time zero.

        """
        a_tilde_init = self.a_tilde_init
        s = _sigma_a_tilde_star(ctypes.byref(self)) 
        g = self.gamma
        z = np.exp(-g*t)
        D = 0.5*s**2
        m0 = self.a_bar

        m = m0 + (a_tilde_init - m0)*z
        v = D*(1 - z**2)/g

        return gaussian(a_tilde, m, v)
        

    def r_star(self, a_tilde, eta):
        """Map given macrostate `(a_tilde, eta)` to the corresponding interest
        rate.

        """

        return _r_star(a_tilde, eta,
                       ctypes.byref(self)) 

    def eta_star(self, a_tilde, r):
        """Given observables `a_tilde` and `r`, compute corresponding eta
        variable. Inverse function of `r_star` for all `a_tilde`.

        """

        return _eta_star(a_tilde, r,
                         ctypes.byref(self)) 



class MCSimu:

    def __init__(self, n_samples=1000, dt=0.01, *args, **kwargs):
        """Initialize the simulation object. Parameter `n_samples` sets the
        number of random samples generated for approximating the
        process probability distributions, parameter `dt` sets the
        length of time-step used in the Euler step.

        """
        self._n_samples = n_samples
        self.pars = Pars(*args, **kwargs)
        self.dt = dt
        self.mcsimu = _mcsimu_make(ctypes.c_int(n_samples))
        if self.mcsimu is None:
            raise MemoryError('Failed allocating mcsimu object')


    def __del__(self):
        _mcsimu_free(self.mcsimu)


    def run_step(self, ts, obs, **kwargs):
        """Run a simulation to generate data for constructing the model
        process' distribution. The parameter `ts` can be either (1)
        time, up to which the simulation is run, or (2) a list of
        times, at which to evaluate the process' distribution. Initial
        values for the `a_tilde` and `r` components of the solution
        are fetched from the parameters passed to the object during
        its initialization; the `delta` component, or equivalently,
        the `a` component, is sampled from its equilibrium distribution
        function.

        If `ts` is a list of length `n`, then `obs` should be a list
        of `n - 1` elements, or `None`. If `obs` is not `None`, it
        should contain a sequence a sequence of `(y, r)` tuples.
        These will be used in the run to reset the values of output
        and interest rate to the observed values.

        Function returns an `MCStepResult` object, that can
        subsequently be used to extract the process' distribution
        function.

        """

        # self.pars.print_values()

        if len(np.shape(ts)) == 0:
            ts = [ts]

        if obs is not None:
            assert np.shape(obs)[0] == np.shape(ts)[0] - 1

        n = np.shape(ts)[0]
            
        status = _mcsimu_init(self.mcsimu,
                              ctypes.byref(self.pars))
        assert status == 0

        t0 = 0.0
        snapshots = []
        
        for istep, t in zip(range(n), ts):
            status = _mcsimu_run(self.mcsimu,
                                 ctypes.c_double(t - t0),
                                 ctypes.c_double(self.dt),
                                 ctypes.byref(self.pars))
            assert status == 0

            t0 = t

            csamp = _mcsimu_get_samples(self.mcsimu)
            data = np.empty(shape=(self._n_samples, xelems))

            for i in range(self._n_samples):
                for j in range(xelems):
                    data[i, j] = csamp[i*xelems + j]

            snapshots.append(data)

            if obs is not None and istep < n - 1:
                # Select the initial value based on whether we
                # interpret y as the firms' output or the total
                # output. The new value for the other y does not
                # matter, it will not be used anywhere.
                if use_total_y:
                    y_new = 1.0
                    ytot_new = obs[istep, 0]
                else:
                    y_new = obs[istep, 0]
                    ytot_new = 1.0

                status = _mcsimu_init_wobs(self.mcsimu,
                                           ctypes.c_double(obs[istep, 1]),
                                           ctypes.c_double(y_new),
                                           ctypes.c_double(ytot_new),
                                           ctypes.byref(self.pars))
                assert status == 0

 
        return MCStepResult(ts, obs, snapshots, self.pars)

    

class MCLikelihoodModel(GenericLikelihoodModel):
    """Maximum likelihood class for doing simple max. likelihood
    estimations for the parameter values.

    """

    def __init__(self, endog, exog, **kwds):
        self.mcsimu = MCSimu(**kwds)
        super(MCLikelihoodModel, self).__init__(endog, exog, **kwds)

        
    def nloglikeobs(self, params):

        self.mcsimu.pars.a_bar = params[0]
        self.mcsimu.pars.gamma = params[1]
        self.mcsimu.pars.sigma_c = params[2]
        self.mcsimu.pars.sigma_c_bar = params[3]
        self.mcsimu.pars.sigma_a = params[4]
        self.mcsimu.pars.kappa = params[5]
        self.mcsimu.pars.r_bar = params[6]
        self.mcsimu.pars.rho = params[7]
        self.mcsimu.pars.rho_bar = params[8]
        self.mcsimu.pars.sin_theta = params[9]
        self.mcsimu.pars.sin_theta_bar = params[10]
        self.mcsimu.pars.a_tilde_init = params[11]

        if use_total_y:
            self.mcsimu.pars.Ntot0_over_ytot0 = params[12]
        else:
            self.mcsimu.pars.N0_over_y0 = params[12]

        (y_init, r_init) = self.endog[0]

        if y_init != 1:
            raise ValueError("Input y data must be normalised so that the initial value is one.")

        self.mcsimu.pars.r_init = r_init
        
        # Each element in the list of exog variables passed to the MLE
        # object is turned into a singleton list.  We're just
        # flattening the list here.
        ts = list(itertools.chain(*self.exog))

        if ts[0] != 0:
            raise ValueError("Input time data should be such that the initial time is zero")

        mcres = self.mcsimu.run_step(ts[1:], self.endog[1:-1])

        return -mcres.loglhood(self.endog[1:])
    

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):

        if start_params is None:
            pars = Pars()
            pars.set_defaults()
            start_params = [
                pars.a_bar,
                pars.gamma,
                pars.sigma_c,
                pars.sigma_c_bar,
                pars.sigma_a,
                pars.kappa,
                pars.r_bar,
                pars.rho,
                pars.rho_bar,
                pars.sin_theta,
                pars.sin_theta_bar,
                pars.a_tilde_init,
                pars.Ntot0_over_ytot0 if use_total_y else pars.N0_over_y0
            ]

        return super(MCLikelihoodModel,
                     self).fit(start_params=start_params,
                               maxiter=maxiter, maxfun=maxfun, **kwds) 
    
class MCStepResult:
    """Instances of `MCStepResult` class encapsulate results from running
    a simulation of the model given parameters and initial values over
    a fixed time-step.

    """
    
    def __init__(self, ts, obs, snapshots, pars):
        """Construct a `MCStepResult` object, containing data evaluated at a
        sequence of times `ts`, with corresponding sets of samples
        contained in the array of arrays `snapshots`.  The used
        parameters are to be passed in as `pars`.

        """

        self.ts = ts
        self.obs = obs
        self.snapshots = snapshots
        self.pars = pars


    def loglhood(self, obs):
        """Compute the log-likelihood of a sequence of observations. Each
        observation in `obs` should be (y, r) pair. The `delta`
        variable is assumed to be at its equilibrium distribution at
        initial point.

        """

        s = 0.0
        for (ob, data) in zip(obs, self.snapshots):
            kde = sm.nonparametric.KDEMultivariate(data[:,[xiyuse,
                                                           xir]],
                                                   'cc')
            s += np.log(kde.pdf(ob))

        return s


    def get_obs_pdf(self):
        """Construct a probability density function for the observed variables
        at final time of observations. Uses a multivariate kernel
        density estimation on the sample data.

        """

        data = self.snapshots[-1]
        kde = sm.nonparametric.KDEMultivariate(data[:,[xiyuse,
                                                       xir]],
                                               'cc')
        return kde.pdf


    def plot_obs_pdf(self, plott='surf'):
        """Plot the probability density function constructed using
        `get_obs_pdf`.

        """

        xlabel = r'$y$'
        ylabel = r'$r$'
        zlabel = r'$f(y, r)$'

        pdf = self.get_obs_pdf()

        data = self.snapshots[-1]

        if plott == 'scatter':
            scatter3d(data[:,xiyuse], data[:,xir], pdf(),
                      xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)

        elif plott == 'surf':
            xmin, xmax = aminmax(data[:,xiyuse])
            ymin, ymax = aminmax(data[:,xir])

            fplot3d(lambda x, y: pdf((x,y)), (xmin, xmax), (ymin, ymax),
                    xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)                    


    def plot_a_tilde_marginal(self):
        """Plots the a_tilde marginal distribution at the final time, from the
        generated sample and its analytic formula. Mainly for testing
        and self-consistency checking purposes.

        """

        data = self.snapshots[-1]
        temp = data[:,iobs_atilde]
        kde = sm.nonparametric.KDEMultivariate(temp, 'c')

        a_tilde_init = self.pars.a_tilde_init

        f = np.vectorize(lambda x: self.pars.a_tilde_pdf(self.ts[-1],
                                                         x, a_tilde_init))
        x0, x1 = aminmax(temp)
        xs = np.linspace(x0, x1, 100)

        plt.plot(xs, kde.pdf(xs), 'o')
        plt.plot(xs, f(xs))
                
        plt.show()


    def plot_hidden_pdf(self, plott='surf'):
        """Plot the (eta, delta) joint probability density function at final
        time of observations. These amount to the 'hidden' variables,
        as eta is observed only indirectly via the interest rate,
        while delta is completely unobserved.

        """

        xlabel = r'$\eta$'
        ylabel = r'$\delta$'
        zlabel = r'$f(\eta, \delta)$'

        data = self.snapshots[-1]

        temp = np.empty(shape=(len(data),2))
        for i in range(len(data)):
            eta = self.pars.eta_star(data[i,xiatilde], data[i,xir])
            temp[i,0] = eta
            temp[i,1] = data[i,xidelta]
        kde = sm.nonparametric.KDEMultivariate(temp, 'cc')
        pdf = kde.pdf
        
        if plott == 'scatter':
            scatter3d(temp[:,0], temp[:,1], kde.pdf(),
                      xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)                    

        elif plott == 'surf':
            xmin, xmax = aminmax(temp[:,0])
            ymin, ymax = aminmax(temp[:,1])

            fplot3d(lambda x, y: pdf((x,y)), (xmin, xmax), (ymin, ymax),
                    xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)                    


def a_tilde_pdf(t, a_tilde, **kwargs):
    """Simple wrapper function to evaluate the a_tilde unconditional pdf
    at time t, point a_tilde, given initial, time zero value of
    a_tilde_init.

    """

    p = Pars(**kwargs)
    return p.a_tilde_pdf(t, a_tilde)

            
def mckde_pdf(t, **kwargs):
    """Generate the `(y, r)` distribution function at time `t`.

    This is a high-level helper function, intended for quick and dirty
    calculations when the observation pdf is only needed for a single
    set of parameters. For somewhat better performance, use `MCSimu`
    and `MCStepResult` objects directly (see examples).

    """
    mcsimu = MCSimu(**kwargs)
    res = mcsimu.run_step(t, None, **kwargs)
    pdf = res.get_obs_pdf()

    return pdf

            
def mckde_eval_pdf(t, y, r, **kwargs):
    """Evaluate the (y, r) distribution at time t and at point (y, r).
    This is a high-level helper function, intended for quick and dirty
    calculations when the observation pdf is only needed to be
    evaluated at a single point. If the pdf needs to be evaluated
    repeatedly for different (y, r), use `mckde_pdf` function instead
    to construct the pdf function for much better performance, or use
    `MCSimu` and `MCStepResult` objects directly (see examples).

    """

    pdf = mckde_pdf(t, **kwargs)
    return pdf((y, r))


def aminmax(arr):
    """Helper function: Find the min and max of an array, and return them
    as a tuple."""

    return (np.amin(arr), np.amax(arr))


def scatter3d(x, y, z,
              xlabel='x', ylabel='y', zlabel='z'):
    """Helper function: Does a simple 3D scatter plot with given labels.

    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    plt.show()


def fplot3d(f, xr, yr, plotpoints=36,
            xlabel='x', ylabel='y', zlabel='z'):
    """Helper function: Plot the function `f` over the range `xr = (x0,
    x1)`, `yr = (y0, y1)`.

    """

    (xmin, xmax) = xr
    (ymin, ymax) = yr

    xs = np.linspace(xmin, xmax, plotpoints)
    ys = np.linspace(ymin, ymax, plotpoints)

    zs = np.array([[f(x, y) for x in xs] for y in ys])

    xs, ys = np.meshgrid(xs, ys)
            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    plt.show()
