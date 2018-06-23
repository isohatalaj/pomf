
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

_r_star = libpomf.pomf_r_star
_r_star.argtypes = [ctypes.c_double, ctypes.c_double,
                    ctypes.c_void_p]
_r_star.restype = ctypes.c_double

_eta_star = libpomf.pomf_eta_star
_eta_star.argtypes = [ctypes.c_double, ctypes.c_double,
                      ctypes.c_void_p]
_eta_star.restype = ctypes.c_double

_delta_equil_stdev = libpomf.pomf_delta_equil_stdev
_delta_equil_stdev.argtypes = [ctypes.c_void_p]
_delta_equil_stdev.restype = ctypes.c_double

_sigma_a_tilde_star = libpomf.pomf_sigma_a_tilde_star
_sigma_a_tilde_star.argtypes = [ctypes.c_void_p]
_sigma_a_tilde_star.restype = ctypes.c_double

_mcsimu_make = libpomf.pomf_mcsimu_make
_mcsimu_make.argtypes = [ctypes.c_int]
_mcsimu_make.restype = ctypes.c_void_p

_mcsimu_free = libpomf.pomf_mcsimu_free
_mcsimu_free.argtypes = [ctypes.c_void_p]
_mcsimu_free.restype = None

_mcsimu_get_samples = libpomf.pomf_mcsimu_get_samples
_mcsimu_get_samples.argtypes = [ctypes.c_void_p]
_mcsimu_get_samples.restype = ctypes.POINTER(ctypes.c_double)

_mcsimu_init = libpomf.pomf_mcsimu_init
_mcsimu_init.argtypes = [ctypes.c_void_p,   # self
                         ctypes.c_double,   # a_tilde
                         ctypes.c_double,   # r
                         ctypes.c_double,   # delta_mean
                         ctypes.c_double,   # delta_stdev
                         ctypes.c_void_p]
_mcsimu_init.restype = ctypes.c_int

_mcsimu_run = libpomf.pomf_mcsimu_run
_mcsimu_run.argtypes = [ctypes.c_void_p,   # self
                        ctypes.c_double,   # t_end
                        ctypes.c_double,   # dt
                        ctypes.c_void_p]
_mcsimu_run.restype = ctypes.c_int


# ******************************************************************************
# PYTHON API BEGINS HERE

dim = ctypes.c_int.in_dll(libpomf, "pomf_dim").value
wdim = ctypes.c_int.in_dll(libpomf, "pomf_wdim").value

# nobs = number of variables in the state vector generated in the
# simulations, integers iobs_atilde, iobs_r, and iobs_delta are
# indices to the corresponding variables in the state vector, ie.
# state[iobs_r] would give the interest rate in that particular datum.
nobs = ctypes.c_int.in_dll(libpomf, "pomf_nobs").value
iobs_atilde = ctypes.c_int.in_dll(libpomf, "pomf_iobs_atilde").value
iobs_r = ctypes.c_int.in_dll(libpomf, "pomf_iobs_r").value
iobs_delta = ctypes.c_int.in_dll(libpomf, "pomf_iobs_delta").value

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
        ("_x_star", dim*ctypes.c_double),
        ("_U", dim*dim*ctypes.c_double)
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
        return [getattr(self, key) for key, _ in
                filter(lambda s: s[0][0] != '_', self._fields_)]
            
    def set_defaults(self):
        """Initialize the object to the 'default' values of the parameters."""

        self.a_bar = 0.070
        self.gamma = 0.25
        self.sigma_a = 0.035
        self.sigma_c = 0.025
        self.sigma_c_bar = 0.025
        self.kappa = -0.5
        self.r_bar = 0.010
        self.rho = 0.070
        self.rho_bar = 0.050        


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

    def a_tilde_pdf(self, t, a_tilde, a_tilde_init):
        """Evaluate the unconditional a_tilde pdf at `a_tilde` at time `t`
        given initial value `a_tilde_init` at time zero.

        """
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


    def run_step(self, ts, obs, delta_mean=None, delta_stdev=None, **kwargs):
        """Run a simulation to generate data for constructing the model
        process' distribution. Input parameter `obs = (a_tilde, r`
        give the initial, known values of TFP process estimate and
        interest rate. Initial TFP estimate error is described by
        parameters `delta_mean` and `delta_stdev`, giving the mean and
        standard deviation of the initial, Gaussian distribution of
        `delta := a - a_tilde`; if not given, these default to their
        equilibrium values. Parameter `ts` can be either (1) time, up
        to which the simulation is run, or (2) list of times, at which
        to sample the process' distribution. 

        If `ts` is a list, then `obs` can also be a list, and the
        entries in `obs` are used to reset the observable distribution
        at each input time (`delta` variables are not reset). The 
        length of `ts` must be one less the length of `obs` (initial element
        in `obs` is always taken as the zero time init. value).

        Function returns a `MCStepResult` object, that can
        subsequently be used to extract the process' distribution
        function.

        """

        # self.pars.print_values()

        if len(np.shape(ts)) == 0:
            ts = [ts]

        if len(np.shape(obs)) == 2:
            assert np.shape(obs)[0] == np.shape(ts)[0] + 1
            assert np.shape(obs)[1] == 2
            a_tilde_init = obs[0, 0]
            r_init = obs[0, 1]
            reset = True
        else:
            a_tilde_init = obs[0]
            r_init = obs[1]
            reset = False

        rmin, rmax = self.pars.r_range(a_tilde_init)
        assert rmax >= r_init and rmin <= r_init

        if delta_mean is None:
            delta_mean = 0.0

        if delta_stdev is None:
            delta_stdev = _delta_equil_stdev(ctypes.byref(self.pars))

        status = _mcsimu_init(self.mcsimu,
                              ctypes.c_double(a_tilde_init),
                              ctypes.c_double(r_init),
                              ctypes.c_double(delta_mean),
                              ctypes.c_double(delta_stdev),
                              ctypes.byref(self.pars))
        assert status == 0

        t0 = 0.0
        snapshots = []
        
        for istep, t in zip(range(len(ts)), ts):
            status = _mcsimu_run(self.mcsimu,
                                 ctypes.c_double(t - t0),
                                 ctypes.c_double(self.dt),
                                 ctypes.byref(self.pars))
            assert status == 0

            t0 = t

            csamp = _mcsimu_get_samples(self.mcsimu)
            data = np.empty(shape=(self._n_samples, nobs))

            for i in range(self._n_samples):
                for j in range(nobs):
                    data[i, j] = csamp[i*nobs + j]

            snapshots.append(data)

            if reset:
                # Negative sigma_delta will prevent resetting the
                # delta components.
                status = _mcsimu_init(self.mcsimu,
                                      ctypes.c_double(obs[istep, 0]),
                                      ctypes.c_double(obs[istep, 1]),
                                      ctypes.c_double(0.0),
                                      ctypes.c_double(-1.0),
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

        self.mcsimu.pars.from_list(params)

        # For some reason, each element in the list of exog variables
        # passed to the MLE object is turned into a singleton list.
        # We're here just flattening the list here.
        ts = list(itertools.chain(*self.exog))

        mcres = self.mcsimu.run_step(ts[1:], self.endog)

        return -mcres.loglhood(self.endog[1:])
    

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):

        if start_params is None:
            pars = Pars()
            pars.set_defaults()
            start_params = pars.to_list()

        return super(MCLikelihoodModel,
                     self).fit(start_params=start_params,
                               maxiter=maxiter, maxfun=maxfun, **kwds) 
    
class MCStepResult:
    """Instances of `MCStepResult` class encapsulate results from running
    a simulation of the model given parameters and initial values over
    a fixed time-step.

    """
    
    def __init__(self, ts, obs, snapshots, pars):
        self.ts = ts
        self.obs = obs
        self.snapshots = snapshots
        self.pars = pars


    def loglhood(self, obs):
        """Compute the log-likelihood of a sequence of observations. Each
        observation in `obs` should (a_tilde, r) pair. The `delta`
        variable is assumed to be at its equilibrium distribution at
        initial point.

        """

        s = 0.0
        for (ob, data) in zip(obs, self.snapshots):
            kde = sm.nonparametric.KDEMultivariate(data[:,[iobs_atilde,
                                                           iobs_r]],
                                                   'cc')
            s += np.log(kde.pdf(ob))

        return s


    def get_obs_pdf(self):
        """Construct a probability density function for the observed variables
        at final time of observations. Uses a multivariate kernel
        density estimation on the sample data.

        """
        data = self.snapshots[-1]
        kde = sm.nonparametric.KDEMultivariate(data[:,[iobs_atilde,
                                                       iobs_r]],
                                               'cc')
        return kde.pdf


    def plot_obs_pdf(self, plott='surf'):
        """Plot the probability density function constructed using
        `get_obs_pdf`.

        """

        xlabel = r'$\tilde{a}$'
        ylabel = r'$r$'
        zlabel = r'$f(\tilde{a}, r)$'

        pdf = self.get_obs_pdf()

        data = self.snapshots[-1]

        if plott == 'scatter':
            scatter3d(data[:,iobs_atilde], data[:,iobs_r], pdf(),
                      xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)

        elif plott == 'surf':
            xmin, xmax = aminmax(data[:,iobs_atilde])
            ymin, ymax = aminmax(data[:,iobs_r])

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

        if len(np.shape(self.obs)) == 2:
            a_tilde_init = self.obs[0,0]
        else:
            a_tilde_init = self.obs[0]

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
            eta = self.pars.eta_star(data[i,iobs_atilde], data[i,iobs_r])
            temp[i,0] = eta
            temp[i,1] = data[i,iobs_delta]
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


def a_tilde_pdf(t, a_tilde, a_tilde_init, **kwargs):
    """Simple wrapper function to evaluate the a_tilde unconditional pdf
    at time t, point a_tilde, given initial, time zero value of
    a_tilde_init.

    """

    p = Pars(**kwargs)
    return p.a_tilde_pdf(t, a_tilde, a_tilde_init)

            
def mckde_pdf(t, a_tilde_init, r_init, **kwargs):
    """Given values of the observables, (a_tilde_init, r_init) at time 0,
    generate the (a_tilde, r) distribution function at time t.

    This is a high-level helper function, intended for quick and dirty
    calculations when the observation pdf is only needed for a single
    set of parameters. For somewhat better performance, use `MCSimu`
    and `MCStepResult` objects directly (see examples).

    """
    mcsimu = MCSimu(**kwargs)
    res = mcsimu.run_step(t, (a_tilde_init, r_init), **kwargs)
    pdf = res.get_obs_pdf()

    return pdf

            
def mckde_eval_pdf(t, a_tilde, r, a_tilde_init, r_init, **kwargs):
    """Given values of the observables, (a_tilde_init, r_init) at time 0,
    evaluate the (a_tilde, r) distribution at time t and at point
    (a_tilde, r).  This is a high-level helper function, intended for
    quick and dirty calculations when the observation pdf is only
    needed to be evaluated at a single point. If the pdf needs to be
    evaluated repeatedly for different (a_tilde, r), use `mckde_pdf`
    function instead to construct the pdf function for much better
    performance, or use `MCSimu` and `MCStepResult` objects directly
    (see examples).

    """

    pdf = mckde_pdf(t, a_tilde_init, r_init, **kwargs)
    return pdf((a_tilde, r))


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
