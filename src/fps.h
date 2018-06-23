
#include <cs.h>

#include "smaux.h"
#include "util.h"
#include "rgrid.h"

/* Set to non-zero to fail whenever cannot guarantee model
   validity. */
#define MAY_FAIL 0

/* If non-zero the algorithm will proceed even when transition rates
 * are negative. When set, convexity of delta coefficients is
 * enforced. */
#define ALLOW_NEGTRANSR 0

/* Testing purposes: If 1, forces delta coefficients to be 1/2 always. 
 * Probably need to set ALLOW_NEGTRANSR to non-zero as well. */
#define FORCE_05_DELTA 0

/* ******************************************************************************** */

/* 
 * Stochastic differential equation:
 *
 *   dx_i_t = mu_i(x_t) dt + s_ij(x_t) dz_j_t
 *
 * Drift vector and diffusion matrix:
 *
 *   D_ij := s_ik s_jk / 2
 *   K_i  := mu_i - d_j D_ij
 *   J_i  := mu_i f - d_j[ D_ij f ]
 *         = [mu_i - d_j D_ij] f - D_ij d_j f
 *         = K_i f - D_ij d_j f
 */

typedef int kfe_sde_eval_t(const double *x,  /*< State vector, dimension `d` */
			   double *b,        /*< Drift vector, dimension `d` */
			   double *s,        /*< Diff matrix,
					       dimension `k*d` where
					       `k` number of driving
					       Brownians. */
			   double *dsdx,     /*< Derivatives of `s`,
					       first `x[0]` matrix of
					       derivatives, then
					       `x[1]` etc., total size
					       `d*d*d`. */
			   void *params);
typedef int kfe_sde_eval_bounds_t(double *bounds,
				  void *params);


typedef struct {
  int n;	/*< State process dimension */
  int k;	/*< Dimension of the driving Brownian process */
  kfe_sde_eval_t *eval;	/*< Drift-diffusion coefficient evaluator */
  kfe_sde_eval_bounds_t *eval_bounds;
			/*< Function to compute boundaries */
} kfe_sde_t;

#define KFE_SDE_NIL {0, 0, NULL, NULL}
const kfe_sde_t kfe_sde_nil;


/* ******************************************************************************** */
/* MONTE-CARLO SOLVER */


/**
 * Normal distributed rng function type.
 */
typedef double kfe_nrng_t(double mu, double sigma, void *);

typedef struct {
  int n_samples;  /*< Number of independent samples of the process to simulate. */
  double dt;      /*< Length of a single time step. */
  double *data;   /*< Sample data. */
  kfe_sde_t sde;  /*< SDE description. */
  double *bounds; /*< Solution domain boundaries. */
  void *params;   /*< SDE parameters. */
  double t;       /*< Total time. */

  kfe_nrng_t *rng;
  void *rng_state;
} kfe_mcsolver_t;

#define KFE_MCSOLVER_NIL {0, 0.0, NULL, KFE_SDE_NIL, NULL, NULL, 0.0, NULL, NULL}
const kfe_mcsolver_t kfe_mcsolver_nil;

/**
 * Initialize Monte-Carlo solver. Must call
 * `kfe_mcsolver_init_state_pars` also before using the solver.
 */
int
kfe_mcsolver_init(kfe_mcsolver_t *self,
		  const kfe_sde_t *sde,
		  int n_samples, double dt,
		  kfe_nrng_t *rng,
		  void *rng_state);

void
kfe_mcsolver_destroy(kfe_mcsolver_t *self);

int
kfe_mcsolver_init_state_pars(kfe_mcsolver_t *self,
			     const double *init_state,
			     void *params);

int
kfe_mcsolver_step(kfe_mcsolver_t *self);


/* ******************************************************************************** */
/* SPARSE MATRIX SOLVER */

/* CHANG-COOPER TYPE DISCRETIZATION
 *
 * ONE AXIS SUPERPOSITIONS
 * d = delta for short
 *
 * f(x + h_i/2) = d(x) f(x) + [1 - d(x)] f(x + h_i)
 *  
 * 
 * TWO AXIS SUPERPOSITIONS
 * g = gamma for short
 *
 * f(x + h_i/2 + h_j/2) = g^0_ij(x) f(x)       + g^1_ij(x) f(x + h_i) 
 *                      + g^2_ij(x) f(x + h_j) + g^3_ij(x) f(x + h_i + h_j)
 *
 * g^0_ij(x) + g^1_ij(x) + g^2_ij(x) + g^3_ij(x) = 1
 *
 * g^0_ij(x) = g^0_ji(x)
 * g^1_ij(x) = g^2_ji(x)   
 * g^3_ij(x) = g^3_ji(x)
 *
 * Symmetry used to compress g coefficient storage, but unit sum rule
 * is not used. In total, there are then 4*n*(n-1)/2 coefficients to
 * store per lattice point. We will only store the upper triangle of
 * each g, i.e. g^l_ij, i = 0,..,n-1, j=i+1,..,n-1 for each l=0,1,2,3.
 *
 */


typedef int kfe_eval_D_t(const double *x,
			 double *D,
			 void *params);
typedef int kfe_eval_K_t(const double *x,
			 double *K,
			 void *params);

typedef struct {
  int n;		/*< Dimension of the state space. */
  kfe_eval_D_t *D;      /*< Diffusion matrix evaluator. */
  kfe_eval_K_t *K;      /*< Drift vector evaluator. */
  void *params;
  double *bounds;	/*< Solution domain boundaries, 2`n` array of
			  lower bound - upper bound pairs. */
} kfe_problem_t;

kfe_problem_t *
kfe_problem_from_sde(kfe_sde_t *sde,
		     void *params);

void
kfe_problem_free(kfe_problem_t *self);


int kfe_gamma_size(int n);
int kfe_gamma_stride(int n);

int
kfe_gamma_index(int n, int l, int i, int j);

int kfe_solver_iwork_size(int n_dim);
int kfe_solver_work_size(int n_dim);

typedef struct {
  kfe_problem_t *prob;
  double *gamma;	/*< Weight factors for two-axis superpositions. */
  double *delta;        /*< Weights for one-axis superpositions. */
  rgrid_t rgrid;        /*< Regular grid descriptor. */

  int *iwork;
  double *work;
} kfe_solver_t;

#define KFE_SOLVER_NIL {NULL, NULL, NULL, RGRID_NIL, NULL, NULL}
const kfe_solver_t kfe_solver_nil;

int
kfe_solver_init(kfe_solver_t *self,
		kfe_problem_t *prob,
		const int *ns);

void
kfe_solver_destroy(kfe_solver_t *self);

cs *
kfe_solver_make_KFE_matrix(kfe_solver_t *self);

int
kfe_solver_eval_D_(kfe_solver_t *self,
		   const double *x,
		   int i, double hi,
		   int j, double hj,
		   double *D,
		   double *x_temp);

int
kfe_solver_eval_K_(kfe_solver_t *self,
		   const double *x,
		   int i, double hi,
		   double *K,
		   double *x_temp);

int
kfe_solver_eval_delta_std(kfe_solver_t *self,
			  int i,
			  const double *x,
			  int h_axis, double h_len,
			  double *delta_i);

int
kfe_solver_eval_gamma_std(kfe_solver_t *self,
			  int i, int j,
			  const double *x,
			  int h_axis, double h_len,
			  int h_axis_2, double h_len_2,
			  double gamma_ij[4]);

int
kfe_solver_eval_beta_plus(kfe_solver_t *self,
			  const double *x,
			  int i,
			  int h_axis, double h_len,
			  double *beta_i);

int
kfe_solver_eval_beta_minus(kfe_solver_t *self,
			   const double *x,
			   int i,
			   int h_axis, double h_len,
			   double *beta_i);

int
kfe_solver_eval_theta_S_plus(kfe_solver_t *self,
			     const double *x,
			     int i, int j,
			     int h_axis, double h_len,
			     int h_axis_2, double h_len_2,
			     double *theta_S);

int
kfe_solver_eval_theta_S_minus(kfe_solver_t *self,
			      const double *x,
			      int i, int j,
			      int h_axis, double h_len,
			      int h_axis_2, double h_len_2,
			      double *theta_S);

int
kfe_solver_eval_theta_A(kfe_solver_t *self,
			const double *x,
			int i, int j,
			int h_axis, double h_len,
			int h_axis_2, double h_len_2,
			double *theta_A);

int
kfe_solver_eval_Gamma(kfe_solver_t *self,
		      const double *x,
		      double *Gamma);
  
int
kfe_solver_set_init_state(kfe_solver_t *self,
			  const double *x_init,
			  double *data);
