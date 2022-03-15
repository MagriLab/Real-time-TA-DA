import sys
import os
os.environ["OMP_NUM_THREADS"] = '1' # imposes only one core
import numpy as np
import matplotlib.pyplot as plt
import h5py
import skopt
from skopt.space import Real
from skopt.learning import GaussianProcessRegressor as GPR
from skopt.learning.gaussian_process.kernels import Matern, WhiteKernel, Product, ConstantKernel
import matplotlib as mpl
from scipy.io import loadmat
import time
plt.style.use('dark_background')
#Latex
mpl.rc('text', usetex = True)
mpl.rc('font', family = 'serif')

folder = folder + '\\'


# In functions there are the ESN functions
exec(open(folder + "Functions.py").read())


folder = folder + 'ESN_bias_data\\'

#### Dataset management

# Wave model
dt = 1/10000  #integration time step
upsample  = 5 #upsample of the ESN from dt



data      = loadmat(folder + 'BIAS_' + filename)
# data      = loadmat('./BIAS_Wave_high_G_eqn_Galerkin_low_sqrt.mat')
string    = 'bias' #'residual'

U        = np.array(data[string].T[::upsample])[:,:]
N_latent = U.shape[1] #number of dimensions

# To normalize in Lyapunov Times when appropriate
t_lyap    = 10.**(-1)       # Lyapunov time
dt        = dt*upsample     # ESN dt
N_lyap    = int(t_lyap/dt)  # number of time steps in one Lyapunov time

# number of time steps for washout, train, validation, test
N_washout = 50
N_train   = 8*N_lyap
N_val     = int(2*N_lyap) #length of the data used is train+val
N_test    = 500*N_lyap

print('')
print('% data used for training:', (N_train+N_val)/(U.shape[0]))
print('Time :',(N_train+N_val)*dt, 's')

#compute norm (normalize inputs by component range)
U_data = U[:N_washout+N_train+N_val]
m = U_data.min(axis=0)
M = U_data.max(axis=0)
norm = M-m

# washout
U_washout = U[:N_washout]
# training + validation data
U_tv  = U[N_washout:N_washout+N_train+N_val-1]
Y_tv  = U[N_washout+1:N_washout+N_train+N_val].reshape(1,N_train+N_val-1,N_latent)

## add noise to inputs and targets during training
noisy = True
noise_level = 0.03 #larger values promote stability in long term, but hinder time accuracy
noises   = np.array([noise_level]) #target noise (if multiple values it's optimized)
data_std = np.std(U,axis=0)
seed = 0                        #to be able to recreate the data
rnd  = np.random.RandomState(seed)

if noisy: #input noise (it is not optimized)
    for i in range(N_latent):
        U_tv[:,i] = U_tv[:,i].copy()                         + rnd.normal(0, noise_level*data_std[i], N_train+N_val-1)

    Y_tv  = np.zeros((len(noises), N_train+N_val-1, N_latent)) #U[N_washout+1:N_washout+N_train+N_val]
    for jj in range(noises.size):
        for i in range(N_latent):
                Y_tv[jj,:,i] = U[N_washout+1:N_washout+N_train+N_val,i].copy()                             + rnd.normal(0, noises[jj]*data_std[i], N_train+N_val-1)


#### ESN Initiliazation Hyperparameters

bias_in = .1    #input bias
bias_out = 1.0  #output bias 
N_units = 100   #units in the reservoir
dim = U.shape[1] # dimension of inputs (and outputs) 
connectivity   = 5 #average connections per row in the state matrix 
sparseness     =  1 - connectivity/(N_units-1) 

tikh = np.array([1e-4,1e-8,1e-12,1e-16])  # Tikhonov factor (optimized in a grid search in these values)


####  Grid Search and Bayesian Optimization 
# Here we set the parameters for Grid Search and Bayesian Optimization.

n_tot = 20   #Total Number of Function Evaluatuions
n_in  = 0    #Number of Initial random points

spec_in     = .5   #range for hyperparameters (spectral radius and input scaling)
spec_end    = 1.2    
in_scal_in  = np.log10(1e-3)
in_scal_end = .3

# In case we want to start from a grid_search, the first n_grid^2 points are from grid search
# if n_grid^2 = n_tot then it is pure grid search
n_grid = 4  # (with n_grid**2 < n_tot you get Bayesian Optimization)

# computing the points in the grid
if n_grid > 0:
    x1    = [[None] * 2 for i in range(n_grid**2)]
    k     = 0
    for i in range(n_grid):
        for j in range(n_grid):
            x1[k] = [spec_in + (spec_end - spec_in)/(n_grid-1)*i,
                     in_scal_in + (in_scal_end - in_scal_in)/(n_grid-1)*j]
            k   += 1

# range for hyperparameters
search_space = [Real(spec_in, spec_end, name='spectral_radius'),
                Real(in_scal_in, in_scal_end, name='input_scaling')]

# ARD 5/2 Matern Kernel with sigma_f in front for the Gaussian Process
kernell = ConstantKernel(constant_value=1.0, constant_value_bounds=(1e-1, 3e0))*                  Matern(length_scale=[0.2,0.2], nu=2.5, length_scale_bounds=(1e-2, 1e1)) 


#Hyperparameter Optimization using either Grid Search or Bayesian Optimization
def g(val):
    
    #Gaussian Process reconstruction
    b_e = GPR(kernel = kernell,
            normalize_y = True, #if true mean assumed to be equal to the average of the obj function data, otherwise =0
            n_restarts_optimizer = 3,  #number of random starts to find the gaussian process hyperparameters
            noise = 1e-10, # only for numerical stability
            random_state = 10) # seed
    
    
    #Bayesian Optimization
    res = skopt.gp_minimize(val,                         # the function to minimize
                      search_space,                      # the bounds on each dimension of x
                      base_estimator       = b_e,        # GP kernel
                      acq_func             = "gp_hedge", # the acquisition function
                      n_calls              = n_tot,      # total number of evaluations of f
                      x0                   = x1,         # Initial grid search points to be evaluated at
                      n_random_starts      = n_in,       # the number of additional random initialization points
                      n_restarts_optimizer = 3,          # number of tries for each acquisition
                      random_state         = 10,         # seed
                           )   
    return res


#### Train and Validate the network

from skopt.plots import plot_convergence

#Number of Networks in the ensemble
ensemble = 1
# Which validation strategy
val  = RVC_Noise
N_fo = 15  # number of folds
N_in = 0   # interval before the first fold
N_fw = N_lyap//2 # how many steps forward the validation interval is shifted (N_fw*N_fo has to be smaller than N_train)

#Quantities to be saved
par      = np.zeros((ensemble, 4))      # GP parameters
x_iters  = np.zeros((ensemble,n_tot,2)) # coordinates in hp space where f has been evaluated
f_iters  = np.zeros((ensemble,n_tot))   # values of f at those coordinates
minimum  = np.zeros((ensemble, 5))      # minima found per each member of the ensemble

tikh_opt = np.zeros(n_tot) #optimal tikhonov
noise_opt= np.zeros(n_tot) #optimal noise  (both then saved in minimum)
Woutt    = np.zeros(((ensemble, N_units+1,dim))) #to save output matrix
Winn     = np.zeros((ensemble,dim+1, N_units))   #to save input matrix
Ws       = np.zeros((ensemble,N_units, N_units)) #to save state matrix

# save the final gp reconstruction for each network
gps      = [None]*ensemble

# to check time
ti       = time.time()

print('')
print('HYPERPARAMETER SEARCH:')
print(str(n_grid) + 'x' +str(n_grid) + ' grid points plus ' + str(n_tot-n_grid**2) + ' points with Bayesian Optimization')
print('')

for i in range(ensemble):
    
    print('Realization    :',i+1)
    
    k        = 0
    
    # Win and W generation
    seed= i+1
    rnd = np.random.RandomState(seed)

    Win = np.zeros((dim+1, N_units))
    for j in range(N_units):
        Win[rnd.randint(0, dim+1),j] = rnd.uniform(-1, 1) #only one element different from zero per row
    
    # practical way to set the sparseness
    W = rnd.uniform(-1, 1, (N_units, N_units)) * (rnd.rand(N_units, N_units) < (1-sparseness))
    spectral_radius = np.max(np.abs(np.linalg.eigvals(W)))
    W /= spectral_radius #scaled to have unitary spec radius
    
    # Bayesian Optimization
    res        = g(val)
    
    
    #Saving Quantities for post_processing anmd future use of the network
    gps[i]     = res.models[-1]    
    gp         = gps[i]
    x_iters[i] = np.array(res.x_iters)
    f_iters[i] = np.array(res.func_vals)
    minimum[i] = np.append(res.x,[tikh_opt[np.argmin(f_iters[i])],
                                  noise_opt[np.argmin(f_iters[i])],res.fun])
    params     = gp.kernel_.get_params()
    key        = sorted(params)
    par[i]     = np.array([params[key[2]],params[key[5]][0], params[key[5]][1], gp.noise_])
    
    Woutt[i]   = train_save_n(U_washout, U_tv, U[N_washout+1:N_washout+N_train+N_val],
                              minimum[i,2],10**minimum[i,1], minimum[i,0], minimum[i,3])
    
    Winn[i]    = Win.copy()
    Ws[i]      = W.copy()
    
    #Plotting Optimization Convergence for each network
    print('')
    print('Time per hyperparameter evaluation:', -(ti - time.time())/n_tot)
    print('Best Results: x', minimum[i,0], 10**minimum[i,1], minimum[i,2], minimum[i,3],
          ', f', -minimum[i,-1])
    print('')


#### Quick test
# Running the networks in the test set.

subplots = 10 #number of plotted intervals
plt.rcParams["figure.figsize"] = (15,2*subplots)
plt.subplots(subplots,1)

N_test   = 100 #number of test set intervals
N_tstart = N_washout + N_train + N_val + 95 #when the first interval starts
NNN      = 10    #length of the test set interval
subloops = 5*NNN #number of updates of the input with correct data inside the test interval
N_intt   = NNN*N_val//subloops #length of each subinterval before correct input data is given

for k in range(ensemble):
    
    # load matrices and hyperparameters
    Win      = Winn[k].copy()
    W        = Ws[k].copy()
    Wout     = Woutt[k].copy()
    rho      = minimum[k,0].copy()
    sigma_in = 10**minimum[k,1].copy()
    
    
    errors   = np.zeros(N_test)
    
    #Different intervals in the test set
    for i in range(N_test):
                
        # data for washout and target in each interval
        U_wash    = U[N_tstart - N_washout +i*N_val : N_tstart + i*N_val]
        Y_t       = U[N_tstart + i*N_val            : N_tstart + i*N_val + N_intt*subloops] 
        
        #washout for each interval
        xa1        = open_loop(U_wash, np.zeros(N_units), sigma_in, rho)[-1]
        
        Yh_t, xa1  = closed_loop_test(N_intt-1, xa1, Y_t[0], Wout, sigma_in, rho)
        
        # Do the multiple subloops inside each test interval
        if subloops > 1:
            for j in range(subloops-1):
                Y_start    = Y_t[(j+1)*N_intt-1].copy() #
#                 Y_start    = Yh_t[-1].copy()# #uncomment this to not update input
                Y1, xa1    = closed_loop_test(N_intt, xa1, Y_start, Wout, sigma_in, rho)
                Yh_t       = np.concatenate((Yh_t,Y1[1:]))
        
        errors[i] = np.log10(np.mean((Yh_t-Y_t)**2)/np.mean(norm**2))
        
#         print(np.log10(np.mean((Yh_t-Y_t)**2)/np.mean(norm**2)))
        if i <subplots:
            plt.subplot(subplots,1,i+1)
            plt.plot(np.arange(N_intt*subloops)*dt,Y_t[:,:2]/norm[:2], 'w')
            plt.plot(np.arange(N_intt*subloops)*dt,Yh_t[:,:2]/norm[:2], '--r')
        
print('Median and max error in test:', np.median(errors), errors.max())  
plt.tight_layout()
plt.savefig(folder + 'Test_run_' + filename + '.pdf')
plt.close()



#### Visualize hyperparameter search


# Plot Gaussian Process reconstruction for each network in the ensemble afte n_tot evaluations
# The GP reconstruction is based on the n_tot function evaluations decided in the search

# points to evaluate the GP at
n_length    = 100
xx, yy      = np.meshgrid(np.linspace(spec_in, spec_end,n_length),
                          np.linspace(in_scal_in, in_scal_end,n_length))
x_x         = np.column_stack((xx.flatten(),yy.flatten()))
x_gp        = res.space.transform(x_x.tolist())  ##gp prediction needs this normalized format 
y_pred      = np.zeros((ensemble,n_length,n_length))

plt.rcParams["figure.figsize"] = (15,5*ensemble)
plt.rcParams["font.size"] = 20


fig = plt.figure()

for i in range(ensemble):
    # retrieve the gp reconstruction
    gp         = gps[i]
    
    plt.subplot(ensemble, 1, 1+i)
    
    amin = np.amin([10,f_iters.max()])
    
    y_pred[i] = np.clip(-gp.predict(x_gp), a_min=-amin,
                        a_max=-f_iters.min()).reshape(n_length,n_length) 
                        # Final GP reconstruction for each realization at the evaluation points
        
    plt.title('Mean GP of realization \#'+ str(i+1))
    
    #Plot GP Mean
    plt.xlabel('Spectral Radius')
    plt.ylabel('Input Scaling (log-scale)')
    CS      = plt.contourf(xx, yy, y_pred[i],levels=20,cmap='Blues')
    cbar = plt.colorbar()
    cbar.set_label('-$\log_{10}$(MSE)',labelpad=15)
    CSa     = plt.contour(xx, yy, y_pred[i],levels=20,colors='black',
                          linewidths=1, linestyles='solid', alpha=0.3)
    
    #   Plot the n_tot search points
    plt.scatter(x_iters[i,:n_grid**2,0],x_iters[i,:n_grid**2,1], c='gray', marker='^', alpha=0.8)
    plt.scatter(x_iters[i,n_grid**2:,0],x_iters[i,n_grid**2:,1],
                c=np.arange(n_tot-n_grid**2), marker='s', cmap='Reds') #bayesian Optimization ones are plotted with sequential color indicating their orderind
    
fig.tight_layout()
plt.savefig(folder + 'Hyperparameter_search_' + filename + '.pdf')
plt.close()


#### Save Results

from scipy.io import savemat

fln = folder + 'ESN_' + filename + '.mat'
with open(fln,'wb') as f:  # need 'wb' in Python3
    savemat(f, {"norm": norm})
    savemat(f, {'hyperparameters':[minimum[0,0], 10**minimum[0,1],bias_in]})
    savemat(f, {"Win": Winn})
    savemat(f, {'W': Ws})
    savemat(f, {"Wout": Woutt})
    savemat(f, {"dt": dt})
    savemat(f, {"N_washout": N_washout*1.})
    savemat(f, {"N_units": N_units+1.})
    savemat(f, {"N_dim": dim*1.})
    savemat(f, {"training_time": (N_train+N_val)*dt*1.})
