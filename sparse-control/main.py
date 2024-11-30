import datetime
import os.path
import pickle
from scipy.io import loadmat
import numpy as np
import math

from numpy.linalg import matrix_rank
from copy import deepcopy

from control_design.control_design import Designer
from control_design.cost_function import CostFunction

# to save results
save_vars = True

# set up control design problem
## import model matrices A and B
exp_id = 8
if exp_id == 1:
    from examples.ex1 import *
elif exp_id == 2:
    from examples.ex2 import *
elif exp_id == 3:
    from examples.ex3 import *
elif exp_id == 4:
    from examples.ex4 import *
elif exp_id == 5:
    from examples.ex5 import *
elif exp_id == 6:
    from examples.ex6 import *
elif exp_id == 8: # Counter Example
    Data = loadmat('./Ipexp/CounterExample.mat')
    A = Data.get('A')
    B = Data.get('B')

## set sparsity constraint
sparsity = max(len(A) - matrix_rank(A), 1)
#sparsity = 5
## set time horizon
h = math.ceil(len(A)/sparsity)
print(f'The control horizon is h = {h}')

## set cost function
cost = 'tr-inv'
cost_func = CostFunction(h, cost)

## run design algorithms
### fully actuated
cost_fully_actuated = cost_func.compute(A, B)
print(f'cost fully actuated: {cost_fully_actuated} \n')
print(f'sparsity: {sparsity}')

### s-sparse greedy
algo = 's-greedy'
designer = Designer(A, B, sparsity, cost_func, algo)
schedule_s_greedy, cost_s_greedy = designer.design()
schedule_s_greedy = [schedule_k for schedule_k in schedule_s_greedy if len(schedule_k) > 0]

print('s-sparse greedy:')
print('input schedule:', schedule_s_greedy)
print(f'cost: {cost_s_greedy} \n')
breakpoint()

### s-sparse greedy + MCMC
if cost_s_greedy == np.inf:
    e0 = 1e-10
designer.set_algo('mcmc')
schedule_s_greedy_mcmc, cost_s_greedy_mcmc = designer.design(schedule=deepcopy(schedule_s_greedy),eps=e0)

print('s-sparse greedy + MCMC:')
print('input schedule:', schedule_s_greedy_mcmc)
print(f'cost: {cost_s_greedy_mcmc} \n')
breakpoint()

### naive greedy
designer.set_algo('greedy')
schedule_greedy, cost_greedy = designer.design(eps=1e-10)
schedule_greedy = [schedule_k for schedule_k in schedule_greedy if len(schedule_k) > 0]

print('greedy:')
print('input schedule:', schedule_greedy)
print(f'cost: {cost_greedy} \n')

### naive greedy + MCMC
designer.set_algo('mcmc')
schedule_greedy_mcmc, cost_greedy_mcmc = designer.design(schedule=deepcopy(schedule_greedy))

print('greedy + MCMC:')
print('input schedule:', schedule_greedy_mcmc)
print(f'cost: {cost_greedy_mcmc} \n')

### naive greedy + MCMC checking rank
schedule_greedy_mcmc_rk, cost_greedy_mcmc_rk = designer.design(schedule=schedule_greedy, check_rank=True)

print('greedy + MCMC with rank check:')
print('input schedule:', schedule_greedy_mcmc_rk)
print(f'cost: {cost_greedy_mcmc_rk}')

# log out results
if save_vars:
    dir_name = 'exp'
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)

    file_name = 'ex' + str(exp_id) + '_' + datetime.datetime.now().strftime('%Y%m%d%I%M')
    with open(dir_name + '/' + file_name + '.pickle', 'wb') as file:
        pickle.dump({
            'A': A,
            'B': B,
            's': sparsity,
            'cost': cost,
            'cost_fully_actuated': cost_fully_actuated,
            's_greedy': (schedule_s_greedy, cost_s_greedy),
            's_greedy_mcmc': (schedule_s_greedy_mcmc, cost_s_greedy_mcmc),
            'greedy': (schedule_greedy, cost_greedy),
            'greedy_mcmc': (schedule_greedy_mcmc, cost_greedy_mcmc),
            'greedy_mcmc_rk': (schedule_greedy_mcmc_rk, cost_greedy_mcmc_rk)
        }, file)
