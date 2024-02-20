# Sparse-Controllability
This Repository contains the source codes for the Sparse Controllability Project. The files include Greedy, Modified Deterministic and Random Actuator Schedulers for Noiseless LDS. MPC based piecewise sparse control for Noisy LDS

### x(k+1) = Ax(k) + Bu(k) k = 0,1,..,K
#### $`A \in R^{n x n}, B \in R^{n x m}`$
#### Helper Files
CtrlMatrix.m - Generate a Controllability Matrix given A, B, K.  
bck_lwr_mtx.m - Generate 2 Large Matrices, One is a Large Block Lower traingle Matrix, Second is column stack of $`A^i`$ i= 1,2,...,K. Used in MPC Formulation.  
Erdos_Renyi.m - Generate a Erdos Renyi Random Graph where each edge is present with probability $`2*ln(n)/n`$.  
KF.m - Kalman Filter.  
KF_prd.m - Kalman Filter return the prediction for the next state..  
OMP.m - Orthogonal Matching Pursuit.  
POMP.m - Piecewise OMP.  
SpaIpDsg.m - Sparse Input Design using POMP given $`x_0, x_f`$ (Algorithm 1 from the Journal).  
SpaIpDsg_1.m - Sparse Input Design using POMP given $`x_0, x_f`$ (Residue Level is Parameterized).  
SparseScheduling.m - Wrapper file for Deterministic Actuator Scheduling (Setup, CleanUp).  
Plotting.m - Organize some plots.  
TestScript.m - Test the new Scripts.  

## Schedulers Routines
DualSet.m - Deterministic Actuator Scheduling (Adapted to Piecewise Sparsity).  
DualSet_2.m - Deterministic Actuator Scheduling for Unweighted Scheduling (Adapted to Piecewise Sparsity Case) -- *Not used in the current Plots*.  
GreedyScheduling_Aopt_1.m - Propsed Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$).  
GreedyScheduling_Aopt_2.m - Reverse Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$) -- *Not used in current Plots*.  
GreedyScheduling_Eopt.m - Greedy Scheduling with E-Optimality (min. Eigenvalue metric) -- *Not used in the current Plots*.  
GreedyScheduling_Static_Aopt_1.m - Propsed Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$) for Fixed Support Case -- *Not used in the current Plots*.  
RandomSamp_Aopt.m - Weighted Random Sampling adapted to Piecewise Sparsity.  
RandomSamp_Aopt_1.m - Original Random Sampling (Average Sparsity) -- *Not used in the current Plots*.  
RandomSamp_Aopt_2.m - UnWeighted Random Sampling (Sample actuators without Replacement) adapted to Piecewise Sparsity.  
ImpSamp.m - Derived from the principles of Random Sampling used to PMF to pick the s actuators with highest probability -- *Not used in the current Plots*.  

## Main Files - Execute to Generate Plots, Data
ActSch_R1_Comparison.m - Sparse Actuator Scheduling Result 1 : Comparison of Greedy, Deterministic, Random Scheudlers -- **Used in CDC Conference Draft**.  
ActSch_R2_PerfVsTime.m - Sparse Actuator Scheduling Result 2 : A-Optimality Performance over Time (K>n) -- Not *used for the current Plots*.  
ActSch_R3_NPerfVsFSparsity.m - Sparse Actuator Scheduling Result 3 : $`\frac{Tr(W_S^{-1})}{Tr(W^{-1})}`$ vs $`\frac{s}{m}`$ for various n and lower and upper bounds -- **Used in CDC Conference Draft**.  
ActSch_R4_TStepsVsSparsity.m - Sparse Actuator Scheduling Result 4 : Bounds for the Number of iterations of Greedy Algorithm to Stop over various Sparsity levels -- *Not used for the current Plots*.  
ActSch_R5_ApprxVsSparsity.m - Sparse Actuator Scheduling Result 5 : Ratio of smallest and largest Eigenvalues for some sparsity level compared  to Fully actuated case -- *Not used for the current Plots*.  
ActSch_R6_CDFPlot.m - Sparse Actuator Scheduling Result 6 : CDF plot of the proposed Greedy algorithm for fixed support and dynamic support **Used in CDC Conference Draft**.  
NsysLS_KF_P_OMP.m - MPC based controller for Noisy LDS, Comparision for both OMP and POMP case.  
POMP_R1_Performance.m - Results for POMP Algorithm for different Sparsity, Trails, Systems, residue levels.  
POMP_R2_Performance.m - Time Steps (K) taken for controllability with Bounds.  
