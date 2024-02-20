# Sparse-Controllability
This Repository contains the source codes for the Sparse Controllability Project. The files include Greedy, Modified Deterministic and Random Actuator Schedulers for Noiseless LDS. MPC based piecewise sparse control for Noisy LDS

### $`x(k+1) = Ax(k) + Bu(k) k = 0,1,..,K`$
#### $`A \in R^{n \times n}, B \in R^{n \times m}`$
### Helper Files
<ul>
  <li>CtrlMatrix.m - Generate a Controllability Matrix given A, B, K.</li>  
  <li>bck_lwr_mtx.m - Generate 2 Large Matrices, One is a Large Block Lower traingle Matrix, Second is column stack of $`A^i, i= 1,2,...,K`$. Used in MPC Formulation.  </li>
  <li>Erdos_Renyi.m - Generate a Erdos Renyi Random Graph where each edge is present with probability $`2\frac{\ln(n)}{n}`$.  </li>
  <li>KF.m - Kalman Filter.  </li>
  <li>KF_prd.m - Kalman Filter return the prediction for the next state..  </li>
  <li>OMP.m - Orthogonal Matching Pursuit.  </li>
  <li>POMP.m - Piecewise OMP.  </li>
  <li>SpaIpDsg.m - Sparse Input Design using POMP given $`x_0, x_f`$ (Algorithm 1 from the Journal).  </li>
  <li>SpaIpDsg_1.m - Sparse Input Design using POMP given $`x_0, x_f`$ (Residue Level is Parameterized).  </li>
  <li>SparseScheduling.m - Wrapper file for Deterministic Actuator Scheduling (Setup, CleanUp).  </li>
  <li>Plotting.m - Organize some plots.  </li>
  <li>TestScript.m - Test the new Scripts.  </li>
</ul>

## Schedulers Routines
<ul>
  <li>DualSet.m - Deterministic Actuator Scheduling (Adapted to Piecewise Sparsity).  </li>
  <li>DualSet_2.m - Deterministic Actuator Scheduling for Unweighted Scheduling (Adapted to Piecewise Sparsity Case) -- *Not used in the current Plots*.  </li>
  <li>GreedyScheduling_Aopt_1.m - Propsed Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$).  </li>
  <li>GreedyScheduling_Aopt_2.m - Reverse Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$) -- *Not used in current Plots*.  </li>
  <li>GreedyScheduling_Eopt.m - Greedy Scheduling with E-Optimality (min. Eigenvalue metric) -- *Not used in the current Plots*.  </li>
  <li>GreedyScheduling_Static_Aopt_1.m - Propsed Greedy Scheduling with A-Optimality ($`Tr(W_S^{-1})`$) for Fixed Support Case -- *Not used in the current Plots*.  </li>
  <li>RandomSamp_Aopt.m - Weighted Random Sampling adapted to Piecewise Sparsity.  </li>
  <li>RandomSamp_Aopt_1.m - Original Random Sampling (Average Sparsity) -- *Not used in the current Plots*.  </li>
  <li>RandomSamp_Aopt_2.m - UnWeighted Random Sampling (Sample actuators without Replacement) adapted to Piecewise Sparsity.  </li>
  <li>ImpSamp.m - Derived from the principles of Random Sampling used to PMF to pick the s actuators with highest probability -- *Not used in the current Plots*. </li>
</ul>

## Main Files - Execute to Generate Plots, Data
<ul>
  <li>ActSch_R1_Comparison.m - Sparse Actuator Scheduling Result 1 : Comparison of Greedy, Deterministic, Random Scheudlers -- **Used in CDC Conference Draft**.  </li>
  <li>ActSch_R2_PerfVsTime.m - Sparse Actuator Scheduling Result 2 : A-Optimality Performance over Time (K>n) -- Not *used for the current Plots*.  </li>
  <li>ActSch_R3_NPerfVsFSparsity.m - Sparse Actuator Scheduling Result 3 : $`\frac{Tr(W_S^{-1})}{Tr(W^{-1})}`$ vs $`\frac{s}{m}`$ for various n and lower and upper bounds -- **Used in CDC Conference Draft**.  </li>
  <li>ActSch_R4_TStepsVsSparsity.m - Sparse Actuator Scheduling Result 4 : Bounds for the Number of iterations of Greedy Algorithm to Stop over various Sparsity levels -- *Not used for the current Plots*.  </li>
  <li>ActSch_R5_ApprxVsSparsity.m - Sparse Actuator Scheduling Result 5 : Ratio of smallest and largest Eigenvalues for some sparsity level compared  to Fully actuated case -- *Not used for the current Plots*.  </li>
  <li>ActSch_R6_CDFPlot.m - Sparse Actuator Scheduling Result 6 : CDF plot of the proposed Greedy algorithm for fixed support and dynamic support -- **Used in CDC Conference Draft**.  </li>
  <li>NsysLS_KF_P_OMP.m - MPC based controller for Noisy LDS, Comparision for both OMP and POMP case.  </li>
  <li>POMP_R1_Performance.m - Results for POMP Algorithm for different Sparsity, Trails, Systems, residue levels.  </li>
  <li>POMP_R2_Performance.m - Time Steps (K) taken for controllability with Bounds.  </li>
</ul>
