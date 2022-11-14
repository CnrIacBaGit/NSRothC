Added features:
1) Solver for fractional order model;
2) Classical scheme for RothC model;
3) Solution time measurement.
Minor corrections:
1) References to files changed to local, case-sensitive and with `/` slash. These correction were needed to make the code run on Linux.
Main changes in code:
1) new main2.py was created modifying main.py
-- import of NS module commented - it is unused further in the code and I had problems importing it on my system;
-- explicit import of scipy.linalg was added - importing scipy as a whole didn't lead to importing linalg on my system;
-- block of code to read additional parameters from xlsx file was added (lines 47-67);
-- call to fractional order solver (in external module) was added in line 171;
-- classical scheme for RothC model was implemented in lines 173-180;
2) fr_rothC.py file was added. It contains the implementation of a fractional order solver
-- function fr_rothC implements the solution starting from time point t0 (in years). Ending time of simulation is passed to it in tend parameter in months. Array of monthly time and C is returned;
-- one step of Crank-Nicholson finite-difference scheme is implemented in scheme2_one_step. The most part of time is here spent on computing the "memory" (integral approximation) in lines 67-74. This procedure needs all already computed solutions. They are stored in old_c array that is ammended in fr_rothC function;
Changes in configuration files:
-- file hoosfield_scenario1a.xlsx was created on the base of hoosfield_scenario1.xlsx. It was ammended with four lines that contain parameters of the fractional order solver. Their description was added in column D of this file. Default behaviour (for config files without these lines) of main2.py remained the same as in the main.py. Slashes in changed config file were reversed to '/'.
