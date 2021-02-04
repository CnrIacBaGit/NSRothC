NSRothC - NonStandard RothC models in Matlab

Angela Martiradonna
Department of Mathematics, University of Bari, Italy
Institute for Applied Mathematics (IAC), CNR, Bari, Italy
angela.martiradonna@uniba.it

NSRothC is a MATLAB routine for implementig non standard RothC type methods \cite{coleman1995rothc} for simaluting the dynamics of soil organic carbon (SOC)  contained in the organic matter. Five compartments of the soil are considered:  DPM (decomposable plant material), RPM (resistant plant material),  BIO (microbial biomass),  HUM  (humified organic matter) and IOM  (inert organic matter). 

NSRothC implemets three numerical schemes for the SOC dynamics:
•	the RothC discrete classical scheme in [1];
•	the Exponential Rosenbrock Euler scheme (ERE) used in [3];
•	a novel NonStandard scheme (NS) proposed in [2].

The code is implemented in Matlab ( version R2017b).

The development and the implementation of the model and the routine have been made possible thanks to the Innonetwork Project COHECO, funded by Regione Puglia, Italy.

The routine can be used under the conditions of CC-BY-NC 2.0. 

A full description of the model is available in:
Diele, F.,  Marangi, C.,  and Martiradonna, A.  Non-standarddiscrete RothC models for soil carbon dynamics. Axioms, submitted.
