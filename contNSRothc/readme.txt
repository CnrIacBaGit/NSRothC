contNSRothC – continuous NonStandard RothC models in Matlab

Angela Martiradonna  - Department of Mathematics, University of Bari, Italy -
Institute for Applied Mathematics (IAC), CNR, Bari, Italy. 
Mail: angela.martiradonna@uniba.it


contNSRothC is a MATLAB routine (version R2017b) for implementig non standard RothC type methods for simaluting the dynamics of soil organic carbon (SOC) contained in the organic matter. 

Five compartments of the soil are considered: DPM (decomposable plant material), RPM (resistant plant material), BIO (microbial biomass), HUM (humified organic matter) and IOM (inert organic matter).

contNSRothC implements three numerical schemes for the SOC dynamics: 
•	the RothC discrete classical scheme in [1];
•	the Exponential Rosenbrock Euler scheme (ERE) used in [3]; 
•	a novel NonStandard scheme (NS) developed in [2].

The development and the implementation of the model and the routine have been made possible thanks to the Innonetwork Project COHECO, funded by Regione Puglia, Italy.
The routine can be used under the conditions of CC-BY-NC 2.0.

A full description of the model is available in: 
Diele, F., Marangi, C., Martiradonna, A. (2021). Non-standard discrete RothC models for soil carbon dynamics. Axioms, 10(2), 56.
