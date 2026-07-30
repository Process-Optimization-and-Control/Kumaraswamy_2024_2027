# FOCAPO - CPC Conference

This folder contains code, results, and figures submitted for the FOCAPO - CPC Conference 2027. Here run a closed-loop simulation for an Nonlinear Model Predictive Controller (NMPC) used to regulate a $CO_2$ pipeline-injection network. The manipulated variables in the syste include the choke valves and pump speed. The NMPC internally uses a first-principles dynamic model for mass and momentum balances. However, it is linked to a dynamic multiphase flow simulator which represents the "plant" model. LedaFlow is triggered to run using javascript at user-defined intervals. Simplified state estimation and offset-free methods are employed here. The NMPC is tested with changes in reference tracking profiles for the well flow rates. Two cases are run here, with and without the offset free method. 

Code for extracting LedaFlow results, and running LedaFlow via javascript are based off Marie Sunde's master thesis. This article further improves the code and runs a simulation including several step-changes in the wellhead flow rate setpoints. 


