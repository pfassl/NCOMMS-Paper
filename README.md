# NCOMMS-Paper
Contains the MATLAB code used for the publication in Nature Communications

The function Kardio_getdata.m is used to extract the data from the raw data .mat files. The extracted data contains unfiltered, AC filtered and the voltage trace.

The function Kardio_ResCap.m calculates access resistance (Rs), membrane resistance (Rm), cell capacitance (Cap), holding current (Ih) as well as the elicited HCN current (I_HCN). From these data, the current density is calculated by calculating I_HCN/Cap.

The function Kardio_AP.m analyzes the characteristics of elicited action potentials in cardiac myocytes. It returns the resting membrane potential (RMP), action potential amplitude (AP_Amp) and the action potential durations at 20%, 50%, and 90% repolarization (APD20, APD50, APD90).
