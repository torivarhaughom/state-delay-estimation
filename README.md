# Extension of Observer Design for Joint State and Input Delay Estimation to MIMO LTI Systems – MATLAB/Simulink Files

This repository contains the MATLAB scripts and Simulink model used in the master's thesis:

**"Extension of Observer Design for Joint State and Input Delay Estimation to MIMO LTI Systems"**

---

## Files Included

### Simulation

- **observer_sim.m**  
  Main script for simulating the switched LPV observer and evaluating performance across different scenarios.

### Experiment

- **Experiments.slx**  
  Simulink file for experiments with the DC motor–propeller system.

- **setup_experiments.m**  
  Initializes variables used in `Experiments.slx`.
  
- **plotting_experiments.m**  
  Generates plots of experimental results. Must be run after `Experiments.slx`.

## Notes

The simulation script and experiment files correspond to the results presented in Chapters 4 and 5 of the thesis.  


## Known Typos and Errors in the Thesis
The following typos and errors were discovered in the thesis after submission:
- The identity matrix in the bottom-right block of the LMI in Equation (5.11) was incorrectly written as `I_{2×2}`. It should be `I_{1×1}`.
- A plotting error in the file `setup_experiments.m` caused the same state estimate to be plotted multiple times under different labels in Chapter 5. This only affected the subplots showing the state estimates. The delay estimate subplots, performance indices, and all conclusions remain correct and unaffected.



