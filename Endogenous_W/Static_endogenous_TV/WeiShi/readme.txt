Attached are the Matlab codes for the Monte Carlo simulations. There is one z variable, and it is used to construct
the weights in the y equation.

    Q_nT: the objective function
    bias_correction:
    [bc,sd] = bias_correction(theta,0)

This script provides the bias corrected estimator and its standard error based on the asymptotic variance formula, 
given the QML estimate theta.

    MonteCarlo, MonteCarlo_run: utility scripts for running the Monte Carlo simulations.


As you can see, depending on the exact specifications of the model, the code requires some tweaking. 
I have a plan to develop a Stata program that automates the process. Not yet started, maybe next year.

