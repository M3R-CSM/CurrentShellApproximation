# CurrentShellApproximation
Calculates the Current-Shell approximation for a real solenoid.

The attached code calculates a current-shell approximation for a real solenoid and plots the percent error for the approximation as compared to a exact solution.  The current-shell approximation is exact for axially magnetized cylindrical magnets and the provided BFieldShell can be used for that case as well.

Please Cite: 
Husa P.L., Saunders B.D., Petruska A.J., Optimal Current Shell Approximation for Solenoids of Rectangular
Cross-Section, IEEE Trans Magnetics, In Review

The params_CC function calculates the optimal approximation using two current shells.  Simillarly, the other params functions calculate teh geometries for combinations of current shells (C) and current rings (R). The plot_approximations makes a plot of the error.  Also provided are BFieldRing, BFieldCurrentShell, and BFieldSolenoid that calculate the exact fields for a ring, shell, and real solenoid. These also provide gradient terms.

Note: The BFieldSolenoid funciton is based on the paper: Conway, J.T., "Trigonometric Integrals for the Magnetic Field of the Coil of Rectangular Cross Section", IEEE Trans Magn, Vol 42/5, 2006.  It provides an explicit solution to the field of a rectangular-crossection solenoid, but requires 10x as much computation as the approximation proposed. 
