#################################
#Test case: Cavity

#################################
#Description: Simulate the evloution of a small droplet in a 3d lid-driven cavity with both DNS and LPP method. 

#################################
#  Input Parameters: (Standard Unit)
#  Ambient fluid: rho1 = 1    mu1 = 1e-5
#  Droplet     :  rho2 = 100  mu2 = 1e-3
#  Surface Tenstion:   sigma = 0.0001
#  Lid velocity:       U0 = 0.5
#  Domain dimension:   L  = 3.2e-4
#  Droplet size:       d  = 2.0e-5

#################################
#  Dimensionless Parameters: 
#  Reynolds number of the cavity flow: ReL = U0*L*rho1/mu1 = 16
#  Stokes number of the droplet:       St  = taup/tauf=(rho2*d^2/(18mu1))/(L/U0) = 0.34
#  Estimate of relative velocity:      |Urel| = taup*FluidAcceleration = St*U0 = 0.17
#  Note: Equilibrium Eulerian Approximation is used to estimate relative velocity (Ferry & Balachandar, IJMF 2001)
#  Weber number of droplet:            We  = rho1*|Urel|^2*d/sigma = 6e-3
#  Reynolds number of droplet:         Re  = rho1*|Urel|*d/mu1 = 0.34

#################################
#  Simulation
#  Command: ./run.sh
#  Compuation time: 2.5hrs with 8 processors
