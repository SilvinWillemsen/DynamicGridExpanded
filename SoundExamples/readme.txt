This file shows the parameters used to create the sound files in this folder. 
Some parameters have been made time varying and "morph" the physics from one state to the next.
All simulations ran at a sample rate of fs = 44100 Hz, and the models are excited every half-second.

The number of states is provided for each file as well as the (fractional) number of intervals over the course of the simulation. 
Arrows indicate parameter values between states. 
Parameters influencing pitch, change in a way that the resulting change in fundamental frequency is as linear as possible. 
Parameters with only one value are unchanged during the simulation.

1D:

- stringInharmonicityChange.wav (2 states, N: 106.22 -> 60.42)
	Length (L): 1
	Density (rho): 7850
	Radius (r): 2e-4
	Tension (T): 157.83 -> 157.58
	Youngs Modulus (E): 2e11
	Freq. Indep. Damping (sigma0): 1
	Freq. Dep. Damping (sigma1): 0.005

- stringOctaveVibraMarimba.wav (4 states, N: 108.35 -> 54.17 -> 15.28 -> 21.56)
	Length (L): 1 -> 0.5 -> 0.5 -> 0.5
	Density (rho): 7850
	Radius (r): 5e-4
	Tension (T): 555 -> 555 -> 0 -> 0
	Youngs Modulus (E): 2e11 -> 2e11 -> 7e13 -> 1.75e13
	Freq. Indep. Damping (sigma0): 1 -> 1 -> 1 -> 1.5
	Freq. Dep. Damping (sigma1): 0.005 -> 0.005 -> 0.005 -> 0.05

2D: 

- thicknessChange.wav (2 states, Nx: 18.99 -> 6.01, Ny: 37.98 -> 12.01)
	Horizontal Length (Lx): 0.5	
	Vertical Length (Ly): 0.5
	Density (rho): 7850
	Thickness (H): 5e-3 -> 5e-2
	Youngs Modulus (E): 2e11
	Freq. Indep. Damping (sigma0): 1 
	Freq. Dep. Damping (sigma1): 0.005 -> 0.01

- lengthChange.wav (3 states, Nx: 10.68, Ny: 5.34 -> 10.68 -> 5.34 )
	Horizontal Length (Lx): 0.5
	Vertical Length (Ly): 0.25 -> 0.5 -> 0.25	
	Density (rho): 7850
	Thickness (H): 5e-3
	Youngs Modulus (E): 2e12
	Freq. Indep. Damping (sigma0): 0.5
	Freq. Dep. Damping (sigma1): 0.005
