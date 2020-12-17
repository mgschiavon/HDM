# Heterodimer competition model

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dY  = (mY * ((aY * (kY^nY)) + (AD^nY))/((kY^nY) + (AD^nY) + (D^nY))) - (g * Y)
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
	end g mA mD mY aY nY kY eP eM;
end