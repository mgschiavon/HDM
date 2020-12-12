# Heterodimer competition model
#   with DN (N) binding AD (A) part

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dY  = (mY * ((aY * (kY^nY)) + (AD^nY))/((kY^nY) + (AD^nY) + (D^nY))) - (g * Y)
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D) + (bM * AN) - (bP * A * N)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D)
		dN  = mN - (g * N)                            + (bM * AN) - (bP * A * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dAN =   - (g * AN)                            - (bM * AN) + (bP * A * N)
	end g mA mD mN mY aY nY kY eP eM bP bM;
end 