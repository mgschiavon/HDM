# Heterodimer competition model
#   with DN (N) binding DBD (D) part

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dY  = (mY * ((aY * (kY^nY)) + (AD^nY))/((kY^nY) + (AD^nY) + (D^nY) + (DN^nY))) - (g * Y)
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D) + (bM * DN) - (bP * D * N)
		dN  = mN - (g * N)                            + (bM * DN) - (bP * D * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dDN =   - (g * DN)                            - (bM * DN) + (bP * D * N)
	end g mA mD mN mY aY nY kY eP eM bP bM;
end 