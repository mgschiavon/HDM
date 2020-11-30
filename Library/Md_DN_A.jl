# Heterodimer competition model
#   with DN (N) binding AD (A) part

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeHD = @ode_def begin
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D) + (bM * AN) - (bP * A * N)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D)
		dN  = mN - (g * N)                            + (bM * AN) - (bP * A * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dAN =   - (g * AN)                            - (bM * AN) + (bP * A * N)
		dY  = (mY * ((a * k) + AD)/(k + AD + D)) - (g * Y)
	end g mA mD mN mY a n k eP eM bP bM;
end 