# Heterodimer competition model
#   with DN (N) binding DBD (D) part

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeHD = @ode_def begin
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D) + (bM * DN) - (bP * D * N)
		dN  = mN - (g * N)                            + (bM * DN) - (bP * D * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dDN =   - (g * DN)                            - (bM * DN) + (bP * D * N)
		dY  = (mY * ((a * k) + AD)/(k + AD + D + DN)) - (g * Y)
	end g mA mD mN mY a n k eP eM bP bM;
end 