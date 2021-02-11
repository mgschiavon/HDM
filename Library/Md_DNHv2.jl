# Heterodimer competition model
#   with DN (N) binding AD (A) part
#   and with DN (N) binding DBD (D) part
#   v2: DN:DBD (DN) cannot bind the promoter

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system -- N:A vs A:D
	odeA = @ode_def begin
		dY  = (mY * ((aY * (kY^nY)) + (AD^nY))/((kY^nY) + (AD^nY) + (D^nY))) - (g * Y)
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D) + (bM * AN) - (bP * A * N)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D)
		dN  = mN - (g * N)                            + (bM * AN) - (bP * A * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dAN =   - (g * AN)                            - (bM * AN) + (bP * A * N)
	end g mA mD mN mY aY nY kY eP eM bP bM;

	# ODE system -- N:D vs A:D
	odeD = @ode_def begin
		dY  = (mY * ((aY * (kY^nY)) + (AD^nY))/((kY^nY) + (AD^nY) + (D^nY))) - (g * Y)
		dA  = mA - (g * A) + (eM * AD) - (eP * A * D)
		dD  = mD - (g * D) + (eM * AD) - (eP * A * D) + (bM * DN) - (bP * D * N)
		dN  = mN - (g * N)                            + (bM * DN) - (bP * D * N)
		dAD =   - (g * AD) - (eM * AD) + (eP * A * D)
		dDN =   - (g * DN)                            - (bM * DN) + (bP * D * N)
	end g mA mD mN mY aY nY kY eP eM bP bM;
end
