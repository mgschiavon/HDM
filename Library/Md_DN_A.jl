# Heterodimer competition model
#   with DN (N) binding AD (A) part

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeHD = @ode_def begin
		dA  =  mC                          - (g * A) + (eM * AD) - (eP * A * D) + (bM * AN) - (bP * A * N)
		dD  = fM(xE,hE,kXE,bE,mE,aE,nE,kE) - (g * D) + (eM * AD) - (eP * A * D)
		dN  = fM(xP,hP,kXP,bP,mP,aP,nP,kP) - (g * N)                            + (bM * AN) - (bP * A * N)
		dAD =                              - (g * AD) - (eM * AD) + (eP * A * D)
		dAN =                              - (g * AN)                            - (bM * AN) + (bP * A * N)
		dY  = (mY * ((a * k) + AD)/(k + AD + D)) - (g * Y)
	end g mC xE hE kXE bE mE aE nE kE xP hP kXP bP mP aP nP kP mY a n k eP eM bP bM;
end 