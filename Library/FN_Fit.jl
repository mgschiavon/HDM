## ODE FITTING PIPELINE - Heterodimer model
#  ODE simulations & other functions
#	Mariana Gómez-Schiavon
#	November, 2020
#		Julia v.1.1.1
#		Required libraries:
#			DifferentialEquations

module fn
	# Required libraries
	using DifferentialEquations
	using Distributions
	
	# iSynTF transcription regulation function
	# INPUT: x  - Total iSynTF concentration
	#        h  - Hormone concentration
	#        kX - Dissociation constant hormone-iSynTF interaction
	#        b  - Fraction of the inactive iSynTF in the nucleus
	#        m  - Maximum synthesis rate
	#        a  - Basal expression of the output gene (in the absence of iSynTF)
	#        n  - Hill coefficie
	#        k  - Dissociation constant
	# OUPUT:      - Synthesis rate value
	function iSynTF_mu(x,h,kX,b,m,a,n,k)
		a2 = 1;
		a1 = -(h + x + kX);
		a0 = h * x;
		xa = (-a1 - sqrt((a1^2) - (4 * a2 * a0)))/(2 * a2);
		xo = xa + (b * (x - xa));
		return (m * (a + ((1 - a) * (xo^n)/((xo^n) + (k^n)))));
	end;

	# Steady state function for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        rtol - Tolerance value for ODE solver
	# OUPUT: ss   - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol)
		pV = [p[i] for i in syst.params];
		ss = try
				solve(ODEProblem(syst,x0,1e6,pV),alg_hint=[:stiff],reltol=rtol,callback=TerminateSteadyState());
			 catch
				println("WARNING: Error in steady state calculation. NaN values assigned instead.")
			 end
		 if(typeof(ss)==Nothing)
			return x0.+NaN
		 else
			return last(ss.u)
		end
	end;

	# ODE dynamics for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        tspan- Time to simulate
	# OUPUT: xD   - Vector of steady state of the ODE system
	function Dyn(syst, p, x0, tspan)
		pV = [p[i] for i in syst.params];
		xD = solve(ODEProblem(syst,x0,tspan,pV),AutoTsit5(Rosenbrock23()),reltol=1e-6);
		return xD
	end;

	# Mean Squared Error (MSE) between model and data
	# INPUT: Y - Vector of model values (i.e. steady state predictions)
	#        D - Vector of data values (i.e. steady state measurements)
	# OUPUT:   - Mean Squared Error
	function MSE(Y,D)
		return (sum((log10.(D) - log10.(Y)).^2)/std(log10.(D)))/length(Y);
	end;
end