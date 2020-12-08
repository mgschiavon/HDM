# Heterodimer models
# ODE simulations & other functions

# Julia v.1.1.1

module fn
	# Required libraries
	using DifferentialEquations
	
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
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT: ss   - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol, uns)
		pV = [p[i] for i in syst.params];
		if (uns == 1)
			s0 = x0;
			t = 0.
			while (any([abs(x) for x in syst(s0, pV, 0.)] .> (rtol * s0)) && t<1e6)
				ss = solve(ODEProblem(syst,s0,10000.,pV),AutoTsit5(Rosenbrock23()),reltol=rtol);
				t += 100.;
				s0 = last(ss.u);
				#println("Time: ",t,", ",s0)
				end;
			return s0
		else
			ss = solve(SteadyStateProblem(syst, x0, pV), SSRootfind());
			# Verify this is the steady state:
			if any(ss.u .< 0)
				ss = solve(SteadyStateProblem(syst, x0, pV), DynamicSS(Rodas5(); reltol=rtol));
				return ss.u
			end
			if any([abs(x) for x in syst(ss.u, pV, 0.)] .> (rtol * ss.u))
				ss = solve(SteadyStateProblem(syst, ss.u, pV), DynamicSS(Rodas5(); reltol=rtol));
				return ss.u
			end
			return ss.u
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
end