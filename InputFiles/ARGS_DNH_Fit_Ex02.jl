mrw  = (pOp  = [:bC1,:bC2],# Parameters to optimize
		pMin = [-5,-5],# Minimum parameter value to explore (log10)
		pMax = [5,5],# Maximum parameter value to explore (log10)
		runs = 4,	# Number of optimization runs
		iter = 20000,# Number of iterations per optimization run
		cov  = [0.1,0.1],# Covariance to calculate parameter random walk
		M    = 2,	# "Mutation step size" for multiplicative random walk
		rnP0 = 0,	# Flag for random initial values of parameters to optimize
		temp = 0,	# Flag for simulated annealing (if 0, MRW)
		prtW = 1);	# Flag for printing each walk step

# Load data to compare:
using CSV
using DataFrames
x = CSV.File("DATA_Fig2B_Mean.csv") |> Tables.matrix;
d = (Hi = x[1,2:end],		# Hormone (Pg) concentrations tested
	 Xl = x[2:end,1],		# Experiment labels
	 Xe = x[2:end,2:end]);	# YFP steady state measurements
# Adjust measurement units (1 a.u.= 0.4 nM):
d.Xe[:,:] *= 0.4;

# RULES:
function mySS(fn,mm,p,d)
	Y = zeros(size(d.Xe));
	# Updating parameters & data according to the used construct:
	for i in 1:3
		p[:bP] = p[:bC2];
		p[:eP] = [p[:eC4],p[:eC5],p[:eC6]][i];
		# Calculate steady state for each hormone concentration:
		for h in 1:length(d.Hi)
			p[:hP] = d.Hi[h];
			pSynth(p,fn.iSynTF_mu);
			# Calculate steady states:
			ss = fn.SS(mm.odeD, p, ones(length(mm.odeD.syms)), 1e-4);
			Y[i,h] = ss[1];
		end
	end
	for i in 1:3
		p[:bP] = p[:bC1];
		p[:eP] = [p[:eC1],p[:eC2],p[:eC3]][i];
		# Calculate steady state for each hormone concentration:
		for h in 1:length(d.Hi)
			p[:hP] = d.Hi[h];
			pSynth(p,fn.iSynTF_mu);
			# Calculate steady states:
			ss = fn.SS(mm.odeA, p, ones(length(mm.odeA.syms)), 1e-4);
			Y[i+3,h] = ss[1];
		end
	end
	return Y;
end
