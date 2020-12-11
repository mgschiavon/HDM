mrw  = (pOp  = [:mU,:mW,:eP],	# Parameters to optimize
		pMin = [-3,-3,-3],		# Minimum parameter value to explore (log10)
		pMax = [3,3,3],			# Maximum parameter value to explore (log10)
		runs = 1000,			# Number of optimization runs
		iter = 1000,			# Number of iterations per optimization run
		cov  = [0.1,0.1,0.1],	# Covariance to calculate parameter random walk
		M    = 10,				# "Mutation step size" for multiplicative random walk
		rnP0 = 0,				# Flag for random initial values of parameters to optimize
		temp = 0,				# Flag for simulated annealing (if 0, MRW)
		prtW = 0);				# Flag for printing each walk step
		
# Load data to compare:
using CSV
using DataFrames
d.Xe = CSV.File("DATA_Fig2B_Mean.csv") |> Tables.matrix;
d.Hi = d.Xe[1,2:end];
d.Xe = d.Xe[5:end,:];

# RULES:
function myMSE(fn,mm,p,d)
	mse = 0;
	# Updating parameters & data according to the used construct:
	for i in 1:length(p[:eC])
		p[:eP] = p[:eC][i];
		D = d.Xe[i,2:end];
		# Calculate steady state for each hormone concentration:
		Y = zeros(length(D));
		for h in 1:length(d.Hi)
			p[:hP] = d.Hi[h];
			pSynth(p,iSynTF_mu);
			# Calculate steady states:
			ss = fn.SS(mm.myODE, p, ones(length(mm.myODE.syms)), 1e-4, 0);
			Y[h] = ss[1];
		end
		# Compare to data & calculate MSE:
		mse += (fn.MSE(X,D)/length(p[:eC]));
	end
	return mse
end