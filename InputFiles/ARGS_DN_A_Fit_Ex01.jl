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
Xe = CSV.File("DATA_Fig2B_Mean.csv") |> Tables.matrix;
Hi = Xe[1,2:end];
Xe = Xe[5:end,:];

# RULES:
function myMSE(fn,mm,p,Xe,Hi)
	# Updating parameters according to the used construct:
	for i in 1:length(p[:eC])
		p[:eP] = p[:eC][i];
		pSynth(p,iSynTF_mu);
		D = Xe[i,2:end];
		# TO DO: Calculate steady states (using function in FN_HDM)
		# TO DO: Compare to data & calculate MSE (using function in FN_HDM)
	end
	return mse
end