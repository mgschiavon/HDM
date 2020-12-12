mrw  = (pOp  = [:mE,:mP,:bP],	# Parameters to optimize
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
x = CSV.File("DATA_Fig2B_Mean.csv") |> Tables.matrix;
d = (Hi = x[1,2:end],		# Hormone (Pg) concentrations tested
	 Xl = x[5:end,1],		# Experiment labels
	 Xe = x[5:end,2:end]);	# YFP steady state measurements
# Adjust measurement units (1 a.u.= 0.4 nM):
d.Xe[:,:] *= 0.4;

# RULES:
function myMSE(fn,mm,p,d)
	Y = zeros(size(d.Xe));
	mse = 0;
	# Updating parameters & data according to the used construct:
	for i in 1:length(p[:eC])
		p[:eP] = p[:eC][i];
		# Calculate steady state for each hormone concentration:
		for h in 1:length(d.Hi)
			p[:hP] = d.Hi[h];
			pSynth(p,iSynTF_mu);
			# Calculate steady states:
			ss = fn.SS(mm.myODE, p, ones(length(mm.myODE.syms)), 1e-4, 0);
			Y[i,h] = ss[1];
		end
		# Compare to data & calculate MSE:
		mse += (fn.MSE(Y[i,:],d.Xe[i,:])/length(p[:eC]));
	end
	return [Y, mse]
end