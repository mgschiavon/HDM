mrw  = (pOp  = [:eC1,:eC2,:mY,:kY,:nY,:aY],# Parameters to optimize
		pMin = [-3,-3,-3,-3,0,-3],# Minimum parameter value to explore (log10)
		pMax = [3,3,3,3,2,0],# Maximum parameter value to explore (log10)
		runs = 4,	# Number of optimization runs
		iter = 20000,# Number of iterations per optimization run
		cov  = [0.1,0.1,0.1,0.1,0.05,0.01],# Covariance to calculate parameter random walk
		M    = 2,	# "Mutation step size" for multiplicative random walk
		rnP0 = 0,	# Flag for random initial values of parameters to optimize
		temp = 0,	# Flag for simulated annealing (if 0, MRW)
		prtW = 1);	# Flag for printing each walk step
		
# Load data to compare:
# NOTE: d.Xe is always the data matrix used to calculate the MSE.
using CSV
using DataFrames
x = CSV.File("DATA_Fig1C_Mean.csv") |> Tables.matrix;
d = (hP = x[1,3:end],				# Hormone (Pg) concentrations tested
	 hE = x[[2,3,4,6,7,8],2],		# Hormone (E2) concentrations tested
	 Xl = x[[2,3,4,6,7,8],1],		# Experiment labels
	 Xe = x[[2,3,4,6,7,8],3:end]);	# YFP steady state measurements
# Adjust measurement units (1 a.u.= 0.4 nM):
d.Xe[:,:] *= 0.4;

# RULES:
function mySS(fn,mm,p,d)
	Y = zeros(size(d.Xe));
	# Updating parameters & data according to the used construct:
	for i in 1:6
		p[:eP] = [p[:eC1],p[:eC1],p[:eC1],p[:eC2],p[:eC2],p[:eC2]][i];
		p[:hE] = d.hE[i];
		# Calculate steady state for each hormone concentration:
		for h in 1:length(d.hP)
			p[:hP] = d.hP[h];
			pSynth(p,fn.iSynTF_mu);
			# Calculate steady states:
			ss = fn.SS(mm.myODE, p, ones(length(mm.myODE.syms)), 1e-4);
			Y[i,h] = ss[1];
		end
	end
	return Y;
end