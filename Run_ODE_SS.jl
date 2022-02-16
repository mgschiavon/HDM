## Running in julia terminal
	cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Heterodimer feedback\\HDM\\")
	iARG = (mm = "HDM",	# Label for motif file
			ex = "Ex01");	# Label for parameters file

# Load functions & system:
using DelimitedFiles
using DifferentialEquations
using Distributions
mm = include(string("Library\\Md_",iARG.mm,".jl"));	# ODE system
fn = include(string("Library\\FN_Fit.jl"));			# Functions

## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file);
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# System parameters (i.e. p structure -list of kinetic parameters)
include(string("InputFiles\\ARGS_",iARG.mm,"_Fit_",iARG.ex,".jl"))	# Fitting rules (i.e. mrw structure -list of metropolis random walk parameters; d structure -experimental data/conditions; mySS -specific function to calculate steady states)
pO = copy(p);

## Update parameters (if needed) and get steady states:
p[:eC1]=0.0049906;p[:eC2]=0.028351;p[:mY]=0.039917;p[:kY]=0.018864;p[:nY]=1;p[:aY]=0.0082908;
mySS(fn,mm,p,d)
