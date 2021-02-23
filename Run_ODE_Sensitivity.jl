## Running in julia terminal
	cd("C:\\Users\\mgsch\\Dropbox (MGS-UCSF)\\MODEL - Heterodimer feedback\\HDM\\")
	iARG = (mm = "HDM",	# Label for motif file
			ex = "Ex01");	# Label for parameters file
			#include("Run_ODE_Sensitivity.jl")

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

## Update parameters (if needed):
#p[:eC1]=5.4378e-05;p[:eC2]=6.0801e-05;p[:eC3]=0.0001204;p[:bP]=0.1309;p[:kY]=0.3556;
pO = copy(p);

## Perturbe each parameter in the list and print output steady states:
for i in [:eC1,:mY,:aY,:kY]
    open(string("OUT_Sensitivity_",iARG.mm,"_",iARG.ex,"_",i,".txt"), "w") do io
		p = copy(pO);
        p[i] *=0.1;
        Y = mySS(fn,mm,p,d);
        writedlm(io, Y,'\t')
		p = copy(pO);
        Y = mySS(fn,mm,p,d);
        writedlm(io, Y,'\t')
        p = copy(pO);
        p[i] *=10;
        Y = mySS(fn,mm,p,d);
        writedlm(io, Y,'\t')
    end
end
