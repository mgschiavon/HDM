## ODE FITTING PIPELINE - Heterodimer model
#	Mariana Gómez-Schiavon
#	November, 2020
#		Julia v.1.1.1
#		Required libraries:
#			DifferentialEquations
#			ParameterizedFunctions
#			DelimitedFiles

## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles\\ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl"))	# Perturbation details
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters

# Load functions & parameters:
using DelimitedFiles
mm = include(string("Library\\Md_",iARG.mm,".jl"));
fn = include(string("Library\\FN_DYs.jl"));
pO = copy(p);

# Run analysis
p = copy(pO);
open(string("OUT_ExDyn_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
	writedlm(io, [vcat("FB","rho","time",[string(i) for i in mm.odeNF.syms])],'\t');
	rtol = 1e-12;
	uns = 0;
	ssR = ones(length(mm.odeFB.syms));
	soR = ones(length(mm.odeNF.syms));
	while(rtol >= 1e-24)
		# Reference steady state:
		ssR = fn.SS(mm.odeFB, p, ssR, rtol, uns);
		# Locally analogous system reference steady state:
		mm.localNF(p,ssR);
		soR = fn.SS(mm.odeNF, p, soR, rtol, uns);
		if(abs(mm.outFB(ssR) - mm.outNF(soR)) > 1e-4)
			rtol *= 1e-3;
			if(rtol < 1e-24)
				println("ERROR: Check NF system (reltol=",rtol*1e3,").")
				println(vcat(pert.p,i,[p[i] for i in mm.odeFB.params],mm.outFB(ssR),mm.outNF(soR)))
				#throw(DomainError("x-("))
				if(abs(mm.outFB(ssR) - mm.outNF(soR))/mm.outFB(ssR) > 0.01)
					flg1 = 0;
					println("SS results excluded!")
				end
			end
		else
			break
		end
	end
	# Feedback system:
	x = fn.Dyn(mm.odeFB, p, ssR, 500.0);
	for i in 1:length(x.t)
		writedlm(io, [vcat(1,p[iARG.pp],x.t[i],x.u[i],"NaN")],'\t');
	end
	p[pert.p] *= pert.d;
	x = fn.Dyn(mm.odeFB, p, last(x.u), 9500.0);
	for i in 1:length(x.t)
		writedlm(io, [vcat(1,p[iARG.pp],x.t[i]+500.0,x.u[i],"NaN")],'\t');
	end
	ssD = fn.SS(mm.odeFB, p, ssR, rtol, uns);
	writedlm(io, [vcat(1,p[iARG.pp],"Inf",ssD,"NaN")],'\t');
	p[pert.p] /= pert.d;
	# No-Feedback system:
	x = fn.Dyn(mm.odeNF, p, soR, 500.0);
	for i in 1:length(x.t)
		writedlm(io, [vcat(0,p[iARG.pp],x.t[i],x.u[i])],'\t');
	end
	p[pert.p] *= pert.d;
	x = fn.Dyn(mm.odeNF, p, last(x.u), 9500.0);
	for i in 1:length(x.t)
		writedlm(io, [vcat(0,p[iARG.pp],x.t[i]+500.0,x.u[i])],'\t');
	end
	soD = fn.SS(mm.odeNF, p, soR, rtol, uns);
	writedlm(io, [vcat(0,p[iARG.pp],"Inf",soD)],'\t');
	p[pert.p] /= pert.d;
end