# HDM
Heterodimer models

Code associated to:

Nguyen, Taylor H., Galen Dods, Mariana Gómez-Schiavon, Muziyue Wu, Zibo Chen, Ryan Kibler, David Baker, Hana El-Samad, and Andrew H. Ng. “Competitive Displacement of De Novo Designed HeteroDimers Can Reversibly Control Protein–Protein Interactions and Implement Feedback in Synthetic Circuits.” GEN Biotechnology 1, no. 1 (February 2022): 91–100. https://doi.org/10.1089/genbio.2021.0011.

## General fitting pipeline

1. [Run_ODE_Fit.txt](https://github.com/mgschiavon/HDM/blob/main/Run_ODE_Fit.txt) : Running instructions for the particular system (labels for ODE system file, `myODElabel`, and for parameters files, `myEXlabel`).

1. Library/Md_myODElabel.jl (e.g. [Library/Md_DN_A.jl](https://github.com/mgschiavon/HDM/blob/main/Library/Md_DN_A.jl)) : File with the ODE system to work with. It is provided by the user, and always named `Md_[LABEL].jl` and located in the `Library\` folder. Then, the identifying label must be specified in the running instructions.

1. InputFiles/ARGS_myODElabel_Par_myEXlabel.jl (e.g. [InputFiles/ARGS_DN_A_Par_Ex01.jl](https://github.com/mgschiavon/HDM/blob/main/InputFiles/ARGS_DN_A_Par_Ex01.jl)) : File with the kinetic parameters associated to the system, and if needed functions required to calculate such parameters. It is provided by the user, and always named `ARGS_[LABEL]_Par_[LABEL].jl` and located in the `InputFiles\` folder. Then, the identifying labels must be specified in the running instructions, and the match the labels for the other associated files.

1. InputFiles/ARGS_myODElabel_Fit_myEXlabel.jl (e.g. [InputFiles/ARGS_DN_A_Fit_Ex01.jl](https://github.com/mgschiavon/HDM/blob/main/InputFiles/ARGS_DN_A_Fit_Ex01.jl)) : File with the parameters and conditions for the fitting process, including the experimental data to compare and the relationship between data and parameters. It is provided by the user, and always named `ARGS_[LABEL]_Fit_[LABEL].jl` and located in the `InputFiles\` folder. Then, the identifying labels must be specified in the running instructions, and the match the labels for the other associated files.

2. [ODE_Fit.jl](https://github.com/mgschiavon/HDM/blob/main/ODE_Fit.jl) : Main pipeline file. *It doesn't need to be modified by the user.*

2. [Library/FN_Fit.jl](https://github.com/mgschiavon/HDM/blob/main/Library/FN_Fit.jl) : Functions library file. *It doesn't need to be modified by the user.*

3. OUT_Fit_myODElabel_myEXlabel.txt (e.g. [OUT_Fit_DN_A_Ex01.txt](https://github.com/mgschiavon/HDM/blob/main/OUT_Fit_DN_A_Ex01.txt)) : Output file produced by the pipeline. First column corresponds to the run number, second column the iteration number for each run, third column the Mean Squared Error (MSE) of the run, and then each column one of the fit parameters. Always named `OUT_Fit_[LABEL]_[LABEL].jl` with the identifying labels specified in the running instructions by the user.
