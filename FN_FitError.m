%%   [HDM] Heterodimer models                %%
% ------------------------------------------- %
% FUNCTION: Calculate fitting error           %

% Created by Mariana Gómez-Schiavon
% November 2020 
% Modified from https://github.com/mgschiavon/GZM_InducibleTF/blob/master/FN_FitError.m

% FN_FitError : Calculate the sum of square errors between model's 
%               steady state (given a set of biophysical paramers and 
%               inducer concentration) and observed data.
%
%   Ef = FN_FitError(p,M,D)
%   p : Structure with the kinetic parameters & conditions
%   M : Transcriptional model to consider
%       ('SimpleModel')
%   D : Measured output (data) matrix
%
%   OUTPUT Ef : Sum of square errors
%
%   See also FN_SS_SimpleModel.m
%   See also RUN_FitMRW.m
%   See also FN_FitMRW.m

function Ef = FN_FitError(p,M,D)
    if(strcmp(M,'SimpleModel'))
        Ye = FN_SS_SimpleModel(p);
    else
        'ERROR: Transcriptional model not defined. Options: SimpleModel.'
    end
    Ef = sum(sum((log10(D)-log10(Ye)).^2)./(2*var(log10(D))));
end