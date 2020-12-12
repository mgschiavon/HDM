%%   [HDM] Heterodimer models                %%
% ------------------------------------------- %
% FUNCTION: Finding steady state solution     %
%           for the simple competition model  %

% Created by Mariana Gómez-Schiavon
% November 2020

%   See also RUN_FitMRW.m

function Y = FN_SS_SimpleModel(p)
    OHFn = @(x,n,k,a) [a+((1-a)*(x.^n)./((x.^n)+(k^n)))];
    % A synthesis function using expanded Hill model:
    fA = zeros(length(p.D_H),length(p.A_H));
    for h = 1:length(p.A_H)
        Xa = roots([1,-(p.A_H(h)+p.A_XT+p.A_KX),p.A_H(h)*p.A_XT]);
        Xa = Xa([Xa<p.A_XT]);
        if(length(Xa)~=1)
            'error -- multiple solutions'
            Xa = NaN;
        end
        fA(:,h) = p.A_m * OHFn(Xa+(p.A_b*(p.A_XT-Xa)),p.A_n,p.A_K,p.A_a);
    end    
    % D synthesis function using expanded Hill model:
    fD = zeros(length(p.D_H),length(p.A_H));
    for h = 1:length(p.D_H)
        Xa = roots([1,-(p.D_H(h)+p.D_XT+p.D_KX),p.D_H(h)*p.D_XT]);
        Xa = Xa([Xa<p.D_XT]);
        if(length(Xa)~=1)
            'error -- multiple solutions'
            Xa = NaN;
        end
        fD(h,:) = p.D_m * OHFn(Xa+(p.D_b*(p.D_XT-Xa)),p.D_n,p.D_K,p.D_a);
    end
    clear h Xa
    % Steady states
    A  = ((-p.eM*p.g)-(p.g^2)+(p.eP*fA)-(p.eP*fD)+sqrt((4*p.eP*p.g*((p.eM*fA)+(p.g*fA)))+(((-p.eM*p.g)-(p.g^2)+(p.eP*fA)-(p.eP*fD)).^2)))/(2*p.eP*p.g);
    D  = (p.eM+p.g)*fD./(p.g*(p.eM+(A*p.eP)+p.g));
    AD = A.*D*p.eP/(p.eM+p.g);
    Y  = (p.Y_m*((p.Y_a*p.Y_k)+AD)./(p.Y_k+AD+D))/p.g;
end
