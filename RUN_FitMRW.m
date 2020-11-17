%%   [HDM] Heterodimer models                %%
% ------------------------------------------- %
% RUN: Find parameters that best fit the data %

% Created by Mariana Gómez-Schiavon
% November 2020

%   See also FN_SS_SimpleModel.m

clear

%% Inputs & conditions
% Transcriptional model:
M = 'SimpleModel_v02';      % Options: 'SimpleModel'
% Data
ExID = 'h37AB';	% Experiment (TF) to consider
    load('DATA_HDM.mat','X');
    p.A_H = X.(ExID).Pg;
    p.D_H = X.(ExID).E2;
    D   = X.(ExID).mu;
    Dsd = X.(ExID).sd;
    clear X
S = [1:1000];   % Random number seed(s) (1 per run)
I = 20000;      % Iterations per fitting run
printAll = 0;   % Flag for printing full random walk
% Kinetic parameters:
    p.A_XT = 0.3042;    % [nM] for pPAB1~pRPL18b (with 1 a.u. = 0.4 nM)
    p.A_KX = 32.5;      % [nM]
    p.A_b  = 0.0317;
    p.A_n  = 2.64;
    p.A_K  = 0.0683;    % [nM]
    p.A_m  = 0.00944;   % [nM/min]
    p.A_a  = 0.0153;
    
    p.D_XT = 0.0607;    % [nM] for pPOP6~pRNR2 (with 1 a.u. = 0.4 nM)
    p.D_KX = 139;       % [nM]
    p.D_b  = 0.000102;
    p.D_n  = 1.4;
    p.D_K  = 0.00731;   % [nM]
    p.D_m  = 0.0159;    % [nM/min]
    p.D_a  = 0.00483;
    
    p.g  = 0.01;        % [1/min]
    p.nM = 0.4;         % Assume the measured fluorescence arbitrary units ([a.u.]) are proportional to the molecule concentration ([nM]), with 1 a.u. = 0.4 nM.
    
    p.eP = 1;           % [1/(nM min)]
    p.eM = 1;           % [1/min]
    
    p.Y_m = 0.1;        % [nM/min]
    p.Y_a  = 0.01;
    p.Y_n  = 1;
    p.Y_k  = 1;         % [nM]
    
% Parameters to fit:
    i = 0;
    i = i + 1;
    f(i).par = 'eP';
    f(i).cov = 0.1;
    f(i).lim = [1e-2,10000];
    i = i + 1;
    f(i).par = 'Y_m';
    f(i).cov = 0.1;
    f(i).lim = [2e-6,2];
    i = i + 1;
    f(i).par = 'Y_a';
    f(i).cov = 0.1;
    f(i).lim = [2e-7,0.2];
    i = i + 1;
    f(i).par = 'Y_n';
    f(i).cov = 0.1;
    f(i).lim = [1e-4,100];
    i = i + 1;
    f(i).par = 'Y_k';
    f(i).cov = 0.1;
    f(i).lim = [1e-4,100];
    i = i + 1;
    f(i).par = 'A_m';
    f(i).cov = 0.1;
    f(i).lim = [2e-6,2]*10;
    clear i

%% Run fitting:
bestP = zeros(length(S),length(f));
minE  = zeros(length(S),1);
for s = S
    cat(2,'Running seed #',num2str(s))
    [bP,mE] = FN_FitMRW(p,M,D,s,f,I,ExID,printAll)
    bestP(s,:) = bP;
    minE(s) = mE;
    if(mod(s,10)==0)
        save(cat(2,'TEMP_MRW_',M,'_',ExID,'.mat'));
    end
end
clear s bP mE ans
save(cat(2,'MRW_',M,'_',ExID,'.mat'));
delete(cat(2,'TEMP_MRW_',M,'_',ExID,'.mat'));
% load(cat(2,'MRW_',M,'_',ExID,'.mat'));

%% Figures
if(printAll)
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 0 18 10];
    fig.Position = fig.PaperPosition;
    C = colormap('parula');
    for s = S
        load(cat(2,'MRW_',M,'_',ExID,'_s',num2str(s),'.mat'),'mrw');
        for i = 1:size(mrw.P,2)
            subplot(2,4,i)
            hold on;
            plot(mrw.P(:,i),'LineWidth',2,'Color',C(s*6,:))
                xlabel('Iterations')
                ylabel(f(i).par)
                xlim([0,I])
                set(gca,'YScale','log')
                box on
        end
        subplot(2,4,8)
        hold on;
        plot(mrw.e,'LineWidth',2,'Color',C(s*6,:))
            xlabel('Iterations')
            ylabel('Error')
            xlim([0,I])
            set(gca,'YScale','log')
            box on
    end
    clear s i a b
    print(gcf,cat(2,'MRW_',M,'_',ExID,'_Runs.png'),'-dpng','-r300')
else
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 5 8];
    fig.Position = fig.PaperPosition;
    
    hist(log10(minE),20)
        xlabel('log_{10}(min(error))')
        ylabel('Count')
        title(cat(2,'Experiment: ',ExID))
        box on
        
        axes('Position',[0.25 0.5 0.3 0.3])
            hist(log10(minE([minE<10])),20)
            xlabel('log_{10}(min(error))')
            ylabel('Count')
            title(cat(2,'min(min(error)) = ',num2str(min(minE))))
            box on
        print(gcf,cat(2,'MRW_',M,'_',ExID,'_minE.png'),'-dpng','-r300')
        
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 8];
    fig.Position = fig.PaperPosition;
    [a b] = sort(minE);
    C = [0 0 0;
    1 0.6 0;
    0.85 0.85 0.85];
    for i = 1:length(f)
        subplot(2,ceil(length(f)/2),i)
        hold on;
        scatter(minE(b([101:1000]))/(size(D,1)*size(D,2)),...
            bestP(b([101:1000]),i),25,...
            'MarkerFaceColor',C(3,:),...
            'MarkerEdgeColor',C(3,:)-0.1)
        scatter(minE(b([1:100]))/(size(D,1)*size(D,2)),...
            bestP(b([1:100]),i),25,...
            'MarkerFaceColor',C(2,:),...
            'MarkerEdgeColor',C(2,:)-[0.1 0.1 0])
        for ii = 1
            scatter(a(ii)/(size(D,1)*size(D,2)),...
                bestP(b(ii),i),25,[0 0 0],'filled')
        end
                ylabel(f(i).par)
                xlabel('MSE')
%                 title(cat(2,'Med(',f(i).par,')=',...
%                     num2str(median(bestP(b([1:1000]),i)),4)))
                title(cat(2,'Med_{Q10}(',f(i).par,')=',...
                    num2str(median(bestP(b([1:100]),i)),4)))
                set(gca,'XScale','log','YScale','log',...
                    'XMinorGrid','off','YMinorGrid','off',...
                    'XTick',10.^[-6:1:6])
                box on
                grid on
    end
        print(gcf,cat(2,'MRW_',M,'_',ExID,'_minExPar.png'),'-dpng','-r300')
end
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 8];
    fig.Position = fig.PaperPosition;
    C = colormap('jet');
    Hi = log2(p.A_H); Hi(length(p.A_H)) = Hi(length(Hi)-1) - 1;
    [a b] = sort(minE);
    for ii = 1:10
        for i = 1:length(f)
            p.(f(i).par) = bestP(b(ii),i);
        end
        if(strcmp(M,'SimpleModel'))
            Ye = FN_SS_SimpleModel(p);
        elseif(strcmp(M,'SimpleModel_v02'))
            Ye = FN_SS_SimpleModel_v02(p);
        else
            'ERROR: Transcriptional model not defined. Options: SimpleHill, HillxBasal, Mechanistic, Allosteric.'
        end
        subplot(2,5,ii)
        hold on
        for i = 1:length(p.D_H)
            plot(Hi,Ye(i,:),'Color',C(i*60,:),'LineWidth',2)
            plot(Hi,D(i,:)*p.nM,'Color',C(i*60,:),'LineStyle','none','Marker','o')
        end
                xlabel('Pg [nM]')
                ylabel('YFP [nM]')
                title(cat(2,'Error = ',num2str(FN_FitError(p,M,D*p.nM))))
                xlim([min(Hi)-0.5 max(Hi)+0.5])
                set(gca,'YScale','log','YGrid','on',...
                    'XTick',Hi([length(p.A_H):-3:1]),'XTickLabel',p.A_H([length(p.A_H):-3:1]),...
                    'XTickLabelRotation',45)
                box on
    end
    clear Hi s i 
    print(gcf,cat(2,'MRW_',M,'_',ExID,'_BestFits.png'),'-dpng','-r300')

%% END