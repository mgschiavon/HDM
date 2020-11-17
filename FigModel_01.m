clear
load('C:\Users\mgsch\Desktop\MRW_SimpleModel_v02_h37AB.mat')

%% Parameters & model
[a b] = sort(minE);
for i = 1:length(f)
    p.(f(i).par) = bestP(b(4),i);
end
y = FN_SS_SimpleModel_v02(p);
clear a b i

%% Figure -- Data fitting
fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [0 0 3 2.5];
fig.Position = fig.PaperPosition;
Hi = log2(p.A_H); Hi(length(p.A_H)) = Hi(length(Hi)-1) - 1;
C = colormap('parula');
hold on;
for i = 1:length(p.D_H)
    plot(Hi,y(i,:)/p.nM,'LineWidth',2,...
        'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
    errorbar(Hi,D(i,:),Dsd(i,:),'LineWidth',1,'LineStyle','none','Marker','o',...
        'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
end
        xlabel('Pg [nM]','FontSize',12)
        ylabel('YFP [a.u.]','FontSize',12)
        title(cat(2,'Model fitting'),'FontSize',12)
%         legend('location','southeast','FontSize',12)
        xlim([min(Hi)-0.5 max(Hi)+0.5])
%         ylim([0.005 10])
        set(gca,'YScale','log','YGrid','on',...
            'YTick','',...
            'XTick',Hi([length(p.A_H):-2:1]),...
            'XTickLabel',p.A_H([length(p.A_H):-2:1]),...
            'XTickLabelRotation',45)
        box on
        annotation(fig,'textbox',[0.15 0.725 1 0.15],...
            'String',{cat(2,'MSE=',...
            num2str(FN_FitError(p,M,D*p.nM)/(size(D,1)*size(D,2)),2))},...
            'FitBoxToText','on','LineWidth',1,'BackgroundColor',[1 1 1],'FontSize',10);
        print(gcf,'RAW_FigModel01','-dpng','-r300')
