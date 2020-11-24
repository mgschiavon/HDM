clear
load('MRW_SimpleModel_v02_h37AB.mat')

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
        ylim([0.005 2])
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
clear i fig C Hi

%% Figure -- Parameter effect
fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [0 0 16 8];
fig.Position = fig.PaperPosition;
Hi = log2(p.A_H); Hi(length(p.A_H)) = Hi(length(Hi)-1) - 1;
C = colormap('parula');
delta = 2.^[-1,1];
for j1 = 1:2
    for j2 = 1:length(f)
        p.(f(j2).par) = p.(f(j2).par)*delta(j1);
        y = FN_SS_SimpleModel_v02(p);
        
        subplot(2,length(f),j2+(6*(j1-1)))
        hold on;
        for i = 1:length(p.D_H)
            plot(Hi,y(i,:)/p.nM,'LineWidth',2,...
                'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
            errorbar(Hi,D(i,:),Dsd(i,:),'LineWidth',1,'LineStyle','none','Marker','o',...
                'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
        end
                xlabel('Pg [nM]','FontSize',12)
                ylabel('YFP [a.u.]','FontSize',12)
                title(cat(2,f(j2).par,' x ',num2str(delta(j1),2),' (MSE=',...
                    num2str(FN_FitError(p,M,D*p.nM)/(size(D,1)*size(D,2)),2),')'),...
                    'FontSize',12)
        %         legend('location','southeast','FontSize',12)
                xlim([min(Hi)-0.5 max(Hi)+0.5])
                ylim([0.005 2])
                set(gca,'YScale','log','YGrid','on',...
                    'YTick','',...
                    'XTick',Hi([length(p.A_H):-2:1]),...
                    'XTickLabel',p.A_H([length(p.A_H):-2:1]),...
                    'XTickLabelRotation',45)
                box on
                
        p.(f(j2).par) = p.(f(j2).par)/delta(j1);
    end
end
print(gcf,'RAW_FigModel01_ParEffect','-dpng','-r300')
clear i j1 j2 fig C Hi delta
        
%%
fig = figure();
fig.Units = 'inches';
fig.PaperPosition = [0 0 3 2.5];
fig.Position = fig.PaperPosition;
Hi = log2(p.A_H); Hi(length(p.A_H)) = Hi(length(Hi)-1) - 1;
C = colormap('parula');

j = 5; delta = 0.25;
p.(f(j).par) = p.(f(j).par)*delta;
y = FN_SS_SimpleModel_v02(p);

hold on;
for i = 1:length(p.D_H)
    plot(Hi,y(i,:)/p.nM,'LineWidth',2,...
        'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
    errorbar(Hi,D(i,:),Dsd(i,:),'LineWidth',1,'LineStyle','none','Marker','o',...
        'DisplayName',cat(2,'E2=',num2str(p.D_H(i))),'Color',C(i*50,:))
end
        xlabel('Pg [nM]','FontSize',12)
        ylabel('YFP [a.u.]','FontSize',12)
        title(cat(2,f(j).par,' x ',num2str(delta,2),' (MSE=',...
            num2str(FN_FitError(p,M,D*p.nM)/(size(D,1)*size(D,2)),2),')'),...
            'FontSize',12)
%         legend('location','southeast','FontSize',12)
        xlim([min(Hi)-0.5 max(Hi)+0.5])
        ylim([0.005 2])
        set(gca,'YScale','log','YGrid','on',...
            'YTick','',...
            'XTick',Hi([length(p.A_H):-2:1]),...
            'XTickLabel',p.A_H([length(p.A_H):-2:1]),...
            'XTickLabelRotation',45)
        box on

p.(f(j).par) = p.(f(j).par)/delta;