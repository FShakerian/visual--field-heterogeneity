clc
clear
close all


%% after Loadind periphery and fovea data
clear perf tim n_trial

%loc = 0; xdata = data_train1; % for fovea
loc = (0:7).*45; xdata = data_ses1; % for periphery

vis_sig1 = [0 0;1 39;40 59;60 79;80 100];

n_s = 0
for n_subject = 1:length(xdata)
    n_s = n_s+1;
    Res = xdata{n_subject};
    vis_sig = unique(abs(sort(Res(:,9))));
    total_t (n_s) = length(Res);
    
    for LL = 1:length(loc)
        
        for vs = 1:size(vis_sig1,1)
            
            ix_4leg = ismember(Res(:,9),-vis_sig1(vs,2):-vis_sig1(vs,1)) & Res(:,8)== loc(LL);%  & Res(:,11)< 3;
            ix_Ch   = ismember(Res(:,9),vis_sig1(vs,1):vis_sig1(vs,2))   & Res(:,8)== loc(LL);%  & Res(:,11)< 3;
            
            n_trial(n_s,LL,vs,:) = [sum(ix_4leg) sum(ix_Ch)];
            perf(n_s,LL,vs,1)    = sum(Res(ix_4leg,12))./sum(ix_4leg);
            perf(n_s,LL,vs,2)    = sum(Res(ix_Ch,12))./sum(ix_Ch);
            std_perf(n_s,LL,vs)  = std(Res(ix_Ch,12));
            
            tim(n_s,LL,vs,1) = nanmean(Res(ix_4leg,11));
            tim(n_s,LL,vs,2) = nanmean(Res(ix_Ch,11));
        end
    end
end

if length(loc) == 1
    perft = perf; timt = tim; ntrialt = n_trial;
    
else
    perfp = perf;
    timp = tim;
    ntrialp = n_trial;
end

clear perf tim n_trial

%%%%% noise is chair report%%%%%%%%%%%%%


%% fig 2.a fovea and lower and upper
clear var_L var_C
vs = 1:5;
var_F = squeeze(perft(:,1,vs,:)); 
var_L = (squeeze(perfp(:,:,vs,1)) - (1-repmat(perfp(:,:,1,1),[1 1 5]))); 
var_C = (squeeze(perfp(:,:,vs,2)) - repmat(perfp(:,:,1,1),[1 1 5])); 
locU = [2 3 4];
locL = [6 7 8];
colC = [0 0.355 0.7410];
colL = [0.8500 0.3250 0.0980];

%% figure 2B barplot fovea-upper-lower
vs = 5;
n_var1 = sum(~isnan(perft(:,1,5,1)));
n_var2 = sum(~isnan(perfp(:,1,1,1)));

mean_x(1,:) = squeeze(nanmean(perft(:,:,vs,:),1));
mean_x(2,:) = squeeze(nanmean(nanmean(perfp(4:end,locU,vs,:),2),1));
mean_x(3,:) = squeeze(nanmean(nanmean(perfp(4:end,locL,vs,:),2),1));

sem_x(1,:) = squeeze(nanstd(perft(:,:,vs,:),[],1))./sqrt(n_var1);
sem_x(2,:) = squeeze(nanstd(nanmean(perfp(4:end,locU,vs,:),2),[],1))./sqrt(n_var2);
sem_x(3,:) = squeeze(nanstd(nanmean(perfp(4:end,locL,vs,:),2),[],1))./sqrt(n_var2);

f = figure('Units','Centimeters','Position',[5 5 15 15]);

h(1) = bar(1:4:9,squeeze(mean_x(:,1))',0.3,'FaceColor',colL,'EdgeColor',colL);
hold on
errorbar(1:4:9,squeeze(mean_x(:,1)),squeeze(sem_x(:,1)),'LineWidth',1.5,'LineStyle', 'none','color',[0 0 0])
hold on
h(2) =  bar(2.25:4:10.5,squeeze(mean_x(:,2)),0.3,'FaceColor',colC,'EdgeColor',colC);
hold on
errorbar(2.25:4:10.5,squeeze(mean_x(:,2)),squeeze(sem_x(:,2)),'LineWidth',1.5,'LineStyle', 'none','color',[0 0 0])

legend(h,{'4leg','chair'});
ylim([0.5 1])
set(gca,'fontsize',16,'XTick',1.5:4:9.5,...
    'XTickLabel',{'Fovea','Upper','Lower'},'XTickLabelRotation',45);
xlabel('Visual field');
ylabel('Performance');

[pf,~]= signrank(perft(1:end,1,vs,1),perft(1:end,1,vs,2));
[pu,~] = signrank(squeeze(nanmean(perfp(1:end,locU,vs,1),2)),squeeze(nanmean(perfp(1:end,locU,vs,2),2)));
[pl,~] = signrank(squeeze(nanmean(perfp(1:end,locL,vs,1),2)),squeeze(nanmean(perfp(1:end,locL,vs,2),2)));
disp (num2str([pf, pu, pl]));

f_diff = perft(:,1,vs,2)-perft(:,1,vs,1);
nanmean(f_diff)
nanstd(f_diff)./sqrt(n_var1)

pu_diff = nanmean(perfp(1:end,locL,vs,2),2)-nanmean(perfp(1:end,locU,vs,2),2);
nanmean(pu_diff)
nanstd(pu_diff)./sqrt(n_var2)

pl_diff = nanmean(perfp(1:end,locL,vs,2),2)-nanmean(perfp(1:end,locL,vs,1),2);
nanmean(pl_diff)
nanstd(pl_diff)./sqrt(n_var2)
%% figure 2-c/ polar plot with min-max range
clear ypolar1 ypolar2
vv=0;
for vs= 2:5
    vv=vv+1;
    ypolarL(vv,:) = nanmean(perfp(1:end,:,vs,1),1);
    ypolarC(vv,:) = nanmean(perfp(1:end,:,vs,2),1);
end
% 
x_boudL = [nanmean(ypolarL,1);min(ypolarL,[],1);max(ypolarL,[],1)];
x_boudC = [nanmean(ypolarC,1);min(ypolarC,[],1);max(ypolarC,[],1)];

% x_boudL = [nanmean(ypolarL,1);nanmean(ypolarL,1)-std(ypolarL,[],1);nanmean(ypolarL,1)+std(ypolarL,[],1)];
% x_boudC = [nanmean(ypolarC,1);nanmean(ypolarC,1)-std(ypolarC,[],1);nanmean(ypolarC,1)+std(ypolarC,[],1)];

figure;
opt_area.err        = 'std';
opt_area.FaceAlpha  = 0.5;
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = '.';
opt_lines.Labels    = true;
opt_lines.Legend    = {'4leg','Chair'};
opt_area.Color      = colC;
opt_lines.Color     = colC;

opt_axes.Background = 'w';
opt_axes.Labels     = {'0','45','90','135','180','225','270','315'};

polygonplot_v3(x_boudC1,opt_axes,opt_lines,opt_area);
hold on
opt_area.Color      = colL;
opt_lines.Color     = colL;

polygonplot_v3(x_boudL1,opt_axes,opt_lines,opt_area);

%% statistic figure 2-c
for ll = 1:8
    
    p_L = squeeze(perfp(:,ll,2:5,1));
    p_L1 = max(p_L,[],2)-min(p_L,[],2);
    %p_L = p_L(:);p_L(isnan(p_L))=[];
    
    p_C = squeeze(perfp(:,ll,2:5,2));
    p_C1 = max(p_C,[],2)-min(p_C,[],2);
    %p_C = p_C(:);p_C(isnan(p_C))=[];
    
    min_max_L(ll) = nanmean(p_L1);
    min_max_C(ll) = nanmean(p_C1);

    std_L(ll) = std(p_L1)/sqrt(16);
    std_C(ll) = std(p_C1)/sqrt(16);

    
    [p(ll),h] = signrank(p_L1,p_C1);
end
%pp = h.*8; % bonfrroni corrected
[hh, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05,'pdep','no');

%% figure 3/6/7-a,b/scatter plot
clear fig fig1 h

loc_field = [6 7 8;2 3 4;4 5 6;8 1 2];
loc_name = { 'Lower','Upper', 'Left', 'Right'};
%close all
vs =5
ll =1
xvar1 = perfp(:,:,vs,1);
xvar2 = perfp(:,:,vs,2);

xL = nanmean(xvar1(:,loc_field(ll,:)),2);
xL(isnan(xL))=[];


xL1 = nanmean(xvar1(:,loc_field(ll+1,:)),2);
xL1(isnan(xL1)) = [];

xC = nanmean(xvar2(:,loc_field(ll,:)),2);
xC(isnan(xC))=[];

xC1 = nanmean(xvar2(:,loc_field(ll+1,:)),2);
xC1(isnan(xC1)) = [];

disp([nanmean(xL) nanstd(xL)./sqrt(length(xL))]);
disp([nanmean(xL1) nanstd(xL1)./sqrt(length(xL1))]);
disp([nanmean(xC) nanstd(xC)./sqrt(length(xC))]);
disp([nanmean(xC1) nanstd(xC1)./sqrt(length(xC1))]);

[pvaL,h,stats] = signrank(xL,xL1);
[pvaC,~] = signrank(xC,xC1);
[pvaLo,~] = signrank(xL,xC);
[pvaU,~] = signrank(xL1,xC1);

disp(strcat('pval=',num2str([pvaL pvaC])));

f = figure('Units','Centimeters','Position',[5 5 15.5 15]);
 %fig(1) = scatter(xL, xL1,[],'fill','SizeData', 100 ,'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3]);
fig(1) = scatter(xL,xL1,[],'fill','SizeData', 100 ,'MarkerEdgeColor',colL,'MarkerFaceColor',colL);
hold on;
fig(2) = scatter(xC,xC1,[],'d','SizeData', 100,'MarkerEdgeColor',colC,'MarkerFaceColor',colC);

xmin = min([xL;xC;xL1;xC1])-0.1;
xmax = max([xL;xC;xL1;xC1])+0.1;
hold on;

set(gca,'fontsize',18,'xtick',0:0.2:2.2,'ytick',0:0.2:2.2)
xlabel([loc_name{ll} ' performance']);
ylabel([loc_name{ll+1} ' performance']);

line([0 1],[0 1],'Color','black','LineStyle','--','linewidth',1);

xlim([0.4 1]);
ylim([0.4 1]);
title({['visual signal ' num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2)) '%'];...
['p(A) = ' num2str(pvaL,'%.4f') ' -p(C) = ' num2str(pvaC,'%.4f')]});
%  title({['vs = ' num2str(vs) '/ r(4leg) = ' num2str(r1) '- pval= ' num2str(round(pval1*100)/100)];...
%      ['r(chair) = ' num2str(r2) '- pval= ' num2str(round(pval2*100)/100)]});
legend(fig,{'Body','Chair'})
%% figure 3-a,b/distribution
%xlab = {'L-U Perf.','L-U Perf.','Le-R Perf.'};
xlab = {'L-U RT.','L-U RT.','Le-R RT.'};

figure;
xhistL = xL- xL1;
xhistC = xC- xC1;

% fig1(1) = histogram(xhistL,'BinWidth',0.1,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3]);
fig1(1) = histogram(xhistL,'BinWidth',0.1,'FaceColor',colL,'EdgeColor',colL);
hold on
fig1(2) = histogram(xhistC,'BinWidth',0.1,'FaceColor',colC,'EdgeColor',colC);
%
xlabel(xlab{ll});
ylabel('Pop')
set(gca,'fontsize',16,'Box','off',...
    'YAxisLocation', 'right','YTick',...
    [4 8 12],...
    'YTickLabel',[4 8 12]);%max([fig1(1).Values,fig1(2).Values])

xlim([-0.6 0.65])

 ylim([0 max([fig1(1).Values,fig1(2).Values])+1])
% ylim([0 max(fig1(1).Values)+1])

disp([nanmean(xhistL) nanstd(xhistL)./sqrt(length(xhistL))]);
disp([nanmean(xhistC) nanstd(xhistC)./sqrt(length(xhistC))]);

hold on
% line([nanmean(xhistL) nanmean(xhistL)],[1 2], 'Color', [0.3 0.3 0.3]);
line([nanmean(xhistL) nanmean(xhistL)],[9 10], 'Color', colL);
hold on
line([nanmean(xhistC) nanmean(xhistC)],[7 8], 'Color', colC);

[~,p_ks] = kstest2(xhistL,xhistC)
% title(['p=' num2str(p_ks)]);
% title({['vs = ' num2str(vis_sig1(vs,:)) '/ 4leg= ' num2str(nanmean(xhist1)) '- pval= ' num2str(round(pval1*100)/100)];...
%     ['chair= ' num2str(nanmean(xhist2)) '- pval= ' num2str(round(pval2*100)/100)]});
% legend({'4leg','chair'},'Location','northwest');
title({['visual signal ' num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2)) '%'];...
['p(A) = ' num2str(p_ks,'%.4f')]});

%% figure 3-a, b/ cumulative figure
figure
[f1,h1]=ecdf(xhistL);
plot(h1,f1,'Color', colL,'LineWidth',2)
hold on
[f1,h1]=ecdf(xhistC);
plot(h1,f1,'Color', colC,'LineWidth',2)
xlabel(xlab{ll});
ylabel('F(x)')
set(gca,'fontsize',24,'xtick',[-0.3,0,0.3]);
xlim([-.3 .3]);
%% figure 3 c - violin plot

ll_L = [1 2 3;3 4 5;5 6 7; 7 8 1];

%ll_L = [1;2;3;4;5;6;7;8];

x_data_m =[];
for ll1= 1:4
    x_data_m(:,ll1,:,1) = nanmean(perfp(:,ll_L(ll1,:),:,1),2);
    x_data_m(:,ll1,:,2) = nanmean(perfp(:,ll_L(ll1,:),:,2),2);

end

xtic_loc = {'Upper-right','Upper-left','Lower-left','Lower-right'};
%xtic_loc = {'0', '45', '90', '135', '180', '215', '270', '315'};
C = {colL,colC};

for vs = 5
    %figure
    for catt = 1
        x_viol=[];
        for ll = 1:4
            x_viol(:,ll) = squeeze(x_data_m(:,ll,vs,catt));
            hold on
            Violin({x_viol(:,ll)},ll,'ViolinColor',{C{catt}});
        end
        
        psig = signrank_test(1:ll,x_viol);
        [groups,pvalu] = sigstar_group(1:ll,psig);
        if ~isempty(groups); H = sigstar(groups,pvalu,0,C{catt}); end
    end
    
    ylim([0.4 1.5]);
    set(gca,'fontsize',14,'XTick',1:8,...
        'XTickLabel',xtic_loc,'XTickLabelRotation',45,...
        'YTick',0.2:0.2:1,'YTickLabel',0.2:0.2:1);
    xlabel('Location');
    ylabel('mean Perf.');
    
    %title(['visual signal= ' num2str(vis_sig1(vs,:))]);
end

x_v = squeeze(x_data_m(:,:,vs,2));
[nanmean(x_v,1); nanstd(x_v,[],1)./sqrt(16)]

%% figure 4-a/ polar plot 2 category in different range of visual signal
close all
for vs= 1
    figure;
%     ypolarL = nanmean(var_L(:,:,vs),1);
%     ypolarC = nanmean(var_C(:,:,vs),1);
    ypolarL = nanmean(1-perfp(:,:,vs,1),1);
    ypolarC = nanmean(perfp(:,:,vs,2),1);
%     ypolarL = nanmean(perecision_an(:,:,vs),1);
% ypolarC = nanmean(perecision_in(:,:,vs),1);

rlim1 = min(min([ypolarL,ypolarC])); rlim2 = max(max([ypolarL,ypolarC])); 
     
    %figure
    %polarplot(deg2rad(0:45:360),[ypolarL,ypolarl(1)],'.-','MarkerSize',30,'Color',[0.3 0.3 0.3]);
    polarplot(deg2rad(0:45:360),[ypolarL,ypolarL(1)],'.-','MarkerSize',30,'LineWidth',1.5,'Color',colL);
    hold on
    polarplot(deg2rad(0:45:360),[ypolarC,ypolarC(1)],'.-','MarkerSize',30,'LineWidth',1.5,'Color',colC);
    hold on
    ax = gca;
    ax.ThetaTick = 0:45:315;
    ax.ThetaTickLabel = [];
    ax.RLim = [rlim1-0.05 rlim2+0.05];
    %ax.RLim = [-0.25 0.5];
    %ax.RTickLabel = [];
    title(['visual signal ' num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2)) '%']);
    set(gca,'fontsize',16);
end

%anova1(1-perfp(:,:,vs,1))

%% figure 4-b / animacy bias in different visual signal

x_diff = ((perfp(:,:,:,1)-perfp(:,:,:,2))./(perfp(:,:,:,1)+perfp(:,:,:,2))); % ABI with recall

x_dp = squeeze((nanmean(perfp(:,:,:,1),1)-nanmean(perfp(:,:,:,2),1))./ sqrt((nanvar(perfp(:,:,:,1),1)+nanvar(perfp(:,:,:,2),1))./2));


close all
xval =[];
semval =[];
for vs = 2:5
    xx = squeeze(x_diff(:,:,vs));
    xval   = nanmean(xx,1);
    semval   = nanstd(xx,[],1)./sqrt(sum(~isnan(xx(:,1))));
    
    xval = [xval,xval(1)];
    
    semval = nanstd(xx,[],1)./sqrt(size(xx,1));
    semval = [semval,semval(1)];
    
    figure
    polarplot(deg2rad(0:45:360),xval,'.-','MarkerSize',30,'LineWidth',1.5,'Color',[0.2 0.2 0.2]);
    ax = gca;
    ax.ThetaTick = 0:45:315;
    ax.ThetaTickLabel = [];
    ax.RLim = [min(xval)-0.03 max(xval)+0.03 ];
    set(gca,'fontsize',16);
    title(['visual signal ' num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2)) '%']);
end
%% statistical test figure 4-b

x_anova = [];
for vs = 5
    x_anova =  squeeze(x_diff(:,:,vs));
    
    x_anova1 = reshape(x_anova,16*8,1);
    
    [p,tb1,stats] = anova1(x_anova);
    %     figure
    %     c = multcompare(stats,'Estimate','row');
end

%% figure 5 colormap  of ABI

i=0;
for vs = 5:-1:2
    i=i+1;
    ytic_map{i} = [num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2))];
end

i=0;
for ll = 8:-1:1
    i=i+1;
    xtic_map{i} = num2str(loc(ll));
end

figure;
colormap('parula');
xmap = (squeeze(nanmean(x_diff,1)))';
imagesc(xmap(end:-1:2,end:-1:1));
c = colorbar; c.Label.String = 'ABI';
xlabel('Location(degree)');
ylabel('Visual signal(%)');
set(gca,'fontsize',16,'YTick',1:4,'YTickLabel',ytic_map,'XTick',1:8,'XTickLabel',xtic_map);
%% figure 5-b

figure;
xval   = squeeze(nanmean(nanmean(x_diff(:,:,2:5),3),1));
semval = squeeze(nanstd(nanmean(x_diff(:,:,2:5),3),[],1)./sqrt(size(x_diff,1)));

bar(1:length(xval),xval,0.8,'FaceColor',[0.5 0.5 0.5],...
    'EdgeColor',[0.5 0.5 0.5],'LineWidth',1);
hold on
errorbar(1:length(xval),xval,semval,...
    'LineWidth',1.5,'LineStyle', 'none','color',[0 0 0])

set(gca,'fontsize',16,'XTick',1:length(xval),'XTickLabel',0:45:315,'XTickLabelRotation',45);
xlabel('Location');

%%xlim([0 5])
ylabel('ABI');

%% figure 5-c

for vs = 2:5
    xtic_nam{vs-1} = [num2str(vis_sig1(vs,1)) '-' num2str(vis_sig1(vs,2))];
end

figure;
xval   = squeeze(nanmean(nanmean(x_diff(:,:,2:5),2),1));
semval = squeeze(nanstd(nanmean(x_diff(:,:,2:5),2),[],1)./sqrt(size(x_diff,1)));

bar(1:length(xval),xval,0.8,'FaceColor',[0.5 0.5 0.5],...
    'EdgeColor',[0.5 0.5 0.5],'LineWidth',1);
hold on
errorbar(1:length(xval),xval,semval,...
    'LineWidth',1.5,'LineStyle', 'none','color',[0 0 0])

set(gca,'fontsize',16,'XTick',1:length(xval),'XTickLabel',xtic_nam,'XTickLabelRotation',45);
xlabel('Visual signal(%)');
ylabel('ABI');
ylim([-0.2 0.2])

%% anova figure 5

x_anova = [];
for ll = 1:8
    x_anova = [x_anova; squeeze(x_diff(4:16,ll,2:5))];
end

[p,tb1,stats] = anova2(x_anova,13);
figure
c = multcompare(stats,'Estimate','column');

set(gca,'fontsize',14,'fontweight','bold','YTick',1:8,'YTickLabel',xtic_loc);
ylabel('Location');

%% anova category visual field

vs = 5;

if sum(ismember(vs,2))==1
pp = perfp(4:end,:,1:end,:);
else
    pp = perfp(1:end,:,1:end,:);
end

pp1=[];
pp2=[];
if sum(ismember(vs,1))==1; vs(vs==1)=[]; end

for v = vs
    pp1 = [pp1;pp(:,:,v,1)];
    pp2 = [pp2;pp(:,:,v,2)];
end

 pp1 = [nanmean(pp1(:,[2 3 4]),2),nanmean(pp1(:,[6 7 8]),2)];
 pp2 = [nanmean(pp2(:,[2 3 4]),2),nanmean(pp2(:,[6 7 8]),2)];

%ranksum(pp1(:,1),pp1(:,2))

pp3 = [pp1;pp2];

pp1 = reshape(pp1,[],1);
pp2 = reshape(pp2,[],1);
pp3 = [pp1, pp2];

[p,tb1,stats] = anova1(pp3);
[p,tb1,stats] = anova2(pp3,size(pp,1)*length(vs));


