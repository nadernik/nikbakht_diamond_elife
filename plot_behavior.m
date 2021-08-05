% MATLAB script for analysis of the behavioral results shown in figures 1-2
% together with figure supplements in Nikbakht, N., & Diamond, M. E. (2021). 
% Conserved visual capacity of rats under red light. eLife, 10, e66429.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
load behavior_data.mat;
LEDs = {'626 nm','652 nm','729 nm','854 nm','930 nm','white LED'};
nrats = 4;
angles = 0:10:90;
x = 0:.1:90;
cols = hot(11);                        % Setting colormap...
cols = cols(1:5,:);                    % Building 6 column colormap...
cols = flipud(cols);                   % Arranging colors upside down...
cols = [cols; lines(1);];              % Adding blue for the white light...
alph = .5;
for f = 1:numel(LEDs)
    for r = 1:nrats+1
        if f<6
            st = 4; % consider only the control condition
        else
            st = 3;
        end
        figure(1); set(gcf,'Name','Figure-2(A-E) and Figure-1(E)','units','normalized','position',[.1 0.1 .7 .6]);
        subplot(2,round(numel(LEDs)/2),f);
        hold on;
        y = cumulativegaussianLapse(DATA.results.fit.percent_v(f,r,:),x);
        linefit = fit(angles', squeeze(DATA.results.percent_v(f,r,:)),'poly1'); % fit a line to control data

        if f >= 4 && f<6
            if r == nrats+1
                psi(f)=plot(x,linefit(x),'color',cols(f,:),'linewidth',3); % plot linear fit
            else
                psi(f)=plot(x,linefit(x),'color',[cols(f,:) alph],'linewidth',1.5);
            end
        elseif f<4
            if r == nrats+1
                psi(f)=plot(x,y,'color',cols(f,:),'linewidth',3);
            else
                psi(f)=plot(x,y,'color',[cols(f,:) alph],'linewidth',1.5);
                
            end
        else
            if r == nrats+1
                psi(f)=plot(x,y,'color',cols(f,:),'linewidth',3);
            else
                psi(f)=plot(x,y,'color',[cols(f,:) alph],'linewidth',1.5);
            end
        end
        if r == nrats+1
            errorbar(angles,squeeze(DATA.results.percent_v(f,r,:)),squeeze(2*DATA.results.stats.binomErrors(f,r,:)),'ko','markerfacecolor',cols(f,:));
        end
        if r == 1
            ylabel('proportion called vertical'); xlabel('angle (deg)');
        end
        set(gca,'xtick',[0,45,90]);
        xlim([-2 92]);
        ylim([-.01 1.01]);axis square;
        beautify(gca,1);
    end
end
subplot(2,round(numel(LEDs)/2),1);
legend(psi,LEDs,'location','southeast');
%% plot daily performance
figure(3);set(gcf,'Name','Figure2-S2','units','normalized','position',[.1 0.1 .5 .4]);
ndays = 10;
% days  = randperm(numel(LEDs)*ndays); %  random days for each condition
days = [7    45    40    55    30    38    43    20     3    41    29     1 ...
    17    28    46     6    10    23     8    12    14    49     2    15    19 ...
    50    27    11    35    32     4    53    21    24    31    59    16    48 ...
    51    44    39    22    42    34    57    60    37     9    58    13    26 ...
    5    47    52    18    54    36    33    56    25]; % to reproduce exactly the same results as shown in FIgure2-S2
fx = 1:ndays:numel(days);
ax = [];
for f = 1:numel(LEDs)
    perfinday = squeeze(DATA.perfinday(f,1:nrats,:,1));
    for i = 1:nrats
        pxd = perfinday(i,:);
        pxd = pxd(~isnan(pxd));
        pxd = pxd(pxd>0);
        if f>3 && f<6
            pxd = pxd(pxd>.3 & pxd<.6);
        end
        perfinday(i,1:numel(pxd)) = pxd;
    end
    if f == 4
        perfinday = perfinday(:,95:end);
    end
    meanperfinday(f,1:length(perfinday)) = nanmedian(perfinday);
    boxplot(perfinday(:,1:ndays),'symbol','.','colors',cols(f,:),'PlotStyle','compact','boxstyle','filled',...
        'factorseparator','auto','factorgap','auto','jitter',0,'positions',days(fx(f):fx(f)+ndays-1)); hold on;
    ax(f) = plot(days(fx(f):fx(f)+ndays-1),meanperfinday(f,1:ndays),'ko','markerfacecolor',cols(f,:),'markersize',8); hold on;
end

xlim([0 numel(days)+1]);
ylim([0 1]);
xticks([1:4:numel(days)]);
xticklabels(1:4:numel(days));
legend('626 nm','652 nm','729 nm','854 nm','930 nm','white LED','location','southeast');legend boxoff;
xlabel('test session');
ylabel('cumulative performance');
beautify(gca,1);
%% plot average performance distributions
perf_avg_626 = DATA.stats.perf_avg_626;
perf_avg_652 = DATA.stats.perf_avg_652;
perf_avg_729 = DATA.stats.perf_avg_729;
perf_avg_854 = DATA.stats.perf_avg_854;
perf_avg_930 = DATA.stats.perf_avg_930;
perf_avg_w   = DATA.stats.perf_avg_w;

figure(4);set(gcf,'Name','Figure2-S1');
bin = 0.001;
[n, xc]=hist(perf_avg_626,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2); hold all;
[n, xc]=hist(perf_avg_652,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2);
[n, xc]=hist(perf_avg_729,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2);
[n, xc]=hist(perf_avg_854,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2);
[n, xc]=hist(perf_avg_930,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2);
[n, xc]=hist(perf_avg_w,0:bin:1);n = n/sum(n);
plot(xc,smooth(n),'linewidth',2);
legend('626 nm','652 nm','729 nm','854 nm','930 nm','white LED','location','northwest');legend boxoff;
xlim([.4 .8]); beautify(gca ,1);
xlabel('average performance'); ylabel('normalized count');
%%
cumperf = DATA.results.cumperf;
cols2 = lines(4);
figure(2);set(gcf,'Name','Figure2 (F)','units','normalized','position',[.1 0.1 .3 .4]);
for r = 1:nrats
    p(r) = plot(1:5,cumperf(r,1:end-1),'-s','color',cols2(r,:),'linewidth',1.5);hold on;
    plot(6,cumperf(r,end),'s','color',cols2(r,:),'linewidth',1.5);
end
p(5) = errorbar(1:5,mean(cumperf(:,1:end-1)),std(cumperf(:,1:end-1))/2,'-ks','linewidth', 2.5);
errorbar(6,mean(cumperf(:,end)),std(cumperf(:,end))/2,'-ks','linewidth', 2.5);
xticks([1:6]);
xticklabels([{'626 nm'};{'652 nm'};{'729 nm'};{'854 nm'};{'930 nm'};{'white LED'}]);xtickangle(45)
xlim([0 7]);
ylabel('cumulative performance');
legend(p,'rat 1','rat 2','rat 3','rat 4','average','location','southwest');legend boxoff;
beautify(gca,1);
