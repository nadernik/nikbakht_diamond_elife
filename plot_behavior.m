% MATLAB script for analysis of the behavioral results shown in figures 1-2
% together with figure supplements in Nikbakht, N., & Diamond, M. E. (2021). 
% Conserved visual capacity of rats under red light. eLife, 10, e66429.
%
% This script analyzes rat behavioral data from visual discrimination tasks
% under different LED light conditions and generates publication-ready figures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and close all figures for clean start
clear all;
close all;
clc;

% Load behavioral data from MAT file
load behavior_data.mat;

% Define LED light conditions used in the experiment
LEDs = {'626 nm','652 nm','729 nm','854 nm','930 nm','white LED'};

% Number of rats in the experiment
nrats = 4;

% Define angle range for orientation discrimination task
angles = 0:10:90;  % Test angles from 0 to 90 degrees in 10-degree increments
x = 0:.1:90;      % Fine resolution for smooth curve plotting

% Set up color scheme for different LED conditions
cols = hot(11);                        % Start with hot colormap (11 colors)
cols = cols(1:5,:);                    % Take first 5 colors for LED conditions
cols = flipud(cols);                   % Flip colors upside down for better contrast
cols = [cols; lines(1);];              % Add blue color for white LED condition

% Transparency level for individual rat data
alph = .5;

% ============================================================================
% PLOT 1: Psychometric curves for each LED condition (Figure 2A-E and Figure 1E)
% ============================================================================
for f = 1:numel(LEDs)  % Loop through each LED condition
    for r = 1:nrats+1  % Loop through rats + 1 (last iteration is for average)
        
        % Determine starting condition based on LED type
        if f<6
            st = 4; % Consider only the control condition for colored LEDs
        else
            st = 3; % Different condition for white LED
        end
        
        % Create main figure for psychometric curves
        figure(1); 
        set(gcf,'Name','Figure-2(A-E) and Figure-1(E)','units','normalized','position',[.1 0.1 .7 .6]);
        subplot(2,round(numel(LEDs)/2),f);  % Create subplot for current LED
        hold on;
        
        % Fit psychometric function to data
        y = cumulativegaussianLapse(DATA.results.fit.percent_v(f,r,:),x);
        
        % Fit linear function to control data for certain LED conditions
        linefit = fit(angles', squeeze(DATA.results.percent_v(f,r,:)),'poly1');
        
        % Plot different curve types based on LED condition
        if f >= 4 && f<6
            % For 854 nm and 930 nm LEDs, plot linear fit
            if r == nrats+1
                psi(f)=plot(x,linefit(x),'color',cols(f,:),'linewidth',3); % Average data
            else
                psi(f)=plot(x,linefit(x),'color',[cols(f,:) alph],'linewidth',1.5); % Individual rat data
            end
        elseif f<4
            % For 626 nm, 652 nm, and 729 nm LEDs, plot psychometric curve
            if r == nrats+1
                psi(f)=plot(x,y,'color',cols(f,:),'linewidth',3); % Average data
            else
                psi(f)=plot(x,y,'color',[cols(f,:) alph],'linewidth',1.5); % Individual rat data
            end
        else
            % For white LED, plot psychometric curve
            if r == nrats+1
                psi(f)=plot(x,y,'color',cols(f,:),'linewidth',3); % Average data
            else
                psi(f)=plot(x,y,'color',[cols(f,:) alph],'linewidth',1.5); % Individual rat data
            end
        end
        
        % Add error bars for average data only
        if r == nrats+1
            errorbar(angles,squeeze(DATA.results.percent_v(f,r,:)),squeeze(2*DATA.results.stats.binomErrors(f,r,:)),'ko','markerfacecolor',cols(f,:));
        end
        
        % Add labels only to first subplot
        if r == 1
            ylabel('proportion called vertical'); 
            xlabel('angle (deg)');
        end
        
        % Configure axis properties
        set(gca,'xtick',[0,45,90]);
        xlim([-2 92]);
        ylim([-.01 1.01]);
        axis square;
        beautify(gca,1);  % Apply professional formatting
    end
end

% Add legend to first subplot
subplot(2,round(numel(LEDs)/2),1);
legend(psi,LEDs,'location','southeast');

% ============================================================================
% PLOT 2: Daily performance across sessions (Figure 2-S2)
% ============================================================================
figure(3);
set(gcf,'Name','Figure2-S2','units','normalized','position',[.1 0.1 .5 .4]);

ndays = 10;  % Number of days to display per condition

% Define specific day order to reproduce exact results from Figure 2-S2
days = [7    45    40    55    30    38    43    20     3    41    29     1 ...
    17    28    46     6    10    23     8    12    14    49     2    15    19 ...
    50    27    11    35    32     4    53    21    24    31    59    16    48 ...
    51    44    39    22    42    34    57    60    37     9    58    13    26 ...
    5    47    52    18    54    36    33    56    25];

% Calculate x-axis positions for each LED condition
fx = 1:ndays:numel(days);
ax = [];

% Process and plot daily performance for each LED condition
for f = 1:numel(LEDs)
    % Extract performance data for current LED condition
    perfinday = squeeze(DATA.perfinday(f,1:nrats,:,1));
    
    % Clean data for each rat
    for i = 1:nrats
        pxd = perfinday(i,:);
        pxd = pxd(~isnan(pxd));  % Remove NaN values
        pxd = pxd(pxd>0);        % Remove zero/negative values
        
        % Special filtering for 854 nm and 930 nm conditions
        if f>3 && f<6
            pxd = pxd(pxd>.3 & pxd<.6);
        end
        
        % Store cleaned data
        perfinday(i,1:numel(pxd)) = pxd;
    end
    
    % Special handling for 854 nm condition
    if f == 4
        perfinday = perfinday(:,95:end);
    end
    
    % Calculate median performance across rats
    meanperfinday(f,1:length(perfinday)) = nanmedian(perfinday);
    
    % Create boxplot for daily performance
    boxplot(perfinday(:,1:ndays),'symbol','.','colors',cols(f,:),'PlotStyle','compact','boxstyle','filled',...
        'factorseparator','auto','factorgap','auto','jitter',0,'positions',days(fx(f):fx(f)+ndays-1)); 
    hold on;
    
    % Add median line
    ax(f) = plot(days(fx(f):fx(f)+ndays-1),meanperfinday(f,1:ndays),'ko','markerfacecolor',cols(f,:),'markersize',8); 
    hold on;
end

% Configure plot appearance
xlim([0 numel(days)+1]);
ylim([0 1]);
xticks([1:4:numel(days)]);
xticklabels(1:4:numel(days));
legend('626 nm','652 nm','729 nm','854 nm','930 nm','white LED','location','southeast');
legend boxoff;
xlabel('test session');
ylabel('cumulative performance');
beautify(gca,1);

% ============================================================================
% PLOT 3: Average performance distributions (Figure 2-S1)
% ============================================================================
% Extract average performance data for each LED condition
perf_avg_626 = DATA.stats.perf_avg_626;
perf_avg_652 = DATA.stats.perf_avg_652;
perf_avg_729 = DATA.stats.perf_avg_729;
perf_avg_854 = DATA.stats.perf_avg_854;
perf_avg_930 = DATA.stats.perf_avg_930;
perf_avg_w   = DATA.stats.perf_avg_w;

% Create figure for performance distributions
figure(4);
set(gcf,'Name','Figure2-S1');

% Set histogram bin size
bin = 0.001;

% Plot smoothed histograms for each LED condition
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

% Add legend and configure plot
legend('626 nm','652 nm','729 nm','854 nm','930 nm','white LED','location','northwest');
legend boxoff;
xlim([.4 .8]); 
beautify(gca ,1);
xlabel('average performance'); 
ylabel('normalized count');

% ============================================================================
% PLOT 4: Cumulative performance comparison (Figure 2F)
% ============================================================================
% Extract cumulative performance data
cumperf = DATA.results.cumperf;

% Set up colors for individual rats
cols2 = lines(4);

% Create figure for cumulative performance
figure(2);
set(gcf,'Name','Figure2 (F)','units','normalized','position',[.1 0.1 .3 .4]);

% Plot individual rat performance
for r = 1:nrats
    p(r) = plot(1:5,cumperf(r,1:end-1),'-s','color',cols2(r,:),'linewidth',1.5);
    hold on;
    plot(6,cumperf(r,end),'s','color',cols2(r,:),'linewidth',1.5);
end

% Plot average performance with error bars
p(5) = errorbar(1:5,mean(cumperf(:,1:end-1)),std(cumperf(:,1:end-1))/2,'-ks','linewidth', 2.5);
errorbar(6,mean(cumperf(:,end)),std(cumperf(:,end))/2,'-ks','linewidth', 2.5);

% Configure axis labels and appearance
xticks([1:6]);
xticklabels([{'626 nm'};{'652 nm'};{'729 nm'};{'854 nm'};{'930 nm'};{'white LED'}]);
xtickangle(45)
xlim([0 7]);
ylabel('cumulative performance');
legend(p,'rat 1','rat 2','rat 3','rat 4','average','location','southwest');
legend boxoff;
beautify(gca,1);
