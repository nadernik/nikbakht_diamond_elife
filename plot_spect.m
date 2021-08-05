% MATLAB script for analysis of the spectral analysis shown in figure 1 (D)
% together with figure supplements in Nikbakht, N., & Diamond, M. E. (2021). 
% Conserved visual capacity of rats under red light. eLife, 10, e66429.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('spectral_data.mat');              % LED spectral measurements
load('plos_s_labda.mat');               % data retinal spectral sensitivity function S(?) using the electroretinogram
                                        % (ERG) from Rocha et al - PloS one, 2016 (DOI: https://doi.org/10.1371/journal.pone.0147318)
c = 299792458;                          % Speed of light [m/s].
h = 6.62607004*10^-34;                  % Planck's constant [m^2*kg/s].
cos_corr_diam = 0.39;                   % Cosine corrector diameter [cm].
pupil_diam = 0.2;                       % Fully dilated rat's eye pupil diameter [cm].
d = 50;                                 % Pupil to light source distance [cm].
cols = hot(11);                         % Setting colormap...
cols = cols(1:5,:);                     % Building 6 column colormap...
cols = flipud(cols);                    % Arranging colors upside down...
cols = [lines(1); cols;[.9 0 .5; .9 0 .2]];    % Adding blue for the white light...

for k = 1:numel(led)
    x = spectrum_raw{1,k};
    y = spectrum_raw{2,k};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cos_corr_area = pi*((cos_corr_diam/2)^2);     % Cosine corrector area [cm^2].
    pupil_area = pi*(pupil_diam/2)^2;             % Fully dilated rat's eye pupil area [cm^2].
    radiation_area = 4*pi*d^2/2;                  % Half the Area of radiation bc of LED shape at d [cm^2].
    pupil_frac = pupil_area/radiation_area;       % Photon fraction passing through the pupil [].
    cos_corr_frac = cos_corr_area/radiation_area; % Photon fraction passing through the cosine corrector [].
    PD = trapz(x, y);                             % Computing Power Density [uW/cm^2]...
    figure(1);
    set(gcf,'Name','Figure1 (D)','position',[100   100   950   550]); % Setting gcf parameters...
    SPD_plot(k) = plot(x, (10^4)*cos_corr_frac*y, '-', 'color', [cols(k, :) 1], 'LineWidth', 2); hold on; % SPD received by the cosine corrector at "d" [uW/m^2/nm].
    xlabel('wavelength (nm)');
    ylabel('Intensity (\muW m^{-2} nm^{-1})');
    fprintf('=================================\rPower density at LED: %s light = %2.2e (uW/cm^2)\r', ...
        led{k}, PD);
    fprintf('Intensity at LED: %s light = %2.2f (mW/cm^2)\r', led{k}, (10^-3)*PD);
    fprintf('Power at cosine corrector at LED: %s light = %2.2f (mW)\r', led{k}, (10^-3)*cos_corr_area*PD);
    axis tight; beautify(gca);
    ymax = ylim;
    axis tight;
    xlim([355 1030]);
    set(gca, 'YScale', 'log');
    beautify(gca, 1);
    figure(5);     % normalized spectra for comparison
    set(gcf,'Name','Figure1-S3','position',[100   100   950   550]);                                % Setting gcf parameters...
    SPD_plot_norm(k) = plot(x, y/max(y), '-', 'color', [cols(k, :) 1], 'LineWidth', 2); hold on;    % SPD received by the cosine corrector at "d" [uW/m^2/nm].
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PHOTON COUNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SPD = (10^-6)*(10^4)*y;                                                % Spectral Power Density [W/m^2/nm].
    SPC = (10^-12)/(h*c).*SPD.*x;                                          % Spectral Photon Count [#photons/um^2/s/nm].
    C = trapz(x, SPC);                                                     % Photon flux density [#photons/um^2/s].
    pupil_C = pupil_frac*C;                                                % Photon Count at pupil [#photons/um^2/s].
    retina_area = 0.8;                                                     % Retina Area of rat eye in cm^2
    pupil_retina_frac = pupil_area/retina_area;                            % ratio of pupil area to retinal surface area []
    retina_C = pupil_retina_frac * pupil_C;                                % Photon Count at retina [#photons/um^2/s].
    disk_outer_segment_area = 4*10^-8;                                     % Area of the disk outer segment cm^2 (https://link.springer.com/content/pdf/10.1023/A:1018563409196.pdf)
    disk_retina_frac = disk_outer_segment_area/retina_area;                % ratio of outer segment area to retinal surface area []
    disk_C = retina_C * disk_retina_frac;                                  % Photon flux density at disk outer segment [#photons/um^2/s].
    disk_photon_flux = disk_C * disk_outer_segment_area * 10^8;            % Photon flux for each photoreceptor [photon/s]
    a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,...
        -6675.007923501414, 1813.525992411163, -215.177888526334,...
        12.487558618387, -0.289541500599]; % from Stockman and Sharpe (2000) eq. (8)
    eq8=@(x)10.^(a(1)+a(2)*x.^2+a(3)*x.^4+a(4)*x.^6+a(5)*x.^8+a(6)*x.^10+a(7)*x.^12+a(8)*x.^14); % Sensitivity spectra of mouse opsin
    S_opsin_template = eq8(log10(x));
    S_activation = S_opsin_template.*(y/max(y)); % effective activation of photoreceptors by the LEDs (peak-normalised spectra)
    R_isomerization(k) = sum(S_activation * disk_photon_flux);             % photoisomerisation rate
    
    figure(2);
    set(gcf,'Name','Figure1-S2','position',[100   100   950   550]);
    SPC_plot(k) = plot(x, SPC*cos_corr_frac, '-', 'color', [cols(k, :) 1], 'LineWidth',2);% Spectral Photon Count [#photons/um^2/s/nm].
    
    hold on;
    xlabel('wavelength (nm)');
    ylabel('Spectral Photon Count [#photons \mum^{-2} s^{-1} nm^{-1}]');
    fprintf('Photon count received by cosine corrector at LED for %s light = %1.2e (#photons/um^2/s)\r', led{k}, C);
    fprintf('Photon count received by rat pupil for %s light = %1.2e (#photons/um^2/s)\r', led{k}, pupil_C);
    fprintf('Photon count received by rat pupil for %s light = %f (log #photons/cm^2/s)\r', led{k}, log10(pupil_C*10^8));
    fprintf('Power density received by rat pupil for %s light= %2.2e (uW/cm^2)\r', led{k}, PD*(pupil_area/radiation_area));
    fprintf('Photon flux density by each photoreceptor for %s light = %1.2e (#photons/um^2/s)\r', led{k}, disk_C);
    fprintf('Photon flux by each photoreceptor for %s light = %1.2e (#photons/s)\r', led{k}, disk_photon_flux);
    fprintf('Photoisomerization rate %s light = %f (P*/photoreceptor/s.10^3)\r', led{k}, R_isomerization(k)/1000);
    
    axis tight;
    xlim([355 1030]);
    set(gca, 'YScale', 'log');
    beautify(gca, 2);
end
led_names = {'White LED','626 nm','652 nm','729 nm','854 nm','930 nm','652 nm-UV-VIS','729 nm UV-VIS'};
figure(1);
legend(SPD_plot(1:k),led_names{1:k},'location','southeast');%legend boxoff;

figure(2);
legend(SPD_plot(1:k),led_names{1:k},'location','southeast');%legend boxoff;
%% generate rat opsin curves - S-, M- opsin as well as rhodopsin
figure(5);
col2 = winter(3);
wl = x;
a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,...
    -6675.007923501414, 1813.525992411163, -215.177888526334,...
    12.487558618387, -0.289541500599]; % from Stockman and Sharpe (2000) eq. (8)
xx = log10(x);
S_opsin_template=10.^(a(1)+a(2)*xx.^2+a(3)*xx.^4+a(4)*xx.^6+a(5)*xx.^8+a(6)*xx.^10+a(7)*xx.^12+a(8)*xx.^14); % Sensitivity spectra of mouse opsin
S_activation = S_opsin_template.*(y/max(y)); % effective activation of photoreceptors by the LEDs (peak-normalised spectra)

peak_nm = [360, 511, 500];
opsins  = {'S-opsin', 'M-opsin', 'rhodopsin'};
wl = 1:1300;
for p=1:numel(peak_nm)
    xx = log10(wl - peak_nm(p));
    spect_y = eq8(xx);
    plot((700:1300)-558,spect_y(700:1300),'linewidth',2,'color',col2(p,:)); hold on;axis tight;
end
plot(s_lambda.wl, s_lambda.fit_fft/max(s_lambda.fit_fft), 'k--', 'LineWidth',2);
ylim([0.0005 1]);
xlim([355 1030]);
xlabel('wavelength (nm)');
ylabel('normalized intensity/sensitivity');
legend([led_names(1:k)  opsins 'retinal sensitivity']); legend boxoff;
beautify(gca,1);
set(gca, 'YScale', 'log');
%% plot the retinal spectral sensitivity function S(lambda) using the electroretinogram (ERG) for Rats
figure(6);
plot(s_lambda.wl, s_lambda.measured, 'wo','markerfacecolor','k'); hold on;
plot(s_lambda.wl, s_lambda.fit_fft, '-', 'LineWidth',2);
plot(s_lambda.wl, s_lambda.uv_peak, 'm--', 'LineWidth',1);
plot(s_lambda.wl, s_lambda.green_peak, 'g--', 'LineWidth',1);
beautify(gca,1); ylim([0 100]);
ylabel('relative sensitivity');
xlabel('wavelength (nm)');
legend('observed', 'best fit','S-cone gaussian fit','M-cones gaussian fit'); legend boxoff;
