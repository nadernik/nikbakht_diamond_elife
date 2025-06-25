% MATLAB script for analysis of the spectral analysis shown in figure 1 (D)
% together with figure supplements in Nikbakht, N., & Diamond, M. E. (2021). 
% Conserved visual capacity of rats under red light. eLife, 10, e66429.
%
% This script analyzes LED spectral data and calculates photon flux, 
% photoisomerization rates, and retinal sensitivity for rat vision studies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and close all figures for clean start
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load spectral measurement data
load('spectral_data.mat');              % LED spectral measurements
load('plos_s_labda.mat');               % Data retinal spectral sensitivity function S(λ) using the electroretinogram
                                        % (ERG) from Rocha et al - PloS one, 2016 (DOI: https://doi.org/10.1371/journal.pone.0147318)

% Physical constants for calculations
c = 299792458;                          % Speed of light [m/s]
h = 6.62607004*10^-34;                  % Planck's constant [m^2*kg/s]

% Experimental setup parameters
cos_corr_diam = 0.39;                   % Cosine corrector diameter [cm]
pupil_diam = 0.2;                       % Fully dilated rat's eye pupil diameter [cm]
d = 50;                                 % Pupil to light source distance [cm]

% Set up color scheme for different LED conditions
cols = hot(11);                         % Start with hot colormap (11 colors)
cols = cols(1:5,:);                     % Take first 5 colors for LED conditions
cols = flipud(cols);                    % Flip colors upside down for better contrast
cols = [lines(1); cols;[.9 0 .5; .9 0 .2]];    % Add blue for white light and additional colors

% ============================================================================
% MAIN ANALYSIS LOOP: Process each LED condition
% ============================================================================
for k = 1:numel(led)
    % Extract wavelength and intensity data for current LED
    x = spectrum_raw{1,k};  % Wavelength values [nm]
    y = spectrum_raw{2,k};  % Intensity values [arbitrary units]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate areas for photon flux calculations
    cos_corr_area = pi*((cos_corr_diam/2)^2);     % Cosine corrector area [cm^2]
    pupil_area = pi*(pupil_diam/2)^2;             % Fully dilated rat's eye pupil area [cm^2]
    radiation_area = 4*pi*d^2/2;                  % Half the area of radiation due to LED shape at distance d [cm^2]
    
    % Calculate photon fractions
    pupil_frac = pupil_area/radiation_area;       % Photon fraction passing through the pupil []
    cos_corr_frac = cos_corr_area/radiation_area; % Photon fraction passing through the cosine corrector []
    
    % Calculate power density
    PD = trapz(x, y);                             % Computing Power Density [μW/cm^2]
    
    % Create main spectral power density plot (Figure 1D)
    figure(1);
    set(gcf,'Name','Figure1 (D)','position',[100   100   950   550]); % Set figure parameters
    
    % Plot SPD received by cosine corrector at distance d [μW/m^2/nm]
    SPD_plot(k) = plot(x, (10^4)*cos_corr_frac*y, '-', 'color', [cols(k, :) 1], 'LineWidth', 2); 
    hold on;
    
    % Add labels
    xlabel('wavelength (nm)');
    ylabel('Intensity (μW m^{-2} nm^{-1})');
    
    % Print power density information
    fprintf('=================================\rPower density at LED: %s light = %2.2e (μW/cm^2)\r', ...
        led{k}, PD);
    fprintf('Intensity at LED: %s light = %2.2f (mW/cm^2)\r', led{k}, (10^-3)*PD);
    fprintf('Power at cosine corrector at LED: %s light = %2.2f (mW)\r', led{k}, (10^-3)*cos_corr_area*PD);
    
    % Configure plot appearance
    axis tight; 
    beautify(gca);
    ymax = ylim;
    axis tight;
    xlim([355 1030]);
    set(gca, 'YScale', 'log');
    beautify(gca, 1);
    
    % Create normalized spectra plot for comparison (Figure 1-S3)
    figure(5);
    set(gcf,'Name','Figure1-S3','position',[100   100   950   550]);
    
    % Plot normalized SPD for comparison
    SPD_plot_norm(k) = plot(x, y/max(y), '-', 'color', [cols(k, :) 1], 'LineWidth', 2); 
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PHOTON COUNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert to spectral power density in SI units
    SPD = (10^-6)*(10^4)*y;                                                % Spectral Power Density [W/m^2/nm]
    
    % Calculate spectral photon count using E = hc/λ
    SPC = (10^-12)/(h*c).*SPD.*x;                                          % Spectral Photon Count [#photons/μm^2/s/nm]
    
    % Calculate total photon flux density
    C = trapz(x, SPC);                                                     % Photon flux density [#photons/μm^2/s]
    
    % Calculate photon count at rat pupil
    pupil_C = pupil_frac*C;                                                % Photon Count at pupil [#photons/μm^2/s]
    
    % Calculate photon count at retina
    retina_area = 0.8;                                                     % Retina area of rat eye in cm^2
    pupil_retina_frac = pupil_area/retina_area;                            % Ratio of pupil area to retinal surface area []
    retina_C = pupil_retina_frac * pupil_C;                                % Photon Count at retina [#photons/μm^2/s]
    
    % Calculate photon flux at individual photoreceptor disks
    disk_outer_segment_area = 4*10^-8;                                     % Area of the disk outer segment cm^2 (from literature)
    disk_retina_frac = disk_outer_segment_area/retina_area;                % Ratio of outer segment area to retinal surface area []
    disk_C = retina_C * disk_retina_frac;                                  % Photon flux density at disk outer segment [#photons/μm^2/s]
    disk_photon_flux = disk_C * disk_outer_segment_area * 10^8;            % Photon flux for each photoreceptor [photon/s]
    
    % Calculate photoisomerization rate using opsin sensitivity
    % Coefficients from Stockman and Sharpe (2000) equation (8) for mouse opsin
    a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,...
        -6675.007923501414, 1813.525992411163, -215.177888526334,...
        12.487558618387, -0.289541500599];
    
    % Define opsin sensitivity function
    eq8=@(x)10.^(a(1)+a(2)*x.^2+a(3)*x.^4+a(4)*x.^6+a(5)*x.^8+a(6)*x.^10+a(7)*x.^12+a(8)*x.^14);
    
    % Calculate opsin template sensitivity
    S_opsin_template = eq8(log10(x));
    
    % Calculate effective activation of photoreceptors by the LEDs (peak-normalized spectra)
    S_activation = S_opsin_template.*(y/max(y));
    
    % Calculate photoisomerization rate
    R_isomerization(k) = sum(S_activation * disk_photon_flux);             % Photoisomerization rate [P*/photoreceptor/s]
    
    % Create spectral photon count plot (Figure 1-S2)
    figure(2);
    set(gcf,'Name','Figure1-S2','position',[100   100   950   550]);
    
    % Plot spectral photon count
    SPC_plot(k) = plot(x, SPC*cos_corr_frac, '-', 'color', [cols(k, :) 1], 'LineWidth',2);
    hold on;
    
    % Add labels
    xlabel('wavelength (nm)');
    ylabel('Spectral Photon Count [#photons μm^{-2} s^{-1} nm^{-1}]');
    
    % Print detailed photon count information
    fprintf('Photon count received by cosine corrector at LED for %s light = %1.2e (#photons/μm^2/s)\r', led{k}, C);
    fprintf('Photon count received by rat pupil for %s light = %1.2e (#photons/μm^2/s)\r', led{k}, pupil_C);
    fprintf('Photon count received by rat pupil for %s light = %f (log #photons/cm^2/s)\r', led{k}, log10(pupil_C*10^8));
    fprintf('Power density received by rat pupil for %s light= %2.2e (μW/cm^2)\r', led{k}, PD*(pupil_area/radiation_area));
    fprintf('Photon flux density by each photoreceptor for %s light = %1.2e (#photons/μm^2/s)\r', led{k}, disk_C);
    fprintf('Photon flux by each photoreceptor for %s light = %1.2e (#photons/s)\r', led{k}, disk_photon_flux);
    fprintf('Photoisomerization rate %s light = %f (P*/photoreceptor/s.10^3)\r', led{k}, R_isomerization(k)/1000);
    
    % Configure plot appearance
    axis tight;
    xlim([355 1030]);
    set(gca, 'YScale', 'log');
    beautify(gca, 2);
end

% ============================================================================
% ADD LEGENDS TO PLOTS
% ============================================================================
% Define LED names for legend
led_names = {'White LED','626 nm','652 nm','729 nm','854 nm','930 nm','652 nm-UV-VIS','729 nm UV-VIS'};

% Add legend to SPD plot
figure(1);
legend(SPD_plot(1:k),led_names{1:k},'location','southeast');

% Add legend to SPC plot
figure(2);
legend(SPD_plot(1:k),led_names{1:k},'location','southeast');

% ============================================================================
% GENERATE RAT OPSIN SENSITIVITY CURVES (Figure 1-S3)
% ============================================================================
figure(5);
col2 = winter(3);  % Color scheme for different opsin types
wl = x;

% Coefficients from Stockman and Sharpe (2000) equation (8) for mouse opsin
a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,...
    -6675.007923501414, 1813.525992411163, -215.177888526334,...
    12.487558618387, -0.289541500599];

% Calculate opsin template sensitivity
xx = log10(x);
S_opsin_template=10.^(a(1)+a(2)*xx.^2+a(3)*xx.^4+a(4)*xx.^6+a(5)*xx.^8+a(6)*xx.^10+a(7)*xx.^12+a(8)*xx.^14);

% Calculate effective activation of photoreceptors by the LEDs (peak-normalized spectra)
S_activation = S_opsin_template.*(y/max(y));

% Define peak wavelengths for different opsin types
peak_nm = [360, 511, 500];  % S-opsin, M-opsin, rhodopsin peak wavelengths
opsins  = {'S-opsin', 'M-opsin', 'rhodopsin'};
wl = 1:1300;

% Plot opsin sensitivity curves
for p=1:numel(peak_nm)
    xx = log10(wl - peak_nm(p));
    spect_y = eq8(xx);
    plot((700:1300)-558,spect_y(700:1300),'linewidth',2,'color',col2(p,:)); 
    hold on;
    axis tight;
end

% Add retinal sensitivity curve from ERG data
plot(s_lambda.wl, s_lambda.fit_fft/max(s_lambda.fit_fft), 'k--', 'LineWidth',2);

% Configure plot appearance
ylim([0.0005 1]);
xlim([355 1030]);
xlabel('wavelength (nm)');
ylabel('normalized intensity/sensitivity');
legend([led_names(1:k)  opsins 'retinal sensitivity']); 
legend boxoff;
beautify(gca,1);
set(gca, 'YScale', 'log');

% ============================================================================
% PLOT RETINAL SPECTRAL SENSITIVITY FUNCTION (Figure 6)
% ============================================================================
figure(6);

% Plot measured ERG data and fitted curves
plot(s_lambda.wl, s_lambda.measured, 'wo','markerfacecolor','k'); 
hold on;
plot(s_lambda.wl, s_lambda.fit_fft, '-', 'LineWidth',2);
plot(s_lambda.wl, s_lambda.uv_peak, 'm--', 'LineWidth',1);
plot(s_lambda.wl, s_lambda.green_peak, 'g--', 'LineWidth',1);

% Configure plot appearance
beautify(gca,1); 
ylim([0 100]);
ylabel('relative sensitivity');
xlabel('wavelength (nm)');
legend('observed', 'best fit','S-cone gaussian fit','M-cones gaussian fit'); 
legend boxoff;
