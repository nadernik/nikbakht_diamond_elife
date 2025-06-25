"""
Analysis of spectral analysis shown in figure 1 (D) together with figure supplements 
in Nikbakht, N., & Diamond, M. E. (2021). Conserved visual capacity of rats under red light. 
eLife, 10, e66429.

This script analyzes LED spectral data and calculates photon flux, 
photoisomerization rates, and retinal sensitivity for rat vision studies.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.integrate import trapz
from scipy.ndimage import gaussian_filter1d

from beautify import beautify
from utils import load_matlab_data

# Physical constants for calculations
C = constants.c                    # Speed of light [m/s]
H = constants.h                    # Planck's constant [m^2*kg/s]

def eq8(x, a):
    """
    Stockman and Sharpe (2000) equation (8) for mouse opsin sensitivity.
    
    Parameters:
    -----------
    x : array-like
        Log10 of wavelength values
    a : array-like
        Coefficients from Stockman and Sharpe (2000)
    
    Returns:
    --------
    array-like: Opsin sensitivity values
    """
    return 10**(a[0] + a[1]*x**2 + a[2]*x**4 + a[3]*x**6 + 
                a[4]*x**8 + a[5]*x**10 + a[6]*x**12 + a[7]*x**14)

def main():
    """Main analysis function."""
    
    # Load spectral measurement data
    spectral_data = load_matlab_data('../spectral_data.mat')
    plos_data = load_matlab_data('../plos_s_labda.mat')
    
    if spectral_data is None or plos_data is None:
        print("Error: Could not load spectral data files")
        return
    
    # Extract data from MATLAB structures
    spectrum_raw = spectral_data['spectrum_raw'][0, 0]
    led = spectral_data['led'][0, 0]
    s_lambda = plos_data['s_lambda'][0, 0]
    
    # Experimental setup parameters
    cos_corr_diam = 0.39                   # Cosine corrector diameter [cm]
    pupil_diam = 0.2                       # Fully dilated rat's eye pupil diameter [cm]
    d = 50                                 # Pupil to light source distance [cm]
    
    # Set up color scheme for different LED conditions
    colors = plt.cm.hot(np.linspace(0, 1, 11))[:5]  # Take first 5 colors
    colors = np.flipud(colors)                       # Flip colors upside down
    colors = np.vstack([plt.cm.Set1(0), colors])     # Add blue for white light
    colors = np.vstack([colors, [[0.9, 0, 0.5], [0.9, 0, 0.2]]])  # Additional colors
    
    # ============================================================================
    # MAIN ANALYSIS LOOP: Process each LED condition
    # ============================================================================
    spd_plots = []
    spc_plots = []
    r_isomerization = []
    
    for k in range(len(led)):
        # Extract wavelength and intensity data for current LED
        x = spectrum_raw[0, k].flatten()  # Wavelength values [nm]
        y = spectrum_raw[1, k].flatten()  # Intensity values [arbitrary units]
        
        # ============================================================================
        # SPD: Spectral Power Density calculations
        # ============================================================================
        # Calculate areas for photon flux calculations
        cos_corr_area = np.pi * ((cos_corr_diam/2)**2)     # Cosine corrector area [cm^2]
        pupil_area = np.pi * (pupil_diam/2)**2             # Fully dilated rat's eye pupil area [cm^2]
        radiation_area = 4 * np.pi * d**2 / 2              # Half the area of radiation due to LED shape at distance d [cm^2]
        
        # Calculate photon fractions
        pupil_frac = pupil_area / radiation_area           # Photon fraction passing through the pupil []
        cos_corr_frac = cos_corr_area / radiation_area     # Photon fraction passing through the cosine corrector []
        
        # Calculate power density
        pd = trapz(y, x)                                   # Computing Power Density [μW/cm^2]
        
        # Create main spectral power density plot (Figure 1D)
        if k == 0:
            fig1, ax1 = plt.subplots(1, 1, figsize=(12, 8))
            ax1.set_title('Figure1 (D)')
        
        # Plot SPD received by cosine corrector at distance d [μW/m^2/nm]
        spd_plot, = ax1.plot(x, (10**4) * cos_corr_frac * y, 
                            color=colors[k], linewidth=2, label=led[k][0])
        spd_plots.append(spd_plot)
        
        # Print power density information
        print(f"=================================")
        print(f"Power density at LED: {led[k][0]} light = {pd:.2e} (μW/cm^2)")
        print(f"Intensity at LED: {led[k][0]} light = {(10**-3)*pd:.2f} (mW/cm^2)")
        print(f"Power at cosine corrector at LED: {led[k][0]} light = {(10**-3)*cos_corr_area*pd:.2f} (mW)")
        
        # Create normalized spectra plot for comparison (Figure 1-S3)
        if k == 0:
            fig5, ax5 = plt.subplots(1, 1, figsize=(12, 8))
            ax5.set_title('Figure1-S3')
        
        # Plot normalized SPD for comparison
        ax5.plot(x, y/np.max(y), color=colors[k], linewidth=2, label=led[k][0])
        
        # ============================================================================
        # PHOTON COUNT: Photon flux and photoisomerization calculations
        # ============================================================================
        # Convert to spectral power density in SI units
        spd = (10**-6) * (10**4) * y                        # Spectral Power Density [W/m^2/nm]
        
        # Calculate spectral photon count using E = hc/λ
        spc = (10**-12) / (H * C) * spd * x                 # Spectral Photon Count [#photons/μm^2/s/nm]
        
        # Calculate total photon flux density
        c_total = trapz(spc, x)                             # Photon flux density [#photons/μm^2/s]
        
        # Calculate photon count at rat pupil
        pupil_c = pupil_frac * c_total                      # Photon Count at pupil [#photons/μm^2/s]
        
        # Calculate photon count at retina
        retina_area = 0.8                                   # Retina area of rat eye in cm^2
        pupil_retina_frac = pupil_area / retina_area        # Ratio of pupil area to retinal surface area []
        retina_c = pupil_retina_frac * pupil_c              # Photon Count at retina [#photons/μm^2/s]
        
        # Calculate photon flux at individual photoreceptor disks
        disk_outer_segment_area = 4e-8                      # Area of the disk outer segment cm^2 (from literature)
        disk_retina_frac = disk_outer_segment_area / retina_area  # Ratio of outer segment area to retinal surface area []
        disk_c = retina_c * disk_retina_frac                # Photon flux density at disk outer segment [#photons/μm^2/s]
        disk_photon_flux = disk_c * disk_outer_segment_area * 10**8  # Photon flux for each photoreceptor [photon/s]
        
        # Calculate photoisomerization rate using opsin sensitivity
        # Coefficients from Stockman and Sharpe (2000) equation (8) for mouse opsin
        a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,
             -6675.007923501414, 1813.525992411163, -215.177888526334,
             12.487558618387, -0.289541500599]
        
        # Calculate opsin template sensitivity
        s_opsin_template = eq8(np.log10(x), a)
        
        # Calculate effective activation of photoreceptors by the LEDs (peak-normalized spectra)
        s_activation = s_opsin_template * (y / np.max(y))
        
        # Calculate photoisomerization rate
        r_iso = np.sum(s_activation * disk_photon_flux)     # Photoisomerization rate [P*/photoreceptor/s]
        r_isomerization.append(r_iso)
        
        # Create spectral photon count plot (Figure 1-S2)
        if k == 0:
            fig2, ax2 = plt.subplots(1, 1, figsize=(12, 8))
            ax2.set_title('Figure1-S2')
        
        # Plot spectral photon count
        spc_plot, = ax2.plot(x, spc * cos_corr_frac, 
                            color=colors[k], linewidth=2, label=led[k][0])
        spc_plots.append(spc_plot)
        
        # Print detailed photon count information
        print(f"Photon count received by cosine corrector at LED for {led[k][0]} light = {c_total:.2e} (#photons/μm^2/s)")
        print(f"Photon count received by rat pupil for {led[k][0]} light = {pupil_c:.2e} (#photons/μm^2/s)")
        print(f"Photon count received by rat pupil for {led[k][0]} light = {np.log10(pupil_c*10**8):f} (log #photons/cm^2/s)")
        print(f"Power density received by rat pupil for {led[k][0]} light= {pd*(pupil_area/radiation_area):.2e} (μW/cm^2)")
        print(f"Photon flux density by each photoreceptor for {led[k][0]} light = {disk_c:.2e} (#photons/μm^2/s)")
        print(f"Photon flux by each photoreceptor for {led[k][0]} light = {disk_photon_flux:.2e} (#photons/s)")
        print(f"Photoisomerization rate {led[k][0]} light = {r_iso/1000:f} (P*/photoreceptor/s.10^3)")
    
    # ============================================================================
    # ADD LEGENDS TO PLOTS
    # ============================================================================
    # Define LED names for legend
    led_names = ['White LED', '626 nm', '652 nm', '729 nm', '854 nm', '930 nm', '652 nm-UV-VIS', '729 nm UV-VIS']
    
    # Configure SPD plot
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('Intensity (μW m^{-2} nm^{-1})')
    ax1.set_xlim(355, 1030)
    ax1.set_yscale('log')
    ax1.legend(led_names[:len(spd_plots)], loc='southeast')
    beautify(ax1, print_ready=True)
    
    # Configure SPC plot
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('Spectral Photon Count [#photons μm^{-2} s^{-1} nm^{-1}]')
    ax2.set_xlim(355, 1030)
    ax2.set_yscale('log')
    ax2.legend(led_names[:len(spc_plots)], loc='southeast')
    beautify(ax2, print_ready=True)
    
    # ============================================================================
    # GENERATE RAT OPSIN SENSITIVITY CURVES (Figure 1-S3)
    # ============================================================================
    # Coefficients from Stockman and Sharpe (2000) equation (8) for mouse opsin
    a = [-188862.970810906644, 90228.966712600282, -2483.531554344362,
         -6675.007923501414, 1813.525992411163, -215.177888526334,
         12.487558618387, -0.289541500599]
    
    # Calculate opsin template sensitivity
    xx = np.log10(x)
    s_opsin_template = eq8(xx, a)
    
    # Calculate effective activation of photoreceptors by the LEDs (peak-normalized spectra)
    s_activation = s_opsin_template * (y / np.max(y))
    
    # Define peak wavelengths for different opsin types
    peak_nm = [360, 511, 500]  # S-opsin, M-opsin, rhodopsin peak wavelengths
    opsins = ['S-opsin', 'M-opsin', 'rhodopsin']
    wl = np.arange(1, 1301)
    
    # Set up colors for different opsin types
    opsin_colors = plt.cm.winter(np.linspace(0, 1, 3))
    
    # Plot opsin sensitivity curves
    for p, (peak, opsin_name, color) in enumerate(zip(peak_nm, opsins, opsin_colors)):
        xx_opsin = np.log10(wl - peak)
        spect_y = eq8(xx_opsin, a)
        ax5.plot((np.arange(700, 1301) - 558), spect_y[699:1300], 
                linewidth=2, color=color, label=opsin_name)
    
    # Add retinal sensitivity curve from ERG data
    s_lambda_wl = s_lambda['wl'][0, 0].flatten()
    s_lambda_fit = s_lambda['fit_fft'][0, 0].flatten()
    ax5.plot(s_lambda_wl, s_lambda_fit / np.max(s_lambda_fit), 
            'k--', linewidth=2, label='retinal sensitivity')
    
    # Configure plot appearance
    ax5.set_ylim(0.0005, 1)
    ax5.set_xlim(355, 1030)
    ax5.set_xlabel('wavelength (nm)')
    ax5.set_ylabel('normalized intensity/sensitivity')
    ax5.set_yscale('log')
    ax5.legend(loc='best')
    beautify(ax5, print_ready=True)
    
    # ============================================================================
    # PLOT RETINAL SPECTRAL SENSITIVITY FUNCTION (Figure 6)
    # ============================================================================
    fig6, ax6 = plt.subplots(1, 1, figsize=(12, 8))
    ax6.set_title('Retinal Spectral Sensitivity Function')
    
    # Extract ERG data
    s_lambda_measured = s_lambda['measured'][0, 0].flatten()
    s_lambda_uv_peak = s_lambda['uv_peak'][0, 0].flatten()
    s_lambda_green_peak = s_lambda['green_peak'][0, 0].flatten()
    
    # Plot measured ERG data and fitted curves
    ax6.plot(s_lambda_wl, s_lambda_measured, 'wo', markerfacecolor='k', 
            markersize=8, label='observed')
    ax6.plot(s_lambda_wl, s_lambda_fit, '-', linewidth=2, label='best fit')
    ax6.plot(s_lambda_wl, s_lambda_uv_peak, 'm--', linewidth=1, label='S-cone gaussian fit')
    ax6.plot(s_lambda_wl, s_lambda_green_peak, 'g--', linewidth=1, label='M-cones gaussian fit')
    
    # Configure plot appearance
    ax6.set_ylim(0, 100)
    ax6.set_xlabel('wavelength (nm)')
    ax6.set_ylabel('relative sensitivity')
    ax6.legend(loc='best')
    beautify(ax6, print_ready=True)
    
    plt.show()

if __name__ == "__main__":
    main() 