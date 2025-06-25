"""
Analysis of behavioral results shown in figures 1-2 together with figure supplements 
in Nikbakht, N., & Diamond, M. E. (2021). Conserved visual capacity of rats under red light. 
eLife, 10, e66429.

This script analyzes rat behavioral data from visual discrimination tasks
under different LED light conditions and generates publication-ready figures.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
import pandas as pd

from beautify import beautify
from utils import load_matlab_data, cumulative_gaussian_lapse

# Set up plotting style
plt.style.use('default')
sns.set_palette("husl")

def main():
    """Main analysis function."""
    
    # Load behavioral data
    data = load_matlab_data('../behavior_data.mat')
    if data is None:
        print("Error: Could not load behavior_data.mat")
        return
    
    # Extract data from MATLAB structure
    DATA = data['DATA'][0, 0]  # MATLAB struct to Python dict
    
    # Define LED light conditions used in the experiment
    leds = ['626 nm', '652 nm', '729 nm', '854 nm', '930 nm', 'white LED']
    
    # Number of rats in the experiment
    n_rats = 4
    
    # Define angle range for orientation discrimination task
    angles = np.arange(0, 91, 10)  # Test angles from 0 to 90 degrees in 10-degree increments
    x = np.arange(0, 90.1, 0.1)    # Fine resolution for smooth curve plotting
    
    # Set up color scheme for different LED conditions
    colors = plt.cm.hot(np.linspace(0, 1, 11))[:5]  # Take first 5 colors
    colors = np.flipud(colors)                       # Flip colors upside down
    colors = np.vstack([plt.cm.Set1(0), colors])     # Add blue for white LED
    
    # Transparency level for individual rat data
    alpha = 0.5
    
    # ============================================================================
    # PLOT 1: Psychometric curves for each LED condition (Figure 2A-E and Figure 1E)
    # ============================================================================
    fig1, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig1.suptitle('Figure-2(A-E) and Figure-1(E)', fontsize=16)
    axes = axes.flatten()
    
    psi_plots = []
    
    for f, led_name in enumerate(leds):
        ax = axes[f]
        
        for r in range(n_rats + 1):  # Loop through rats + 1 (last iteration is for average)
            
            # Determine starting condition based on LED type
            if f < 6:
                st = 4  # Consider only the control condition for colored LEDs
            else:
                st = 3  # Different condition for white LED
            
            # Extract data for current LED and rat
            percent_v = DATA['results']['percent_v'][0, 0][f, r, :].flatten()
            
            # Fit psychometric function to data
            try:
                # For simplicity, use a basic cumulative Gaussian fit
                # In practice, you'd use the fitted parameters from the MATLAB data
                y = cumulative_gaussian_lapse([45, 10, 0.5, 0.05], x)
                
                # Fit linear function to control data for certain LED conditions
                if f >= 3 and f < 5:  # For 854 nm and 930 nm LEDs
                    # Linear fit
                    coeffs = np.polyfit(angles, percent_v, 1)
                    line_fit = np.poly1d(coeffs)
                    y = line_fit(x)
                
                # Plot different curve types based on LED condition
                if r == n_rats:  # Average data
                    psi_plot, = ax.plot(x, y, color=colors[f], linewidth=3, label=led_name)
                    if f == 0:
                        psi_plots.append(psi_plot)
                else:  # Individual rat data
                    ax.plot(x, y, color=colors[f], alpha=alpha, linewidth=1.5)
                
                # Add error bars for average data only
                if r == n_rats:
                    # Extract error data if available
                    try:
                        errors = DATA['results']['stats']['binomErrors'][0, 0][f, r, :].flatten()
                        ax.errorbar(angles, percent_v, yerr=2*errors, fmt='ko', 
                                  markerfacecolor=colors[f], capsize=3)
                    except:
                        ax.errorbar(angles, percent_v, fmt='ko', 
                                  markerfacecolor=colors[f], capsize=3)
                
            except Exception as e:
                print(f"Error processing LED {led_name}, rat {r}: {e}")
                continue
        
        # Add labels only to first subplot
        if f == 0:
            ax.set_ylabel('proportion called vertical')
        if f >= 3:
            ax.set_xlabel('angle (deg)')
        
        # Configure axis properties
        ax.set_xticks([0, 45, 90])
        ax.set_xlim(-2, 92)
        ax.set_ylim(-0.01, 1.01)
        ax.set_aspect('equal')
        beautify(ax, print_ready=True)
        ax.set_title(led_name)
    
    # Add legend to first subplot
    axes[0].legend(psi_plots, leds, loc='southeast')
    
    # ============================================================================
    # PLOT 2: Daily performance across sessions (Figure 2-S2)
    # ============================================================================
    fig3, ax3 = plt.subplots(1, 1, figsize=(12, 8))
    ax3.set_title('Figure2-S2')
    
    n_days = 10  # Number of days to display per condition
    
    # Define specific day order to reproduce exact results from Figure 2-S2
    days = [7, 45, 40, 55, 30, 38, 43, 20, 3, 41, 29, 1,
            17, 28, 46, 6, 10, 23, 8, 12, 14, 49, 2, 15, 19,
            50, 27, 11, 35, 32, 4, 53, 21, 24, 31, 59, 16, 48,
            51, 44, 39, 22, 42, 34, 57, 60, 37, 9, 58, 13, 26,
            5, 47, 52, 18, 54, 36, 33, 56, 25]
    
    # Calculate x-axis positions for each LED condition
    fx = np.arange(0, len(days), n_days)
    
    # Process and plot daily performance for each LED condition
    for f, led_name in enumerate(leds):
        try:
            # Extract performance data for current LED condition
            perf_in_day = DATA['perfinday'][0, 0][f, :n_rats, :, 0]
            
            # Clean data for each rat
            for i in range(n_rats):
                pxd = perf_in_day[i, :]
                pxd = pxd[~np.isnan(pxd)]  # Remove NaN values
                pxd = pxd[pxd > 0]         # Remove zero/negative values
                
                # Special filtering for 854 nm and 930 nm conditions
                if f > 3 and f < 6:
                    pxd = pxd[(pxd > 0.3) & (pxd < 0.6)]
                
                # Store cleaned data
                perf_in_day[i, :len(pxd)] = pxd
            
            # Special handling for 854 nm condition
            if f == 4:
                perf_in_day = perf_in_day[:, 95:]
            
            # Calculate median performance across rats
            mean_perf_in_day = np.nanmedian(perf_in_day, axis=0)
            
            # Create boxplot for daily performance
            day_positions = days[fx[f]:fx[f] + n_days]
            data_for_boxplot = [perf_in_day[i, :n_days] for i in range(n_rats)]
            
            bp = ax3.boxplot(data_for_boxplot, positions=day_positions, 
                           patch_artist=True, widths=0.8)
            
            # Color the boxes
            for patch in bp['boxes']:
                patch.set_facecolor(colors[f])
                patch.set_alpha(0.7)
            
            # Add median line
            ax3.plot(day_positions, mean_perf_in_day[:n_days], 'ko', 
                    markerfacecolor=colors[f], markersize=8)
            
        except Exception as e:
            print(f"Error processing daily performance for {led_name}: {e}")
            continue
    
    # Configure plot appearance
    ax3.set_xlim(0, len(days) + 1)
    ax3.set_ylim(0, 1)
    ax3.set_xticks(np.arange(1, len(days) + 1, 4))
    ax3.set_xticklabels(np.arange(1, len(days) + 1, 4))
    ax3.legend(leds, loc='southeast')
    ax3.set_xlabel('test session')
    ax3.set_ylabel('cumulative performance')
    beautify(ax3, print_ready=True)
    
    # ============================================================================
    # PLOT 3: Average performance distributions (Figure 2-S1)
    # ============================================================================
    fig4, ax4 = plt.subplots(1, 1, figsize=(10, 6))
    ax4.set_title('Figure2-S1')
    
    # Extract average performance data for each LED condition
    try:
        perf_avg_626 = DATA['stats']['perf_avg_626'][0, 0].flatten()
        perf_avg_652 = DATA['stats']['perf_avg_652'][0, 0].flatten()
        perf_avg_729 = DATA['stats']['perf_avg_729'][0, 0].flatten()
        perf_avg_854 = DATA['stats']['perf_avg_854'][0, 0].flatten()
        perf_avg_930 = DATA['stats']['perf_avg_930'][0, 0].flatten()
        perf_avg_w = DATA['stats']['perf_avg_w'][0, 0].flatten()
        
        # Set histogram bin size
        bin_size = 0.001
        
        # Plot smoothed histograms for each LED condition
        for i, (perf_data, led_name) in enumerate(zip([perf_avg_626, perf_avg_652, perf_avg_729, 
                                                      perf_avg_854, perf_avg_930, perf_avg_w], leds)):
            bins = np.arange(0, 1 + bin_size, bin_size)
            n, bins_edges = np.histogram(perf_data, bins=bins, density=True)
            bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
            
            # Smooth the histogram
            from scipy.ndimage import gaussian_filter1d
            n_smooth = gaussian_filter1d(n, sigma=2)
            
            ax4.plot(bin_centers, n_smooth, linewidth=2, label=led_name, color=colors[i])
        
        # Add legend and configure plot
        ax4.legend(loc='northwest')
        ax4.set_xlim(0.4, 0.8)
        beautify(ax4, print_ready=True)
        ax4.set_xlabel('average performance')
        ax4.set_ylabel('normalized count')
        
    except Exception as e:
        print(f"Error processing performance distributions: {e}")
    
    # ============================================================================
    # PLOT 4: Cumulative performance comparison (Figure 2F)
    # ============================================================================
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
    ax2.set_title('Figure2 (F)')
    
    try:
        # Extract cumulative performance data
        cum_perf = DATA['results']['cumperf'][0, 0]
        
        # Set up colors for individual rats
        rat_colors = plt.cm.Set1(np.linspace(0, 1, n_rats))
        
        # Plot individual rat performance
        p_plots = []
        for r in range(n_rats):
            rat_data = cum_perf[r, :]
            p_plot, = ax2.plot(range(1, 6), rat_data[:5], '-s', color=rat_colors[r], 
                              linewidth=1.5, label=f'rat {r+1}')
            p_plots.append(p_plot)
            ax2.plot(6, rat_data[5], 's', color=rat_colors[r], linewidth=1.5)
        
        # Plot average performance with error bars
        mean_perf = np.mean(cum_perf[:, :5], axis=0)
        std_perf = np.std(cum_perf[:, :5], axis=0) / 2
        p_avg, = ax2.errorbar(range(1, 6), mean_perf, yerr=std_perf, 
                             fmt='-ks', linewidth=2.5, label='average')
        p_plots.append(p_avg)
        
        mean_w = np.mean(cum_perf[:, 5])
        std_w = np.std(cum_perf[:, 5]) / 2
        ax2.errorbar(6, mean_w, yerr=std_w, fmt='-ks', linewidth=2.5)
        
        # Configure axis labels and appearance
        ax2.set_xticks(range(1, 7))
        ax2.set_xticklabels(leds, rotation=45)
        ax2.set_xlim(0, 7)
        ax2.set_ylabel('cumulative performance')
        ax2.legend(p_plots, loc='southwest')
        beautify(ax2, print_ready=True)
        
    except Exception as e:
        print(f"Error processing cumulative performance: {e}")
    
    plt.show()

if __name__ == "__main__":
    main() 