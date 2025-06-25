"""
Function beautify renders plots and figures using matplotlib for high-quality graphics
suitable for publication or presentation. It also modifies axis properties for
consistent, professional appearance.

This is a Python translation of the original MATLAB beautify.m function.

INPUTS:
    ax - matplotlib axes object (optional, defaults to current axes)
    print_ready - boolean flag for print-ready settings (optional, defaults to False)

OUTPUTS:
    None - modifies the current figure and axes properties

AUTHOR: NADER NIKBAKHT, SISSA 2013
TRANSLATED TO PYTHON: 2024
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

def beautify(ax=None, print_ready=False):
    """
    Apply professional formatting to matplotlib plots.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes, optional
        Axes object to format. If None, uses current axes.
    print_ready : bool, optional
        If True, applies print-ready settings (larger font size).
    """
    
    # Get current axes if none provided
    if ax is None:
        ax = plt.gca()
    
    # Set font properties
    ax.set_fontfamily('Arial')
    
    # Configure axis appearance properties for professional look
    ax.set_box_on(True)                    # Display box around plot
    ax.tick_params(direction='out')        # Tick marks point outward
    ax.tick_params(length=6)               # Length of tick marks (in points)
    ax.tick_params(which='minor', length=0)  # Disable minor ticks
    ax.grid(axis='y', visible=False)       # Disable y-axis grid
    
    # Set line width for axes
    for spine in ax.spines.values():
        spine.set_linewidth(1)
    
    # Set figure background to white for clean appearance
    ax.figure.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Apply print-ready settings if requested
    if print_ready:
        # Increase font size for better readability in publications
        ax.tick_params(labelsize=14)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        if ax.get_title():
            ax.title.set_size(14)
        if ax.get_legend():
            ax.legend(fontsize=14) 