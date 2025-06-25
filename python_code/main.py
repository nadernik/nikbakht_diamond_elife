"""
Main script for rat vision analysis.

This script provides a command-line interface to run the behavioral and spectral
analyses from the Nikbakht & Diamond (2021) study.
"""

import argparse
import sys
import os

def main():
    """Main function with command-line interface."""
    
    parser = argparse.ArgumentParser(
        description='Rat vision analysis from Nikbakht & Diamond (2021)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --behavior                    # Run behavioral analysis only
  python main.py --spectral                    # Run spectral analysis only
  python main.py --all                         # Run all analyses
  python main.py --behavior --save-figs        # Run behavioral analysis and save figures
        """
    )
    
    parser.add_argument('--behavior', action='store_true',
                       help='Run behavioral analysis (plot_behavior.py)')
    parser.add_argument('--spectral', action='store_true',
                       help='Run spectral analysis (plot_spect.py)')
    parser.add_argument('--all', action='store_true',
                       help='Run both behavioral and spectral analyses')
    parser.add_argument('--save-figs', action='store_true',
                       help='Save figures to files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='figures',
                       help='Output directory for saved figures (default: figures)')
    
    args = parser.parse_args()
    
    # Check if any analysis was specified
    if not any([args.behavior, args.spectral, args.all]):
        parser.print_help()
        return
    
    # Create output directory if saving figures
    if args.save_figs:
        os.makedirs(args.output_dir, exist_ok=True)
        # Set matplotlib to non-interactive backend
        import matplotlib
        matplotlib.use('Agg')
    
    # Run behavioral analysis
    if args.behavior or args.all:
        print("Running behavioral analysis...")
        try:
            from plot_behavior import main as run_behavior
            run_behavior()
            print("Behavioral analysis completed successfully.")
        except Exception as e:
            print(f"Error in behavioral analysis: {e}")
            return
    
    # Run spectral analysis
    if args.spectral or args.all:
        print("Running spectral analysis...")
        try:
            from plot_spect import main as run_spectral
            run_spectral()
            print("Spectral analysis completed successfully.")
        except Exception as e:
            print(f"Error in spectral analysis: {e}")
            return
    
    print("All analyses completed successfully!")

if __name__ == "__main__":
    main() 