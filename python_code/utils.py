"""
Utility functions for rat vision analysis.
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm

def cumulative_gaussian_lapse(params, x):
    """
    Cumulative Gaussian function with lapse rate.
    
    Parameters:
    -----------
    params : array-like
        [mu, sigma, gamma, lambda] where:
        - mu: mean (threshold)
        - sigma: standard deviation (slope)
        - gamma: guess rate (lower asymptote)
        - lambda: lapse rate (upper asymptote)
    x : array-like
        Input values (e.g., angles)
    
    Returns:
    --------
    array-like: Psychometric function values
    """
    mu, sigma, gamma, lambda_param = params
    
    # Cumulative Gaussian
    z = (x - mu) / sigma
    p = gamma + (1 - gamma - lambda_param) * norm.cdf(z)
    
    return p

def fit_psychometric_function(angles, responses, n_trials):
    """
    Fit psychometric function to behavioral data.
    
    Parameters:
    -----------
    angles : array-like
        Stimulus angles
    responses : array-like
        Number of "vertical" responses
    n_trials : array-like
        Total number of trials per angle
    
    Returns:
    --------
    tuple: (fitted_params, success_flag)
    """
    # Initial parameter guesses
    p0 = [45, 10, 0.5, 0.05]  # mu, sigma, gamma, lambda
    
    # Bounds for parameters
    bounds = [(0, 90), (0.1, 50), (0, 0.5), (0, 0.1)]
    
    # Define objective function (negative log likelihood)
    def objective(params):
        p_pred = cumulative_gaussian_lapse(params, angles)
        # Avoid log(0) issues
        p_pred = np.clip(p_pred, 1e-6, 1-1e-6)
        ll = responses * np.log(p_pred) + (n_trials - responses) * np.log(1 - p_pred)
        return -np.sum(ll)
    
    # Fit using scipy.optimize
    try:
        result = optimize.minimize(objective, p0, bounds=bounds, method='L-BFGS-B')
        return result.x, result.success
    except:
        return p0, False

def load_matlab_data(filename):
    """
    Load MATLAB .mat file data.
    
    Parameters:
    -----------
    filename : str
        Path to .mat file
    
    Returns:
    --------
    dict: Loaded data
    """
    try:
        from scipy import io
        return io.loadmat(filename)
    except:
        try:
            import h5py
            data = {}
            with h5py.File(filename, 'r') as f:
                for key in f.keys():
                    data[key] = np.array(f[key])
            return data
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None 