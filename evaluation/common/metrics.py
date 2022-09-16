# Metrics are meant to be calculated using xarray
# These are looking like a duck
import numpy as np
import xarray as xr


def mean_err(mod = None, obs = None, skipna=False, dim=None): 
    assert len(mod) == len(obs)
    return (mod - obs).mean(skipna=skipna, dim=dim)

def rmse(mod = None, obs = None, skipna=False, dim=None):
    assert len(mod) == len(obs)
    err = mod - obs
    return np.sqrt((err * err).mean(skipna=skipna, dim=dim))

def nrmse(mod = None, obs = None, skipna=False, dim=None):
    _rmse = rmse(mod=mod, obs=obs, skipna=skipna, dim=dim)
    return _rmse / obs.mean(skipna=skipna, dim=dim)

def max_abs_err(mod = None, obs = None, skipna=False, dim=None):
    assert len(mod) == len(obs)
    err = mod - obs
    return abs(err).max(skipna=skipna, dim=dim)

def max_abs_rel_err(mod = None, obs = None, skipna=True, dim=None):
    assert len(mod) == len(obs)
    abs_rel_err = abs((mod - obs) / obs)
    abs_rel_err = xr.where(np.isinf(abs_rel_err), np.nan, abs_rel_err)
    return abs_rel_err.max(skipna=skipna, dim=dim)

# taken from spotpy
def nse(mod = None, obs = None, skipna=False, dim=None):
    assert len(mod) == len(obs)
    mean_obs = obs.mean(skipna=skipna, dim=dim)
    err = mod - obs
    numerator = (err * err).sum(skipna=skipna, dim=dim)
    obs_mean_resid = obs - mean_obs
    denominator = (obs_mean_resid * obs_mean_resid).sum(skipna=skipna, dim=dim)
    return (1 - (numerator / denominator))

stat_dict = {
    'mean_err': mean_err, 
    'rmse': rmse, 
    'nrmse': nrmse,
    'max_abs_err': max_abs_err,
    'max_abs_rel_err': max_abs_rel_err,
    'nse': nse,
}

stat_lims = {
    'mean_error': (None, None),
    'rmse': (0.000, None),
    'nrmse': (0.000, None),
    'max_abs_err': (0.000, None),
    'max_abs_rel_err': (0.000, None),
    'nse': (None, 1),
}
