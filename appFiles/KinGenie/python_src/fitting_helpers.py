"""
This module contains helper functions for fitting kinetic binding models to surface-based data
Author: Osvaldo Burastero
LICENSE:

MIT License

Copyright (c) [2025] [Osvaldo Burastero]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import  numpy as np
import pandas as pd

from scipy           import stats

from scipy.optimize  import curve_fit
from scipy.optimize  import minimize_scalar
from scipy.integrate import solve_ivp
from scipy.linalg    import expm
from scipy.linalg    import solve

def concat_signal_lst(signal_lst):

    """
    Concatenate a list of signals into a single array

    Args:
        signal_lst (list): List of signals to concatenate, each signal is a numpy array

    Returns:
        allSignal (np.ndarray): Concatenated signal
    """

    try:
        allSignal = np.concatenate(signal_lst)
    except:
        allSignal = signal_lst

    return (allSignal)

def rss_p(rrs0, n, p, alfa):

    """
    Given the residuals of the best fitted model,
    compute the desired residual sum of squares for a 1-alpha confidence interval
    This is used to compute asymmetric confidence intervals for the fitted parameters

    Args:
        rrs0 (float): residual sum of squares of the model with the best fit
        n (int): number of data points
        p (int): number of parameters
        alfa (float): desired confidence interval

    Returns:
        rss (float): residual sum of squares for the desired confidence interval
    """

    critical_value = stats.f.ppf(q=1 - alfa, dfn=1, dfd=n - p)

    return rrs0 * (1 + critical_value / (n - p))

def get_rss(y, y_fit):

    """
    Compute the residual sum of squares

    Args:

        y (np.ndarray): observed values
        y_fit (np.ndarray): fitted values

    Returns:

        rss (np.ndarray): residual sum of squares

    """

    residuals = y - y_fit
    rss       = np.sum(residuals ** 2)

    return rss

def get_desired_rss(y, y_fit, n, p,alpha=0.05):

    """
    Given the observed and fitted data,
    find the residual sum of squares required for a 1-alpha confidence interval

    Args:

        y (np.ndarray): observed values
        y_fit (np.ndarray): fitted values
        n (int): number of data points
        p (int): number of parameters
        alpha (float): desired confidence interval

    Returns:

        rss (np.ndarray): residual sum of squares

    """

    rss = get_rss(y, y_fit)

    return rss_p(rss, n, p, alpha)

def detect_time_list_continuos(assoc_time_lst,disso_time_lst,tolerance=3):

    """
    Detect which association steps come directly after a dissociation step
    Useful for single-cylce kinetics

    Args:
        assoc_time_lst (list): List of association time arrays
        disso_time_lst (list): List of dissociation time arrays
        tolerance (float):     Tolerance for the time difference (in seconds)
    Returns:
        continuos (list): List of booleans indicating if the association phase have a dissociation phase just before
    """

    continuos = []

    for i,element in enumerate(assoc_time_lst):

        if i == 0:

            continuos.append(element[0] < tolerance)

        else:

            prev_time = disso_time_lst[i-1][-1]
            continuos.append(element[0] < prev_time+tolerance)

    return continuos

def steady_state_one_site(C,Rmax,Kd):

    """
    Calculate the steady state signal for a given concentration of ligand
    If the concentration is zero, the signal is zero
    If the concentration is infinite, the signal is Rmax
    The units of Kd must match the units of C (ligand concentration)
    Rmax depends on the amount of loaded receptor

    Args:
        C (np.ndarray): Concentration of the analyte
        Rmax (float): Maximum response of the analyte
        Kd (float): Equilibrium dissociation constant
    Returns:
        signal (np.ndarray): Steady state signal, according to the given parameters
    """

    C    = np.array(C)
    Rmax = np.array(Rmax)
    signal = Rmax*C/(Kd + C)
    return signal

def fit_steady_state_one_site(signal_lst, ligand_lst,initial_parameters,
                              low_bounds, high_bounds,fixed_Kd = False,Kd_value = None):

    """
    Fits a one-site binding model to a set of steady state signals
    Args:
        signal_lst (list): List of signals to fit, each signal is a numpy array
        ligand_lst (list): List of ligand concentrations, each concentration is a numpy array
        initial_parameters (list): Initial guess for the parameters
        low_bounds (list): Lower bounds for the parameters
        high_bounds (list): Upper bounds for the parameters
        fixed_Kd (bool): If True, Kd is fixed to Kd_value
        Kd_value (float): Value of Kd to use if fixed_Kd is True
    Returns:
        global_fit_params (list): Fitted parameters
        cov (np.ndarray): Covariance matrix of the fitted parameters
        fitted_values (list): Fitted values for each signal, same dimensions as signal_lst
    """

    all_signal = concat_signal_lst(signal_lst)

    start = 1
    if fixed_Kd:

        # relax bounds
        low_bounds  = [x / 5 for x in low_bounds]
        high_bounds = [x * 5 for x in high_bounds]
        start       = 0

    def fit_fx(dummyVariable, *args):

        Kd = Kd_value if fixed_Kd else args[0]

        Rmax_all = args[start:]

        signal = [steady_state_one_site(C,Rmax,Kd) for C,Rmax in zip(ligand_lst,Rmax_all)]

        return np.concatenate(signal, axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    Kd = Kd_value if fixed_Kd else global_fit_params[0]

    Rmax_all      = global_fit_params[start:]
    fitted_values = [steady_state_one_site(C, Rmax, Kd) for C, Rmax in zip(ligand_lst, Rmax_all)]

    return global_fit_params, cov, fitted_values

def steady_state_one_site_asymmetric_ci95(kd_estimated,signal_lst, ligand_lst,initial_parameters,
                                          low_bounds, high_bounds,rss_desired):

    """
    Calculate the asymmetric confidence interval for the steady-state signal
    Args:
        kd_estimated (float): Estimated Kd value
        signal_lst (list): List of signals to fit, each signal is a numpy array
        ligand_lst (list): List of ligand concentrations, each concentration is a numpy array
        initial_parameters (list): Initial guess for the parameters (without the Kd!)
        low_bounds (list): Lower bounds for the parameters (without the Kd!)
        high_bounds (list): Upper bounds for the parameters (without the Kd!)
        rss_desired (float): Maximum residual sum of squares

    Returns:
        ci95 (np.ndarray): 95 % asymmetric confidence interval for the Kd value

    """

    def f_to_optimize(Kd):

        fit_params, _, fit_vals = fit_steady_state_one_site(signal_lst, ligand_lst,
                                                  initial_parameters,
                                                  low_bounds,
                                                  high_bounds,
                                                  fixed_Kd = True,
                                                  Kd_value = Kd)

        rss = get_rss(concat_signal_lst(signal_lst), concat_signal_lst(fit_vals))

        return np.abs(rss - rss_desired)

    boundsMin = np.array([kd_estimated/5e4,kd_estimated])
    boundsMax = np.array([kd_estimated,kd_estimated*5e4])

    kd_min95 = minimize_scalar(f_to_optimize, bounds=boundsMin,method='bounded')
    kd_max95 = minimize_scalar(f_to_optimize, bounds=boundsMax,method='bounded')

    kd_min95, kd_max95 = kd_min95.x, kd_max95.x

    ci95 = np.array([kd_min95, kd_max95])

    return ci95

def ode_one_site_association(t,s,s_max,koff,Kd,analyte_conc):

    """
    Args:
        t (float): Time variable.
        s (float): Signal at time t.
        s_max (float): Maximum signal.
        koff (float): Dissociation rate constant.
        Kd (float): Equilibrium Dissociation constant.
        analyte_conc (float): Concentration of the analyte.

    Returns:

          dsdt (float): signal over time
    """

    dsdt = koff / Kd * analyte_conc * (s_max - s) - koff * s

    return  dsdt

def ode_one_site_dissociation(t,s,koff):

    """
    Args:
        t (float): Time variable.
        s (float): Signal at time t.
        koff (float): Dissociation rate constant.

    Returns:

          dsdt (float): signal over time

    """

    dsdt =  - koff * s

    return  dsdt

def ode_one_site_mass_transport_association(t,y,params,analyte_conc):

    """
    ODE for Mass Transport-limited binding

    See https://www.cell.com/AJHG/fulltext/S0006-3495(07)70982-7

    Args:

        t (float): Time variable.
        y (list): List of the state variables [s1,cs]
        params (list): List of the parameters [K_d, k_off, k_tr, s_max]

    Returns:

        ds1_dt (float): signal over time
        dcs_dt (float): analyte surface concentration over time

    """

    s1, cs = y
    K_d, k_off, k_tr, s_max = params

    k_on = k_off / K_d

    c0 = analyte_conc - cs

    ds1_dt = k_on * cs * (s_max - s1) - k_off * s1
    dcs_dt = k_tr * (c0 - cs) - ds1_dt

    return [ds1_dt, dcs_dt]

def ode_one_site_mass_transport_dissociation(t,s1,params):

    """
    ODE for Mass Transport-limited binding

    See Equation 7 from https://pmc.ncbi.nlm.nih.gov/articles/PMC4134667/

    Args:

        t (float): Time variable.
        y (list): List of the state variables [s1]
        params (list): List of the parameters [K_d, k_off, k_tr, s_max]

    Returns:

        ds_dt (float): signal over time

    """

    K_d, k_off, k_tr, s_max = params

    k_on = k_off / K_d

    ds_dt = (-k_off * s1) / (1 + (k_on / k_tr) * (s_max - s1))

    return ds_dt

def solve_ode_one_site_mass_transport_association(t,s1_0,cs_0,analyte_conc,K_d,k_off,k_tr,s_max,t0=0):

    """
    Solves the ODE for the one site mass transport model - association phase

    Args:
        t (np.ndarray): Time variable
        s1_0 (float):   Initial signal
        cs_0 (float):   Initial analyte surface concentration
        analyte_conc (float): Concentration of the analyte
        K_d (float):    Equilibrium dissociation constant
        k_off (float):  Dissociation rate constant
        k_tr (float):   Mass transport rate constant
        s_max (float): Maximum signal
        t0 (float):     Initial time offset
    Returns:
        signal (np.ndarray): Signal overtime
    """

    t = t + t0

    out = solve_ivp(ode_one_site_mass_transport_association,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=[s1_0,cs_0],args=([K_d,k_off,k_tr,s_max],analyte_conc),method="LSODA")

    signal = out.y[0]
    return signal

def solve_ode_one_site_mass_transport_dissociation(t,s1_0,K_d,k_off,k_tr,s_max,t0=0):

    """
    Solves the ODE for the one site mass transport model - dissociation phase

    Args:
        t (np.ndarray): Time variable
        s1_0 (float):   Initial signal
        K_d (float):    Equilibrium dissociation constant
        k_off (float):  Dissociation rate constant
        k_tr (float):   Mass transport rate constant
        s_max (float): Maximum signal
        t0 (float):     Initial time offset
    Returns:
        signal (np.ndarray): Signal overtime
    """

    t = t + t0

    out = solve_ivp(ode_one_site_mass_transport_dissociation,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=[s1_0],args=([K_d,k_off,k_tr,s_max],),method="LSODA")

    signal = out.y[0]
    return signal

def ode_mixture_analyte_association(t,Ris,C_TOT,Fis,Ris_max,koffs,Kds):

    """
    We assume a Langmuir 1:1 interaction between each analyte (Ai) and the ligand (L):

    Ai + L <-> AiL, ∀i in [1,N]

    Based on https://doi.org/10.1038/s41598-022-18450-y

    This model is equivalent to the heterogenous analyte model

    The parameters to fit are the factor weights (Fis), the maximum response of each analyte (Ris_max),
    the dissociation rate constants (koffs), and the equilibrium dissociation constants (Kds)

    Args:

        t (np.ndarray):         Time variable
        Ris (list):             Response of each analyte
        C_TOT (float):          Total concentration of the analyte mixture
        Fi (np.ndarray):        Factor weights
        Ris_max (np.ndarray):   Maximum response of each analyte
        koffs (np.ndarray):     Dissociation rate constants
        Kds (np.ndarray):       Equilibrium dissociation constants

    Returns:
        dRis (np.ndarray):     Rate of change of the response of each analyte

    """

    kons = koffs / Kds

    dRis = np.zeros(len(Ris))

    choclo = (1 - np.sum(Ris / Ris_max))

    for i in range(len(Ris)):

        r_i    = Ris[i]
        kon_i  = kons[i]
        F_i    = Fis[i]
        Rmax_i = Ris_max[i]
        koff_i = koffs[i]

        dRi_dt = kon_i * F_i * C_TOT * Rmax_i * choclo - koff_i * r_i
        dRis[i] = dRi_dt

    return dRis

def solve_ode_mixture_analyte_association(t,Ris0,C_TOT,Fis,Ris_max,koffs,Kds,t0=0):

    """
    Solves the ODE for the mixture analyte model - association phase
    Args:
        t (np.ndarray): Time variable
        Ris0 (list):    Initial response of each analyte
        C_TOT (float):  Total concentration of the analyte mixture
        Fis (np.ndarray): Factor weights
        Ris_max (np.ndarray): Maximum response of each analyte
        koffs (np.ndarray): Dissociation rate constants
        Kds (np.ndarray): Equilibrium dissociation constants
        t0 (float):     Initial time offset

    Returns:
        out.y (np.ndarray): Response of EACH ANALYTE over time
    """

    t = t + t0

    out = solve_ivp(ode_mixture_analyte_association,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=Ris0,args=(C_TOT,Fis,Ris_max,koffs,Kds),method="LSODA")

    return out.y

def ode_mixture_analyte_dissociation(t,Ris,koffs):

    """
    We assume a Langmuir 1:1 interaction between each analyte (Ai) and the ligand (L):

    Ai + L <-> AiL, ∀i in [1,N]

    Based on https://doi.org/10.1038/s41598-022-18450-y

    This model is equivalent to the heterogenous analyte model

    The parameters to fit are the dissociation rate constants (koffs)

    Args:

        t (np.ndarray):         Time variable
        Ris (list):             Response of each analyte
        koffs (np.ndarray):     Dissociation rate constants

    Returns:
        dRis_dt (np.ndarray):     Rate of change of the response of each analyte
    """

    dRis_dt = - koffs * Ris

    return dRis_dt

def solve_ode_mixture_analyte_dissociation(t,Ris0,koffs,t0=0):

    """
    Solves the ODE for the mixture analyte model - dissociation phase
    Args:
        t (np.ndarray): Time variable
        Ris0 (list):    Initial response of each analyte
        koffs (np.ndarray): Dissociation rate constants
        t0 (float):     Initial time offset

    Returns:
        out.y (np.ndarray): Response of EACH ANALYTE over time
    """

    t = t + t0

    out = solve_ivp(ode_mixture_analyte_dissociation,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=Ris0,args=(koffs,),method="LSODA")

    return out.y

def solve_ode_one_site_association(t,s0,s_max,koff,Kd,analyte_conc,t0=0):

    """
    Solves the ODE for the one site model - association phase

    Args:
        t (np.ndarray): Time variable
        s0 (float):     Initial signal
        s_max (float):  Maximum signal
        koff (float):   Dissociation rate constant
        Kd (float):     Equilibrium dissociation constant
        analyte_conc (float): Concentration of the analyte
        t0 (float):     Initial time offset
    Returns:
        out.y[0] (np.ndarray): Signal over time
    """

    t = t + t0

    out = solve_ivp(ode_one_site_association,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=[s0],args=(s_max,koff,Kd,analyte_conc),method="LSODA")

    return out.y[0]

def one_site_association_analytical(t,s0,s_max,k_off,Kd,analyte_conc,t0=0):

    """
    Analytical solution for the one site association model

    Args:
        t (np.ndarray): Time variable
        s0 (float):     Initial signal
        s_max (float):  Maximum signal
        k_off (float):  Dissociation rate constant
        Kd (float):     Equilibrium dissociation constant
        analyte_conc (float): Concentration of the analyte
        t0 (float):     Initial time offset

    Returns:
        s_t (np.ndarray): Signal over time
    """

    A = analyte_conc
    t = t - t0
    # Precompute constants
    rate = k_off * (A + Kd) / Kd
    s_eq = (A / (A + Kd)) * s_max

    # Solution for s(t)
    s_t = s_eq + (s0 - s_eq) * np.exp(-rate * t)

    return s_t

def one_site_dissociation_analytical(t,s0,k_off,t0=0):

    """
    Analytical solution for the one site model - dissociation phase

    Args:
        t (np.ndarray): Time variable
        s0 (float):     Initial signal
        k_off (float):  Dissociation rate constant
        t0 (float):     Initial time offset

    Returns:
        s_t (np.ndarray): Signal over time
    """

    t   = t + t0
    s_t = s0 * np.exp(-k_off*t)

    return s_t

def solve_ode_one_site_dissociation(t,s0,koff,t0=0):

    """
    Solves the ODE for the one site model - dissociation phase
    Args:
        t (np.ndarray): Time variable
        s0 (float):     Initial signal
        k_off (float):  Dissociation rate constant
        t0 (float):     Initial time offset
    Returns:
        out.y[0] (np.ndarray):  Signal over time
    """

    t = t + t0

    out = solve_ivp(ode_one_site_dissociation,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=[s0],args=(koff,),method="LSODA")

    return out.y[0]

def expand_parameter_list(parameter_lst,id_list):

    """

    Given a list of n-parameters, such as [1,3] and another list
    containing IDs, such as [0,0,0,0,0,1,1]
    we will create a new list where the elements are repeated
    according to the IDs. In this case: [1,1,1,1,1,3,3]

    Args:
        parameter_lst (list):   n-Parameters
        id_list (list):         m-IDs
    Returns:
        expanded_parameters (list): m-Parameters according to the IDs
    """

    expanded_parameters = [parameter_lst[i] for i in id_list]
    return expanded_parameters

def fit_one_site_association(signal_lst, time_lst, analyte_conc_lst,
                             initial_parameters,low_bounds, high_bounds,
                             smax_idx=None,
                             shared_smax = False,
                             fixed_t0 = True,
                             fixed_Kd = False, Kd_value = None,
                             fixed_koff = False, koff_value = None):

    """
    Global fit to a list of association traces - one-to-one binding model

    Args:
        signal_lst (list): List of signals to fit, each signal is a numpy array
        time_lst (list): List of time arrays
        analyte_conc_lst (list): List of analyte concentrations, each element is a numpy array
        initial_parameters (list): Initial guess for the parameters
        low_bounds (list): Lower bounds for the parameters
        high_bounds (list): Upper bounds for the parameters
        smax_idx (list): List of indices for the s_max parameters, used if shared_smax is TRUE
        shared_smax (bool): If True, the s_max parameters are shared between the signals
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_Kd (bool): If True, Kd is fixed to Kd_value
        Kd_value (float): Value of Kd to use if fixed_Kd is True
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True

    Returns:
        global_fit_params (list): Fitted parameters
        cov (np.ndarray): Covariance matrix of the fitted parameters
        fitted_values (list): Fitted values for each signal, same dimensions as signal_lst

    """

    all_signal = concat_signal_lst(signal_lst)

    time_lst = [np.array(t) for t in time_lst]

    start = 3 - sum([fixed_t0, fixed_Kd, fixed_koff])

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
            Kd, Koff, T0, Smax1, Smax2, Smax3, ...
        """

        Kd = Kd_value if fixed_Kd else args[0]

        Koff = koff_value if fixed_koff else args[1 - sum([fixed_Kd])]

        t0 = args[2 - sum([fixed_Kd, fixed_koff])] if not fixed_t0 else 0

        signal = []

        for i in range(len(time_lst)):

            t            = time_lst[i]
            analyte_conc = analyte_conc_lst[i]

            s_max        = args[start+smax_idx[i]] if shared_smax else args[start + i]

            res = one_site_association_analytical(t,0,s_max,Koff,Kd,analyte_conc,t0)

            signal.append(res)

        return np.concatenate(signal, axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    predicted_curve = fit_fx(1, *global_fit_params)

    fitted_values_assoc = []

    init = 0
    for t in time_lst:

        end  = init + len(t)

        predicted_curve_i = predicted_curve[init:end]

        fitted_values_assoc.append(predicted_curve_i)

        init = end

    return global_fit_params, cov, fitted_values_assoc

def fit_one_site_dissociation(signal_lst, time_lst,
                             initial_parameters,low_bounds, high_bounds,
                             fixed_t0   = True,
                             fixed_koff = False, koff_value = None,
                             fit_s0=True):

    """
    Global fit to a list of dissociation traces - one-to-one binding model

    Args:
        signal_lst (list): List of signals to fit, each signal is a numpy array
        time_lst (list): List of time arrays
        initial_parameters (list): Initial guess for the parameters
        low_bounds (list): Lower bounds for the parameters
        high_bounds (list): Upper bounds for the parameters
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True
        fit_s0 (bool): If True, s0 is fitted, otherwise it is fixed to the first value of the signal

    Returns:
        global_fit_params (list): Fitted parameters
        cov (np.ndarray): Covariance matrix of the fitted parameters
        fitted_values (list): Fitted values for each signal, same dimensions as signal_lst

    """

    all_signal = concat_signal_lst(signal_lst)

    time_lst = [np.array(t) for t in time_lst]

    # Start at time zero
    time_lst = [t - t[0] for t in time_lst]

    if not fit_s0:

        s0_all = [s[0] for s in signal_lst]

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
           Koff, T0, S01, S02, S03, ...
        """

        Koff = koff_value if fixed_koff else args[0]

        t0 = args[1 - sum([fixed_koff])] if not fixed_t0 else 0

        if fit_s0:

            s0_all = args[( 2 - sum([fixed_koff,fixed_t0]) ):]

        signal  = [one_site_dissociation_analytical(t,s0,Koff,t0) for t,s0 in zip(time_lst,s0_all)]

        return np.concatenate(signal, axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    Koff = koff_value if fixed_koff else global_fit_params[0]

    if fit_s0:

        s0_all = global_fit_params[( 2 - sum([fixed_koff,fixed_t0]) ):]

    fitted_values = [solve_ode_one_site_dissociation(t,s0,Koff) for t,s0 in zip(time_lst,s0_all)]

    return global_fit_params, cov, fitted_values

def guess_initial_signal(assoc_time_lst, assoc_signal_lst):

    """
    Guess the initial signal for each signal in the list, by trying to fit a single exponential

    Args:
        assoc_time_lst (list): List of association time arrays
        assoc_signal_lst (list): List of association signal arrays

    Returns:
        s0s (list): List of initial signals for each association signal
    """
    s0s = []

    for t,y in zip(assoc_time_lst,assoc_signal_lst):

        t = t - t[0]
        y = y[t < 10]
        t = t[t < 10]

        try:

            a0,a1,kobs = fit_single_exponential(y,t)
            s0s.append(a0 + a1)

        except:

            s0s.append(y[0])

    return s0s

def fit_one_site_assoc_and_disso(assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                 disso_signal_lst, disso_time_lst,
                                 initial_parameters,low_bounds, high_bounds,
                                 smax_idx=None,
                                 shared_smax = False,
                                 fixed_t0 = True,
                                 fixed_Kd = False,   Kd_value = None,
                                 fixed_koff = False, koff_value = None):

    """
    Global fit to a set of association and dissociation traces - one-to-one binding model
    Args:
        assoc_signal_lst (list): List of association signals to fit, each signal is a numpy array
        assoc_time_lst (list): List of association time arrays
        analyte_conc_lst (list): List of analyte concentrations, each element is a numpy array
        disso_signal_lst (list): List of dissociation signals to fit, each signal is a numpy array
        disso_time_lst (list): List of dissociation time arrays
        initial_parameters (list): Initial guess for the parameters
        low_bounds (list): Lower bounds for the parameters
        high_bounds (list): Upper bounds for the parameters
        smax_idx (list): List of indices for the s_max parameters, used if shared_smax is TRUE
        shared_smax (bool): If True, the s_max parameters are shared between traces
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_Kd (bool): If True, Kd is fixed to Kd_value
        Kd_value (float): Value of Kd to use if fixed_Kd is True
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True
    Returns:
        global_fit_params (list): Fitted parameters
        cov (np.ndarray): Covariance matrix of the fitted parameters
        fitted_values_assoc (list): Fitted values for each association signal, same dimensions as assoc_signal_lst
        fitted_values_disso (list): Fitted values for each dissociation signal, same dimensions as disso_signal_lst
    """

    # Set a flag for the association that were done after a dissociation step
    # We do this empirically by detecting if the initial time is greater than 2 seconds
    initial_signal_at_zero = [time[0] < 2 for time in assoc_time_lst]

    all_signal_assoc = concat_signal_lst(assoc_signal_lst)
    all_signal_disso = concat_signal_lst(disso_signal_lst)

    time_lst_assoc = [np.array(t) for t in assoc_time_lst]
    time_lst_disso = [np.array(t) for t in disso_time_lst]

    time_lst_disso = [t - t[0] for t in time_lst_disso]

    start = 2 - sum([fixed_Kd, fixed_koff])

    n_t0s = len(np.unique(smax_idx))*(not fixed_t0)

    # The user may have deleted one of the association-dissociation steps
    # In that case, we can't use for the next step the previous association signal
    # We'll use the first signal value in that cases
    continuos_time = detect_time_list_continuos(time_lst_assoc,disso_time_lst)
    s0s            = guess_initial_signal(assoc_time_lst, assoc_signal_lst)

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
            Kd, Koff, T0, Smax1, Smax2, Smax3, ...
        """

        Kd   = Kd_value if fixed_Kd else args[0]

        Koff = koff_value if fixed_koff else args[1 - sum([fixed_Kd])]

        signal_a = []
        signal_d = []

        i = 0
        for t_assoc,t_dissoc in zip(time_lst_assoc,time_lst_disso):

            analyte_conc = analyte_conc_lst[i]

            t0     = args[start + smax_idx[i]] if not fixed_t0 else 0

            s_max  = args[start+smax_idx[i]+n_t0s ]   if shared_smax  else args[start + i + n_t0s]

            # Set s0 to zero if initial_signal_at_zero, otherwise set it to the value of the previous dissociation
            if np.logical_or(i == 0,initial_signal_at_zero[i]) and continuos_time[i]:

                s0 = 0

            elif continuos_time[i]:

                s0 = signal_d[i-1][-1]

            else:

                s0 = s0s[i] # Use the first data point from a single-exponential fit

            y_pred = one_site_association_analytical(t_assoc-t_assoc[0],s0,s_max,Koff,Kd,analyte_conc,t0)

            signal_a.append(y_pred)

            s0 = y_pred[-1]

            y_pred = one_site_dissociation_analytical(t_dissoc, s0, Koff)
            signal_d.append(y_pred)

            i += 1

        signal_assoc_fit  = np.concatenate(signal_a, axis=0)
        signal_dissoc_fit = np.concatenate(signal_d, axis=0)

        return np.concatenate([signal_assoc_fit,signal_dissoc_fit], axis=0)

    all_signal = np.concatenate([all_signal_assoc,all_signal_disso], axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    predicted_curve = fit_fx(1, *global_fit_params)

    fitted_values_assoc = []

    init = 0
    for t in time_lst_assoc:

        end  = init + len(t)

        predicted_curve_i = predicted_curve[init:end]

        fitted_values_assoc.append(predicted_curve_i)

        init = end

    fitted_values_disso = []

    for t in time_lst_disso:

        end  = init + len(t)

        predicted_curve_i = predicted_curve[init:end]

        fitted_values_disso.append(predicted_curve_i)

        init = end

    return global_fit_params, cov, fitted_values_assoc, fitted_values_disso

def fit_one_site_assoc_and_disso_ktr(assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                 disso_signal_lst, disso_time_lst,
                                 initial_parameters,low_bounds, high_bounds,
                                 smax_idx=None,
                                 shared_smax = False,
                                 fixed_t0 = True,
                                 fixed_Kd = False,   Kd_value = None,
                                 fixed_koff = False, koff_value = None):

    """
    Global fit to a set of association and dissociation traces - one-to-one with mass transport limitation binding model
    Args:
        assoc_signal_lst (list): List of association signals to fit, each signal is a numpy array
        assoc_time_lst (list): List of association time arrays
        analyte_conc_lst (list): List of analyte concentrations, each element is a numpy array
        disso_signal_lst (list): List of dissociation signals to fit, each signal is a numpy array
        disso_time_lst (list): List of dissociation time arrays
        initial_parameters (list): Initial guess for the parameters
        low_bounds (list): Lower bounds for the parameters
        high_bounds (list): Upper bounds for the parameters
        smax_idx (list): List of indices for the s_max parameters, used if shared_smax is TRUE
        shared_smax (bool): If True, the s_max parameters are shared between traces
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_Kd (bool): If True, Kd is fixed to Kd_value
        Kd_value (float): Value of Kd to use if fixed_Kd is True
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True
    Returns:
        global_fit_params (list): Fitted parameters
        cov (np.ndarray): Covariance matrix of the fitted parameters
        fitted_values_assoc (list): Fitted values for each association signal, same dimensions as assoc_signal_lst
        fitted_values_disso (list): Fitted values for each dissociation signal, same dimensions as disso_signal_lst
    """

    all_signal_assoc = concat_signal_lst(assoc_signal_lst)
    all_signal_disso = concat_signal_lst(disso_signal_lst)

    time_lst_assoc = [np.array(t) for t in assoc_time_lst]
    time_lst_disso = [np.array(t) for t in disso_time_lst]

    initial_signal_at_zero = [time[0] < 2 for time in assoc_time_lst]

    time_lst_disso = [t - t[0] for t in time_lst_disso]

    start = 2 - sum([fixed_Kd, fixed_koff])

    n_t0s = len(np.unique(smax_idx))*(not fixed_t0)
    n_ktr = len(np.unique(smax_idx))

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
            Kd, Koff, T0, Smax1, Smax2, Smax3, ...
        """

        Kd   = Kd_value if fixed_Kd else args[0]

        Koff = koff_value if fixed_koff else args[1 - sum([fixed_Kd])]

        signal_a = []
        signal_d = []

        i = 0
        for t_assoc,t_dissoc in zip(time_lst_assoc,time_lst_disso):

            analyte_conc = analyte_conc_lst[i]

            t0     = args[start + smax_idx[i]] if not fixed_t0 else 0

            k_tr   = args[start + smax_idx[i]] if fixed_t0 else args[start + n_t0s + smax_idx[i]]

            s_max  = args[start+smax_idx[i]+n_t0s+n_ktr ]   if shared_smax  else args[start + i + n_t0s + n_ktr]

            # Set s0 to zero if initial_signal_at_zero, otherwise set it to the value of the previous dissociation
            if np.logical_or(i == 0,initial_signal_at_zero[i]):

                s0 = 0

            else:

                s0 = signal_d[i-1][-1]

            y_pred    = solve_ode_one_site_mass_transport_association(t_assoc-t_assoc[0],s0,analyte_conc/2,analyte_conc,Kd,Koff,k_tr,s_max,t0)

            signal_a.append(y_pred)

            s0 = y_pred[-1]

            y_pred = solve_ode_one_site_mass_transport_dissociation(t_dissoc,s0,Kd,Koff,k_tr,s_max)

            signal_d.append(y_pred)

            i += 1

        signal_assoc_fit  = np.concatenate(signal_a, axis=0)
        signal_dissoc_fit = np.concatenate(signal_d, axis=0)

        return np.concatenate([signal_assoc_fit,signal_dissoc_fit], axis=0)

    all_signal = np.concatenate([all_signal_assoc,all_signal_disso], axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    predicted_curve = fit_fx(1, *global_fit_params)

    fitted_values_assoc = []

    init = 0
    for t in time_lst_assoc:

        end  = init + len(t)

        predicted_curve_i = predicted_curve[init:end]

        fitted_values_assoc.append(predicted_curve_i)

        init = end

    fitted_values_disso = []

    for t in time_lst_disso:

        end  = init + len(t)

        predicted_curve_i = predicted_curve[init:end]

        fitted_values_disso.append(predicted_curve_i)

        init = end

    return global_fit_params, cov, fitted_values_assoc, fitted_values_disso


def one_site_assoc_and_disso_asymmetric_ci95(kd_estimated,rss_desired,
                                             assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                             disso_signal_lst, disso_time_lst,
                                             initial_parameters,low_bounds, high_bounds,
                                             smax_idx=None,
                                             shared_smax = False,
                                             fixed_t0 = True,
                                             fixed_koff = False, koff_value = None):

    """
    Calculate the asymmetric confidence interval for the Kd value, given a desired RSS value
    Args:
        kd_estimated (float): Estimated Kd value
        rss_desired (float): Desired RSS value
        assoc_signal_lst (list): List of association signals to fit, each signal is a numpy array
        assoc_time_lst (list): List of association time arrays
        analyte_conc_lst (list): List of analyte concentrations, each element is a numpy array
        disso_signal_lst (list): List of dissociation signals to fit, each signal is a numpy array
        disso_time_lst (list): List of dissociation time arrays
        initial_parameters (list): Initial guess for the parameters, without the Kd value!
        low_bounds (list): Lower bounds for the parameters, without the Kd value!
        high_bounds (list): Upper bounds for the parameters, without the Kd value!
        smax_idx (list): List of indices for the s_max parameters, used if shared_smax is TRUE
        shared_smax (bool): If True, the s_max parameters are shared between traces
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True

    Returns:
        kd_min95 (float): Minimum Kd value for the 95% confidence interval
        kd_max95 (float): Maximum Kd value for the 95% confidence interval
    """

    boundsMax = np.array([kd_estimated*1e2, kd_estimated * 1e4]) * 1e3

    # Guess starting point for the upper bound
    test_factors = [1,1.1,1.5, 2, 5, 10, 25,50,100]
    for i,test_factor in enumerate(test_factors):

        if test_factor == 1:
            continue

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
            disso_signal_lst, disso_time_lst,
            initial_parameters, low_bounds, high_bounds,
            smax_idx=smax_idx,
            shared_smax=shared_smax,
            fixed_t0=fixed_t0,
            fixed_Kd=True, Kd_value=kd_estimated*test_factor,
            fixed_koff=fixed_koff, koff_value=koff_value)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        if rss1 + rss2 > rss_desired:
            boundsMax = np.array([kd_estimated*test_factors[i-1],
                                  kd_estimated *test_factor]) * 1e3  # Scale to nanomolar units
            break

    boundsMin = np.array([kd_estimated/1e4, kd_estimated / 1e2]) * 1e3

    # Guess starting point for the lower bound
    for i,test_factor in enumerate(test_factors):

        if test_factor == 1:
            continue

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
            disso_signal_lst, disso_time_lst,
            initial_parameters, low_bounds, high_bounds,
            smax_idx=smax_idx,
            shared_smax=shared_smax,
            fixed_t0=fixed_t0,
            fixed_Kd=True, Kd_value=kd_estimated/test_factor,
            fixed_koff=fixed_koff, koff_value=koff_value)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        if rss1 + rss2 > rss_desired:
            boundsMin = np.array([kd_estimated/test_factor,
                                  kd_estimated *test_factors[i-1]]) * 1e3  # Scale to nanomolar units
            break

    def f_to_optimize(Kd):

        Kd = Kd / 1e3 # Input Kd is given in nanomolar, so we scale back to micromolar units

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
                                                            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                                            disso_signal_lst, disso_time_lst,
                                                            initial_parameters,low_bounds, high_bounds,
                                                            smax_idx=smax_idx,
                                                            shared_smax=shared_smax,
                                                            fixed_t0=fixed_t0,
                                                            fixed_Kd=True,Kd_value=Kd,
                                                            fixed_koff=fixed_koff,koff_value=koff_value)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        return np.abs(rss1 + rss2 - rss_desired)

    kd_min95 = minimize_scalar(f_to_optimize, bounds=boundsMin,method='bounded')
    kd_max95 = minimize_scalar(f_to_optimize, bounds=boundsMax,method='bounded')

    kd_min95, kd_max95 = kd_min95.x, kd_max95.x

    # Rescale back the Kd to micromolar
    kd_min95, kd_max95 = kd_min95 / 1e3, kd_max95 / 1e3

    return kd_min95, kd_max95

def one_site_assoc_and_disso_asymmetric_ci95_koff(koff_estimated,rss_desired,
                                             assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                             disso_signal_lst, disso_time_lst,
                                             initial_parameters,low_bounds, high_bounds,
                                             smax_idx=None,
                                             shared_smax = False,
                                             fixed_t0 = True):
    """
    Calculate the asymmetric confidence interval for the koff value, given a desired RSS value
    Args:
        koff_estimated (float): Estimated koff value
        rss_desired (float): Desired RSS value
        assoc_signal_lst (list): List of association signals to fit, each signal is a numpy array
        assoc_time_lst (list): List of association time arrays
        analyte_conc_lst (list): List of analyte concentrations, each element is a numpy array
        disso_signal_lst (list): List of dissociation signals to fit, each signal is a numpy array
        disso_time_lst (list): List of dissociation time arrays
        initial_parameters (list): Initial guess for the parameters, without the Kd value!
        low_bounds (list): Lower bounds for the parameters, without the Kd value!
        high_bounds (list): Upper bounds for the parameters, without the Kd value!
        smax_idx (list): List of indices for the s_max parameters, used if shared_smax is TRUE
        shared_smax (bool): If True, the s_max parameters are shared between traces
        fixed_t0 (bool): If True, t0 is fixed to 0, otherwise we fit it
        fixed_koff (bool): If True, koff is fixed to koff_value
        koff_value (float): Value of koff to use if fixed_koff is True

    Returns:
        k_min95 (float): Minimum koff value for the 95% confidence interval
        k_max95 (float): Maximum koff value for the 95% confidence interval
    """

    boundsMax = np.array([koff_estimated*1e2, koff_estimated * 1e4]) * 1e3

    # Guess starting point for the upper bound
    test_factors = [1,1.02,1.1,1.5, 2, 5, 10, 25,50,100]
    for i,test_factor in enumerate(test_factors):

        if test_factor == 1:
            continue

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
            disso_signal_lst, disso_time_lst,
            initial_parameters, low_bounds, high_bounds,
            smax_idx=smax_idx,
            shared_smax=shared_smax,
            fixed_t0=fixed_t0,
            fixed_koff=True, koff_value=koff_estimated*test_factor)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        if rss1 + rss2 > rss_desired:
            boundsMax = np.array([koff_estimated*test_factors[i-1],
                                  koff_estimated *test_factor]) * 1e3
            break

    boundsMin = np.array([koff_estimated/1e4, koff_estimated / 1e2]) * 1e3

    # Guess starting point for the lower bound
    for i,test_factor in enumerate(test_factors):

        if test_factor == 1:
            continue

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
            disso_signal_lst, disso_time_lst,
            initial_parameters, low_bounds, high_bounds,
            smax_idx=smax_idx,
            shared_smax=shared_smax,
            fixed_t0=fixed_t0,
            fixed_koff=True, koff_value=koff_estimated*test_factor)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        if rss1 + rss2 > rss_desired:
            boundsMin = np.array([koff_estimated/test_factor,
                                  koff_estimated *test_factors[i-1]]) * 1e3
            break

    def f_to_optimize(Koff):

        Koff = Koff / 1e3

        _, _, fitted_values_assoc, fitted_values_disso = fit_one_site_assoc_and_disso(
                                                            assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                                            disso_signal_lst, disso_time_lst,
                                                            initial_parameters,low_bounds, high_bounds,
                                                            smax_idx=smax_idx,
                                                            shared_smax=shared_smax,
                                                            fixed_t0=fixed_t0,
                                                            fixed_koff=True,
                                                            koff_value=Koff)

        rss1 = get_rss(concat_signal_lst(assoc_signal_lst), concat_signal_lst(fitted_values_assoc))
        rss2 = get_rss(concat_signal_lst(disso_signal_lst), concat_signal_lst(fitted_values_disso))

        return np.abs(rss1 + rss2 - rss_desired)

    k_min95 = minimize_scalar(f_to_optimize, bounds=boundsMin,method='bounded')
    k_max95 = minimize_scalar(f_to_optimize, bounds=boundsMax,method='bounded')

    k_min95, k_max95 = k_min95.x, k_max95.x

    # Rescale back the k_off
    k_min95, k_max95 = k_min95 / 1e3, k_max95 / 1e3

    return k_min95, k_max95

def differential_matrix_association_induced_fit(koff,kon,kc,krev,a):

    """
    Association step
    set of differential equations for
    -dS/dt and dS1/dt

    S    = R([E2S] + [E1S]) ; Total signal
    S1   = R[E2S]           ; Signal produced by E2S
    Smax = R([Et])          ; Maximum signal

    Note: E is here the molecule attached to the surface, S is the subtrate
    The signal produced by E2S and E1S is the same.

    Args:
        koff (float): rate constant for E2S -> E2 + S
        kon (float): rate constant for E2 + S -> E2S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1
        a (float): concentration of the analyte (molecule being flown)

    Returns:
         matrix (np.ndarray): differential matrix
    """

    row1 = [-koff - kon * a, -koff]
    row2 = [-kc, -kc - krev]

    matrix = np.array([row1, row2])

    return matrix

def differential_matrix_association_conformational_selection(koff,kon,kc,krev,a):

    """
    Conformational selection model
    E1     <-> E2
    E2 + S <-> E2S

    Association step
    set of differential equations for
    (d[E1]*R)/dt and d[E2S]*R/dt. The signal is only given by E2S

    Args:
        koff (float): rate constant for E2S -> E2 + S
        kon (float): rate constant for E2 + S -> E2S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1
        a (float): concentration of the analyte (molecule being flown)
    Returns:
        matrix (np.ndarray): differential matrix
    """

    row1 = [-kc - krev, -krev]
    row2 = [-kon*a,-kon*a - koff]

    matrix = np.array([row1, row2])

    return matrix

def differential_matrix_dissociation_induced_fit(koff,kc,krev):

    """
    Dissociation step
    set of differential equations for
    -ds/dt and ds1/dt

    Args:
        koff (float): rate constant for E2S -> E2 + S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1

    Returns:
        matrix (np.ndarray): differential matrix
    """

    row1 = [-koff, -koff]
    row2 = [-kc, -kc - krev]

    matrix = np.array([row1, row2])

    return matrix

def differential_matrix_dissociation_conformational_selection(koff,kc,krev):

    """
    Dissociation step
    set of differential equations for
    (d[E1]*R)/dt and d[E2S]*R/dt.

    The signal is only proportional to E2S

    Args:
        koff (float): rate constant for E2S -> E2 + S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1

    Returns:
        matrix (np.ndarray): differential matrix
    """

    row1 = [-kc -krev, -krev]
    row2 = [0, -koff]

    matrix = np.array([row1, row2])

    return matrix

def constant_vector_induced_fit(koff, a, kc):

    """
    Calculate the constant vector for the induced fit model
    Args:
        koff (float): rate constant for E·S -> E + S
        a (float): concentration of the analyte (molecule being flown)
        kc (float): rate constant for E·S -> ES, induced fit step
    Returns:
        b_col (np.ndarray): constant vector
    """

    return [koff * a, kc * a]

def constant_vector_conformational_selection(kon,smax,krev,a):

    """
    Calculate the constant vector for the conformational selection model
    Args:
        kon (float): rate constant for E2 + S -> E2S
        smax (float): signal proportional to the complex (E2S)
        krev (float): rate constant for E2 -> E1
        a (float): concentration of the analyte (molecule being flown)
    Returns:
        b_col (np.ndarray): constant vector
    """

    return   [krev*smax, kon*a*smax]

def solve_state_t(time, A, initial_conditions, steady_state_value):

    """
    solve x(t) = x* + exp(At)(x0 - x*)
    x(t) is the system state at time t
    x*   is the steady state solution
    x0   is the initial state

    Args:
        time (float): time
        A (np.ndarray): differential matrix
        initial_conditions (np.ndarray): initial conditions
        steady_state_value (np.ndarray): steady state values
    Returns:
        state (np.ndarray): state at time t
    """

    v = initial_conditions - steady_state_value
    eAtv = expm(A * time).dot(v)
    state = steady_state_value + eAtv

    return state

def solve_steady_state(A, b):
    """
    Calculate the steady state solution of the system of equations
    Args:
        A (np.ndarray): differential matrix
        b (np.ndarray): constant vector
    Returns:
        steady_state (np.ndarray): steady state solution
    """
    return -solve(A, b)

def solve_induced_fit_association(time, a_conc, kon, koff, kc, krev,sP1L=0,sP2L=0,smax=0):

    """
    Obtain the signal for the induced fit model (surface-based)

    We assume that the signal is given by the complex E·S and ES is the same
    In other words, the complex before and after the conformational change produces the same signal

    Args:

        t (float): time
        a_conc (float): concentration of the analyte (molecule being flown)
        kon (float): rate constant for E + S -> E·S
        koff (float): rate constant for E·S -> E + S
        kc (float): rate constant for E·S -> ES
        krev (float): rate constant for ES -> E·S
        t0 (float): initial time
        smax (float): signal proportional to the complex (E·S or ES)

    Returns:

       signal (np.ndarray): signal over time

    """
    time = time - np.min(time) # Start at zero

    b_col = constant_vector_induced_fit(koff, smax, kc)

    # Initial conditions for St - R(E·S) - R(ES), and R(ES)
    # where St is the max signal, R(E·S) is the signal produced by E·S, and R(ES) is the signal produced by ES

    initial_conditions = np.array([smax-sP1L-sP2L, sP2L])

    m = differential_matrix_association_induced_fit(koff, kon, kc, krev, a_conc)

    res_steady_state = solve_steady_state(m, b_col)

    states = np.array([solve_state_t(t, m, initial_conditions, res_steady_state) for t in time])

    df = pd.DataFrame({
        'signal': smax - states[:, 0],
        'sP1L': smax - states[:, 0] - states[:, 1],
        'sP2L': states[:, 1]
    })

    return df

def solve_conformational_selection_association(time, a_conc, kon, koff, kc, krev,smax=0,sP1=0,sP2L=0):

    """
    Obtain the signal for the conformational selection model (surface-based)

    We assume that the signal is given by the complex E2S only

    Args:

        time (float): time
        a_conc (float): concentration of the analyte (molecule being flown)
        kon (float): rate constant for E2 + S -> E2S
        koff (float): rate constant for E2S -> E2 + S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1
        smax (float): signal proportional to the complex (E2S)

    Returns:

       signal (np.ndarray): signal over time

    """

    time = time - np.min(time) # Start at zero

    b_col   = constant_vector_conformational_selection(kon, smax, krev, a_conc)

    if sP2L == 0:

        ka_conf            = kc / krev
        initial_conditions = np.array([smax / (ka_conf + 1), 0]) # Assume pre-equilibrium of E1 and E2

    else:

        initial_conditions = np.array([sP1, sP2L])

    m = differential_matrix_association_conformational_selection(koff, kon, kc, krev, a_conc)

    res_steady_state = solve_steady_state(m, b_col)

    states = np.array([solve_state_t(t, m, initial_conditions, res_steady_state) for t in time])

    df = pd.DataFrame({
        'signal': states[:, 1], # Signal produced by E2S / P2L, depending on the notation
        'sP1': states[:, 0],
        'sP2': smax - states[:, 1] - states[:, 0]
    })

    return df

def solve_induced_fit_dissociation(time, koff, kc, krev,s0=0,sP2L=0,smax=0):

    """
    Obtain the dissociation signal for the induced fit model (surface-based)

    We assume that the signal given by the complex E·S and ES is the same
    In other words, the complex before and after the conformational change produces the same signal

    Args:

        t (float): time
        koff (float): rate constant for E·S -> E + S
        kc (float): rate constant for E·S -> ES
        krev (float): rate constant for ES -> E·S
        t0 (float): initial time
        s0 (float): initial signal
        sP2L (float): signal proportional to the complex ES (NOT the sum of E·S and ES)
        smax (float): signal proportional to the complex (E·S or ES)

    Returns:

       signal (np.ndarray): signal over time

    """
    time = time - np.min(time) # Start at zero

    b_col = constant_vector_induced_fit(koff, smax, kc)
    m     = differential_matrix_dissociation_induced_fit(koff, kc, krev)

    res_steady_state = solve_steady_state(m, b_col)

    initial_conditions = np.array([smax-s0, sP2L])

    states = np.array([solve_state_t(t, m, initial_conditions, res_steady_state) for t in time])

    df = pd.DataFrame({
        'signal': smax - states[:, 0],
        'sP1L': smax - states[:, 0] - states[:, 1],
        'sP2L': states[:, 1]
    })

    return df

def solve_conformational_selection_dissociation(time, koff, kc, krev,smax=0,sP1=0,sP2L=0):

    """
    Obtain the signal for the conformational selection model (surface-based)

    We assume that the signal is given by the complex E2S only

    Args:

        time (float): time
        koff (float): rate constant for E2S -> E2 + S
        kc (float): rate constant for E1 -> E2
        krev (float): rate constant for E2 -> E1
        smax (float): signal proportional to the complex (E2S)

    Returns:

       signal (np.ndarray): signal over time

    """

    time = time - np.min(time) # Start at zero

    b_col   = constant_vector_conformational_selection(0, smax, krev, 0)

    initial_conditions = np.array([sP1, sP2L])

    m = differential_matrix_dissociation_conformational_selection(koff, kc, krev)

    res_steady_state = solve_steady_state(m, b_col)

    states = np.array([solve_state_t(t, m, initial_conditions, res_steady_state) for t in time])

    df = pd.DataFrame({
        'signal': states[:, 1],
        'sP1': states[:, 0],
        'sP2': smax - states[:, 1] - states[:, 0]
    })

    return df

def fit_single_exponential(y,t):

    """

    Fit a single exponential to a signal

    Args:

        y (np.ndarray): signal
        t (np.ndarray): time

    Returns:

        fit_params (np.ndarray): fitted parameters

    """
    # Convert to numpy array, if needed,
    if not isinstance(y, np.ndarray):
        y = np.array(y)
    if not isinstance(t, np.ndarray):
        t = np.array(t)

    p0 = [0,np.max(y)*5, 0.1]

    t = t - np.min(t)  # Start at zero

    def fit_fx(t,a0,a1,kobs):

        return a0+a1*np.exp(-kobs*t)

    fit_params, cov = curve_fit(fit_fx, t, y,p0=p0)

    return fit_params

def fit_double_exponential(y,t):

    """

    Fit a double exponential to a signal

    Args:

        y (np.ndarray): signal
        t (np.ndarray): time

    Returns:

        fit_params (np.ndarray): fitted parameters

    """

    # Estimate a0 as the long-time limit
    a0 = np.mean(y[-10:])  # Take the last few points

    # Estimate kobs1 and kobs2 from log-transformed data
    log_y = np.log(np.abs(y - a0))
    slope1, _ = np.polyfit(t[:10], log_y[:10], 1)  # First few points
    slope2, _ = np.polyfit(t[-10:], log_y[-10:], 1)  # Last few points
    kobs1, kobs2 = -slope1, -slope2

    # Estimate a1 and a2
    a1_a2_sum = y[0] - a0
    a1 = a1_a2_sum / 2  # Roughly split contributions
    a2 = a1_a2_sum / 2

    # Initial parameter estimates
    p0 = [a0, a1, a2, kobs1, kobs2]

    def fit_fx(t,a0,a1,a2,kobs,kobs2):

        return a0 + a1*np.exp(-kobs*t) + a2*np.exp(-kobs2*t)

    fit_params, cov = curve_fit(fit_fx, t, y,p0=p0)

    return fit_params

##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them
##### CAUTION: The functions below still need to be tested and validated - do not use them

def assoc_one_site_asymmetric_ci95(kd_estimated, rss_desired,
                                   assoc_list, time_assoc_lst, lig_conc_lst,
                                   p0, low_bounds, high_bounds,
                                   smax_idx, shared_smax):
    '''

    CAUTION: This function still need to be tested and validated - do not use it
    '''

    def f_to_optimize(Kd):
        fit_params, _, fit_vals = fit_one_site_association(assoc_list, time_assoc_lst, lig_conc_lst,
                                                           p0, low_bounds, high_bounds,
                                                           smax_idx=smax_idx,
                                                           shared_smax=shared_smax,
                                                           fixed_Kd=True, Kd_value=Kd
                                                           )

        rss = get_rss(concat_signal_lst(assoc_list), concat_signal_lst(fit_vals))

        return np.abs(rss - rss_desired)

    boundsMin = np.array([kd_estimated / 10, kd_estimated])
    boundsMax = np.array([kd_estimated, kd_estimated * 100])

    kd_min95 = minimize_scalar(f_to_optimize, bounds=boundsMin, method='bounded')
    kd_max95 = minimize_scalar(f_to_optimize, bounds=boundsMax, method='bounded')

    kd_min95, kd_max95 = kd_min95.x, kd_max95.x

    return np.array([kd_min95, kd_max95])

def ode_bivalent_analyte_association(t,y,A, k_off1, k_off2, Kd1, Kd2):

    """
    Caution - THIS FUNCTION STILL NEEDS TO BE TESTED AND VALIDATED!

    This model assumes that because of limited distance between two adjacent binding sites on the surface,
    the bivalent analyte can form a bridged complex.
    This interaction is linked, meaning that the formation of AB2 complex cannot occur before the formation of AB,
    and AB cannot dissociate before the dissociation of AB2. This avidity effect results in a slower apparent dissociation rate
    than would be expected if the interaction followed a 1:1 binding profile.

    Reference: https://www.sartorius.com/resource/blob/742330/05671fe2de45d16bd72b8078ac28980d/octet-biomolecular-binding-kinetics-application-note-4014-en-1--data.pdf

    Args:

        t (np.ndarray): Time variable
        y (list): List of the state variables [B, AB, AB2]
        A (float): Concentration of the analyte
        k_off1 (float): Dissociation rate constant for the first binding event
        k_off2 (float): Dissociation rate constant for the second binding event
        Kd1 (float): Equilibrium dissociation constant for the first binding event
        Kd2 (float): Equilibrium dissociation constant for the second binding event

    Returns:

        dB_dt (float): signal over time
        dAB_dt (float): signal over time
        dAB2_dt (float): signal over time

    """

    B, AB, AB2 = y

    k_on1 = k_off1 / Kd1
    k_on2 = k_off2 / Kd2

    dB_dt = - (2 * k_on1 * A * B - k_off1 * AB) - (k_on2 * AB * B - 2 * k_off2 * AB2)
    dAB_dt = (2 * k_on1 * A * B - k_off1 * AB) - (k_on2 * AB * B - 2 * k_off2 * AB2)
    dAB2_dt = (k_on2 * AB * B - 2 * k_off2 * AB2)

    return [dB_dt, dAB_dt, dAB2_dt]

def solve_ode_bivalent_analyte_association(t,B0,AB0,AB20,A,k_off1,k_off2,Kd1,Kd2,t0=0):

    """
    Caution - THIS FUNCTION STILL NEEDS TO BE TESTED AND VALIDATED!

    Solves the ODE for the bivalent analyte association model
    """

    t = t + t0

    out = solve_ivp(ode_bivalent_analyte_association,t_span=[np.min(t), np.max(t)],
                    t_eval=t,y0=[B0,AB0,AB20],args=(A,k_off1,k_off2,Kd1,Kd2),method="LSODA")

    return out.y


def fit_two_sites_assoc_and_disso(assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                  disso_signal_lst, disso_time_lst,
                                  initial_parameters,low_bounds, high_bounds,
                                  smax_idx    = None,
                                  shared_smax = False,
                                  fixed_t0 = True,
                                  fixed_Kd1 = False,   Kd1_value = None,
                                  fixed_koff1 = False, koff1_value = None,
                                  fixed_Kd2 = False,   Kd2_value = None,
                                  fixed_koff2 = False, koff2_value = None,
                                  same_affinity=False):

    """
    CAUTION - THIS FUNCTION IS STILL UNDER DEVELOPMENT! IT NEEDS TO BE TESTED AND VALIDATED!
    Fits a two-site binding model to a set of association and dissociation signals

    The order of the parameters is:

        For shared S-max, two datasets with different loading levels:

        global_Kd1, global_Koff1, global_Kd2, global_Koff2,
        subglobal_T01, subglobal_T02,
        subglobal_S1max1, subglobal_S1max2, (Smax of the ligand #1)
        subglobal_S2max1,subglobal_S2max2   (Smax of the ligand #2)

        For non-shared S-max, two datasets with different loading levels and three traces each:

        global_Kd1, global_Koff1, global_Kd2, global_Koff2,
        subglobal_T01, subglobal_T02,
        subglobal_S1maxTrace1, subglobal_S1maxTrace2,subglobal_S1maxTrace3, (Smax of the ligand #1)
        subglobal_S1maxTrace4, subglobal_S1maxTrace5,subglobal_S1maxTrace6, (Smax of the ligand #1)
        subglobal_S2maxTrace1, subglobal_S2maxTrace2,subglobal_S2maxTrace3, (Smax of the ligand #2)
        subglobal_S2maxTrace4, subglobal_S2maxTrace5,subglobal_S2maxTrace6, (Smax of the ligand #2)

    """

    all_signal_assoc = concat_signal_lst(assoc_signal_lst)
    all_signal_disso = concat_signal_lst(disso_signal_lst)

    time_lst_assoc = [np.array(t) for t in assoc_time_lst]
    time_lst_disso = [np.array(t) for t in disso_time_lst]

    time_lst_disso = [t - t[0] for t in time_lst_disso]

    # Remove the first parameters if we have the same affinity
    if same_affinity:

        initial_parameters = initial_parameters[2:]
        low_bounds         = low_bounds[2:]
        high_bounds        = high_bounds[2:]

        start = 2 - sum([fixed_Kd1, fixed_koff1])

    else:

        start = 4 - sum([fixed_Kd1, fixed_koff1,fixed_Kd2, fixed_koff2])

    n_unq_smax = len(np.unique(smax_idx))

    n_t0s = n_unq_smax*(not fixed_t0)

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
            Kd, Koff, T0, Smax1, Smax2, Smax3, ...
        """

        Kd1   = Kd1_value   if fixed_Kd1   else args[0]
        Koff1 = koff1_value if fixed_koff1 else args[1 - sum([fixed_Kd1])]

        if same_affinity:

            Kd2   = Kd1
            Koff2 = Koff1
        else:

            Kd2   = Kd2_value   if fixed_Kd2   else args[2 - sum([fixed_Kd1, fixed_koff1])]
            Koff2 = koff2_value if fixed_koff2 else args[3 - sum([fixed_Kd1, fixed_koff1, fixed_Kd2])]

        signal_a = []

        last_y1 = []
        last_y2 = []

        for i,t in enumerate(time_lst_assoc):

            analyte_conc = analyte_conc_lst[i]

            t0     = args[start + smax_idx[i]] if not fixed_t0 else 0

            s_max1  = args[start+smax_idx[i]+n_t0s              ]   if shared_smax  else args[start + i + n_t0s]
            s_max2  = args[start+smax_idx[i]+n_t0s + n_unq_smax ]   if shared_smax  else args[start + i + n_t0s + + n_unq_smax ]

            y1 = one_site_association_analytical(t,0,s_max1,Koff1,Kd1,analyte_conc,t0)
            y2 = one_site_association_analytical(t,0,s_max2,Koff2,Kd2,analyte_conc,t0)

            last_y1.append(y1[-1])
            last_y2.append(y2[-1])

            signal_a.append(y1+y2)

        signal_d = []

        for i,t in enumerate(time_lst_disso):

            s01        = last_y1[i]
            s02        = last_y2[i]

            y1 = one_site_dissociation_analytical(t,s01,Koff1)
            y2 = one_site_dissociation_analytical(t,s02,Koff2)

            signal_d.append(y1+y2)

        signal_assoc_fit  = np.concatenate(signal_a, axis=0)
        signal_dissoc_fit = np.concatenate(signal_d, axis=0)

        return np.concatenate([signal_assoc_fit,signal_dissoc_fit], axis=0)

    all_signal = np.concatenate([all_signal_assoc,all_signal_disso], axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    Kd1   = Kd1_value if fixed_Kd1 else global_fit_params[0]
    Koff1 = koff1_value if fixed_koff1 else global_fit_params[1 - sum([fixed_Kd1])]

    if same_affinity:
        Kd2   = Kd1
        Koff2 = Koff1
    else:
        # If we have different affinities, we need to get the Kd2 and Koff2
        # from the global_fit_params
        Kd2   = Kd2_value if fixed_Kd2 else global_fit_params[2 - sum([fixed_Kd1, fixed_koff1])]
        Koff2 = koff2_value if fixed_koff2 else global_fit_params[3 - sum([fixed_Kd1, fixed_koff1, fixed_Kd2])]

    if shared_smax:

        unique_ids = np.unique(smax_idx)
        n_unq      = len(unique_ids)

        s1_max_unq = global_fit_params[(start+n_t0s):(start+n_t0s+n_unq)]
        s2_max_unq = global_fit_params[(start+n_t0s+n_unq):]

        s1_max_all = expand_parameter_list(s1_max_unq,smax_idx)
        s2_max_all = expand_parameter_list(s2_max_unq,smax_idx)

    else:

        s1_max_all     = global_fit_params[start+n_t0s:start+n_t0s+n_unq_smax]
        s2_max_all     = global_fit_params[start+n_t0s+n_unq_smax:]

    if not fixed_t0:

        t0_unq = global_fit_params[(start):(start+n_t0s)]
        t0_all = expand_parameter_list(t0_unq,smax_idx)

    else:

        t0_all = np.repeat(0,len(time_lst_assoc))

    fitted_values_assoc = []

    last_y1 = []
    last_y2 = []

    for t,s_max1,s_max2,analyte_conc,t0 in zip(time_lst_assoc,s1_max_all,s2_max_all,analyte_conc_lst,t0_all):

        y1 = one_site_association_analytical(t,0,s_max1,Koff1,Kd1,analyte_conc,t0)
        y2 = one_site_association_analytical(t,0,s_max2,Koff2,Kd2,analyte_conc,t0)

        last_y1.append(y1[-1])
        last_y2.append(y2[-1])

        fitted_values_assoc.append(y1+y2)

    fitted_values_disso = []

    for i,t in enumerate(time_lst_disso):

        s01   = last_y1[i]
        s02   = last_y2[i]

        y1 = one_site_dissociation_analytical(t,s01,Koff1)
        y2 = one_site_dissociation_analytical(t,s02,Koff2)

        fitted_values_disso.append(y1+y2)

    return global_fit_params, cov, fitted_values_assoc, fitted_values_disso

def fit_induced_fit_sites_assoc_and_disso(assoc_signal_lst, assoc_time_lst, analyte_conc_lst,
                                  disso_signal_lst, disso_time_lst,
                                  initial_parameters,low_bounds, high_bounds,
                                  smax_idx    = None,
                                  shared_smax = False,
                                  fixed_t0 = True,
                                  fixed_kon1 = False,  kon1_value = None,
                                  fixed_koff1 = False, koff1_value = None,
                                  fixed_kon2 = False,  kon2_value = None,
                                  fixed_koff2 = False, koff2_value = None,
                                  same_affinity=False):

    """
    CAUTION - THIS FUNCTION IS STILL UNDER DEVELOPMENT! IT NEEDS TO BE TESTED AND VALIDATED!
    Fits a two-site binding model to a set of association and dissociation signals

    The order of the parameters is:

        For shared S-max, two datasets with different loading levels:

        global_Kd1, global_Koff1, global_Kd2, global_Koff2,
        subglobal_T01, subglobal_T02,
        subglobal_S1max1, subglobal_S1max2, (Smax of the ligand #1)
        subglobal_S2max1,subglobal_S2max2   (Smax of the ligand #2)

        For non-shared S-max, two datasets with different loading levels and three traces each:

        global_Kd1, global_Koff1, global_Kd2, global_Koff2,
        subglobal_T01, subglobal_T02,
        subglobal_S1maxTrace1, subglobal_S1maxTrace2,subglobal_S1maxTrace3, (Smax of the ligand #1)
        subglobal_S1maxTrace4, subglobal_S1maxTrace5,subglobal_S1maxTrace6, (Smax of the ligand #1)
        subglobal_S2maxTrace1, subglobal_S2maxTrace2,subglobal_S2maxTrace3, (Smax of the ligand #2)
        subglobal_S2maxTrace4, subglobal_S2maxTrace5,subglobal_S2maxTrace6, (Smax of the ligand #2)

    """

    all_signal_assoc = concat_signal_lst(assoc_signal_lst)
    all_signal_disso = concat_signal_lst(disso_signal_lst)

    time_lst_assoc = [np.array(t) for t in assoc_time_lst]
    time_lst_disso = [np.array(t) for t in disso_time_lst]

    time_lst_disso = [t - t[0] for t in time_lst_disso]

    start = 4 - sum([fixed_kon1, fixed_koff1,fixed_kon2, fixed_koff2])

    n_unq_smax = len(np.unique(smax_idx))

    n_t0s = n_unq_smax*(not fixed_t0)

    def fit_fx(dummyVariable, *args):

        """
        Arguments order:
            Kd, Koff, T0, Smax1, Smax2, Smax3, ...
        """

        kon1  = kon1_value  if fixed_kon1  else args[0]
        koff1 = koff1_value if fixed_koff1 else args[1 - sum([fixed_kon1])]
        kon2  = kon2_value  if fixed_kon2  else args[2 - sum([fixed_kon1, fixed_koff1])]
        koff2 = koff2_value if fixed_koff2 else args[3 - sum([fixed_kon1, fixed_koff1, fixed_kon2])]

        signal_a = []

        last_y1 = []

        for i,t in enumerate(time_lst_assoc):

            analyte_conc = analyte_conc_lst[i]

            #t0     = args[start + smax_idx[i]] if not fixed_t0 else 0

            s_max  = args[start+smax_idx[i]+n_t0s              ]   if shared_smax  else args[start + i + n_t0s]

            df = solve_induced_fit_association(t,analyte_conc,s_max,kon1,koff1,kon2,koff2)
            last_row = df.iloc[-1]
            last_y1.append(last_row)

            signal_a.append(df[['signal']].to_numpy())

        signal_d = []

        for i,t in enumerate(time_lst_disso):

            last_row   = last_y1[i]
            sP2L       = last_row[['sP2L']].to_numpy()[0]
            s0         = last_row[['signal']].to_numpy()[0]

            df = solve_induced_fit_dissociation(t, koff1, kon2, koff2,s0=s0,sP2L=sP2L,smax=s_max)

            signal_d.append(df[['signal']].to_numpy())

        signal_assoc_fit  = np.concatenate(signal_a, axis=0)
        signal_dissoc_fit = np.concatenate(signal_d, axis=0)

        return np.concatenate([signal_assoc_fit,signal_dissoc_fit], axis=0)

    all_signal = np.concatenate([all_signal_assoc,all_signal_disso], axis=0)

    global_fit_params, cov = curve_fit(fit_fx, 1, all_signal,
                                       p0=initial_parameters,
                                       bounds=(low_bounds, high_bounds))

    kon1   = kon1_value if fixed_kon1 else global_fit_params[0]
    koff1  = koff1_value if fixed_koff1 else global_fit_params[1 - sum([fixed_kon1])]
    kon2   = kon2_value if fixed_kon2 else global_fit_params[2 - sum([fixed_kon1, fixed_koff1])]
    koff2  = koff2_value if fixed_koff2 else global_fit_params[3 - sum([fixed_kon1, fixed_koff1, fixed_kon2])]

    if shared_smax:

        unique_ids = np.unique(smax_idx)
        n_unq      = len(unique_ids)

        s_max_unique = global_fit_params[start:(start+n_unq)]

        s_max_all = expand_parameter_list(s_max_unique,smax_idx)

    else:

        s_max_all     = global_fit_params[start:]

    if not fixed_t0:

        t0_unq = global_fit_params[(start):(start+n_t0s)]
        t0_all = expand_parameter_list(t0_unq,smax_idx)

    else:

        t0_all = np.repeat(0,len(time_lst_assoc))

    fitted_values_assoc = []

    last_y1 = []

    for t,s_max,analyte_conc,t0 in zip(time_lst_assoc,s_max_all,analyte_conc_lst,t0_all):

        df = solve_induced_fit_association(t,analyte_conc,s_max,kon1,koff1,kon2,koff2)

        last_row = df.iloc[-1]
        last_y1.append(last_row)

        fitted_values_assoc.append(df[['signal']].to_numpy())

    fitted_values_disso = []

    for i,t in enumerate(time_lst_disso):

        sP2L = last_y1[i][['sP2L']].to_numpy()[0]
        s0   = last_y1[i][['signal']].to_numpy()[0]
        s_max = s_max_all[i]

        df = solve_induced_fit_dissociation(t, koff1, kon2, koff2,s0=s0,sP2L=sP2L,smax=s_max)

        fitted_values_disso.append(df[['signal']].to_numpy())

    return global_fit_params, cov, fitted_values_assoc, fitted_values_disso
