import numpy as np
import pandas as pd
import os
import re

from helpers2        import *
from fitting_helpers_in_solution import *

"""
Classes for the analysis of kinetics (solution-based) data
Code written by Osvaldo Burastero

No warranty whatsoever
If you have questions please contact me:    oburastero@gmail.com
"""

class KinGenie_csv_in_solution:

    """
    A class used to represent a KinGenie csv file, which can be exported from the Simulation panel

    For this kind of experiments, we only expect one step per trace
    In other words, we only have the interaction phase (e.g., no baseline or dissociation)

    Attributes

        name (str):                    name of the experiment
        fn (str):                      file name
        xs (np.array):                 list of x values (time, length n, one per sensor)
        ys (np.array):                 list of y values (length n, one per sensor)
        no_traces (int):               number of traces
        traces_names (list):           list of traces names (length n, one per trace)
        traces_names_unique (list):    list of unique traces (length n, one per sensor)
        conc_df (pd.DataFrame):        dataframe with the ligand and protein concentration information
    """

    def __init__(self,name):

        self.name                = name
        self.xs                  = None
        self.ys                  = None
        self.no_sensors          = None
        self.traces_names        = None
        self.traces_names_unique = None
        self.traces_loaded       = False
        self.type                = 'kingenie_csv_in_solution'
        self.conc_df             = None

    def create_unique_traces_names(self):

        self.traces_names_unique = [self.name + ' ' + trace_name for trace_name in self.traces_names]

        return None

    def read_csv(self,file):

        """
        Read the KinGenie csv file

        Example:

            Time	Signal	Protein_concentration_micromolar Ligand_concentration_micromolar
            0	    0	    5	                                0.1
            0.5	    1	    5	                                0.1

        Results:

            It creates the attributes

                self.xs
                self.ys
                self.no_sensors
                self.sensor_names
                self.sensor_names_unique

        """

        df = pd.read_csv(file)

        # Find the ligand concentrations
        ligand_conc = df['Ligand_concentration_micromolar']

        # Find the protein concentrations
        protein_conc = df['Protein_concentration_micromolar']

        # Combine protein and ligand concentration into one array
        combined_concs_array = [str(x) + ' and ' + str(y) for x,y in zip(protein_conc,ligand_conc)]

        # Add a new column to the dataframe
        df['Combined_concentration'] = combined_concs_array

        # Find the unique combined concentrations
        combined_concs_unq = np.unique(combined_concs_array)

        self.no_traces   = len(combined_concs_unq)

        # Use fake names
        self.traces_names = ['sim. trace ' + str(i+1) for i in range(self.no_traces)]

        # Initiate self.xs and self.ys
        self.xs = []
        self.ys = []

        protein_conc_unqs = []
        ligand_conc_unqs  = []

        # Now populate, for each sensor self.xs and self.ys
        for i,cc in enumerate(combined_concs_unq):

            # Extract the rows with the same combined concentrations
            df_temp = df[df['Combined_concentration'] == cc]

            protein_conc_unqs.append(df_temp['Protein_concentration_micromolar'].to_numpy()[0])
            ligand_conc_unqs.append(df_temp['Ligand_concentration_micromolar'].to_numpy()[0])

            # Extract the time values of the association/dissociation phase
            time_int = df_temp['Time'].to_numpy()

            # Populate self.xs, we include time signal inside a list to make it compatible with the surface-based code
            self.xs += [time_int]

            # Extract the signal values
            signal = df_temp['Signal'].to_numpy()

            # Populate self.ys, we include time signal inside a list to make it compatible with the surface-based code
            self.ys += [signal]

        # Now generate a fake ligand concentration df
        df_traces = pd.DataFrame({
            'Trace':self.traces_names,
            'Protein_concentration_micromolar':protein_conc_unqs,
            'Ligand_concentration_micromolar':ligand_conc_unqs,
            'SampleID':'simulation',
            'Experiment':self.name})

        self.conc_df = df_traces

        self.create_unique_traces_names()
        self.traces_loaded = True

        return None

    def get_step_xy(self,trace_name,type='y'):

        """
        Return the x or y values of a certain step

        Args:

            trace_name (str): name of the trace
            type (str):        x or y

        Returns:

            x or y (np.n) values of the step

        """

        trace_id = self.traces_names.index(trace_name)

        if type == 'x':

            return self.xs[trace_id][0]

        else:

            return self.ys[trace_id][0]

class Kinetics_fitter_in_solution:

    """

    A class used to fit kinetics data with shared thermodynamic parameters

    Attributes
    ----------

        names (list):                   list of experiment names
        assoc (list):                   list of interaction signals, one per replicate,
                                        each signal is a numpy matrix of size n*m where n is the number of time points and m
                                        is the number of protein/ligand concentrations combinations
        lig_conc (list):                list of ligand concentrations, one per replicate
                                        each element contains a numpy array with the ligand concentrations
        protein_conc (list):            list of protein concentrations, one per replicate
                                        each element contains a numpy array with the protein concentrations
        time_assoc (list):              list of time points for the association signals, one per replicate
                                        each element contains a numpy array with the time points
        signal_ss (list):               list of steady state signals, one per replicate
                                        each element contains a numpy array with the steady state signals
        signal_ss_fit (list):           list of steady state fitted signals, one per replicate
                                        each element contains a numpy array with the steady state fitted signals
        signal_assoc_fit (list):        list of association kinetics fitted signals, one per replicate
                                        each element contains a numpy array with the association kinetics fitted signals
        fit_params_kinetics (pd.Dataframe):     dataframe with the fitted parameters for the association / dissociation kinetics
        fit_params_ss (pd.Dataframe):           dataframe with the values of the fitted parameters - steady state

    """

    def __init__(self,time_assoc_lst,association_signal_lst,lig_conc_lst,
                 protein_conc_lst=None,name_lst=None):

        self.names            = name_lst                 # String, experiment name
        self.assoc            = association_signal_lst   # Signal to fit,  numpy matrix of size n*m
        self.lig_conc         = lig_conc_lst             # Ligand concentration
        self.protein_conc     = protein_conc_lst         # Protein concentration
        self.time_assoc       = time_assoc_lst           # Time,   numpy matrix of size n*m

        self.signal_ss        = None  # Steady state signal
        self.signal_ss_fit    = None  # Steady state fitted signal
        self.signal_assoc_fit = None  # Association kinetics fitted signal
        self.signal_disso_fit = None  # Association kinetics fitted signal

    def set_experiments_id(self):

        """
        Create a 1D list that contains the ID corresponding to each protein/ligand concentration

        For example, if we have two experiments, the first one with 7 ligand concentrations
        and the second one with 2 ligand concentrations the ID should be:

        id = [0,0,0,0,0,0,0,1,1]

        """

        experiments_id = []

        for i,x in enumerate(self.lig_conc):

            experiments_id += [i for _ in x]

        self.experiments_id = experiments_id

        return None

    def get_steady_state(self):

        self.signal_ss = [np.median(assoc[-10:, :],axis=0) for assoc in self.assoc]

        return None

    def clear_fittings(self):

        self.signal_ss_fit        = None  # Steady state fitted signal
        self.signal_assoc_fit     = None  # Association kinetics fitted signal
        self.signal_disso_fit     = None  # Association kinetics fitted signal
        self.fit_params_kinetics  = None  # Fitted parameters for the association / dissociation kinetics
        self.fit_params_ss        = None  # Values of the fitted parameters - steady state

        return None

    def kf_fit_steady_state(self):

        '''
        This function fits the steady state signal
        '''

        self.clear_fittings()

        if self.signal_ss is None:

            self.get_steady_state()

        Kd_init = np.median(self.lig_conc[0])
        p0      = [Kd_init] + [np.max(signal) for signal in self.signal_ss]

        kd_min  = np.min(self.lig_conc[0]) / 1e2
        kd_max  = np.max(self.lig_conc[0]) * 1e2

        low_bounds  = [kd_min]  + [x*0.5 for x in p0[1:]]
        high_bounds = [kd_max]  + [x*4   for x in p0[1:]]

        fit, cov, fit_vals = fit_steady_state_one_site(self.signal_ss,self.lig_conc,p0,low_bounds,high_bounds)

        self.Kd_ss   = fit[0]
        Rmax         = fit[1:]

        self.signal_ss_fit = fit_vals

        n = np.sum([len(signal) for signal in self.signal_ss])
        p = len(p0)

        rss_desired = get_desired_rss(concat_signal_lst(self.signal_ss), concat_signal_lst(fit_vals), n, p)

        minKd, maxKd = steady_state_one_site_asymmetric_ci95(self.Kd_ss, self.signal_ss, self.lig_conc, p0[1:],
                                                             low_bounds[1:], high_bounds[1:], rss_desired)

        # Create a dataframe with the fitted parameters
        df_fit = pd.DataFrame({'Kd [µM]': self.Kd_ss,
                               'Kd_min95': minKd,
                               'Kd_max95': maxKd,
                               'Rmax': Rmax, 'Name': self.names})

        self.fit_params_ss = df_fit

        # Fit the steady state signal
        return None


    def kf_fit_one_site_association(self,shared_smax=True):

        self.clear_fittings()

        # Try to fit first the dissociation curves to get a better estimate of Koff
        try:

            self.kf_fit_one_site_dissociation()

            p0          = [self.Kd_ss,self.Koff_disso]
            low_bounds  = [self.Kd_ss/5e2,self.Koff_disso/5e2]
            high_bounds = [self.Kd_ss*5e2,self.Koff_disso*5e2]

            self.clear_fittings()

        except:

            p0 = [self.Kd_ss,0.01]

            low_bounds  = [self.Kd_ss/5e2,1e-7]
            high_bounds = [self.Kd_ss*5e2,10]

        smax_guesses = [np.max(assoc) for assoc in self.assoc]

        if shared_smax:

            i = 0
            for smax_guess,assoc in zip(smax_guesses,self.assoc):

                p0.append(smax_guess)
                low_bounds.append(smax_guess/10)
                high_bounds.append(smax_guess*10)

                i         += 1

        else:

            for smax_guess,assoc in zip(smax_guesses,self.assoc):

                p0       += np.repeat(smax_guess,assoc.shape[1]).tolist()

            low_bounds  += [x/10 for x in p0[2:]]
            high_bounds += [x*10 for x in p0[2:]]

        fit, cov, fit_vals = fit_one_site_association(
            self.assoc_list,self.time_assoc_lst,self.lig_conc_lst,
            p0,low_bounds,high_bounds,smax_idx=self.experiments_id,shared_smax=shared_smax
        )

        self.signal_assoc_fit = fit_vals

        self.Kd_asso   = fit[0]
        self.Koff_asso = fit[1]
        self.Smax_asso = fit[2:]

        # Create a dataframe with the fitted parameters
        df_fit = pd.DataFrame({'Kd [µM]':   self.Kd_asso,
                               'k_off [1/s]': self.Koff_asso,
                               'Smax': self.Smax_asso})

        error     = np.sqrt(np.diag(cov))
        rel_error = error/fit * 100

        df_error = pd.DataFrame({'Kd [µM]':   rel_error[0],
                                 'k_off': rel_error[1],
                                 'Smax': rel_error[2:]})

        self.fit_params_kinetics       = df_fit
        self.fit_params_kinetics_error = df_error

        return  None