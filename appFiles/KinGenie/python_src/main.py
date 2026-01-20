import base64
import xml.etree.ElementTree as ET
import json
import numpy as np
import pandas as pd
import os
import re
from collections import Counter
import copy

from helpers2        import *
from fitting_helpers import *
from fitting_helpers_in_solution import *
from main_insolution import *
from scipy.stats import trim_mean

"""
Classes for the analysis of kinetics (biolayer interferometry) data
Code mainly written by Osvaldo Burastero
Functions 'read_sensor_data' and 'convert_to_numbers' were adapted from code provided by Stephan Niebling

No warranty whatsoever
If you have questions please contact me:    oburastero@gmail.com
"""

factor_conc_to_micro = {'nM':1e-3, 'µM':1, 'mM':1e3, 'M':1e6, 'mg/ml':1e3, 'µg/ml':1}

class surface_based_experiment:

    def __init__(self):
        """

        Initialize the instance

        """

        self.traces_loaded = False
        self.sample_plate_loaded = False
        self.xs = None
        self.ys = None
        self.sensor_names = None
        self.ligand_conc_df = None
        self.df_steps = None
        self.steps_performed = None

    def create_unique_sensor_names(self):

        self.sensor_names_unique = [self.name + ' ' + sensor_name for sensor_name in self.sensor_names]

        return None

    def subtraction_one_to_one(self, sensor_name1, sensor_name2,inplace=True):

        """

        Subtract the signal of sensor2 from sensor1

        Args:

            sensor_name1 (str): name of the sensor to subtract from
            sensor_name2 (str): name of the sensor to subtract
            inplace (bool):     if True, the subtraction is done in place, otherwise a new sensor is created

        Results:

            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df

        """

        new_sensor_name = sensor_name1 + ' - ' + sensor_name2

        if new_sensor_name in self.sensor_names:

            new_sensor_name = new_sensor_name + ' rep'

        sensor1 = self.sensor_names.index(sensor_name1)
        sensor2 = self.sensor_names.index(sensor_name2)

        if not self.traces_loaded:
            print("No traces loaded")
            return None
        # Check if sensors are compatible
        if len(self.xs[sensor1]) != len(self.xs[sensor2]):
            print("Sensors have different number of steps")
            return None

        if len(self.xs[sensor1][0]) != len(self.xs[sensor2][0]):
            print("Sensors have different number of points")
            return None

        # Subtract
        if inplace:

            for i in range(len(self.xs[sensor1])):
                self.ys[sensor1][i] -= self.ys[sensor2][i]
                self.sensor_names[sensor1] = new_sensor_name

            self.ligand_conc_df['Sensor'] = self.ligand_conc_df['Sensor'].replace(sensor_name1,new_sensor_name)

        else:
            ys = []
            for i in range(len(self.xs[sensor1])):
                ys.append(self.ys[sensor1][i] - self.ys[sensor2][i])

            # Fill instance
            self.xs.append(self.xs[sensor1])
            self.ys.append(ys)
            self.sensor_names.append(new_sensor_name)

            # Add new sensor name to the ligand conc df
            previous_row        = self.ligand_conc_df[self.ligand_conc_df['Sensor'] == sensor_name1]
            new_row             = previous_row.copy()
            new_row['Sensor']   = new_sensor_name
            new_row['SampleID'] = new_row['SampleID'] + ' bl subtracted'

            self.ligand_conc_df = pd.concat([self.ligand_conc_df,new_row])

        self.create_unique_sensor_names()

        return None

    def subtraction(self,list_of_sensor_names,reference_sensor,inplace=True):

        """
        Apply the subtract operation to a list of sensors

        Args:
            sensor_name1 (str): name of the sensor to subtract from
            list_of_sensor_names (list): list of sensor names to subtract
            inplace (bool):     if True, the subtraction is done in place, otherwise a new sensor is created
        Results:
            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df
        """

        if not isinstance(list_of_sensor_names, list):
            list_of_sensor_names = [list_of_sensor_names]

        for sensor_name in list_of_sensor_names:
            self.subtraction_one_to_one(sensor_name, reference_sensor, inplace=inplace)

        return  None

    def average(self,list_of_sensor_names,new_sensor_name='Average'):

        """

        Average the signals of the sensors in the list

        Args:

            list_of_sensor_names (list): list of sensor names to average
            new_sensor_name (str):       name of the new sensor

        Results:

            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df

        """

        # Check if sensors are loaded
        if not self.traces_loaded:
            print("No traces loaded")
            return None

        if new_sensor_name in self.sensor_names:

            new_sensor_name = new_sensor_name + ' rep'

        ys = []

        num_sensors = len(list_of_sensor_names)
        sensor1 = self.sensor_names.index(list_of_sensor_names[0])

        for i in range(len(self.xs[sensor1])):

            sensors = [self.sensor_names.index(sensor_name) for sensor_name in list_of_sensor_names]

            sum_ys = sum(self.ys[sensor][i] for sensor in sensors)
            ys.append(sum_ys / num_sensors)

        # Fill instance
        self.xs.append(self.xs[sensor1])
        self.ys.append(ys)
        self.sensor_names.append(new_sensor_name)

        # Add new sensor name to the ligand conc df
        previous_row        = self.ligand_conc_df[self.ligand_conc_df['Sensor'] == list_of_sensor_names[0]]
        new_row             = previous_row.copy()
        new_row['Sensor']   = new_sensor_name
        new_row['SampleID'] = new_row['SampleID'] + ' averaged'

        self.ligand_conc_df = pd.concat([self.ligand_conc_df,new_row])

        self.create_unique_sensor_names()

        return None

    def align_association(self,sensor_names,inplace=True,new_names = False):

        """

        Align the BLI traces based on the signal before the association step(s)

        Args:

            sensor_names (str or list): name of the sensor(s) to align
            inplace (bool):             if True, the alignment is done in place, otherwise a new sensor is created


        Results:

            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df

        """

        if not isinstance(sensor_names, list):
            sensor_names = [sensor_names]

        # Find the index of the association steps
        association_steps_indices = self.df_steps.index[self.df_steps['Type'] == 'ASSOC'].to_numpy()

        # Find the dissociation steps indexes
        dissociation_steps_indices = self.df_steps.index[self.df_steps['Type'] == 'DISASSOC'].to_numpy()

        # Remove all association steps that come directly after a dissociation step (useful for single cycle kinetics)
        for idx in association_steps_indices:
            if idx-1 in dissociation_steps_indices:
                association_steps_indices = np.delete(association_steps_indices, np.where(association_steps_indices == idx)[0])

        sensor_indices = [self.sensor_names.index(sensor_name) for sensor_name in sensor_names]

        # Determine the usage of new sensor names
        use_new_names = not inplace or (inplace and new_names)

        for sensor in sensor_indices:

            # Start - Determine the new sensor name
            new_sensor_name = self.sensor_names[sensor] + ' aligned' if use_new_names else self.sensor_names[sensor]

            if new_sensor_name in self.sensor_names and not inplace:
                new_sensor_name += ' rep'
            # End of - Determine the new sensor name

            # Create a copy of the list
            ys = copy.deepcopy(self.ys[sensor])

            for i, association_step_index in enumerate(association_steps_indices):

                #  Subtract the first point of the previous baseline step
                last_point = np.mean(self.ys[sensor][association_step_index-1][-10:])

                if i == 0:

                    for step in range(association_step_index-1):

                        value = self.ys[sensor][step] - last_point

                        if inplace:

                            self.ys[sensor][step] = value

                        else:

                            ys[step]    = value

                for step in range(association_step_index-1,self.no_steps):

                    value = self.ys[sensor][step] - last_point

                    if inplace:

                        self.ys[sensor][step] = value

                    else:

                        ys[step] = value

            if inplace:

                #Replace in the ligand conc df the sensor name
                if use_new_names:

                    self.ligand_conc_df['Sensor'] = self.ligand_conc_df['Sensor'].replace(self.sensor_names[sensor],new_sensor_name)
                    self.sensor_names[sensor]     = new_sensor_name

            else:

                self.xs.append(self.xs[sensor])
                self.ys.append(ys)
                self.sensor_names.append(new_sensor_name)

                # Add the new sensor name to the ligand conc df
                previous_row        = self.ligand_conc_df[self.ligand_conc_df['Sensor'] == self.sensor_names[sensor]]
                new_row             = previous_row.copy()
                new_row['Sensor']   = new_sensor_name
                new_row['SampleID'] = new_row['SampleID'] + ' aligned'

                self.ligand_conc_df = pd.concat([self.ligand_conc_df,new_row])

        self.create_unique_sensor_names()

        return None

    def align_dissociation(self,sensor_names,inplace=True,new_names = False):

        """
        Align the BLI traces based on the signal before the dissociation step(s)

        Args:

            sensor_names (str or list): name of the sensor(s) to align
            inplace (bool):             if True, the alignment is done in place, otherwise a new sensor is created

        Results:

            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df

        """

        if not isinstance(sensor_names, list):
            sensor_names = [sensor_names]

        # Find the index of the dissociation steps
        dissociation_steps_indices = self.df_steps.index[self.df_steps['Type'] == 'DISASSOC'].to_numpy()

        sensor_indices = [self.sensor_names.index(sensor_name) for sensor_name in sensor_names]

        use_new_names = not inplace or (inplace and new_names)

        for sensor in sensor_indices:

            # Determine the new sensor name
            new_sensor_name = self.sensor_names[sensor] + ' diss. aligned' if use_new_names else self.sensor_names[sensor]

            if new_sensor_name in self.sensor_names and not inplace:
                new_sensor_name += ' rep'

            ys = (self.ys).copy()

            for diss_step_index in dissociation_steps_indices:

                #  Subtract the difference between the steps
                last_point = np.mean(self.ys[sensor][diss_step_index-1][-10:])
                next_point = np.mean(self.ys[sensor][diss_step_index][:10])

                diff = next_point - last_point

                value = self.ys[sensor][diss_step_index] - diff

                if inplace:

                    self.ys[sensor][diss_step_index] = value

                else:

                    ys[sensor][diss_step_index] = value

            if inplace:

                #Replace in the ligand conc df the sensor name
                self.ligand_conc_df['Sensor'] = self.ligand_conc_df['Sensor'].replace(self.sensor_names[sensor],new_sensor_name)

                self.sensor_names[sensor] = new_sensor_name

            else:

                self.xs.append(self.xs[sensor])
                self.ys.append(ys)
                self.sensor_names.append(self.sensor_names[sensor] + ' diss. aligned')

                # Add the new sensor name to the ligand conc df
                previous_row        = self.ligand_conc_df[self.ligand_conc_df['Sensor'] == self.sensor_names[sensor]]
                new_row             = previous_row.copy()
                new_row['Sensor']   = new_sensor_name
                new_row['SampleID'] = new_row['SampleID'] + ' diss. aligned'

                self.ligand_conc_df = pd.concat([self.ligand_conc_df,new_row])

        self.create_unique_sensor_names()

        return None

    def discard_steps(self,sensor_names,step_types=['KREGENERATION','LOADING']):

        """

        Discard the steps of the sensors in the list

        Args:

            sensor_names (str or list): name of the sensor(s) to analyse
            step_types (str or list):    type of the steps to discard

        Results:

            It modifies the attributes self.xs, self.ys, self.sensor_names and self.ligand_conc_df

        """
        if not isinstance(sensor_names, list):
            sensor_names = [sensor_names]

        if not isinstance(step_types, list):
            step_types = [step_types]

        sensor_indices = [self.sensor_names.index(sensor_name) for sensor_name in sensor_names]

        for step_type in step_types:

            step_indices = self.df_steps.index[self.df_steps['Type'] == step_type].to_numpy()

            for step_index in step_indices:

                for sensor in sensor_indices:

                    self.ys[sensor][step_index] = np.repeat(np.nan,len(self.ys[sensor][step_index]))

        return None

    def get_step_xy(self,sensor_name,location_loading,
                    location_sample,step_type='ASSOC',
                    replicate=1,type='y'):

        """
        Return the x or y values of a certain step

        Args:

            sensor_name (str): name of the sensor
            location_sample (int):    column location of the sample. If zero, we assume we only have one location
            location_loading (int):    column location of the loading. If zero, we assume we only have one location
            step_type (str):   type of the step, ASSOC or DISASSOC only
            replicate (int):   replicate number
            type (str):        x or y

        Returns:

            x or y (np.n) values of the step

        """

        # Try to convert to integer, the variables location_sample, location_loading and Replicate
        # If it fails, raise an error
        try:
            location_sample  = int(location_sample)
            location_loading = int(location_loading)
            replicate        = int(replicate)
        except ValueError:
            raise ValueError("location_sample, location_loading and replicate must be integers")

        # Verify we have the correct data types
        for var, expected_type, name in [
            (sensor_name, str, "sensor_name"),
            (location_sample, int, "location_sample"),
            (location_loading, int, "location_loading"),
            (step_type, str, "step_type"),
            (replicate, int, "replicate"),
            (type, str, "type"),
        ]:
            if not isinstance(var, expected_type):
                raise TypeError(f"{name} must be a {expected_type.__name__}")

        sensor = self.sensor_names.index(sensor_name)

        cond   = self.df_steps['Type']             == 'ASSOC'

        if location_sample != 0:

            cond = np.logical_and(cond,self.df_steps['Column_location'] == str(location_sample))

        if location_loading != 0:

            cond = np.logical_and(cond,self.df_steps['Loading_location'] == str(location_loading))

        step_index = self.df_steps[cond].index.to_numpy()[replicate-1] + 1*(step_type == 'DISASSOC')

        if type == 'x':

            time = self.xs[sensor][step_index]

            try:

                # Find if we have single-cycle kinetics or multi-cycle kinetics
                previous_type = self.df_steps['Type'][step_index - 2]
                single_cycle = previous_type == step_type

            except:
                single_cycle = False

            if not single_cycle:

                # If the step_type is an association step, subtract the first data point
                if step_type == 'ASSOC':

                    time = time - time[0]

                # If the step_type is a dissociation step, subtract the first data point of the previous step
                else:

                    time = time - self.xs[sensor][step_index-1][0]

            else:

                # Find the index of the first association, from the single cycle
                # Iterate over the previous steps, two at a time, until we find a step that is not a step of the same type
                for i in range(step_index-2,0,-2):
                    if self.df_steps['Type'][i] != step_type:
                        break

                # If the step_type is an association step, subtract the first data point
                if step_type == 'ASSOC':

                    time = time - self.xs[sensor][i+2][0]

                # If the step_type is a dissociation step, subtract the first data point of the previous step
                else:

                    time = time - self.xs[sensor][i+1][0]

            return time

        else:

            return self.ys[sensor][step_index]

class BLI_experiment(surface_based_experiment):

    """
    A class used to represent a BLI experiment

    Attributes
    ----------

        name (str):                 name of the experiment
        fns (list):                 list of file names (length n, one per sensor)
        xs (list):                  list of x values (time, length n, one per sensor)
        ys (list):                  list of y values (length n, one per sensor)
        exp_info (list):            list of dictionaries with experimental information
        step_info (list):           list of dictionaries with step information
        no_steps (int):             number of steps
        no_sensors (int):           number of sensors
        sample_column (np.ndarray): array of sample column information (96 elements, one per well)
        sample_row (np.ndarray):    array of sample row information    (96 elements, one per well)
        sample_type (list):         list of sample types   (96 elements, one per well)
        sample_id (list):           list of sample ids     (96 elements, one per well)
        sensor_names_unique (list): list of unique sensor names (length n, one per sensor)
        sensor_names (list):        list of sensor names (length n, one per sensor)
        df_steps    (pd.DataFrame): dataframe with the steps information
        ligand_conc_df (pd.DataFrame): dataframe with the ligand concentration information

        ligand_conc_df.head(2):
        Sensor Analyte_location  Concentration_micromolar SampleID  Replicate Experiment
        A1        5                    0.1300       wt          1          t
        B1        5                    0.0692       wt          1          t

        traces_loaded (bool):       True if traces are loaded
        sample_plate_loaded (bool): True if sample plate information is loaded
        sample_conc (np.array):     array with the sample concentrations (96 elements, one per well)
        sample_conc_labeled (list): list with the sample concentrations labeled (96 elements, one per well)
        steps_performed (pd.DataFrame): dataframe with the steps performed

        steps_performed.head(2):
        #Step         Type              Column   Time
        1            Regeneration      1        25
        2  BaselineNeutralization      2        90
        3         BaselineLoading      3        50
    """

    def __init__(self,name):

        """

        Initialize the instance

        """

        self.name = name
        self.fns = None
        self.xs = None
        self.ys = None
        self.exp_info = None
        self.step_info = None
        self.no_steps = None
        self.no_sensors = None
        self.sample_column = None
        self.sample_row = None
        self.sample_type = None
        self.sample_id = None
        self.sensor_names_unique = None
        self.type = 'BLI_experiment'

    def read_sensor_data(self, files,names=None):

        """

        Read the sensor data from the .frd files

        Results:

            It creates the attributes

                self.traces_loaded
                self.xs
                self.ys
                self.exp_info
                self.step_info
                self.no_steps
                self.no_sensors
                self.sensor_names
                self.df_steps
                self.ligand_conc_df

        """

        if names is None:
            names = files

        if not isinstance(files, list):
            files = [files]
            names = [names]

        fns = [fn for fn,name in zip(files,names) if '.frd' in name]

        if len(fns) < 1:
            self.traces_loaded = False
            return None
        else:
            self.fns = fns

        # Initialize dictionaries with data
        xs, ys, all_expinfo, all_stepinfo, more_info = [], [], [], [], []
        for fn in fns:
            # Load file
            tree = ET.parse(fn)
            root = tree.getroot()

            # Extract experimental info
            all_expinfo.append(etree_to_dict(root.find('ExperimentInfo')))

            # Initialize lists for each file
            x_values, y_values, step_info = [], [], []
            more_dict = {'FlowRate': [], 'StepType': [], 'StepName':[], 'StepStatus':[], 'ActualTime':[], 'CycleTime':[]}
            for step in root.find('KineticsData'):
                for step_x in step.findall('AssayXData'):
                    # Convert string to binary
                    data_text = bytes(step_x.text, 'utf-8')
                    # Convert to base64
                    decoded = base64.decodebytes(data_text)
                    # And now convert to float32 array
                    data_values = np.array(np.frombuffer(decoded, dtype=np.float32))
                    x_values.append(data_values)
                for step_y in step.findall('AssayYData'):
                    # Convert string to binary
                    data_text = bytes(step_y.text, 'utf-8')
                    # Convert to base64
                    decoded = base64.decodebytes(data_text)
                    # And now convert to float32 array
                    data_values = np.array(np.frombuffer(decoded, dtype=np.float32))
                    y_values.append(data_values)
                for step_data in step.findall('CommonData'):
                    step_info.append(etree_to_dict(step_data))
                for tag in ['FlowRate', 'StepType', 'StepName', 'StepStatus', 'ActualTime', 'CycleTime']:
                    for step_data in step.findall(tag):
                        more_dict[tag].append(step_data.text)

            xs.append(x_values)
            ys.append(y_values)
            all_stepinfo.append(combine_dicts(step_info))
            more_info.append(more_dict)

        # Merge all_stepinfo and more_info
        for i in range(len(all_stepinfo)):
            all_stepinfo[i] = {**all_stepinfo[i], **more_info[i]}

        # Fill instance
        self.xs   = xs
        self.ys   = ys
        self.exp_info  = all_expinfo
        self.step_info = all_stepinfo
        # Convert text to floats
        self.convert_to_numbers()

        self.no_steps = len(self.step_info[0]['ActualTime'])

        self.no_sensors = len(self.fns)

        self.sensor_names = [self.exp_info[i]['SensorName'] for i in range(self.no_sensors)]

        steps_names = self.step_info[0]['StepName']
        steps_types = self.step_info[0]['StepType']
        steps_start = self.step_info[0]['StartTime'] / 1000 # To seconds
        steps_loc   = self.step_info[0]['SampleLocation']

        self.df_steps = pd.DataFrame({'#Step':np.arange(len(steps_names))+1,
                                      'Name':steps_names,
                                      'Type':steps_types,
                                      'Start':steps_start,
                                      'Column_location':steps_loc})

        # We need to include the loading location in self.df_steps
        loading_location = []
        for row in self.df_steps.iterrows():
            step_type = row[1]['Type']
            if step_type == 'ASSOC':
                # Find the previous loading step
                for i in range(row[0],0,-1):
                    if self.df_steps.iloc[i]['Type'] == 'LOADING':
                        loading_location.append(self.df_steps.iloc[i]['Column_location'])
                        break
            else:
                loading_location.append(np.nan)

        self.df_steps['Loading_location'] = loading_location

        sensor_locs_all = np.concatenate([self.step_info[i]['SampleLocation']       for i in range(self.no_sensors)])
        sensor_type_all = np.concatenate([self.step_info[i]['StepType']             for i in range(self.no_sensors)])

        sensor_molar_conc_all = np.concatenate([self.step_info[i]['MolarConcentration']  for i in range(self.no_sensors)])
        sensor_mass_conc_all  = np.concatenate([self.step_info[i]['Concentration']       for i in range(self.no_sensors)])

        sensor_conc_all = []

        for i in range(len(sensor_molar_conc_all)):

            if sensor_molar_conc_all[i] < 0:

                sensor_conc_all.append(sensor_mass_conc_all[i])

            else:

                sensor_conc_all.append(sensor_molar_conc_all[i])

        sensor_conc_all = np.array(sensor_conc_all)

        sample_id_all   = np.concatenate([self.step_info[i]['SampleID']   for i in range(self.no_sensors)])

        sensor_name_rep = np.concatenate([np.repeat(self.exp_info[i]['SensorName'],len(self.step_info[i]['Concentration'])) for i in range(self.no_sensors)])

        conc_units = np.concatenate([self.step_info[i]['MolarConcUnits'] for i in range(self.no_sensors)])

        df_all = pd.DataFrame({'Sensor':sensor_name_rep,
                               'Analyte_location':sensor_locs_all,
                               'Type':sensor_type_all,
                               'Concentration_micromolar':sensor_conc_all,
                               'ConcUnits':conc_units,
                               'SampleID':sample_id_all})

        # For each association step, find the corresponding loading step
        # and add the loading step to the dataframe as a column next to the association

        # Add empty column to the data frame
        loading_location  = []
        loading_sample_id = []

        for i in range(len(df_all)):

            row = df_all.iloc[i]
            # Find if row is association step
            if row['Type'] == 'ASSOC':

                # Find the previous loading step
                for i in range(i,0,-1):

                    if df_all.iloc[i]['Type'] == 'LOADING':

                        loading_location.append(df_all.iloc[i]['Analyte_location'])
                        loading_sample_id.append(df_all.iloc[i]['SampleID'])
                        break


        # Keep only association or dissociation steps
        df = df_all[df_all['Type'] == 'ASSOC'].copy()

        # ADD column loading_location
        df['Loading_location']  = loading_location

        # Replace None with empty string in loading_sample_id
        loading_sample_id = [x if x is not None else '' for x in loading_sample_id]

        # Include the loading_sample_id, if we have more than one unique value
        unq_loading_ids = np.unique(loading_sample_id)

        if len(unq_loading_ids) > 1:

            # Combine the sample id with the loading id
            df['SampleID'] = df['SampleID'] + ' - ' + loading_sample_id

        # Remove the Type column
        df = df.drop(columns=['Type'])

        # Sort by location and sensor name
        df = df.sort_values(by=['Loading_location','Analyte_location','Sensor'])

        # Add rep column
        sizes = df.groupby(['Loading_location','Analyte_location','Sensor']).size().reset_index(name="Repetitions")

        rep_number = []

        for i in range(len(sizes)):
            rep_number.extend(np.arange(sizes['Repetitions'].iloc[i])+1)

        # Group by sensor and location
        df['Replicate'] = rep_number

        # Sort the dataframe first by sensor, second by Loading location,
        # Third by replicate and finally by location

        df = df.sort_values(by=['Analyte_location','Replicate','Loading_location','Sensor'])

        df['Factor'] = df.apply(lambda row: factor_conc_to_micro[row['ConcUnits']], axis=1)

        df['Concentration_micromolar'] = df['Concentration_micromolar'] * df['Factor']

        # Remove the factor column and conc units
        df = df.drop(columns=['Factor','ConcUnits'])

        # Add the experiment name
        df['Experiment'] = self.name

        self.ligand_conc_df = df

        self.create_unique_sensor_names()

        self.traces_loaded = True

        return None

    def read_sample_plate_info(self,files,names=None):

        """
        Read the sample plate information from the .fmf file

        Results:

            It creates the attributes

                self.sample_column
                self.sample_row
                self.sample_type
                self.sample_id
                self.sample_conc
                self.sample_conc_labeled
                self.sample_plate_loaded
                self.steps_performed

        """

        if names is None:
            names = files

        if not isinstance(files, list):
            files = [files]
            names = [names]

        index = next((i for i, s in enumerate(names) if 'ExpMethod.fmf' in s), None)

        if index is None:

            self.sample_plate_loaded = False
            return None

        file = files[index]

        tree = ET.parse(file)
        root = tree.getroot()

        sample_types     = [x.text for x in root.findall(".//SampleType")]
        sample_locations = [x.text for x in root.findall(".//SampleLoc")]
        sample_ids       = [x.text for x in root.findall(".//SampleID")]

        sample_conc_molar    = np.array([float(x.text) for x in root.findall(".//SampleMolarConc")])
        sample_conc_mass     = np.array([float(x.text) for x in root.findall(".//SampleConc")])

        sample_conc          = sample_conc_mass

        sel_ids = ['SAMPLE' in s for s in sample_types]

        counter = 0
        for i in range(len(sample_conc)):

            if sel_ids[i]:

                if sample_conc_molar[counter] > 0:

                    sample_conc[i] = sample_conc_molar[counter]

                counter += 1

        ConcUnits      = [x.text for x in root.findall(".//ConcUnits")][0]
        MolarConcUnits = [x.text for x in root.findall(".//MolarConcUnits")][0]

        factors = [factor_conc_to_micro[MolarConcUnits] if 'SAMPLE' in st else factor_conc_to_micro[ConcUnits] for st in sample_types]

        sample_conc = sample_conc * np.array(factors)
        sample_conc = np.round(sample_conc, 5)

        sample_column = np.array([int(re.sub(r'\D', '', text)) for text in sample_locations])
        sample_row    = np.array([re.sub(r'\d+', '', text)     for text in sample_locations])

        self.sample_column = sample_column
        self.sample_row    = sample_row
        self.sample_type   = sample_types
        self.sample_id     = sample_ids

        self.sample_conc    = sample_conc

        sample_conc_labeled = [f"{x} µM" if t == 'KSAMPLE' and x >= 0 else f"{x} µg/ml" if t != 'KSAMPLE' and x >= 0 else '' for x, t in zip(sample_conc, sample_types)]

        self.sample_conc_labeled = sample_conc_labeled

        data_name     = [x.text for x in root.findall(".//DataName")]
        assay_time    = [x.text for x in root.findall(".//AssayTime")]

        steps_info_df      = pd.DataFrame({'Type':data_name,'Time':assay_time})

        data_name     = [x.text for x in root.findall(".//StepDataName")]
        data_col      = [x.text for x in root.findall(".//SampleCol")]

        steps_performed = pd.DataFrame({'#Step':np.arange(len(data_name))+1,'Type':data_name,'Column':data_col})

        steps_performed = pd.merge(steps_performed, steps_info_df, on='Type', how='left')

        self.steps_performed = steps_performed

        self.sample_plate_loaded = True

        return None

    def convert_to_numbers(self):

        """

        Convert the strings in the step info to numbers

        Results:

            It modifies the attribute self.step_info

        """

        # List of entries in step info
        entries = ['Concentration', 'MolarConcentration', 'MolecularWeight', 'Temperature', 'StartTime',
                   'AssayTime', 'FlowRate', 'ActualTime', 'CycleTime']

        for entry in entries:
            for sensor in range(len(self.fns)):
                # Do sanity check
                try:
                    self.step_info[sensor][entry] = np.array(self.step_info[sensor][entry], dtype=float)
                except:

                    print("Erroneous entry found for %s and sensor %i: %s" % (entry, sensor, self.step_info[sensor][entry]))
                    print("Will set it to -1. Needs to be corrected")
                    # Correct erroneous value
                    for i in range(len(self.step_info[sensor][entry])):
                        try:
                            float(self.step_info[sensor][entry][i])
                        except:
                            self.step_info[sensor][entry][i] = -1
                    self.step_info[sensor][entry] = np.array(self.step_info[sensor][entry], dtype=float)
        return None

class Gator_experiment(surface_based_experiment):

    def __init__(self, name):

        self.name = name
        self.fns = None
        self.xs = None
        self.ys = None
        self.exp_info = None
        self.df_steps = None

        self.no_steps = None
        self.no_sensors = None
        self.sample_column = None
        self.sample_row = None
        self.sample_type = None
        self.sample_id = None
        self.sensor_names_unique = None
        self.ligand_conc_df = None
        self.type = 'Gator_experiment'

    def read_all_gator_data(self, files, names=None):

        self.traces_loaded = False

        if names is None:
            names = [os.path.basename(file) for file in files]

        # Find the file ExperimentStep.ini and read it
        for file, name in zip(files, names):
            if 'ExperimentStep.ini' in name:
                self.read_experiment_ini(file)

        # Find the file Setting.ini and read it
        for file, name in zip(files, names):
            if 'Setting.ini' in name:
                self.read_settings_ini(file)

        # Find the sensor csvs and read them
        self.read_sensor_data(files, names)

        return None

    def read_experiment_ini(self, file):

        """
        Read the experiment ini file

        Example format:

            [Experiment]
            Num=4
            [Experiment1]
            Step1=0,0,60,0
            Step2=1,60,104,1
            Step3=5,104,164,0
            Step4=6,164,464,2
            Step5=7,464,764,3
            Num=5
            ProbeIndex=12
            [ExperimentLog]
            Step1="{"StepType":0,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":0,"Speed":400,"Time":600,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step2="{"StepType":3,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":1,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step3="{"StepType":4,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":2,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step4="{"StepType":3,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":1,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step5="{"StepType":4,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":2,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step6="{"StepType":3,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":1,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step7="{"StepType":4,"AssayNumber":1,"PlateIndex":1,"PlateColumnType":12,"RowLocation":0,"Location":2,"Speed":1000,"Time":5,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":0,"Status":0}"
            Step8="{"StepType":2,"AssayNumber":1,"PlateIndex":0,"PlateColumnType":12,"RowLocation":0,"Location":12,"Speed":1000,"Time":60,"SumTime":0,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":true,"kst":0,"Status":0}"
            Step9="{"StepType":2,"AssayNumber":1,"PlateIndex":0,"PlateColumnType":12,"RowLocation":0,"Location":13,"Speed":400,"Time":44,"SumTime":60,"ProbeIndex":0,"PickerUsed":1,"AllChannelUsed":[255],"ChannelUsed":0,"bAssayFirstStep":false,"kst":1,"Status":1}"

        """

        with open(file, 'r') as f:
            lines = f.read().splitlines()

            # Find the line [ExperimentLog]
            for i, line in enumerate(lines):
                if line.startswith('[ExperimentLog]'):
                    start_index = i + 1
                    break

            # Count the number of items of the next line, when stripping by commas
            nr_items = len(lines[start_index].split(','))

            step_strings = []
            # Find all lines with the same amount of items and store them
            for line in lines[start_index:]:
                if len(line.split(',')) == nr_items:
                    step_strings.append(line)

        # List to hold dictionaries
        step_dicts = []

        # Iterate over the list of strings
        for step_string in step_strings:
            # Extract the JSON part of the string
            json_part = step_string.split('=', 1)[1].strip('"')
            # Convert the JSON string into a dictionary
            step_dict = json.loads(json_part)
            # Append the dictionary to the list
            step_dicts.append(step_dict)

        # Convert the list of dictionaries into a DataFrame
        df = pd.DataFrame(step_dicts)

        # Remove extra columns
        df = df[['kst', 'Location', 'Time', 'SumTime', 'AssayNumber', 'Speed']]
        # Add the column with the stepnumber to the beginning
        df.insert(0, 'StepNumber', np.arange(len(df)))

        # Change column name 'Location' to 'Column'
        df.rename(columns={'Location': 'Column'}, inplace=True)

        # If column location is larger than 12, then subtract 12 from it
        df['Column'] = df['Column'] - 11

        # If kst equals 2, then the step is a association
        # If kst equals 3, then the step is a dissociation

        # Guess if the step was 'ASSOC' or 'DISASSOC',
        # depending on the value of kst
        df['kst'] = df['kst'].replace({2: 'ASSOC', 3: 'DISASSOC', 1: 'LOADING', 0: 'OTHER'})

        # Change column name 'kst' to 'Type'
        df.rename(columns={'kst': 'Type'}, inplace=True)

        # Update the 'Type' column to 'BASELINE' for rows preceding 'ASSOC'
        df.loc[df['Type'].shift(-1) == 'ASSOC', 'Type'] = 'BASELINE'

        self.df_steps = df

        return None

    def read_settings_ini(self, file):

        """
        Read the settings ini file

        We assume a plate format, so wells 1 to 8 correspond to the first column,
        wells 9 to 16 correspond to the second column, etc.

        Example format:

            [BasicInformation]
            PreExperiment="{"AssayDescription":"","AssayUser":"Marianna","CreationTime":"03-31-2025 14:01:50","ModificationTime":"03-31-2025 14:44:10","StartExperimentTime":"03-31-2025 14:44:10","EndExperimentTime":"03-31-2025 15:55:27","AnotherSavePath":"","AssayType":2,"PreAssayShakerASpeed":400,"PreAssayShakerBSpeed":400,"PreAssayTime":600,"GapTime":100,"AssayShakerATemperature":30,"AssayShakerBTemperature":30,"MachineRealType":"Prime","idleShakerATemperature":30,"idleShakerBTemperature":30,"PlateAType":0,"bPlateAFlat":false,"bRegeneration":true,"bRegenerationStart":true,"RegenerationNum":99999,"RegenerationMode":0,"ShakerASpeedDeviation":10,"ShakerBSpeedDeviation":10,"ShakerATempDeviation":2,"ShakerBTempDeviation":2,"SpectrometerTempDeviation":1,"ParentName":"K Result","ResultName":"4x ab vs egfr 03-31-2025 14-44-10","SoftWareVersion":null,"SeriesNo":"GA00092","ErrorChannelList":[0,0,0,0,0,0,0,0],"ErrorPreAssayList":[0,0,0,0,0,0,0,0],"ErrorSampleList":[],"ErrorRegenerationList":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"AnalysisSettingList":[],"ReportSettingList":[],"PlateColumns":[12,12],"threshold":Infinity,"bsingle":false,"bThresh":false,"strFileID":null}"
            ExperimentType=Kinetics
            SaveTime=03/31/2025 15:55:27
            Frequency=10Hz
            DataVersion=2.15.5.1221
            [Regeneration]
            RegenerationEnable=1
            Repeat=3
            RTime=5
            RSpeed=1000
            NTime=5
            NSpeed=1000
            [SampleKinetics]
            Well1=4,-1,-1,-1,,
            Well2=4,-1,-1,-1,,
            Well3=4,-1,-1,-1,,
            Well4=4,-1,-1,-1,,
            Well5=4,-1,-1,-1,,
            Well6=4,-1,-1,-1,,
            Well7=4,-1,-1,-1,,
            Well8=4,-1,-1,-1,,
            Well9=5,30,150,200,,Cetuximab
            Well10=5,30,150,200,,Cetuximab
            Well11=5,30,150,200,,Cetuximab

        """

        with open(file, 'r') as f:
            lines = f.read().splitlines()

            # Find the molar concentration units
            # below the line MolarConcentrationUnit=NM
            for i, line in enumerate(lines):
                if line.startswith('MolarConcentrationUnit='):
                    units = line.split('=')[1].strip().lower()
                    if units == 'nm':
                        factor = 1e-3
                    elif units == 'mm':
                        factor = 1e3
                    elif units == 'pm':
                        factor = 1e-6
                    elif units == 'm':
                        factor = 1e6
                    else:
                        factor = 1
                    break

            # Find the line [SampleKinetics]
            for i, line in enumerate(lines):
                if line.startswith('[SampleKinetics]'):
                    start_index = i + 1
                    break

            # Count the number of items of the next line, when stripping by commas
            nr_items = len(lines[start_index].split(','))

            sample_strings = []
            # Find all lines with the same amount of items and store them
            for line in lines[start_index:]:
                if len(line.split(',')) == nr_items:
                    sample_strings.append(line)

            # Leave lines starting with 'Well'
            sample_strings = [line for line in sample_strings if line.startswith('Well')]

            # Now we need to convert wells that go from 1 to 96 into a combination of letters and numbers
            # They are ordered by rows, each column has 8 rows
            concs = []

            locations = np.concatenate([np.repeat(x, 8) for x in range(1, 13)], axis=0)
            sensors = np.concatenate([(65 + np.arange(8)) for _ in range(1, 13)], axis=0)
            # Apply chr function to sensors
            sensors = [chr(x) for x in sensors]

            labels = []

            for string in sample_strings:
                # We need to extract the concentration : Well49=1,17,170,100,,EGFR
                # We need to split by commas and take the fourth element
                conc = string.split('=')[1].split(',')[3]
                conc = float(conc) if conc != '' else 0

                concs.append(conc)

                # Extract the label, last column
                label = string.split(',')[-1]
                labels.append(label)

            # Create the ligand concentration dataframe
            df_sensor_prev = pd.DataFrame({'Sensor': sensors,
                                           'Location': locations,
                                           'Concentration_micromolar': concs,
                                           'SampleID': labels})

            loading_labels = []
            loading_columns = []
            sample_labels = []
            location = []
            sensor = []
            concs = []

            loading_columns_extended = []

            # Iterate over the associations and find loading column

            for row_index in range(len(self.df_steps)):

                step_type = self.df_steps.iat[row_index, self.df_steps.columns.get_loc('Type')]
                if step_type == 'ASSOC':

                    loading_column = find_loading_column(self.df_steps, row_index)

                    df_temp = df_sensor_prev[df_sensor_prev['Location'].isin([loading_column])]
                    loading_labels += df_temp['SampleID'].to_list()

                    loading_columns.append(loading_column)

                    association_column = self.df_steps.iat[row_index, self.df_steps.columns.get_loc('Column')]

                    location += [association_column for _ in range(8)]
                    sensor += [chr(65 + i) for i in range(8)]

                    loading_columns_extended += [loading_column for _ in range(8)]

                    df_temp = df_sensor_prev[df_sensor_prev['Location'].isin([association_column])]
                    sample_labels += df_temp['SampleID'].to_list()

                    concs += df_temp['Concentration_micromolar'].to_list()
                else:
                    loading_columns.append('NA')

            # Include the loading column in df_steps
            self.df_steps['Loading_location'] = loading_columns

            sample_id = [x + '-' + y for x, y in zip(sample_labels, loading_labels)]

            replicates = [1 for _ in range(len(sample_id))]

            unq_sample_id = list(set(sample_id))

            for unq in unq_sample_id:

                # Find the index with that sample id
                idx = [i for i, x in enumerate(sample_id) if x == unq]

                replicates[idx[0]] = 1

                idx_cnt = 1
                for i in idx[1:]:
                    if sample_id[i] == sample_id[i - 1]:
                        idx_cnt += 1
                    replicates[i] = (idx_cnt - 1) // 8 + 1

            # Create the ligand concentration dataframe
            df_sensor = pd.DataFrame({'Sensor': sensor,
                                      'Analyte_location': location,
                                      'Concentration_micromolar': concs,
                                      'Loading_location': loading_columns_extended,
                                      'SampleID': sample_id,
                                      'Replicate': replicates,
                                      'Experiment': self.name})

            # Multiply concentrations by the factor
            df_sensor['Concentration_micromolar'] = df_sensor['Concentration_micromolar'] * factor

            sensor_names = df_sensor['Sensor'].unique().tolist()
            self.sensor_names = sensor_names

            self.create_unique_sensor_names()

            self.ligand_conc_df = df_sensor

        return None

    def read_sensor_data(self, files, names=None):

        if names is None:
            names = files

        if not isinstance(files, list):
            files = [files]
            names = [names]

        fns   = [fn for fn, name in zip(files, names) if '.csv' in name and 'Channel' in name]
        names = [name for name in names if '.csv' in name and 'Channel' in name]

        if len(fns) < 1:
            return None
        else:
            self.fns = fns

        # Initialize dictionaries with data
        xs = [[] for _ in range(len(self.sensor_names))]
        ys = [[] for _ in range(len(self.sensor_names))]

        # Find steps with data, where the SumTime is zero,
        # or if not zero, it is the same as the previous step

        steps_with_data = []

        for i in range(1, len(self.df_steps)):

            sum_time_i = self.df_steps.iat[i, self.df_steps.columns.get_loc('SumTime')]
            sum_time_prev = self.df_steps.iat[i - 1, self.df_steps.columns.get_loc('SumTime')]

            if sum_time_i == 0 or sum_time_i == sum_time_prev:
                continue
            else:
                steps_with_data.append(i - 1)

        # Create a copy of self.df_steps that will contain all the data and remove steps with no data
        self.df_steps_all = self.df_steps.copy()

        # Select the steps with data only
        self.df_steps = self.df_steps.iloc[steps_with_data]

        # Remove the stepnumber column
        self.df_steps = self.df_steps.drop(columns=['StepNumber'])

        # Reset the stepnumber column
        self.df_steps['StepNumber'] = np.arange(len(self.df_steps))

        # Reset the index of step_info_with_data
        self.df_steps = self.df_steps.reset_index(drop=True)

        # Rename the column 'Column' to 'Column_location' in df_steps
        self.df_steps.rename(columns={'Column': 'Column_location'}, inplace=True)

        # Move the column StepNumber to the first position
        column_to_move = self.df_steps.pop('StepNumber')
        self.df_steps.insert(0, 'StepNumber', column_to_move)

        # Move the column Loading_location to the fifth position
        column_to_move = self.df_steps.pop('Loading_location')
        self.df_steps.insert(4, 'Loading_location', column_to_move)

        # Convert the columns Loading_location and Column_location to string, for compatibility with get_step_xy
        self.df_steps['Loading_location'] = self.df_steps['Loading_location'].astype(str)
        self.df_steps['Column_location']  = self.df_steps['Column_location'].astype(str)

        self.no_steps = len(self.df_steps)

        # xs will have one element per sensor
        # each element will have as many subelements as steps

        # Extract the assay number and channel number from the file names
        # and sort the files by assay number and channel number
        assay_number = [int(name.split('_')[1].split('Channel')[0]) for name in names]
        channel_number = [int(name.split('Channel')[1].split('.')[0]) for name in names]

        # Find the indices that would sort the arrays by assay number and channel number
        sorted_indices = np.lexsort((channel_number, assay_number))
        # Sort the file names and assay numbers
        fns   = [fns[i] for i in sorted_indices]
        names = [names[i] for i in sorted_indices]

        for fn, name in zip(fns, names):

            try:

                assay_number_i = int(name.split('_')[1].split('Channel')[0])
                channel_number_i = int(name.split('Channel')[1].split('.')[0])

                df = pd.read_csv(fn, skiprows=1, header=None)

                x = df.iloc[:, 1].to_numpy(dtype=float)
                y = df.iloc[:, 0].to_numpy(dtype=float)

                # Subset self.df_steps_all according to the assay number
                temp_df = self.df_steps_all[self.df_steps_all['AssayNumber'] == assay_number_i]

                # Filter which step numbers are in step_with_data
                temp_df = temp_df[temp_df['StepNumber'].isin(steps_with_data)]

                # Cumulative time for this assay
                temp_df['SumTime2'] = temp_df['Time'].cumsum()

                # Find starting sumtime for the assay
                start_sumtime = temp_df.iloc[0, temp_df.columns.get_loc('SumTime')]

                # For each step, include the data into xs and ys
                for i in range(len(temp_df)):

                    if i == 0:
                        start_time = 0
                    else:
                        start_time = temp_df.iloc[i - 1, temp_df.columns.get_loc('SumTime2')]

                    step_time = temp_df.iloc[i, temp_df.columns.get_loc('Time')]

                    end_time = start_time + step_time

                    sel_idx = np.logical_and(x >= start_time, x < end_time)

                    y_sel = y[sel_idx]
                    x_sel = x[sel_idx]

                    # Add the overall sumtime to x
                    x_sel += start_sumtime

                    # Append the second column to xs
                    xs[channel_number_i - 1].append(x_sel)

                    # Append the first column to ys
                    ys[channel_number_i - 1].append(y_sel)

                self.traces_loaded = True
                self.xs = xs
                self.ys = ys

            except:

                pass

        return None

class KinGenie_csv:

    """
    A class used to represent a KinGenie csv file, which can be exported from the Simulation panel

    Example:

        Time	Signal	Smax	Analyte_concentration_micromolar_constant	Cycle
        0	0	5	0.1	1
        0.5	0.00498512709952641	5	0.1	1
        1	0.00994051854979453	5	0.1	1
        1.5	0.0148661741864963	5	0.1	1
        2	0.0197623830271752	5	0.1	1
        2.5	0.0246293723316269	5	0.1	1
        3	0.0294671469344429	5	0.1	1
        3.5	0.0342760214300913	5	0.1	1
        4	0.0390561847030298	5	0.1	1
        4.5	0.043807638803537	5	0.1	1
        5	0.0485304884515415	5	0.1	1

    Attributes

        name (str):                 name of the experiment
        fn (str):                   file name
        xs (np.array):              list of x values (time, length n, one per sensor)
        ys (np.array):              list of y values (length n, one per sensor)
        no_sensors (int):           number of sensors
        sensor_names (list):        list of sensor names (length n, one per sensor)
        sensor_names_unique (list): list of unique sensor names (length n, one per sensor)
        ligand_conc_df (pd.DataFrame): dataframe with the ligand concentration information

    """

    def __init__(self,name):

        self.name = name
        self.xs   = None
        self.ys   = None
        self.no_sensors          = None
        self.sensor_names        = None
        self.sensor_names_unique = None
        self.traces_loaded       = False
        self.type = 'kingenie_csv'

    def create_unique_sensor_names(self):

        self.sensor_names_unique = [self.name + ' ' + sensor_name for sensor_name in self.sensor_names]

        return None

    def read_csv(self,file):

        """
        Read the KinGenie csv file

        Results:

            It creates the attributes

                self.xs
                self.ys
                self.no_sensors
                self.sensor_names
                self.sensor_names_unique

        """

        df = pd.read_csv(file)

        # Read only the first Smax value
        smax_unq = df['Smax'].unique()
        df       = df[df['Smax'] == smax_unq[0]]

        # Detect if we have the column cycle
        is_single_cycle = 'Cycle' in df.columns and len(df['Cycle'].unique()) > 1

        # Generate one 'fake' sensor per ligand concentration, if it's not a single cycle
        ligand_conc = df['Analyte_concentration_micromolar_constant'].unique()

        # Remove zero concentration
        ligand_conc = ligand_conc[ligand_conc != 0]

        if is_single_cycle:
            self.no_sensors = 1
        else:
            self.no_sensors = len(ligand_conc)

        # Use fake names
        self.sensor_names = ['sim. sensor ' + str(i+1) for i in range(self.no_sensors)]

        # Initiate self.xs and self.ys, one empty list per sensor
        self.xs = [[] for _ in range(self.no_sensors)]
        self.ys = [[] for _ in range(self.no_sensors)]

        steps_start = []

        if not is_single_cycle:

            # Now populate, for each sensor self.xs and self.ys
            for i in range(self.no_sensors):

                # Find the start index of the association phase
                start_idx = df[df['Analyte_concentration_micromolar_constant'] == ligand_conc[i]].index[0]

                if i == self.no_sensors - 1:
                    end_idx = df.shape[0]
                else:
                    end_idx   = df[df['Analyte_concentration_micromolar_constant'] == ligand_conc[i+1]].index[0]

                df_temp = df.iloc[start_idx:end_idx]

                # Extract the association phase and dissociation phase
                df_temp_asso  = df_temp[df_temp['Analyte_concentration_micromolar_constant'] > 0]
                df_temp_disso = df_temp[df_temp['Analyte_concentration_micromolar_constant'] == 0]

                # Extract the time values of the association/dissociation phase
                time_asso = df_temp_asso['Time'].to_numpy()
                time_diss = df_temp_disso['Time'].to_numpy()

                if i == 1:

                    steps_start.append(0)
                    steps_start.append(time_asso[-1])

                # Populate self.xs
                self.xs[i].append(time_asso)
                self.xs[i].append(time_diss)

                # Extract the signal values
                signal_asso = df_temp_asso['Signal'].to_numpy()
                signal_diss = df_temp_disso['Signal'].to_numpy()

                # Populate self.ys
                self.ys[i].append(signal_asso)
                self.ys[i].append(signal_diss)

            # Now generate a fake ligand concentration df
            df_sensor = pd.DataFrame({'Sensor':self.sensor_names,
                                      'Analyte_location':'1',
                                      'Concentration_micromolar':ligand_conc,
                                      'SampleID':'simulation',
                                      'Replicate':1,
                                      'Loading_location': '2',
                                      'Experiment':self.name})

            self.ligand_conc_df = df_sensor

            # Generate a fake steps df
            steps_names = ['Association','Dissociation']
            steps_types = ['ASSOC','DISASSOC']
            steps_loc   = ['1','3']
            steps_load  = ['2', 'NA']

            self.df_steps = pd.DataFrame({'#Step':np.arange(len(steps_names))+1,
                                          'Name':steps_names,
                                          'Type':steps_types,
                                          'Start':steps_start,
                                          'Column_location':steps_loc,
                                          'Loading_location':steps_load})

        else:

            cycles = df['Cycle'].unique()
            cycles.sort()

            n_cycles = len(cycles)

            analyte_location = []
            steps_loc        = []

            for cycle in cycles:

                analyte_location.append(str(cycle+1))
                steps_loc.append(str(cycle+1))
                steps_loc.append(str(999))

                df_temp = df[df['Cycle'] == cycle]

                # Extract the association phase and dissociation phase
                df_temp_asso = df_temp[df_temp['Analyte_concentration_micromolar_constant'] > 0]
                df_temp_disso = df_temp[df_temp['Analyte_concentration_micromolar_constant'] == 0]

                # Extract the time values of the association/dissociation phase
                time_asso = df_temp_asso['Time'].to_numpy()
                time_diss = df_temp_disso['Time'].to_numpy()

                steps_start.append(time_asso[0])
                steps_start.append(time_asso[-1])

                # Populate self.xs
                self.xs[0].append(time_asso)
                self.xs[0].append(time_diss)

                # Extract the signal values
                signal_asso = df_temp_asso['Signal'].to_numpy()
                signal_diss = df_temp_disso['Signal'].to_numpy()

                # Populate self.ys
                self.ys[0].append(signal_asso)
                self.ys[0].append(signal_diss)

            # Now generate a fake ligand concentration df
            df_sensor = pd.DataFrame({'Sensor': self.sensor_names*n_cycles,
                                      'Analyte_location': analyte_location,
                                      'Concentration_micromolar': ligand_conc,
                                      'SampleID': 'simulation',
                                      'Replicate': 1,
                                      'Loading_location': '1',
                                      'Experiment': self.name})

            self.ligand_conc_df = df_sensor

            # Generate a fake steps df
            steps_names = ['Association', 'Dissociation'] * n_cycles
            steps_types = ['ASSOC', 'DISASSOC']           * n_cycles
            steps_load  = ['1', 'NA']                   * n_cycles

            self.df_steps = pd.DataFrame({'#Step': np.arange(len(steps_names)) + 1,
                                          'Name': steps_names,
                                          'Type': steps_types,
                                          'Start': steps_start,
                                          'Column_location': steps_loc,
                                          'Loading_location': steps_load})

        self.create_unique_sensor_names()

        self.traces_loaded = True

        return None

    def get_step_xy(self,sensor_name,location_loading,
                    location_sample,step_type='ASSOC',
                    replicate=1,type='y'):

        """
        Return the x or y values of a certain step

        Args:

            sensor_name (str): name of the sensor
            location_sample (int):    column location of the sample. If zero, we assume we only have one location
            location_loading (int):    column location of the loading. If zero, we assume we only have one location
            step_type (str):   type of the step, ASSOC or DISASSOC only
            replicate (int):   replicate number
            type (str):        x or y

        Returns:

            x or y (np.n) values of the step

        """

        # Try to convert to integer, the variables location_sample, location_loading and Replicate
        # If it fails, raise an error
        try:
            location_sample  = int(location_sample)
            location_loading = int(location_loading)
            replicate        = int(replicate)
        except ValueError:
            raise ValueError("location_sample, location_loading and replicate must be integers")

        # Verify we have the correct data types
        for var, expected_type, name in [
            (sensor_name, str, "sensor_name"),
            (location_sample, int, "location_sample"),
            (location_loading, int, "location_loading"),
            (step_type, str, "step_type"),
            (replicate, int, "replicate"),
            (type, str, "type"),
        ]:
            if not isinstance(var, expected_type):
                raise TypeError(f"{name} must be a {expected_type.__name__}")

        sensor = self.sensor_names.index(sensor_name)

        cond   = self.df_steps['Type']             == 'ASSOC'

        if location_sample != 0:

            cond = np.logical_and(cond,self.df_steps['Column_location'] == str(location_sample))

        if location_loading != 0:

            cond = np.logical_and(cond,self.df_steps['Loading_location'] == str(location_loading))

        step_index = self.df_steps[cond].index.to_numpy()[replicate-1] + 1*(step_type == 'DISASSOC')

        if type == 'x':

            time = self.xs[sensor][step_index]

            try:

                # Find if we have single-cycle kinetics or multi-cycle kinetics
                previous_type = self.df_steps['Type'][step_index - 2]
                single_cycle  = previous_type == step_type

            except:

                single_cycle = False

            if not single_cycle:

                # If the step_type is an association step, subtract the first data point
                if step_type == 'ASSOC':

                    time = time - time[0]

                # If the step_type is a dissociation step, subtract the first data point of the previous step
                else:

                    time = time - self.xs[sensor][step_index-1][0]

                return time

            else:

                # Find the index of the first association, from the single cycle
                # Iterate over the previous steps, two at a time, until we find a step that is not a step of the same type

                for i in range(step_index-2,0,-2):
                    if self.df_steps['Type'][i] != step_type:
                        time = time - self.xs[sensor][i + 1 + int(step_type == 'ASSOC')][0]
                        return time

                time =  time - self.xs[sensor][0][0]

                return time

        else:

            return self.ys[sensor][step_index]

class Kinetics_fitter:

    """

    A class used to fit kinetics data with shared thermodynamic parameters

    Attributes
    ----------

        names (list):                   list of experiment names
        assoc (list):                   list of association signals, one per replicate,
                                        each signal is a numpy matrix of size n*m where n is the number of time points and m
        assoc (list):                   list containing the association signals
        disso (list):                   list containing the dissociation signals
        smax_id (list):                 list containing the Smax IDs (maximum amplitude identifiers)


        disso (list):                   list of dissociation signals
                                        each signal is a numpy matrix of size n*m where n is the number of time points and m
                                        is the number of ligand concentrations (different from zero)
        lig_conc (list):                list of ligand concentrations, one per replicate
                                        each element contains a numpy array with the ligand concentrations
        time_assoc (list):              list of time points for the association signals, one per replicate
                                        each element contains a numpy array with the time points
        time_disso (list):              list of time points for the dissociation signals, one per replicate
                                        each element contains a numpy array with the time points
        signal_ss (list):               list of steady state signals, one per replicate
                                        each element contains a numpy array with the steady state signals
        signal_ss_fit (list):           list of steady state fitted signals, one per replicate
                                        each element contains a numpy array with the steady state fitted signals
        signal_assoc_fit (list):        list of association kinetics fitted signals, one per replicate
                                        each element contains a numpy array with the association kinetics fitted signals
        signal_disso_fit (list):        list of dissociation kinetics fitted signals, one per replicate
                                        each element contains a numpy array with the dissociation kinetics fitted signals
        fit_params_kinetics (pd.Dataframe):     dataframe with the fitted parameters for the association / dissociation kinetics
        fit_params_ss (pd.Dataframe):           dataframe with the values of the fitted parameters - steady state

    """

    def __init__(self,time_assoc_lst,association_signal_lst,lig_conc_lst,time_diss_lst,
                 dissociation_signal_lst=None,smax_id=None,name_lst=None,is_single_cycle=False):

        self.names            = name_lst
        self.assoc_lst        = association_signal_lst
        self.disso_lst        = dissociation_signal_lst
        self.lig_conc_lst     = lig_conc_lst
        self.time_assoc_lst   = time_assoc_lst
        self.time_disso_lst   = time_diss_lst
        self.is_single_cycle  = is_single_cycle

        self.signal_ss        = None  # Steady state signal
        self.signal_ss_fit    = None  # Steady state fitted signal
        self.signal_assoc_fit = None  # Association kinetics fitted signal
        self.signal_disso_fit = None  # Association kinetics fitted signal

        # We need to rearrange smax_id to start at 0, then 1, then 2, etc.

        smax_id_unq, idx = np.unique(smax_id, return_index=True)
        smax_id_unq = smax_id_unq[np.argsort(idx)]

        smax_id_new = []
        for i,unq in enumerate(smax_id_unq):
            smax_id_new += [i for _ in range(len(np.where(smax_id == unq)[0]))]

        self.smax_id = smax_id_new

    def get_steady_state(self):

        """
        This function calculates the steady state signal and groups it by smax ID
        """

        signals_steady_state = [np.median(assoc[-10:]) for assoc in self.assoc_lst]

        # Create a new list, that will contain one element per Smax ID
        # each element will be a list of steady-state signals

        self.signal_ss            = [] # convert to list of lists, one list per unique smax id
        self.lig_conc_lst_per_id  = [] # convert to list of lists, one list per unique smax id

        # Obtain the smax guesses, the maximum signal
        smax_guesses_unq    = [] # One element per association signal
        smax_guesses_shared = [] # One element per Smax ID

        for i,smax_id in enumerate(np.unique(self.smax_id)):

            idx = np.where(self.smax_id == smax_id)[0]

            self.signal_ss.append([signals_steady_state[i] for i in idx])
            self.lig_conc_lst_per_id.append([self.lig_conc_lst[i] for i in idx])

            smax = np.max(self.signal_ss[i])
            smax_guesses_unq.append(smax*1.5)
            smax_guesses_shared += [smax*1.5 for _ in range(len(idx))]

        self.smax_guesses_unq    = smax_guesses_unq
        self.smax_guesses_shared = smax_guesses_shared

        return None

    def clear_fittings(self):

        self.signal_ss_fit        = None  # Steady state fitted signal
        self.signal_assoc_fit     = None  # Association kinetics fitted signal
        self.signal_disso_fit     = None  # Association kinetics fitted signal
        self.fit_params_kinetics  = None  # Fitted parameters for the association / dissociation kinetics
        self.fit_params_ss        = None  # Values of the fitted parameters - steady state

        return None

    def fit_steady_state(self):

        '''
        This function fits the steady state signal
        '''

        self.clear_fittings()

        if self.signal_ss is None:

            self.get_steady_state()

        try:

            Kd_init = np.median(self.lig_conc_lst_per_id[0])
            p0      = [Kd_init] + [np.max(signal) for signal in self.signal_ss]

            kd_min  = np.min(self.lig_conc_lst_per_id[0]) / 1e3
            kd_max  = np.max(self.lig_conc_lst_per_id[0]) * 1e3

            # Find the upper bound for the Kd
            if Kd_init >= 1:
                upper_bound = 1e3
            else:
                upper_bound = 1e2

            low_bounds  = [kd_min]  + [x*0.5          for x in p0[1:]]
            high_bounds = [kd_max]  + [x*upper_bound  for x in p0[1:]]

            fit, cov, fit_vals = fit_steady_state_one_site(
                self.signal_ss,self.lig_conc_lst_per_id,
                p0,low_bounds,high_bounds)

            self.Kd_ss   = fit[0]
            Smax         = fit[1:]

            # Re fit the data if the Smax is at the upper bound
            difference = np.abs(Smax - high_bounds[1:])

            # Check if any element is less than one
            if any(difference < 1):

                high_bounds[1:] = [x*1e3 for x in high_bounds[1:]]
                fit, cov, fit_vals = fit_steady_state_one_site(
                    self.signal_ss,self.lig_conc_lst_per_id,
                    p0,low_bounds,high_bounds)

                self.Kd_ss   = fit[0]
                Smax         = fit[1:]

            self.p0 = p0
            self.low_bounds = low_bounds
            self.high_bounds = high_bounds

            self.signal_ss_fit = fit_vals

            n = np.sum([len(signal) for signal in self.signal_ss])
            p = len(p0)

            rss_desired = get_desired_rss(concat_signal_lst(self.signal_ss), concat_signal_lst(fit_vals), n, p)

            minKd, maxKd = steady_state_one_site_asymmetric_ci95(
                self.Kd_ss, self.signal_ss, self.lig_conc_lst_per_id, p0[1:],
                low_bounds[1:], high_bounds[1:], rss_desired)

            # Create a dataframe with the fitted parameters
            df_fit = pd.DataFrame({'Kd [µM]': self.Kd_ss,
                                   'Kd_min95': minKd,
                                   'Kd_max95': maxKd,
                                   'Smax': Smax,
                                   'Name': self.names})

            self.fit_params_ss = df_fit

            self.Smax_upper_bound_factor = 50

            # Set factor to limit the Smax upper bound
            if self.Kd_ss >= 10:
                self.Smax_upper_bound_factor = 1e3
            elif self.Kd_ss >= 1:
                self.Smax_upper_bound_factor = 1e2

        except:
            # Generate a df filled with NAs

            df_fit = pd.DataFrame({'Kd [µM]': np.nan,
                                      'Kd_min95': np.nan,
                                      'Kd_max95': np.nan,
                                      'Smax': np.nan,
                                      'Name': self.names})

            self.fit_params_ss = df_fit

        # Fit the steady state signal
        return None

    def fit_one_site_association(self,shared_smax=True):

        # Initial guess for Kd_ss for single_cycle_kinetics
        if self.is_single_cycle:
            self.Kd_ss = np.median(self.lig_conc_lst)

        self.clear_fittings()

        # Try to fit first the dissociation curves to get a better estimate of Koff
        try:

            self.fit_one_site_dissociation()

            p0          = [self.Kd_ss,self.Koff_disso]
            low_bounds  = [self.Kd_ss/7e2,self.Koff_disso/7e2]
            high_bounds = [self.Kd_ss*7e2,self.Koff_disso*7e2]

            self.clear_fittings()

        except:

            p0 = [self.Kd_ss,0.01]

            low_bounds  = [self.Kd_ss/7e2,1e-7]
            high_bounds = [self.Kd_ss*7e2,10]

        if shared_smax:

            smax_guesses = self.smax_guesses_unq
        else:
            smax_guesses = self.smax_guesses_shared

        p0 += smax_guesses

        low_bounds  += [x/10 for x in p0[2:]]
        high_bounds += [x*self.Smax_upper_bound_factor for x in p0[2:]]

        self.p0           = p0
        self.low_bounds  = low_bounds
        self.high_bounds = high_bounds

        fit, cov, fit_vals = fit_one_site_association(
            self.assoc_lst,self.time_assoc_lst,self.lig_conc_lst,
            p0,low_bounds,high_bounds,smax_idx=self.smax_id,shared_smax=shared_smax
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

    def fit_one_site_dissociation(self,time_limit=0):

        self.clear_fittings()

        disso_lst     = self.disso_lst
        time_disso_lst = self.time_disso_lst

        # Fit only some data. If time_limit = 0, fit all data
        if time_limit > 0:

            disso_lst      = [x[t < (np.min(t)+time_limit)] for x,t in zip(disso_lst,time_disso_lst)]
            time_disso_lst = [t[t < (np.min(t)+time_limit)] for t in time_disso_lst]

        p0 = [0.1] + [np.max(signal) for signal in disso_lst]

        low_bounds  = [1e-7] + [x/4 for x in p0[1:]]
        high_bounds = [10]   + [x*4 for x in p0[1:]]

        self.p0           = p0
        self.low_bounds  = low_bounds
        self.high_bounds = high_bounds

        fit, cov, fit_vals = fit_one_site_dissociation(disso_lst,time_disso_lst,
                                                       p0,low_bounds,high_bounds)

        self.signal_disso_fit = fit_vals

        self.Koff_disso = fit[0]

        # generate dataframe with the fitted parameters
        df_fit = pd.DataFrame({'k_off': self.Koff_disso,
                               'S0': fit[1:]})

        self.fit_params_kinetics = df_fit

        # Add the fitted errors to the dataframe
        error     = np.sqrt(np.diag(cov))
        rel_error = error/fit * 100

        df_error = pd.DataFrame({'k_off [1/s]': rel_error[0],
                                 'S0':          rel_error[1:]})

        self.fit_params_kinetics_error = df_error

        return  None

    def single_exponential(self):

        # Fit one exponential to each signal
        k_obs = []

        for y,t in zip(self.assoc_lst,self.time_assoc_lst):

            try:

                k_obs.append(fit_single_exponential(y,t)[2])

            except:

                k_obs.append(np.nan)

        self.k_obs = k_obs

        return None

    def fit_one_site_assoc_and_disso(self,shared_smax=True,fixed_t0=True,fit_ktr=False):

        # Initial guess for Kd_ss for single_cycle_kinetics
        if self.is_single_cycle:
            self.Kd_ss = np.median(self.lig_conc_lst)

        self.clear_fittings()

        kd_low_bound = np.min([self.Kd_ss/1e3, np.min(self.lig_conc_lst)])

        # Try to fit first the dissociation curves to get a better estimate of Koff
        try:

            self.fit_one_site_dissociation()

            p0          = [self.Kd_ss,self.Koff_disso]
            low_bounds  = [kd_low_bound,self.Koff_disso/7.5e2]
            high_bounds = [self.Kd_ss*7.5e2,self.Koff_disso*7.5e2]

        except:

            p0 = [self.Kd_ss,0.01]

            low_bounds  = [kd_low_bound,1e-7]
            high_bounds = [self.Kd_ss*7e2,10]

        if shared_smax:

            smax_guesses = self.smax_guesses_unq

        else:

            smax_guesses = self.smax_guesses_shared

        p0 += smax_guesses

        low_bounds  += [x/2 for x in p0[2:]]
        high_bounds += [x*self.Smax_upper_bound_factor for x in p0[2:]]

        smax_param_start = 2

        if fit_ktr:

            for i, _ in enumerate(self.smax_guesses_unq):
                p0.insert(2, 1e-4)
                low_bounds.insert(2, 1e-7)
                high_bounds.insert(2, 1)
                smax_param_start += 1

        if not fixed_t0:

            for i,_ in enumerate(self.smax_guesses_unq):
                p0.insert(2,0)
                low_bounds.insert(2,-0.01)
                high_bounds.insert(2,0.1)
                smax_param_start += 1

        self.p0          = p0
        self.low_bounds  = low_bounds
        self.high_bounds = high_bounds

        if fit_ktr:

            fit, cov, fit_vals_assoc, fit_vals_disso = fit_one_site_assoc_and_disso_ktr(
                self.assoc_lst,self.time_assoc_lst,self.lig_conc_lst,
                self.disso_lst,self.time_disso_lst,
                p0,low_bounds,high_bounds,
                smax_idx=self.smax_id,
                shared_smax=shared_smax,
                fixed_t0=fixed_t0
            )

            # Refit the data if the Ktr is at the upper bound
            if 1 - fit[2] < 0.01:

                p0[2] = 1
                low_bounds[2]  = 0.01
                high_bounds[2] = 100
                fit, cov, fit_vals_assoc, fit_vals_disso = fit_one_site_assoc_and_disso_ktr(
                    self.assoc_lst,self.time_assoc_lst,self.lig_conc_lst,
                    self.disso_lst,self.time_disso_lst,
                    p0,low_bounds,high_bounds,
                    smax_idx=self.smax_id,
                    shared_smax=shared_smax,
                    fixed_t0=fixed_t0
                )

        else:

            fit, cov, fit_vals_assoc, fit_vals_disso = fit_one_site_assoc_and_disso(
                self.assoc_lst,self.time_assoc_lst,self.lig_conc_lst,
                self.disso_lst,self.time_disso_lst,
                p0,low_bounds,high_bounds,
                smax_idx=self.smax_id,
                shared_smax=shared_smax,
                fixed_t0=fixed_t0)

        self.signal_assoc_fit = fit_vals_assoc
        self.signal_disso_fit = fit_vals_disso

        self.Kd_asso   = fit[0]
        self.Koff_asso = fit[1]
        self.Smax_asso = fit[smax_param_start:]

        # Create a dataframe with the fitted parameters
        df_fit = pd.DataFrame({'Kd [µM]':     self.Kd_asso,
                               'k_off [1/s]': self.Koff_asso,
                               'Smax':        self.Smax_asso})

        # Include the Kon, derived from the Kd and Koff
        df_fit['(Derived) k_on [1/µM/s]'] = df_fit['k_off [1/s]'] / df_fit['Kd [µM]']

        error     = np.sqrt(np.diag(cov))
        rel_error = error/fit * 100

        df_error = pd.DataFrame({'Kd [µM]':   rel_error[0],
                                 'k_off [1/s]': rel_error[1],
                                 'Smax': rel_error[smax_param_start:]})

        if fit_ktr:

            n_ktr      = len(np.unique(self.smax_id))
            idx_start  = 2+(not fixed_t0)*n_ktr
            ktrs       = fit[idx_start:smax_param_start]
            ktrs_error = rel_error[idx_start:smax_param_start]

            if not shared_smax:

                ktr_all       = expand_parameter_list(ktrs, self.smax_id)
                ktr_error_all = expand_parameter_list(ktrs_error, self.smax_id)

                df_fit['Ktr']   = ktr_all
                df_error['Ktr'] = ktr_error_all

            else:

                df_fit['Ktr']   = ktrs
                df_error['Ktr'] = ktrs_error

        self.fit_params_kinetics       = df_fit
        self.fit_params_kinetics_error = df_error

        return  None

    def calculate_ci95(self,shared_smax=True,fixed_t0=True,fit_ktr=False):

        try:

        # Compute asymmetrical 95% confidence intervals
            if not fit_ktr:

                exp_signal_concat = concat_signal_lst(self.assoc_lst+self.disso_lst)
                fit_signal_concat = concat_signal_lst(self.signal_assoc_fit+self.signal_disso_fit)
                n                 = len(exp_signal_concat)
                p                 = len(self.p0)

                rss_desired = get_desired_rss(exp_signal_concat,
                                              fit_signal_concat,
                                              n, p)


                Kd_min, Kd_max = one_site_assoc_and_disso_asymmetric_ci95(
                    self.Kd_asso,rss_desired,
                    self.assoc_lst,self.time_assoc_lst,
                    self.lig_conc_lst,
                    self.disso_lst,self.time_disso_lst,
                    self.p0[1:],self.low_bounds[1:],self.high_bounds[1:],
                    self.smax_id,shared_smax=shared_smax,fixed_t0=fixed_t0)

                p0 = self.p0[:1] + self.p0[2:]
                low_bounds = self.low_bounds[:1] + self.low_bounds[2:]
                high_bounds = self.high_bounds[:1] + self.high_bounds[2:]

                koff_min, koff_max  = one_site_assoc_and_disso_asymmetric_ci95_koff(
                    self.Koff_asso,rss_desired,
                    self.assoc_lst,self.time_assoc_lst,
                    self.lig_conc_lst,
                    self.disso_lst,self.time_disso_lst,
                    p0,low_bounds,high_bounds,
                    self.smax_id,shared_smax=shared_smax,fixed_t0=fixed_t0)

                # Convert ci95_Kd and ci95_koff to a nice Table
                header = ['Parameter','95% CI lower','95% CI upper']
                row1   = ['Kd [µM]',Kd_min,Kd_max]
                row2   = ['k_off [1/s]',koff_min,koff_max]

                df = pd.DataFrame([row1,row2],columns=header)

                self.fit_params_kinetics_ci95 = df

        except:

            # Generate an empty dataframe if the procedure did not work
            self.fit_params_kinetics_ci95 = pd.DataFrame()

        return  None

class Kinetics_analyzer:

    def __init__(self):

        self.experiments      = {}
        self.experiment_names = []

    def delete_experiment(self,experiment_names):

        if not isinstance(experiment_names, list):
            experiment_names = [experiment_names]

        for experiment_name in experiment_names:

            if experiment_name in self.experiment_names:

                del self.experiments[experiment_name]
                self.experiment_names.remove(experiment_name)

        return None

    def add_experiment(self,experiment,experiment_name):

        self.delete_experiment(experiment_name)

        self.experiments[experiment_name] = experiment
        self.experiment_names.append(experiment_name)

        return None

    def init_fittings(self):

        self.fittings       = {}
        self.fittings_names = []

        return None

    def add_fitting(self,fitting,fitting_name):

        if fitting_name in self.fittings_names:

                del self.fittings[fitting_name]
                self.fittings_names.remove(fitting_name)

        self.fittings[fitting_name] = fitting
        self.fittings_names.append(fitting_name)

        return None

    def get_experiment_properties(self, variable,fittings=False):

        if fittings:

            return [getattr(self.fittings[fitting_name], variable) for fitting_name in self.fittings_names]

        else:

            return [getattr(self.experiments[exp_name], variable) for exp_name in self.experiment_names]

    def merge_ligand_conc_df(self):

        # Combine ligand concentration data frames from experiments
        dfs = [exp.ligand_conc_df for exp in self.experiments.values()]
        df  = pd.concat(dfs, ignore_index=True)

        # Add a 'Select' column with all values set to True
        df['Select'] = True

        columns_to_convert     = ["Analyte_location", "Loading_location", "Replicate"]
        df[columns_to_convert] = df[columns_to_convert].apply(lambda col: col.astype(int))

        # For each combination of SampleID, Experiment, and Replicate, create a unique ID
        df['Smax_ID'] = (df['SampleID'].astype(str) +
                         df['Experiment'].astype(str) +
                         df['Replicate'].astype(str))
        df['Smax_ID'] = pd.factorize(df['Smax_ID'])[0].astype(int)

        # Move specified columns to the end
        columns_to_move = ['Analyte_location', 'Loading_location', 'Replicate', 'Experiment']
        for col in columns_to_move:
            if col in df.columns:
                df = df[[c for c in df.columns if c != col] + [col]]

        # Rename the column 'Concentration_micromolar' to '[Analyte] (μM)'
        df.rename(columns={'Concentration_micromolar': '[Analyte] (μM)'}, inplace=True)

        # Delete columns if they have only one unique value
        for col in ['Experiment', 'Replicate', 'Analyte_location', 'Loading_location']:
            if col in df.columns and df[col].nunique() == 1:
                df.drop(columns=[col], inplace=True)

        # Create the vector of columns that cannot be edited
        possible_fixed_columns = ['Sensor', 'Analyte_location', 'Loading_location', 'Experiment', 'Replicate']
        fixed_columns_ids      = [df.columns.get_loc(col) for col in possible_fixed_columns if col in df.columns]

        self.combined_ligand_conc_df = df.copy()

        return None

    def generate_fittings(self,df):

        """
        Given a data frame, extract the kinetics data and create a fitting object for each sample.
        The data frame should contain the following columns:

        """

        # List of messages to be printed to the console
        messages          = []

        time_diss_all  =  []
        diss_all       =  []
        time_assoc_all =  []
        assoc_all      =  []

        # Reset dataframe index
        df.reset_index(drop=True, inplace=True)

        # Find if 'Experiment' is in the data frame column names
        df_colnames    = df.columns
        have_exp_column = 'Experiment'        in df_colnames
        have_rep_column = 'Replicate'         in df_colnames
        have_analyte_loc = 'Analyte_location' in df_colnames
        have_loading_loc = 'Loading_location' in df_colnames

        if not have_exp_column:
            exp = self.experiment_names[0]
        if not have_rep_column:
            replicate = 1
        if not have_analyte_loc:
            analyte_loc = 0
        if not have_loading_loc:
            loading_loc = 0

        # Find which column has the analyte concentration, with the word '[Analyte]' and rename it to 'Concentration_micromolar'
        for column in df.columns:
            if '[Analyte]' in column:
                df.rename(columns={column: 'Concentration_micromolar'}, inplace=True)

        for row in range(len(df)):
            if have_exp_column:
                exp = df.loc[row, 'Experiment']
            if have_rep_column:
                replicate = int(df.loc[row, 'Replicate'])
            if have_analyte_loc:
                analyte_loc = int(df.loc[row, 'Analyte_location'])
            if have_loading_loc:
                loading_loc = int(df.loc[row, 'Loading_location'])

            sensor = df.loc[row, 'Sensor']

            diss_all.append(self.experiments[exp].get_step_xy(sensor, loading_loc, analyte_loc, 'DISASSOC', replicate, 'y'))
            assoc_all.append(self.experiments[exp].get_step_xy(sensor, loading_loc, analyte_loc, 'ASSOC', replicate, 'y'))

            time_diss  = self.experiments[exp].get_step_xy(sensor, loading_loc, analyte_loc, 'DISASSOC', replicate, 'x')
            time_assoc = self.experiments[exp].get_step_xy(sensor, loading_loc, analyte_loc, 'ASSOC', replicate, 'x')

            time_diss_all.append(time_diss)
            time_assoc_all.append(time_assoc)

        # Find same samples with different experiment IDs
        unique_sample = df['SampleID'].unique()

        # Group by sample
        for unq in unique_sample:

            assoc_lst = []
            diss_lst  = []

            time_assoc_lst = []
            time_diss_lst  = []

            lig_conc_vec = []
            smax_id_vec  = []

            ids = df.index[df['SampleID'] == unq].tolist()

            # Iterate over the Smax IDs
            unq_smax_id = df.loc[ids, 'Smax_ID'].unique()

            message_for_sample = False

            for smax_id in unq_smax_id:
                ids = df.index[
                    (df['SampleID'] == unq) &
                    (df['Smax_ID'] == smax_id) &
                    (df['Concentration_micromolar'] > 0) &
                    (df['Select'])
                ].tolist()

                if len(ids) == 0:
                    continue

                if not message_for_sample:
                    messages.append("Creating a new fitting object for sample: " + unq)
                    message_for_sample = True

                assoc_sel = [assoc_all[i] for i in ids]
                diss_sel  = [diss_all[i]  for i in ids]

                time_assoc_sel = [time_assoc_all[i] for i in ids]
                time_diss_sel  = [time_diss_all[i]  for i in ids]
                lig_conc       = df.loc[ids, 'Concentration_micromolar'].tolist()
                smax_ids       = df.loc[ids, 'Smax_ID'].tolist()

                # Append to the lists
                assoc_lst.extend(assoc_sel)
                diss_lst.extend(diss_sel)
                time_assoc_lst.extend(time_assoc_sel)
                time_diss_lst.extend(time_diss_sel)
                lig_conc_vec.extend(lig_conc)
                smax_id_vec.extend(smax_ids)

            if len(assoc_lst) == 0:
                continue

            # Find initial times of the association signal
            time_inits = [t[0] for t in time_assoc_lst]

            # Find the index that sorts them
            sorted_indices = [index for index, _ in sorted(enumerate(time_inits), key=lambda x: x[1])]

            # Sort them
            time_assoc_lst = [time_assoc_lst[i] for i in sorted_indices]
            assoc_lst      = [assoc_lst[i]      for i in sorted_indices]
            time_diss_lst  = [time_diss_lst[i]  for i in sorted_indices]
            diss_lst       = [diss_lst[i]       for i in sorted_indices]
            lig_conc_vec   = [lig_conc_vec[i]   for i in sorted_indices]
            smax_id_vec    = [smax_id_vec[i]    for i in sorted_indices]

            fit = Kinetics_fitter(
                time_assoc_lst=time_assoc_lst,
                association_signal_lst=assoc_lst,
                lig_conc_lst=lig_conc_vec,
                time_diss_lst=time_diss_lst,
                dissociation_signal_lst=diss_lst,
                smax_id=smax_id_vec,
                name_lst=[f"{unq}_id_{smax}" for smax in unq_smax_id],
                is_single_cycle=any([t[0] > 1 for t in time_assoc_lst])
            )

            fit.get_steady_state()

            self.add_fitting(fit, unq)

        return messages

    def submit_steady_state_fitting(self):

        for kf in self.fittings.values():

            if not kf.is_single_cycle:
                kf.fit_steady_state()

            else:

                kf.Smax_upper_bound_factor = 1e2 # Normal values for lower than micromolar affinity

        return None

    def submit_kinetics_fitting(self,fitting_model,fitting_region,linkedSmax=False):

        for kf in self.fittings.values():

            try:

                if fitting_model == 'one_to_one' and fitting_region == 'association_dissociation':
                    kf.fit_one_site_assoc_and_disso(shared_smax=linkedSmax)
                    kf.single_exponential()

                if fitting_model == 'one_to_one_mtl' and fitting_region == 'association_dissociation':
                    kf.fit_one_site_assoc_and_disso(shared_smax=linkedSmax, fit_ktr=True)

                if fitting_model == 'one_to_one' and fitting_region == 'association':
                    kf.fit_one_site_association(shared_smax=linkedSmax)
                    kf.single_exponential()

                if fitting_model == 'one_to_one' and fitting_region == 'dissociation':
                    kf.fit_one_site_dissociation()

            except:

                pass

        return None

    def calculate_asymmetric_error(self,shared_smax=True,fixed_t0=True,fit_ktr=False):

        for kf in self.fittings.values():

            kf.calculate_ci95(
                shared_smax=shared_smax,
                fixed_t0=fixed_t0,
                fit_ktr=fit_ktr
            )

        return None

test0 = False
if test0:
    bli = BLI_experiment('test')
    folder = './www/test_bli_folder/'
    files = os.listdir(folder)
    files = [folder + file for file in files]

    # sort files by name
    files.sort()

    bli.read_sensor_data(files)

    print(bli.sensor_names)

    bli.align_association(bli.sensor_names)
    bli.align_dissociation(bli.sensor_names)
    bli.subtraction(bli.sensor_names[:7],bli.sensor_names[7])

    pyKinetics = Kinetics_analyzer()
    pyKinetics.add_experiment(bli, 'test')
    pyKinetics.merge_ligand_conc_df()
    pyKinetics.init_fittings()

    # Set only the first 7 rows to True
    pyKinetics.combined_ligand_conc_df['Select'] = False
    pyKinetics.combined_ligand_conc_df.loc[:6, 'Select'] = True

    print(pyKinetics.combined_ligand_conc_df)

    pyKinetics.generate_fittings(pyKinetics.combined_ligand_conc_df)

    fit  = pyKinetics.fittings[pyKinetics.fittings_names[0]]

    fit.get_steady_state()
    fit.fit_steady_state()
    fit.fit_one_site_assoc_and_disso()

    # Plot the fittings
    import matplotlib.pyplot as plt

    x1 = np.concatenate(fit.time_assoc_lst,axis=0)
    x2 = np.concatenate(fit.time_disso_lst,axis=0)
    y1 = np.concatenate(fit.assoc_lst,axis=0)
    y2 = np.concatenate(fit.disso_lst,axis=0)

    y1fit = np.concatenate(fit.signal_assoc_fit,axis=0)
    y2fit = np.concatenate(fit.signal_disso_fit,axis=0)

    plt.scatter(x1,y1,s=0.1,color='blue')
    plt.scatter(x2,y2,s=0.1,color='blue')
    plt.scatter(x1,y1fit,s=0.1,color='green')
    plt.scatter(x2,y2fit,s=0.1,color='green')

    plt.show()

