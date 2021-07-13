# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%

import matplotlib.pyplot as plt
import pandas
import numpy as np
import xarray as xr
from matplotlib.dates import num2date, date2num
import matplotlib
import datetime
import os 
import pdb

from zutils import qc_conventions

default_attrs = {
    'title': '', 
    'institution': 'The University of Western Australia', 
    'source': '', 
    'project': '', 
    'history': '', 
    'references': '', 
    'comment': '', 
    'Conventions': 'CF-1.7', 
    'trip_recovered': '', 
    'trip_deployed': '', 
    'site': '', 
    'instrument_make': '',  
    'instrument_model': '',
    'disclaimer': ''} 

class xrwrap():
    
    folder = ''
    file_ = ''
    
    so = None
    eo = None

    first_good = None
    last_good = None

    _qc_conv = qc_conventions.base()
    _default_attrs = default_attrs
    _outname = None
    
    _attrs = default_attrs

    def fullpath(self, i=0):
        
        if type(self.file_) == list:
            return '{folder}/{file_}'.format(folder=self.folder, file_=self.file_[i])
        else:
            return '{folder}/{file_}'.format(folder=self.folder, file_=self.file_)
    
    # def get_outname(self, keep_name='False', fileappend=''):

    #     if not keep_name:

    #         if len(self.file_) == 0:
    #             outname = self.file_[0]
    #         else:
    #             outname = self.file_[0] + ' PlusOther'
    #     else:
    #         outname = '{}_{}_{}{}'.format(self._attrs['project'], self._attrs['trip_recovered'], self._attrs['site'], fileappend)

    #     return outname

    def generate_outname(self, keep_name='False', fileappend=''):
        """
        Function to autogenerate a file name. Should not do this continuously.
        """
        if not keep_name:

            if len(self.file_) == 0:
                outname = self.file_[0]
            else:
                outname = self.file_[0] + ' PlusOther'
        else:
            outname = '{}_{}_{}{}'.format(self._attrs['project'], self._attrs['trip_recovered'], self._attrs['site'], fileappend)

        return outname

    # @property
    def outname(self, fileappend=''):
        """
        generate an outname only if outname is empty.
        """
        if self._outname is None:
            self._outname = self.generate_outname(fileappend='')

        return self._outname

    def load(self, folder, file_):
        """
        Initialise from netcdf rather than from a raw input file. 
        """

        self.folder = folder
        self.file_ = file_

        nc_file = '{folder}//{file_}'.format(folder=self.folder, file_=file_)
        
        ds = xr.open_dataset(nc_file)

        self.ds = ds
        self._attrs = ds.attrs
        self.ds.attrs = {}

        print('Loaded from NC. Class attributes taken from NC. NC attrs cleared.')

        pass

    def wrap(self, ds):
        """
        Initialise from an xarray dataset. Got to be carful here that the dataset is of the exact format or lots of errors will be thrown. Could used this if 2 NC files
        are merged outside of here and then brought back in. That's the only usage I can immagine.  
        """

        self.ds = ds
        self._attrs = ds.attrs
        self.ds.attrs = {}

        print('Wrapped an existing xarray dataset. Class attributes taken from the dataset. Dataset attrs cleared.')

    def export(self, final=False):
        """
        Base export class. Overloading will likely be necessary in many cases.
        """
        outname = self.generate_outname()
        
        self.ds.attrs = self._attrs

        if final:
            self.ds.attrs['Disclaimer'] = self.disclaimer
            outname = outname + 'finalised'

        self.ds.close() # Force close
        self.ds.to_dataframe().to_csv('{folder}//{file_}.csv'.format(folder=self.folder, file_=outname))

        nc_file = '{folder}//{file_}.nc'.format(folder=self.folder, file_=outname)
        self.ds.to_netcdf(path=nc_file)

        return self.ds # It may be useful in many cases to return somethng here. Just return self.sd for consistency with subclasses.  

    def time_trim(self, first_good, last_good):

        if first_good is None:
            first_good = self.ds.time.values[0]

        if last_good is None:
            last_good = self.ds.time.values[-1]

        self.first_good = first_good
        self.last_good = last_good

        print('Good data from {} to {}'.format(first_good, last_good))

        # This sometimes fails if the DataSet time is timezone aware and the trim times
        # are timezone naive. Dolfyn seems to throw both:
            #  type(rr.ds.time.values[1]) = pandas._libs.tslibs.timestamps.Timestamp
            #  type(rr.ds.time.values[1]) = numpy.datetime64
        # and only the former is datetime aware. Ultimately I'll be going away from Dolfyn. 
        # there is no fix for this, user must be aware of this. 
        self.ds = self.ds.sel(time=slice(self.first_good, self.last_good))

        print('Trimmed Time')

    @property  
    def qc_conv(self):
        return self._qc_conv

    @qc_conv.setter
    def qc_conv(self, new_qc_conv):
        """
        Set function for _qc_conv
        """
        if isinstance(new_qc_conv, qc_conventions):
            self._qc_conv = new_qc_conv
        else:
            raise(Exception(("Must ba a QC Convention object")))

    def run_qc_compliance(self, flag_name):
        """
        Wrapper for the qc_conv.run_compliance method.

        Basically it checks that your QC flag [ds data_var with the name flag_name] complies with the QC conventions specified by the _qc_conv attribute.
        """

        self.qc_conv.run_compliance(self.ds[flag_name])

    def associate_qc_flag(self, var_name, flag_name):
        """
        Assign a QC flag to a DataArray, and add the flag to the dataset if necessary.
        Inputs:
                - var_name: The name of the data variable which will have a QC flag associated. Cannot QC dimensions [coords].
                - flag_name: The name of QC flag to be associated/added. The prefix 'qc_' will be added automatically
        
        """

        # Add the qc_ prefix
        flag_name = 'qc_' + flag_name
        
        # Check if var_name exists
        if not var_name in self.ds.data_vars:
            raise(Exception('Variable does not exist'))

        var_shape = self.ds[var_name].values.shape

        # Check if flag_name exists
        if not flag_name in self.ds.data_vars:
            # if not create the flag
            flag_array = self.ds[var_name].copy()
            flag_array.attrs={}
            flag_array.values[:] = -999
            self.ds[flag_name] = flag_array

            # flag = -999*np.ones_like(self.ds[var_name])
            # self.ds[flag_name] = flag
            self.qc_conv.run_compliance(self.ds[flag_name])
        elif not self.ds[flag_name].values.shape == var_shape:
            # If so check the size of the variable against the size of the flag
            raise(Exception("The selected QC flag does not have the right dimensions for this variable."))

        # Check this variable doesn't already have a qc code. Multiples not yet allowed.
        if 'qc_variable' in self.ds[var_name].attrs:
            raise(Exception("A QC variable has already been associated to the variable {}!".format(var_name)))

        for i in self.ds.data_vars:
            if 'is_qc_flag' in self.ds[i].attrs and 'associated data variables' in self.ds[i].attrs: # It might be a QC flag
                if  self.ds[i].attrs['is_qc_flag']: # It's a QC flag
                    if var_name in self.ds[i].attrs['associated_data_variables'].split(';'):
                        raise(Exception("The variable {} has already been associated to the QC flag {}!".format(var_name, i)))

        # Associate this flag to this variable
        self.ds[flag_name].attrs['associated_data_variables'] += ';'+var_name
        self.ds[var_name].attrs['qc_variable'] = flag_name

        # Trim leading ; if necessary
        if self.ds[flag_name].attrs['associated_data_variables'][0] == ';':
            self.ds[flag_name].attrs['associated_data_variables'] = self.ds[flag_name].attrs['associated_data_variables'][1::]
    
    @property  
    def valid_qc_flags(self):
        """
        All valid QC flags. Checks for:
            - 'qc_' prefix
            - presence of 'is_qc_flag' attribute
            - value of 'is_qc_flag' attribute

        The checks could certainly be more thorough.
        """

        out = []
        for i in self.ds.data_vars:
            if len(i) < 4:
                continue
            if not i[0:3] == 'qc_':
                continue
            if not 'is_qc_flag' in self.ds[i].attrs:
                continue
            if self.ds[i].attrs['is_qc_flag']:
                out.append(i)

        return(out)

    def update_qc_flag(self, flag_name, index_name, start, end, flag_value, comment=None, verbose=False):
            """
            This is a base function to update a QAQC flag. Inputs:
                - flag_name: The name of the netcdf variable corrsponding to the QC flag being edited. Can be:
                                - a single flag name
                                - a list of flag names
                                - '*'
                            if '*' is used all QC flags will be updated 

                - index_name: The name of the netcdf variable which is being used to as an index to 
                                identify which points to in the QC flag variable are to be edited.
                - start: The first point in the index to change the value of
                - end: The last point in the index to change the value of
                - flag_value: new value which the QC flag is to become wherever the index is between start and end 
            """
            
            if isinstance(flag_name, str): 
                if flag_name == '*': 
                    flag_names = self.valid_qc_flags
                else:
                    flag_names = [flag_name]
            else: # Assume list
                flag_names = flag_name

            for flag_name in flag_names:

                if not index_name in self.ds[flag_name].coords:
                    continue

                logind1 = self.ds[index_name] >= start
                logind2 = self.ds[index_name] <= end
                logind = ~np.logical_and(logind1, logind2)

                # print(logind)

                self.ds[flag_name] = self.ds[flag_name].where(logind, flag_value)
                
                ind = np.where(logind)[0]
                logind = self.ds[flag_name] > 0

                if not comment is None: # Log comment on the QAQC.
                    string = 'Flagged {} values with code "{}" and user comment "{}"'.format(len(self.ds[index_name].values)-len(ind), flag_value, comment)
                    self.add_comment('O2 Metocean', string, data_var=flag_name)

                if verbose:
                    
                    print(ind)
                    print(len(ind))
                    
                    print('There are {} flagged data points.'.format(np.sum(logind.values)))
                    print('')        

    def flip_qc_value(self, flag_name, value_out=-999, value_in=0, verbose=False):
        """
        Flip one QC value for another. By default it will flip -999 to 0

        flag_name can be:
                                - a single flag name
                                - a list of flag names
                                - '*'
                            if '*' is used all QC flags will be updated 
        
        NOTE: doesn't actually check that the input variable is QC convention compliant. 
        
        """
        
        if isinstance(flag_name, str): 
                if flag_name == '*': 
                    flag_names = self.valid_qc_flags
                else:
                    flag_names = [flag_name]
        else: # Assume list
            flag_names = flag_name

        for flag_name in flag_names:

            logind = self.ds[flag_name] == value_out
            logind = ~logind

            self.ds[flag_name] = self.ds[flag_name].where(logind, value_in)

            if verbose:
                self.get_fig_text()

                logind = self.ds[flag_name] > 0
                print('There are {} positive flag points.'.format(np.sum(logind.values)))
                logind = self.ds[flag_name] < 0
                print('There are {} negative flag points.'.format(np.sum(logind.values)))
                print('')   

    def get_qaqc_var(self, var_name):
        """
        Retrun a QAQC'd copy of the data array by var_name. This just routes to the main function.
        """
        
        da = get_qaqc_var(self.ds, var_name)
        
        return da 

    # def get_qaqc_var_old(self, var_name, flag_name, axis=None):
    #     """
    #     Retrun a QAQC'd copy of the data array by var_name, with the corresponding flag_name

    #     optional input axis allows the flag to be reduced in dimension by a logical any.

    #     FUTURE UPDATE:
    #         - Get the QAQC flag name from the attributes of the variable
    #     """
        
    #     for i in np.arange(0, 5):
    #         print('WARNING: THE FLAG NAME SHOULD BE IN THE ATTRIBUTES OF THE DATA ARRAY.')

    #     da = self.ds[var_name].copy()
    #     ind = self.ds[flag_name].values > 0

    #     if not axis is None:
    #         ind = np.any(ind, axis=axis)

    #     da.values[ind] = np.nan
        
    #     return da 

    def has_dates(self):

        if type(self.so) == type(None):
            return False
        elif type(self.eo) == type(None):
            return False
        else:
            return True

    def get_fig_text(self, fileappend=''):
        """
        Get a text string to describing the dataset.
        """

        delim = ' | '
        s = ''
        k = self._attrs.keys()
        
        if 'project' in k:
            s = s + delim + self._attrs['project']
            
        if 'trip_recovered' in k:
            s = s + delim + self._attrs['trip_recovered']
            
        if 'site' in k:
            s = s + delim + self._attrs['site']
            
        if not fileappend == '':
            s = s + delim + fileappend
            
        s = s[len(delim)::]

        print(s)

        return s

    def add_history(self, author, string, data_var=None):
        """
        Add CF Compliant string to the history attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """

        attr = 'history'
        self.add_string(attr, author, string, data_var=data_var)
            
    def add_comment(self, author, string, data_var=None):
        """
        Add CF Compliant string to the comments attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """

        attr = 'comment'
        self.add_string(attr, author, string, data_var=data_var)
            
    def add_string(self, attr, author, string, data_var=None):
        """
        Add CF Compliant string to an attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """

        if data_var is None:
            obj = self.ds
        else:
            obj = self.ds[data_var]

        new_string = datetime.datetime.now().isoformat() + ': ' + '[{}]'.format(author) + ' ' + string
        
        if not attr in obj.attrs.keys():
            obj.attrs[attr] = new_string
        elif obj.attrs[attr] in ['', '?', 'blank']:
            obj.attrs[attr] = new_string
        else:
            obj.attrs[attr] = obj.attrs[attr] + ';' + new_string
            
        return new_string

    def update_attribute(self, attribute_name, attribute_value, strict=True):
        """
        This function updates the hidden attributes property of the class. The attribute must exist in the default attributes dictionary.
        """

        if (attribute_name in self._default_attrs) or (not strict):
            print('Setting attribute "{}" to "{}"'.format(attribute_name, attribute_value))
            self._attrs[attribute_name] = attribute_value
        else:
            raise(Exception('{} is not a valid attribute.'.format(attribute_name)))

    def update_attributes_with_dict(self, attribute_dict):
        """
        This function updates the hidden attributes property of the class. The attribute must exist in the default attributes dictionary.
        """

        for attribute_name in attribute_dict.keys():
            
            self.update_attribute(attribute_name, attribute_dict[attribute_name])

    def parse_attributes(self, ds_to_check=None):
        """
        Generic function to check consistency of attributes between the object itself and the properties of the class.
        """

        if ds_to_check is None:
            ds_to_check = self.ds

        print("Parsing attributes.")
        for i in ds_to_check.attrs.keys():
            if i in self._attrs.keys():
                print("{} is both a property of the object and an attribute of the dataset".format(i))
                if ds_to_check.attrs[i] == self._attrs[i]:
                    print("     ... and they are equal")
                else:
                    print("     ... and they NOT are equal!!!")

        ds_to_check.attrs = self._attrs

    @property
    def disclaimer(self):
        """

        """

        disc = """These data were prepared for a specific purpose at a specific site, which may or may not be disclosed in the data or metadata. 
        Use of these data for any purpose, without the express written consent of the author is not permitted. 
        The data are provided WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
        The author is not liable in any way for consequences arising from any use of these data or any derivatives in any application including but not limited to design or decision-making processes.
        """

        return disc
  
def get_qaqc_var(ds, var_name):
    """
    Retrun a QAQC'd copy of the data array by var_name
    """
    
    da = ds[var_name].copy() # return a QCd copy

    if not 'qc_variable' in ds[var_name].attrs:
        print('Variable {} has no QAQC flag. Returning raw.'.format(var_name))
    else:    
        flag_name = ds[var_name].attrs['qc_variable'] # Currently assumes only one qc_variable
        ind = ds[flag_name].values > 0 # Currently does not remove suspect data. Should be an option.
        print('Blanking {} values.'.format(np.sum(ind)))

        values=da.values
        values[ind] = np.nan
        da.values=values
    
    return da 

### Generic XR stuff

def select_calendar_month(X, year_month):
    """
    Return a subset of the Dataset or DataArray limited to the specified year and month. Uses rather crude indexing methods. I know 
    there are much more sopistocated inbuilt functions but I'm not familiar with them.
    """

    def calendar_month(year, month):
        """
        For a given year and month return the date of the begining of the month and the date of the beginning of the next month
        """
        start = datetime.datetime(year, month, 1)
        if month == 12:
            end = datetime.datetime(year+1, 1, 1)
        else:
            end = datetime.datetime(year, month+1, 1)
        print(start)
        print(end)
        return start, end
        
    year, month = year_month
    
    start, end = calendar_month(year, month)
    
    X_cm = X.sel(time = slice(start, end))
    
    return X_cm, [start, end]

def split_by_timegap(X, timename='time', hours=1):
    """
    Split a DataArray along the time dimension based on a set time gap.
    
        Returns a list of DataArrays with the gaps removed.
        
    """
        
    time = X[timename].values
    dt = np.diff(time)

    # print(len(time))
        
    i = [int(ii) for ii in np.where(dt>np.timedelta64(hours, 'h'))[0]] + [len(time)-1]
    i = [-1] + list(set(i))
    i.sort()
    
    # print('Split index')
    # print(i)
    
    Xs = [X.isel({timename: np.arange(i[j]+1, i[j+1])}) for j in np.arange(0, len(i)-1)]
    
    return Xs
