#%%
 
# if True:
#     from ..utils import validation
# else:
#     import validation

import datetime
import numpy as np

import pandas as pd

# Matlab datum is 0-Jan-1990. Datetime won't cop that.
matlab_datum  = datetime.datetime.strptime('1-1-1900', '%d-%m-%Y') - datetime.timedelta(days=1)

def linspacetime(start, end, n=None, dt_s=None):
    """
    Create a linearly spaced np array of datetime objects. Currently uses pandas but could be done 
    with linspace and one of the ordinal time tools if I ever start to care about dependencies.
    
    Inputs:
     - start: datetime object for start of array. 
     - end: datetime object for end of array. 
     - n: number of points in array [use None if using dt_s]
     - dt_s: number of seconds between points in the array [use None if using n]

     If specifying dt_s, array will start from start, and be spaced perfectly by dt_s, ending at or before
     end.
    """

    if n==None and dt_s==None:
        raise Exception('Must specify n or dt_s')
    elif (not n==None) and (not dt_s==None):
        raise Exception('Must specify only one of n or dt_s')
    elif not dt_s==None: # must calculate new end and n
        dt_hopeful = end-start 
        n_spaces = np.floor(dt_hopeful.seconds/dt_s)

        dt = datetime.timedelta(seconds=n_spaces*dt_s) # new total dt
        end = start + dt

        n = np.floor(dt.seconds/dt_s)+1
        
        # print('New end is: {}'.format(end))
        # print('New n is: {}'.format(n))
        # print('New t is: {}'.format(t))
        
    datelist = pd.date_range(start, end, periods=n).to_pydatetime()
    # datelist = pd.date_range(start, end, periods=n).to_list()

    return np.array(datelist)

def is_well_spaced(date, fs, epsilon=0.01):
    """
    A function to determine whether a time array is well spaced at fs. This is necessary to determine 
    if we wish to block 'blindly' average over time T = n/fs where n is a constant set numbers of 
    samples and fs is the sample frequency.
    
    Inputs:
     - date: numpy row vector of datetime objects. 
     - fs: sample frequncy in Hz
     
    Optional inputs:
     - epsilon: maximum deviation from the specified fs
    """

    # validation.is_np_date_rowvector(date)

    if type(date[0]) == datetime.datetime:
        max_dt = np.max(np.diff(date))
        min_dt = np.min(np.diff(date))
    elif type(date[0]) == np.datetime64:
        # max_dt = datetime.timedelta(seconds=np.max(np.diff(date)).item().total_seconds())
        # min_dt = datetime.timedelta(seconds=np.min(np.diff(date)).item().total_seconds())
        max_dt = datetime.timedelta(seconds=np.max(np.diff(date))/np.timedelta64(1, 's'))
        min_dt = datetime.timedelta(seconds=np.min(np.diff(date))/np.timedelta64(1, 's'))

    else:
        raise Exception

    max_allow = datetime.timedelta(seconds=1/fs)*(1+epsilon)
    min_allow = datetime.timedelta(seconds=1/fs)*(1-epsilon)

    if max_dt > max_allow:
        raise Exception('Maximum time spacing of {} greater than allowable [{}].'.format(max_dt, max_allow))
    elif min_dt < min_allow:
        raise Exception('Minimum time spacing of {} smaller than allowable [{}].'.format(min_dt, min_allow))
    else:
        print('Time vector is well spaced')

def datetime_2_timestamp(date_in):

    date_out = np.array([datetime.datetime.timestamp(date_in[i]) for i in np.arange(0, len(date_in))])

    return date_out

def timestamp_2_datetime(date_in):

    date_out = np.array([datetime.datetime.fromtimestamp(date_in[i]) for i in np.arange(0, len(date_in))])
    date_out = np.array(date_out)

    return date_out

def datetime64_2_timestamp(date_in):

    date_out = datetime64_2_datetime(date_in)
    date_out = datetime_2_timestamp(date_out)

    return date_out
    
def timestamp_2_datetime64(date_in):

    date_out = timestamp_2_datetime(date_in)
    date_out = datetime_2_datetime64(date_out)

    return date_out
    
def datetime64_2_datetime(date_in):
    
    date_out = [datetime64_2_datetime_single(date_in[i]) for i in np.arange(0, len(date_in))]
    date_out = np.array(date_out)

    return date_out

def datetime64_2_datetime_single(date_in):
    
    timestamp = ((date_in - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's'))
    date_out = datetime.datetime.utcfromtimestamp(timestamp)
    
    return date_out

def datetime_2_datetime64(date_in):
    
    date_out = [np.datetime64(date_in[i]) for i in np.arange(0, len(date_in))]
    date_out = np.array(date_out)

    return date_out

def datetime_2_datetime64_single(date_in):
    
    date_out = np.datetime64(date_in)

    return date_out

def str_to_dt(string, format_lang = 'python', dateformat = 'dd-mmm-yyyy HH:MM:SS:FFF'):

    if format_lang.lower() == 'matlab':
        dateformat = matformat_2_py(dateformat)

    return [datetime.datetime.strptime(string[i], dateformat) for i in np.arange(0, len(string))]
    
def num2date_lk(mpltime):
    """
    Reimplementation of datetime num2date to stop date errors. That occur continually when updating dolfyn or datetime. 
    
    This was copied from the Levi Kilcher's Dolfyn package, hence the lk suffix. 
    """
    
    if isinstance(mpltime, np.ndarray):
        try:
            n = len(mpltime)
        except TypeError:
            pass
        else:
            out = np.empty(n, dtype='O')
            for idx, val in enumerate(mpltime.flat):
                out[idx] = num2date_lk(val)
            out.shape = mpltime.shape
            return out
        
    if np.isnan(mpltime):
        return None
    return datetime.datetime.fromordinal(int(mpltime)) + datetime.timedelta(days=mpltime % 1)


def matplotlibtime_2_matlabdatenum():

    pass

def matlabdatenum_2_matplotlibtime():

    pass

def datenum_2_matlabdatenum():

    pass

def matlabdatenum_2_datenum():

    pass

def matformat_2_py(matformat):
    """
    Simple hack code to allow me to use matlab date formats pn python. At the moment I 
    just use common formats but should ultimately switch to a regext type approach to 
    find and replace. Not going for that as I should just learn the python formats. 
    """

    if matformat == 'dd-mmm-yyyy HH:MM:SS':
        pyformat = '%d-%b-%Y %H:%M:%S'
    elif matformat == 'dd-mmm-yyyy HH:MM:SS:FFF':
        pyformat = '%d-%b-%Y %H:%M:%S:%f'
    else:
        raise FormatException('Have not hard coded this date format yet')

    print('For personal reference, {0} translates to {1}'.format(matformat, pyformat))

    return pyformat

def FormatException(Exception):

    pass

def tests():

    n = 100
    ans = linspacetime(pd.datetime.today(), pd.datetime.today()+datetime.timedelta(1), n=100)
    if not len(ans) == n:
        raise(Exception)

    n = 100
    dt = datetime.timedelta(seconds=10.001)
    dt_s=1/10
    ans = linspacetime(pd.datetime.today(), pd.datetime.today()+dt, dt_s=dt_s)

    is_well_spaced(ans, 1/dt_s)

    ans = datetime_2_datetime64(ans)

    is_well_spaced(ans, 1/dt_s)

    print('There are {} elements in this {} second period spaced at {}.'.format(len(ans), dt.seconds, dt_s))

    print('Time Module OK')

pass

def np_datetime64_strftime(t, fmt):
    """
    Format a np.datetime64 or list of np.datetime64 objects as a string.
    
    returns: 
        list of strings
        
    """
    
    if not type(t) == list:
        t = [t]
        
    d = [pd.to_datetime(str(t_)).strftime(fmt) for t_ in t]
    
    return d
    
tests()

def seconds_since(timein, basetime=datetime.datetime(1990,1,1)):
    """
    Converts a list or array of datetime object into an array of seconds since "basetime"
    
    From MRayson SFODA
    """
    
    # Convert the time to seconds

    # Workout the input type
    isarray = False
    isdatetime = False
    isdatetime64 = False
    if isinstance(timein, list) or isinstance(timein, np.ndarray):
        isarray = True
        if isinstance(timein[0], datetime.datetime):
             isdatetime = True
        elif isinstance(timein[0], np.datetime64):
             isdatetime64 = True

    # Could use the funcitons above but this seems fine for now. 
    if isarray and isdatetime64:
        time0 = np.datetime64(basetime).astype('<M8[us]')
        tsec = ((timein.astype('<M8[us]') - time0)).astype(np.float64)*1e-6

    elif isarray and isdatetime:
        tsec = np.array([(t-basetime).total_seconds() for t in timein])

    else: # Assume its a datetime.datetime object (this is not good..)
        if isinstance(timein, datetime.datetime):
            tsec = (timein-basetime).total_seconds()

        else:
            time0 = np.datetime64(basetime).astype('<M8[us]')
            try:
                #raise Exception, 'Time conversion error'
                # This case shouldn't arise
                tsec = ((timein.astype('<M8[us]') - time0)).astype(np.float64)*1e-6
            except:
                tsec = (timein-time0).total_seconds() # TimeDelta object

    return tsec

def time_interp1(time_out, time_in, data_in):
    """
    np.interp for timeseries data. 
    
    """
    
    # Convert times to an integer
    time_out_sec = seconds_since(time_out)
    time_in_sec = seconds_since(time_in)
    
    return np.interp(time_out_sec, time_in_sec, data_in)

#%%