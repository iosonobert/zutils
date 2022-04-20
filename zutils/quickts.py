

import numpy as np, scipy.signal as signal
import zutils.time as ztime


def quick_butter(time, data, T_cut_seconds, order=4, btype='lowpass'):
    """
    Very quick butterworth filter that works with python time objects. 
    
    args:
        time: numpy array of time objects that can be converted by zutils.time.seconds_since and zutils.time.is_well_spaced
        data: np.array to be filtered. Can be higher dimensional but goes into signal.butter as is.
        T_cut_seconds: cutoff period in seconds
        
    kwargs: 
        btype: 'lowpass' or 'highpass'
        order: order of the butterworth

    """
    
    # Convert to seconds
    t_sec = ztime.seconds_since(time)

    dt = t_sec[1] - t_sec[0]
    f_s_Hz = 1/dt # my sampling freq
    f_N_Hz = 1/dt # Nyqvist frequency in function of my sampling frequency

    # Validation [should be fine for mopdel data]
    ztime.is_well_spaced(time, f_s_Hz)

    # Filter design
    f_cut_Hz = 1/T_cut_seconds

    f_cut_norm = f_cut_Hz/f_N_Hz  # normalized cut_off frequency

    b, a= signal.butter(order, f_cut_norm, btype)

    return signal.filtfilt(b, a, data)

def quick_interp(time, data, **kwargs):
    
    dt_sec = kwargs.pop('dt_sec', 60)
    start = kwargs.pop('start', None)
    end = kwargs.pop('end', None)
    
    dt_file = time[1] - time[0]
    print('File dt initial: {}'.format(dt_file))
    print('dt out: {}'.format(dt_sec))
    
    ref_time = time[0]

    time_in_secs = ztime.seconds_since(time, ref_time)
    
    if start is None:
        start_sec = time_in_secs[0]
    else:
        start = np.datetime64(start)
        start_sec = ztime.seconds_since(start, ref_time)
        
    if end is None:
        end_sec = time_in_secs[-1]
    else:
        end = np.datetime64(end)
        end_sec = ztime.seconds_since(end, ref_time)

    time_out_secs = np.arange(start_sec, end_sec, dt_sec)
    time_out_np64 = np.array([ref_time + np.timedelta64(int(i), 's') for i in time_out_secs])

    if len(data.shape) == 1:
        var_interp = np.interp(time_out_secs, time_in_secs, data) 
    else:
        var_interp = np.array([np.interp(time_out_secs, time_in_secs, var_row) for var_row in data])

    return time_out_secs, time_out_np64, var_interp