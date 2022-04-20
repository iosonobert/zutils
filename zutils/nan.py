import numpy as np

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def fillnans(y, first_gone=None):
    """
    Fill nans down the columns of y. To do it the other way transpose the input.
    """

    if not type(None) == type(first_gone):
        y[first_gone::, :] = np.nan
    
    s = y.shape
    if len(s) == 1:
        flat = True
        y = y[:, None]
    elif len(s) == 2:
        flat = False
    else:
        raise(Exception)
    
    s = y.shape
    nn, mm = s
    
    for m in np.arange(0, mm):
                
        nans, x = nan_helper(y[:, m])

        if sum(nans)>0 and sum(nans)<len(nans):
            y[nans, m]= np.interp(x(nans), x(~nans), y[~nans, m])

    if flat:
        y = y.flatten()

    return y
