import numpy as np

def jft(s, d, s_bins=5, d_bins=8, s_max=None, clean_nans = True, outfile=None):
    """
    Joint frequency table in %. Specify an outfile to save to csv.

    Inputs:
        s - speed array
        d - direction array
        s_bins [optional, default=5] number of speed bins
        d_bins [optional, default=5] number of direction bins
        s_max [optional] maximum speed to use
        clean_nans [optional, default=True - If False and there are nans in the input arrays the output won't add up to 100%. 
        outfile [optional] output file name.
    """

    if clean_nans:
        bad = np.isnan(s)
        s = s[~bad]
        d = d[~bad]
        bad = np.isnan(d)
        s = s[~bad]
        d = d[~bad]

    n = len(s)
    if s_max is None:
        s_max = np.ceil(np.nanmax(s))
    
    table = np.histogram2d(s, d, bins = [s_bins, d_bins], range = [[0, s_max], [0, 360]])
    
    ss = table[1]
    ds = table[2]
    jft_table = 100*table[0]/n
    jft_list = jft_table.tolist()
    
    top = ['{} - {}'.format(ds[i], ds[i+1]) for i in np.arange(0, len(ds)-1)]
    left = ['{} - {}'.format(ss[i], ss[i+1]) for i in np.arange(0, len(ss)-1)]

    import csv

    if not outfile is None:
        
        with open(outfile,"w", newline='') as file:

            out = csv.writer(file, delimiter=',')
            out.writerow(['spd\dir %']+top)

            for row in np.arange(0, len(left)):

                out.writerow(left[row:row+1]+jft_list[row])

    return table