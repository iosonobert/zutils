# zutils
Simple utilities. Mostly things I use to interface with existing python packages. As such it unashamedly uses has many dependencies upon these packages e.g. 
- datetime 
- matplotlib
- numpy 
- pandas 
- xarray 
- jupyter [for examples]

Many outdated items I wrote years ago when I was first transitioning from Matlab to help me do things I could do easily in Matlab, but struggled with in python. I didn't really know what I was doing at that time, and the underlying modules have developed since. 

Many modules are slowly being split off from the turbotools package [one of my main packages for ocean field data analysis] as I change my mind as to where is the most sensible places for these to sit [I'd like those other packages to be self contained where possible, but I wan't to limit duplication also]. 

Going to actually try document this now, but it is expected that the documentation will be poor, incomplete, and will be rapidly outdated. 

There are some example notebooks in the Examples folder, and I think these will be preferred to detailed documentation.

## Modules:

### file
Utils for working with filenames or reading and writing simple data structures to file.
- drop_extension
- change_extension
- write_dict_as_str. Writes a dict as a text file. 
- read_dict_from_str. Reads text file produced by write_dict_as_str

### plotting
Basically an interface for matplotlib. Classes and methods:
- axis_layer: a class to lay out axes exactly such that journal ready figures can be produced. Much more control than I've been able to attain with subplot2grid etc. 
- hide_xticks
- hide_yticks

[Click here](https://github.com/iosonobert/zutils/blob/master/Examples/plotting.ipynb) for examples

### pyMOS
Utils for the pyMOS [python IMOS] toolbox prepped for archiving KISSME and RS19 data.

### qc_conventions
Trying to devise a class for QAQC conventions used to archive KISSME and RS19 datasets.
Base class for QAQC conventions. 
- base: Base class for QAQC conventions.

### stats
I don't like the name of this module
- jft

### time 
Mostly tools for datetime.datetime objects which are my preferred oobject. Many tools for transferring between the myriad of python datetime classes, some for allowing me
to be more Matlab-like in my interaction with datetime objects. I should have actually made a matlabdatetime class but I didn't really know what I was doing when I first wrote this stuff. 
- linspacetime: linearly space a vector of datetime objects.
- is_well_spaced: A function to determine whether a time array is well spaced at fs.
- datetime_2_timestamp
- timestamp_2_datetime
- datetime64_2_timestamp
- timestamp_2_datetime64
- datetime64_2_datetime
- datetime64_2_datetime_single
- datetime_2_datetime64
- datetime_2_datetime64_single
- str_to_dt: Matlab hack 
- matformat_2_py: Simple hack code to allow me to use matlab date formats pn python. At the moment I 
    just use common formats but should ultimately switch to a regext type approach to 
    find and replace. Not going for that as I should just learn the python formats. 
- matplotlibtime_2_matlabdatenum [NOT IMPLEMENTED]
- matlabdatenum_2_matplotlibtime [NOT IMPLEMENTED]
- datenum_2_matlabdatenum [NOT IMPLEMENTED]
- matlabdatenum_2_datenum [NOT IMPLEMENTED]
- np_datetime64_strftime: Format a np.datetime64 or list of np.datetime64 objects as a string.

### xr
Utils for working with xarray. 
- select_calendar_month: Return a subset of the Dataset or DataArray limited to the specified year and month.
- split_by_timegap: Split a DataArray along the time dimension based on a set time gap.

### xrwrap
An actual wrapper for xarray datasets. Should be pulled into the xr module. 
- xrwrap: An actual wrapper for xarray datasets. Should be pulled into the xr module. Does not inherit xarray, but takes the dataset as a property. This is another class with mant methods that controls interaction with the DataSet, but many operations still require pulling the DateSet out. 

