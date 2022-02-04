
import sys
import os
import fileinput
import numpy as np 

from netCDF4 import Dataset
import glob

def append_and_print(dump_list, str, prnt):
    """
    Utility function that prints strings and adds them to a list. 
    """
    dump_list += [str]
    if prnt:
        print(str)
    
    return dump_list

def check_first_4_bytes(fn):
    """
    Checks first 4 bytes of a file. This can be used to identify the file type. 
    Files with user data before the header may affect this.

    On distinguishing files:

        ```
        "The short answer is that under most circumstances, you should not care, if you use version 4.0 or later of the netCDF library to access data in the file. But the difference is indicated in the first four bytes of the file, which are 'C', 'D', 'F', '\001' for the classic netCDF CDF-1 format; 'C', 'D', 'F', '\002' for the 64-bit offset CDF-2 format; 'C', 'D', 'F', '\005' for the 64-bit data CDF-5 format; or '\211', 'H', 'D', 'F' for an HDF5 file, which could be either a netCDF-4 file or a netCDF-4 classic model file. (HDF5 files may also begin with a user-block of 512, 1024, 2048, ... bytes before what is actually an 8-byte signature beginning with the 4 bytes above.)

        With netCDF version 4.0 or later, there is an easy way that will distinguish between netCDF-4 and netCDF-4 classic model files, using the "-k" option to ncdump to determine the kind of file, for example:"
        ```


    from anoher site:

        ```
        "We refer "CDF-1" as the identification string, "magic", occupying the first 4 bytes of a netCDF file. The string can be "CDF1", "CDF2", or "CDF5".

        CDF-1 and CDF-2 are also referred by the ESDS Community Standard as NetCDF Classic and 64-bit Offset File Formats, respectively. See [ESDS-RFC-011v2.0]

        The difference between CDF-1 and CDF-2 is only in the VERSION byte (\x01 vs. \x02) and the OFFSET entity, a 64-bit instead of a 32-bit offset from the beginning of the file. See CDF-2 file format specification for the detailed specifications of both CDF-1 and CDF-2.

        Below is an older version of CDF file format specification"
        ```

    """
    
    print('Tyr first 4 bytes of "' + fn + '" are:')
    with open(fn, "rb") as file:

        bytel = [file.read(1) for i in np.arange(0, 4)]
        print(bytel)
        
        if bytel == [b'C', b'D', b'F', b'\x01']:
            print('This looks like a CDF-1 file [NETCDF3_CLASSIC]')
            format = 'NETCDF3_CLASSIC'
        elif bytel == [b'C', b'D', b'F', b'\x02']:
            print('This looks like a CDF-2 file [NETCDF3_64BIT]')
            format = 'NETCDF3_64BIT'
        elif bytel[1::] ==  [b'H', b'D', b'F']:
            print('This looks like a NetCDF4 file of some type, guessing NETCDF4 but could be NETCDF4_CLASSIC')
            format = 'NETCDF4'
        else:
            print('This byte structure is unknown, might be a file type I''m not aware of. It might not be a netcdf file, or it might have a user data string.')
            format = None   
    
    print('    ')
    return format

    
def nc_compare(lhs_file, rhs_file):
    """
    Compare 2 netcdfs. Note, when readint the file contents, both are opened as 'NETCDF4' format. 
    Not sure this is the best idea.
    """
    
    print('\\ COMPARING FILETYPES:')
    format_lhs = check_first_4_bytes(lhs_file)
    format_rhs = check_first_4_bytes(rhs_file)
    
    print('\n\n\\ COMPARING REMAINING FILE DIFFERENCES (only differences are shown):')

    lhs_dump_list = ncdump(lhs_file, verbose=False, format=format_lhs)
    rhs_dump_list = ncdump(rhs_file, verbose=False, format=format_rhs)
    
    ldifferences = 0
    rdifferences = 0
    for entry in lhs_dump_list:
        if entry in rhs_dump_list:
            pass
        else:
            ldifferences += 1
            print('"{}" \n\t\t is not in the right hand file [{}]!\n'.format(entry, rhs_file))

    print(' ')

    for entry in rhs_dump_list:
        if entry in lhs_dump_list:
            pass
        else:
            rdifferences += 1
            print('"{}" \n\t\t is not in the left hand file [{}]!\n'.format(entry, lhs_file))

    print('\\ SUMMARY:')
    if format_lhs == format_rhs:
        print('  \\ File types are the same!:')
    else:
        print('  \\ WARNING: FILE TYPES NOT THE SAME!:')
        
    print('  \\ There were {} items in the RHS file not seen in the LHS file (more details above)!:'.format(rdifferences))
    print('  \\ There were {} items in the LHS file not seen in the RHS file (more details above)!:'.format(rdifferences))
    
def ncdump(filename, verbose=True, format='NETCDF4'):
    """
    Do an ncdump to screen and return all a list of strings, where each string is a line in the dump. Also checks first 4 bytes for file type. 

    Inputs:
        - filename: name of the netcdf file
        - verb: Bool controlling whether to print to screen. If false just return the list. 
        - format: format of the netcdf file: THE FUNCTION SHOULD DETERMINE THIS.

    Outputs:
        - dump_list: a list of strings, where each string is one line in the dump.  
    """

    check_first_4_bytes(filename)

    DS = Dataset(filename, mode='r', format=format)

    dump_list = ncdump_from_DS(DS, verbose=verbose)
    
    return dump_list

def ncdump_from_DS(DS, verbose=True):
    """
    Do an ncdump to screen and return all a list of strings, where each string is a line in the dump. 
    """
    
    dump_list = []
    
    dump_list = append_and_print(dump_list, 'dimensions:', verbose)

    for dimension_name in DS.dimensions.keys():
        dim = DS.dimensions[dimension_name]
        
        if dim.isunlimited():
            dump_list = append_and_print(dump_list, '\t{} = {} [UNLIMITED]'.format(dimension_name, dim.size), verbose)
        else:
            dump_list = append_and_print(dump_list, '\t{} = {}'.format(dimension_name, dim.size), verbose)

    dump_list = append_and_print(dump_list, 'variables:', verbose)

    for variable_name in DS.variables.keys():
    #     print('        {}     :   {}'.format(attr, DS.getncattr(attr)))
        var = DS.variables[variable_name]
        dimstr = ''
        for dimension in var.dimensions:
            dimstr += dimension+', '

        dimstr = dimstr[0:-2]

        dump_list = append_and_print(dump_list, '\t{} {}({})'.format(var.dtype, variable_name, dimstr), verbose)

        for attr in var.ncattrs():
            if type(var.getncattr(attr)) == str:
                dump_list = append_and_print(dump_list, '\t\t{}:{} = "{}" ;'.format(variable_name, attr, var.getncattr(attr)), verbose)
            else:
                dump_list = append_and_print(dump_list, '\t\t{}:{} = {} ;'.format(variable_name, attr, var.getncattr(attr)), verbose)


    dump_list = append_and_print(dump_list,'// global attributes:', verbose)

    for attr in DS.ncattrs():
    #     print('        {}     :   {}'.format(attr, ifile.getncattr(attr)))
        if type(DS.getncattr(attr)) == str:
            dump_list = append_and_print(dump_list, '\t\t:{} = "{}" ;'.format(attr, DS.getncattr(attr)), verbose)
        else:
            dump_list = append_and_print(dump_list, '\t\t:{} = {} ;'.format(attr, DS.getncattr(attr)), verbose)
    
    return dump_list 

def ncdump_old(filename, verbose=True, format='NETCDF4'):
    
    nc_fid = Dataset(infile, mode='r', format=format)

    ncdump_from_Dataset_old(nc_fid, verbose=True)

def ncdump_from_Dataset_old(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print ("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars



def nccopy(DS, copy_name):
    """
    The idea here is to be able to directly copy a netcdf file. It's under way - not fully implemented. 
    """

    raise(Exception('Function not ready yet.'))

    with Dataset(copy_name,'w',format='NETCDF4') as new_ncfile :

        dump_list = []

        dump_list = append_and_print(dump_list, 'dimensions:')

        for dimension_name in DS.dimensions.keys():
            dim = DS.dimensions[dimension_name]

            if dim.isunlimited():
                dump_list = append_and_print(dump_list, '\t{} = {} [UNLIMITED]'.format(dimension_name, dim.size))
            else:
                dump_list = append_and_print(dump_list, '\t{} = {}'.format(dimension_name, dim.size))

        dump_list = append_and_print(dump_list, 'variables:')

        for variable_name in DS.variables.keys():
        #     print('        {}     :   {}'.format(attr, DS.getncattr(attr)))
            var = DS.variables[variable_name]
            dimstr = ''
            for dimension in var.dimensions:
                dimstr += dimension+', '

            dimstr = dimstr[0:-2]

            dump_list = append_and_print(dump_list, '\t{} {}({})'.format(var.dtype, variable_name, dimstr))
            new_ncfile.createVariable(variable_name, var.dtype)
            
            for attr in var.ncattrs():
                if type(var.getncattr(attr)) == str:
                    dump_list = append_and_print(dump_list, '\t\t{}:{} = "{}" ;'.format(variable_name, attr, var.getncattr(attr)))
                else:
                    dump_list = append_and_print(dump_list, '\t\t{}:{} = {} ;'.format(variable_name, attr, var.getncattr(attr)))


        print('// global attributes:')

        for attr in DS.ncattrs():
        #     print('        {}     :   {}'.format(attr, ifile.getncattr(attr)))
            if type(DS.getncattr(attr)) == str:
                dump_list = append_and_print(dump_list, '\t\t:{} = "{}" ;'.format(attr, DS.getncattr(attr)))
            else:
                dump_list = append_and_print(dump_list, '\t\t:{} = {} ;'.format(attr, DS.getncattr(attr)))

        return dump_list 