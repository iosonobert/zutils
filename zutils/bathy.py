import pandas as pd
import matplotlib.pyplot as plt 
import shapefile


def xyzbathy_2_contourshp(infiles, outfile, levels, field_name='DEPTH', ):
    """
    Function to convert xyz data into shapefiles of contours. Nominally it is for bathymetry, 
    a different field_name can be used for other scalars. 

    Inputs:
        infiles - xyx csv [or iterable of] with headers X, Y, Z [in capitals]
        outfile - name of the shapefile to be produced
        levels - the levels the contours are to be specified at.
        field_name - the field name that the scalar field (Z, column) will be written 
                to in the shp file. 
    """

    if type(infiles) == str:
        infiles = [infiles]
    
    CSs = []
    
    for infile in infiles:
        
        print('Contouring {}'.format(infile))
        
        df = pd.read_csv(infile)

        CS = plt.tricontour(df.X, df.Y, df.Z, 15, levels=levels, linewidths=0.5, colors='k')
        CSs += [CS]
        
    print('Done contouring.')
    plt.show()

    if False:
        for CS in CSs:
            for ii, segs in enumerate(CS.allsegs):

                print('{} m: {}'.format(levels[ii], len(segs)))

                for seg in segs:

                    plt.plot(seg[:, 0], seg[:, 1], 'k:', linewidth=0.5)

    with shapefile.Writer(outfile, shapeType=3) as w:
    
        w.field(field_name, 'N') # N = numbers
        
        for CS in CSs:
            for ii, segs in enumerate(CS.allsegs):
                print('{} m: {}'.format(levels[ii], len(segs)))

                for seg in segs:
                    plt.plot(seg[:, 0], seg[:, 1], 'k:', linewidth=0.5) 
                    w.line([seg])
                    w.record(DEPTH=levels[ii])
                    
    return CSs