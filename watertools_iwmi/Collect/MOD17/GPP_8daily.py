import sys
from watertools_iwmi.Collect.MOD17.DataAccessGPP import DownloadData


def main(Dir, Startdate, Enddate, latlim, lonlim, cores=False, Waitbar = 1, hdf_library = None, remove_hdf = 1):
    """
    This function downloads MOD17 8-daily GPP data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    cores -- amount of cores used
    Waitbar -- 1 (Default) will print a waitbar
    """

    # print('\nDownload 8-daily MODIS GPP data for period %s till %s' %(Startdate, Enddate))
    output = DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar, cores, hdf_library, remove_hdf)
    
    return output