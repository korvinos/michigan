import os
import glob
from nansat import *
from scipy.ndimage.filters import gaussian_filter, median_filter
from matplotlib.colors import LogNorm
from scipy.ndimage.interpolation import rotate
from multiprocessing import Pool
from matplotlib.patches import Polygon

ifile = '/nfs0/data_ocolor/michigan/sentinel2a/S2A_OPER_PRD_MSIL1C_PDMC_20160904T205606_R126_V20160903T164322_20160903T164911.SAFE'

def stich(band):
    global gdirs, n_obj
    # Generate 2d array which will have shape like #d and filled by nan
    bandArray = np.zeros(d.shape()) + np.nan
    for gdir in gdirs:
        bfile = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B%02d.jp2' % band))[0]
        n = Nansat(bfile)
        # Reprojection of data according to domain; eResampleAlg 1 is Bilinear
        n.reproject(d, eResampleAlg=1, addmask=False)
        # Get data array from nansat object 
        bdata = n[1] 
        bandArray[bdata > 0] = bdata[bdata > 0]
    
    n_obj.add_band(bandArray, parameters={'name': 'Rrs_%s' % band})

def create_domain():
    pass

# north sea
gdirs = []
granules = ['16TER', '16TFR', '16TEQ', '16TFQ']

# get directories of relevant granules
for granule in granules:
    gdirs += sorted(glob.glob(os.path.join(ifile, 'GRANULE', '*_T%s_*' % granule)))

# get lon/lat limits
# Lists for accumulation of lon/lat values from each granule
lons = []
lats = []
for gdir in gdirs:
    # Find absolute path to relevant granule
    b0file = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B01.jp2'))[0]
    # Get granule by nansat 
    n = Nansat(b0file)
    lon, lat = n.get_corners()
    # Add min/max values of long and lat to list
    lons += list(lon)
    lats += list(lat)

# Create domain according to max and min values of lon and lat
d = Domain(n.vrt.get_projection(), '-lle %f %f %f %f -tr 60 60' % (
            min(lons), min(lats), max(lons), max(lats)))

# Create base nansat obj with nan array on a board 
n_obj = Nansat(domain=d)
bands = np.arange(1, 11)
stich = np.vectorize(stich)
stich(bands)

n_obj.export('output/S2A_OPER_PRD_MSIL1C_PDMC_20160904T205606_R126_V20160903T164322_20160903T164911.SAFE.nc')