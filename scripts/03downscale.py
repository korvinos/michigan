# Create image from granules
def stich(band):
    global gdirs
    # Generate 2d array which will have shape like #d and filled by nan
    bandArray = np.zeros(d.shape()) + np.nan
    for gdir in gdirs:
        bfile = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B%02d.jp2' % band))[0]
        n = Nansat(bfile)
        # Reprojection of data according to domain; eResampleAlg 1 is Bilinear
        n.reproject(d, eResampleAlg=1, addmask=False)
        # Get data array from nansat object 
        bdata = n[1]
        # 
        bandArray[bdata > 0] = bdata[bdata > 0]

    return bandArray

# north sea
gdirs = []

# Four granules which covering Sandy Bear Dunes Region
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
n = Nansat(domain=d, array=np.zeros(d.shape(), np.int8))

# Save created file on a disk. So save a file which contain a domain
n.export(os.path.split(ifile)[1] + '_domain.nc')

# list of requested bands
# Information about bands: https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/resolutions/radiometric
bands = range(1, 12)

# Number of threadds
# Documentation: https://docs.python.org/2/library/multiprocessing.html#module-multiprocessing.pool
p = Pool(7)

# For each granule get reprojected array on each band from bands
bData = p.map(stich, bands)
# Transformate map.object to numpy array
bData = np.array(bData)

p.terminate()
ofile = os.path.split(ifile)[1] + '.npz'
np.savez_compressed(ofile, bData)