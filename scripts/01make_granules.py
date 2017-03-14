# Get list of granules
# What is granules: https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/product-types
# Add path to sentinel image <ifile>

gdirs = sorted(glob.glob(os.path.join(ifile, 'GRANULE', '*')))

lons = [] # List of longitudes for each granule 
lonVec = [] # ???
lats = [] # List of latitudes for each granule
latVec = [] # ???
labels = [] # List on names for each granule
# Part I: Generation of domain
# For each granule in granules list / gdirs
for gdir in gdirs:
    # Get B01.jp2 image (granule) path 
    # About jp2: https://ru.wikipedia.org/wiki/JPEG_2000
    b01file = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B01.jp2'))[0]
    
    # Get granule image by nansat
    n = Nansat(b01file)
    # Get max/min lat and long of granule 
    lon, lat = n.get_border()
    ## print lon, lat
    lons += list(lon)
    lats += list(lat)
    
    # For what this two rows ??? 
    lonVec.append(lon)
    latVec.append(lat)
    # Add name/label of granule to list
    labels.append(b01file.split('_')[-2][1:])

# Create Domain which is based on max and min lat and long from all (!!!) granules
# 1000 x 1000 is resolution of MODIS-a img ? 
d = Domain(n.vrt.get_projection(), '-lle %f %f %f %f -tr 1000 1000' % (
    min(lons), min(lats), max(lons), max(lats)))

## print(d)

# Create name for export of img which all granules
mapfile = os.path.split(ifile)[1] + '_map.png'
## print(mapfile)
d.write_map(mapfile, resolution='h',
            lonVec=lonVec, latVec=latVec,
            lonBorder=1, latBorder=1,
            dpi=150,
            labels=labels)

# Part II: Reprojection of granules
# Generate 2d array which will have shape like #d annd fill by nan
s2array = np.zeros(d.shape()) + np.nan 
# For each granule in granules list / gdirs
for gdir in gdirs:
    # Get B01.jp2 image (granule) path     
    b01file = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B01.jp2'))[0]
    ## print b01file
    n = Nansat(b01file)
    # Reprojection of granule according #d, 0 is nearest neighbour
    n.reproject(d, eResampleAlg=0)
    # Get data array from nansat obj
    b1 = n[1]
    # All areas in #s2array array will filled by vales from granule
    # As result, after loop we will get one object which contain data from all granules 
    s2array[b1 > 0] = b1[b1 > 0]

# Create nansat obj according #d which contain watter mask
wm = Nansat(domain=d).watermask()[1]
# Why wm == 2 ?? From documentation: Create numpy array with watermask (water=1, land=0)
s2array[wm == 2] = 0

# Create name of file
qlfile = os.path.split(ifile)[1] + '_ql.png'
# Save on disk full generated image: #s2array
plt.imsave(qlfile, s2array, vmin=1050, vmax=1400)