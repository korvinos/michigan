# Create new geolocation for MODIS-a image
sfile = ifile # <ifile> is path to sentinel image
mfile = '/nfs0/data_ocolor/michigan/MODISa/A2016248183000.L2_LAC_OC.nc'
# Generate path to domain file 
domainfile = os.path.split(sfile)[1] +  '_domain.nc'
# get domain file by nansat
n2 = Nansat(domainfile)

ofile = os.path.split(mfile)[1] + '_pro.nc'
# GCP_COUNT is frequncy of rows used for gettin coordinates array ? 
# Get MODIS-a image by nansat
n = Nansat(mfile, GCP_COUNT=40)
# Remove geolocation
n.vrt.remove_geolocationArray()
n.vrt.tps = True
n.reproject_GCPs()
print n.time_coverage_start

# add index of pixels
index = np.arange(0, n.shape()[0] * n.shape()[1]).reshape(n.shape()).astype('int32')
n.add_band(index, parameters={'name': 'index'})

n.reproject(n2, addmask=False)
#bands = ['Rrs_412', 'Rrs_443', 'Rrs_488', 'Rrs_531', 'Rrs_547', 'Rrs_555', 'Rrs_645', 'Rrs_667', 'Rrs_678', 'angstrom', 'aot_869', 'chlor_a']
#bands = ['Rrs_667', 'Rrs_555', 'Rrs_488', 'Rrs_443', 'Rrs_412']
bands = ['index', 'Rrs_412', 'Rrs_443', 'Rrs_488', 'Rrs_531','Rrs_555', 'Rrs_645', 'Rrs_667', 'Rrs_678']

nexp = Nansat(domain=n2)
for band in bands:
    print band,
    bandArray = n[band]
    print bandArray.dtype
    nexp.add_band(bandArray, parameters={'name': band})

nexp.export(ofile)