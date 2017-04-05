from nansat import Nansat, Domain
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

def geolocation(mfile, domain, final_path):
    print mfile
    n = Nansat(mfile, GCP_COUNT=40)
    # Remove geolocation
    n.vrt.remove_geolocationArray()
    n.vrt.tps = True
    n.reproject_GCPs()
    print n.time_coverage_start

    # add index of pixels
    index = np.arange(0, n.shape()[0] * n.shape()[1]).reshape(n.shape()).astype('int32')
    n.add_band(index, parameters={'name': 'index'})

    n.reproject(domain, addmask=False)
    # bands = ['Rrs_412', 'Rrs_443', 'Rrs_488', 'Rrs_531', 'Rrs_547', 'Rrs_555', 'Rrs_645', 'Rrs_667',
    #          'Rrs_678', 'angstrom', 'aot_869', 'chlor_a']
    # bands = ['Rrs_667', 'Rrs_555', 'Rrs_488', 'Rrs_443', 'Rrs_412']
    bands = ['index', 'Rrs_412', 'Rrs_443', 'Rrs_488', 'Rrs_531', 'Rrs_555', 'Rrs_645', 'Rrs_667', 'Rrs_678']

    nexp = Nansat(domain=domain)
    for band in bands:
        print band,
        bandArray = n[band]
        print bandArray.dtype
        nexp.add_band(bandArray, parameters={'name': band})

    nexp.export(os.path.join(final_path, os.path.split(mfile)[1] + postfix))

geolocation = np.vectorize(geolocation)
