from nansat import Nansat, Domain
import numpy as np
import os
import re
import glob


class Data:
    # Different bands sets for each sensors
    wavelengths = {
        'modis': {
            'full': [412, 443, 488, 531, 547, 555, 645, 667, 678],
            '1x1km_bands': [412, 443, 488, 531, 645, 678]
        },

        'sentinel2': [443, 490, 560, 665, 705, 740, 783, 842, 945, 1375, 1610, 2190]
    }

    # Sandy Bear Dunes domain
    # <pixel_size> spatial resolution of each pixel in meters.
    # Based on characteristics of sentinel2a. Min: 10 m; Mean: 20 m; Max: 60 m;
    # To get more about sentinel2a bands see:
    # https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/resolutions/radiometric
    pixel_size = 60

    # <x_resol> and <y_resol> are resolution of image in kilometers converted to chused resolution
    # of sentinel2a image (see upper)
    x_resolution, y_resolution = 122 * (1000 / pixel_size), 78 * (1000 / pixel_size)

    # <domain> is domain which covers Sandy Bear Dunes Region
    # To see region on a map see: https://github.com/korvinos/michigan/blob/master/michigan.geojson
    # To get characteristics of region see michigan.geojson
    domain = Domain('+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs', '-lle -86.3 44.6 -85.2 45.3 -ts %s %s'
                    % (x_resolution, y_resolution))

    BANDS = {
        '01': {'wavelength': 443, 'resolution': 60},
        '02': {'wavelength': 490, 'resolution': 10},
        '03': {'wavelength': 560, 'resolution': 10},
        '04': {'wavelength': 665, 'resolution': 10},
        '05': {'wavelength': 705, 'resolution': 20},
        '06': {'wavelength': 740, 'resolution': 20},
        '07': {'wavelength': 783, 'resolution': 20},
        '08': {'wavelength': 842, 'resolution': 10},
        #    '8A':{ 'wavelength': 865, 'resolution': 20 },
        '09': {'wavelength': 945, 'resolution': 60},
        '10': {'wavelength': 1375, 'resolution': 60},
        '11': {'wavelength': 1610, 'resolution': 20},
        '12': {'wavelength': 2190, 'resolution': 20}
    }

    # <granules> is list of granules which covers Sandy Bear Dunes region.
    # To get more information about granules see:
    # https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/product-types
    granules = ['16TER', '16TFR', '16TEQ', '16TFQ']

    def __init__(self, ifile, domain):
        """
        :param ifile: str, file path  
        :param domain: <nansat.domain.Domain> object
        """

        self.domain = domain
        self.m_file = None
        self.s_file = None

        if re.match(r'A', ifile) is None:
            self.s_file = ifile
        else:
            self.m_file = ifile

    def modis_geo_location(self, wavelengths, save_path='./', gcp_count=70):
        """
        :param wavelengths: list, list of wavelengths 
        :param save_path: str
        :param gcp_count: str
        :return: <nansat.nansat.Nansat> object, an object with a new geo location 
        """
        if self.m_file is None:
            raise IOError

        print self.m_file
        n = Nansat(self.m_file, GCP_COUNT=gcp_count)
        # Remove geo location
        n.vrt.remove_geolocationArray()
        n.vrt.tps = True
        n.reproject_GCPs()
        print n.time_coverage_start

        # add index of pixels
        index = np.arange(0, n.shape()[0] * n.shape()[1]).reshape(n.shape()).astype('int32')
        n.add_band(index, parameters={'name': 'index'})
        n.reproject(self.domain, addmask=False)

        bands = ['Rrs_%s' % band for band in wavelengths]
        bands.insert(0, 'index')

        n_export = Nansat(domain=self.domain)
        for band in bands:
            print band
            band_arr = n[band]
            n_export.add_band(band_arr, parameters={'name': band})
            n_export.export(os.path.join(save_path, os.path.split(self.m_file)[1] + '_reprojected.nc'))

        return n_export

    def stich(self, domain, bands, gdirs):
        for band in bands:
            print('Band:', band)
            # Generate 2d array which will have shape like <d> and filled by nan
            bandArray = np.zeros(domain.shape()) + np.nan
            for gdir in gdirs:
                bfile = glob.glob(os.path.join(gdir, 'IMG_DATA', '*_B%s.jp2' % band))[0]
                n = Nansat(bfile)
                # Reprojection of data according to domain; eResampleAlg 1 is Bilinear
                n.reproject(d, eResampleAlg=1, addmask=False)
                # Get data array from nansat object
                bdata = n[1]
                bandArray[bdata > 0] = bdata[bdata > 0]

            n_obj.add_band(bandArray, parameters={'name': 'Rrs_%s' % BANDS[band]['wavelength']})

        return n_obj

    def s2_downscale(self):

        if re.match(r'S2A', self.s_file):
            raise IOError

        gdirs = []
        for granule in self.granules:
            gdirs += sorted(glob.glob(os.path.join(self.s_file, 'GRANULE', '*_T%s_*' % granule)))

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
        print('Domain created')

        # Create base nansat object which domain covers all four <granules>
        bands = self.wavelengths['sentinel2']
        n_obj = self.stich(d, bands, gdirs)
        n_obj.reproject(self.domain)
        export_name = os.path.split(self.s_file)[1] + '_full.nc'
        n_obj.export("output/%s" % export_name)
