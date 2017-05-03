from nansat import Nansat, Domain
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import glob


class Data:
    # Different bands sets for each sensors
    wavelengths = {
        'modis': {
            'full': [412, 443, 488, 531, 547, 555, 645, 667, 678],  # All MODISa channels
            '1x1km_bands': [412, 443, 488, 531, 645, 678],  # Only 1x1 km spatial resolution bands
            'blue_off': [443, 469, 488, 531, 547, 555, 645, 667, 678],  # Without 412 nm
            'red_off': [412, 443, 469, 488, 531, 547, 555, 645, 667],  # Without 678 nm
            'blue_and_red_off': [443, 469, 488, 531, 547, 555, 645, 667],  # Without 412 nm and 678 nm
            'red_off_full': [412, 443, 469, 488, 531, 547, 555],
            '1x1km_bands_2': [443, 488, 531, 645, 678]  # Only 1x1 km spatial resolution bands
        },

        'sentinel2': {
            '01': 443,
            '02': 490,
            '03': 560,
            '04': 665,
            '05': 705,
            '06': 740,
            '07': 783,
            '08': 842,
            '8A': 865,
            '09': 945,
            '10': 1375,
            '11': 1610,
            '12': 2190
        }
    }

    #         '01': {'wavelength': 443, 'resolution': 60},
    #         '02': {'wavelength': 490, 'resolution': 10},
    #         '03': {'wavelength': 560, 'resolution': 10},
    #         '04': {'wavelength': 665, 'resolution': 10},
    #         '05': {'wavelength': 705, 'resolution': 20},
    #         '06': {'wavelength': 740, 'resolution': 20},
    #         '07': {'wavelength': 783, 'resolution': 20},
    #         '08': {'wavelength': 842, 'resolution': 10},
    # #        '8A':{ 'wavelength': 865, 'resolution': 20},
    #         '09': {'wavelength': 945, 'resolution': 60},
    #         '10': {'wavelength': 1375, 'resolution': 60},
    #         '11': {'wavelength': 1610, 'resolution': 20},
    #         '12': {'wavelength': 2190, 'resolution': 20}
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
    sbd_dom = Domain('+proj=latlong +datum=WGS84 +ellps=WGS84 +no_defs', '-lle -86.3 44.6 -85.2 45.3 -ts %s %s'
                    % (x_resolution, y_resolution))

    # <granules> is list of granules which covers Sandy Bear Dunes region.
    # To get more information about granules see:
    # https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/product-types
    granules = ['16TER', '16TFR', '16TEQ', '16TFQ']

    def __init__(self, ifile, domain=None):
        """
        :param ifile: str, file path  
        :param domain: <nansat.domain.Domain> object
        """

        self.domain = self.sbd_dom

        if domain is not None:
            self.domain = domain

        self.ifile = ifile

    def modis_geo_location(self, wavelengths_set='1x1km_bands', save_path='./', gcp_count=40):
        """
        :param wavelengths: list, list of wavelengths 
        :param save_path: str
        :param gcp_count: str
        :return: <nansat.nansat.Nansat> object, an object with a new geo location 
        """
        file_name = os.path.split(self.ifile)[-1]
        print file_name

        if re.match(r'A', file_name) is None:
            raise IOError

        n = Nansat(self.ifile, GCP_COUNT=gcp_count)
        # Remove geo location
        n.vrt.remove_geolocationArray()
        n.vrt.tps = True
        n.reproject_GCPs()
        print n.time_coverage_start

        # add index of pixels
        index = np.arange(0, n.shape()[0] * n.shape()[1]).reshape(n.shape()).astype('int32')
        n.add_band(index, parameters={'name': 'index'})
        n.reproject(self.domain, addmask=False)

        bands = ['Rrs_%s' % band for band in self.wavelengths['modis'][wavelengths_set]]
        bands.insert(0, 'index')

        n_export = Nansat(domain=self.domain)
        for band in bands:
            print band
            band_arr = n[band]
            n_export.add_band(band_arr, parameters={'name': band})
            n_export.export(os.path.join(save_path, os.path.split(self.ifile)[1] + '_reprojected.nc'))

        return n_export

    def stich(self, domain, bands, gdirs, s2_user=False):
        """
        :param domain: 
        :param bands: list
        :param gdirs: list, list of directories paths 
        :param s2_user: bool, is it corrected file
        :return: 
        """
        # Create base nansat object which domain covers all four <granules>
        n_obj = Nansat(domain=domain)
        for band in bands:
            print('Band:', band, self.wavelengths['sentinel2'][band])
            # Generate 2d array which will have shape like <d> and filled by nan
            band_arr = np.zeros(domain.shape()) + np.nan
            for gdir in gdirs:
                data_path = os.path.join(gdir, 'IMG_DATA', '*_B%s.jp2' % band)

                if s2_user:
                    data_path = os.path.join(gdir, 'IMG_DATA', 'R60m', '*_B%s_60m.jp2' % band)
                try:
                    bfile = glob.glob(data_path)[0]
                    n = Nansat(bfile)
                    # Reprojection of data according to domain; eResampleAlg 1 is Bilinear
                    n.reproject(domain, eResampleAlg=1, addmask=False)
                    # Get data array from nansat object
                    bdata = n[1]
                    band_arr[bdata > 0] = bdata[bdata > 0]
                except IndexError:
                    print "Data path <%s> doesn't exist" % data_path
                    pass

            n_obj.add_band(band_arr, parameters={'name': 'Rrs_%s' % self.wavelengths['sentinel2'][band]})

        return n_obj

    def s2_downscale(self, save_path='./'):
        file_name = os.path.split(self.ifile)[-1]
        print file_name
        s2_user = False

        if re.match(r'S2A', file_name) is None:
            raise IOError

        gdirs = []
        for granule in self.granules:

            gdirs += sorted(glob.glob(os.path.join(self.ifile, 'GRANULE', '*_T%s_*' % granule)))

        # get lon/lat limits
        # Lists for accumulation of lon/lat values from each granule
        lons = []
        lats = []

        for gdir in gdirs:
            data_path = os.path.join(gdir, 'IMG_DATA', '*_B01.jp2')

            # We need other path pattern for corrected S2 data
            if re.match(r'S2A_USER', file_name) is not None:
                data_path = os.path.join(gdir, 'IMG_DATA', 'R60m', '*_B01_60m.jp2')
                s2_user = True

            # Find absolute path to relevant granule
            b0file = glob.glob(data_path)[0]
            # Get granule by nansat
            n = Nansat(b0file)
            lon, lat = n.get_corners()
            # Add min/max values of long and lat to list
            lons += list(lon)
            lats += list(lat)

        # Create domain according to max and min values of lon and lat
        d = Domain(n.vrt.get_projection(), '-lle %f %f %f %f -tr 60 60' % (min(lons), min(lats), max(lons), max(lats)))
        print('Domain created')

        n_obj = self.stich(d, sorted(self.wavelengths['sentinel2'].keys()), gdirs, s2_user=s2_user)
        n_obj.reproject(self.domain)
        export_name = os.path.split(self.ifile)[1] + '_reprojected.nc'
        n_obj.export(os.path.join(save_path, export_name))

    def s2_make_granules(self, save_path='./'):

        if re.match(r'S2A', os.path.split(self.ifile)[-1]) is None:
            raise IOError

        # Get list of granules
        # What is granules: https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/product-types
        # Add path to sentinel image <ifile>

        gdirs = sorted(glob.glob(os.path.join(self.ifile, 'GRANULE', '*')))

        lons = []  # List of longitudes for each granule
        lonVec = []  # ???
        lats = []  # List of latitudes for each granule
        latVec = []  # ???
        labels = []  # List on names for each granule
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
        mapfile = os.path.split(self.ifile)[1] + '_map.png'
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
        qlfile = os.path.split(self.ifile)[1] + '_ql.png'
        # Save on disk full generated image: #s2array
        plt.imsave(os.path.join(save_path, qlfile), s2array, vmin=1050, vmax=1400)
