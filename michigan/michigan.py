from fusion import Fusion
import numpy as np
from nansat import Nansat
from boreali import Boreali, lm
from sklearn.cluster import KMeans


class MichiganProcessing(Fusion):

    BATHYMETRY_PATH = './requirements/michigan_lld.grd'

    def __init__(self, m_file, s_file=None, fuse=False, domain=False, reproject=False):
        """
        :param m_file: 
        :param s_file: 
        :param fusion: 
        :param wavelengths: 
        :param domain: 
        """
        # TODO: Check out info about inheritance of <__init__> method

        if domain:
            self.domain = domain
        else:
            self.domain = self.sbd_dom

        if fuse:
            try:
                Fusion.__init__(self, m_file, s_file)
                self.ifile = self.fusion()
            except TypeError:
                print 'Requested Sentinel-2 file was not found'

        else:
            self.ifile = Nansat(m_file)
            if reproject:
                self.ifile.reproject(self.domain)

    def get_bottom(self, bathymetry_path=BATHYMETRY_PATH):
        bathymetry = Nansat(bathymetry_path)
        bathymetry.reproject(self.domain)
        # preparing of bottom field
        h = bathymetry[1]
        # all points there h >= 0 will marked as np.nan
        h = np.where(h >= 0, np.nan, np.float32(h) * -1)
        # the mask of land
        h_mask = np.where(np.isfinite(h), np.nan, np.array(1))

        # the mask which shows 10 meters depth zone
        h_10m = np.where(h <= 10, np.array(1), np.nan)
        return h

    def boreali_processing(self, wavelengths_set='1x1km_bands', bottom_type=0, osw_mod='on', bathymetry_path=BATHYMETRY_PATH,
                           hydro_optic='michigan'):

        wavelengths = self.wavelengths['modis'][wavelengths_set]

        cpa_limits = [0.01, 3,
                      0.01, 1,
                      0.01, 1, 10]

        b = Boreali(hydro_optic, wavelengths)
        theta = np.zeros_like(self.ifile[2])
        custom_n = Nansat(domain=self.ifile)
        band_rrs_numbers = list(map(lambda x: self.ifile._get_band_number('Rrs_' + str(x)), wavelengths))

        for index in range(0, len(wavelengths)):
            rrsw = self.ifile[band_rrs_numbers[index]] / (0.52 + 1.7 * self.ifile[band_rrs_numbers[index]])
            custom_n.add_band(rrsw, parameters={'name': 'Rrsw_' + str(wavelengths[index]),
                                                'units': 'sr-1',
                                                'wavelength': wavelengths[index]})

            # If we want to use OSW mod, we will need to add Rrs data in custom_n obj
            if osw_mod == 'on':
                custom_n.add_band(self.ifile[band_rrs_numbers[index]], parameters={'name': 'Rrs_' + str(wavelengths[index]),
                                                                              'units': 'sr-1',
                                                                              'wavelength': wavelengths[index]})
        # Creating of the mask
        # All pixels marked as -0.015534 in img will marked as 0.0 in the mask
        mask = np.where(self.ifile[2] != np.float(-0.015534), np.array(64.0), np.array(0.0))
        # Validation of mask according to bathymetry data.
        # If in the bathymetry pixel was marked as np.nan, in mask he will marked as 0.0
        # else nothing
        h = self.get_bottom(bathymetry_path=bathymetry_path)
        mask = np.where(np.isnan(h), np.array(0.), mask)
        # Adding of mask into custom_n obj
        custom_n.add_band(mask, parameters={'name': 'mask'})

        # h is trigger for OWS processing mod.
        # If want to star OSW we will need to add h to Boreali.process
        # else marked it as None
        if osw_mod == 'on':
            depth = h
        else:
            depth = None

        cpa = b.process(custom_n, cpa_limits, mask=custom_n['mask'], depth=depth, theta=theta, threads=4)

        custom_n.add_band(array=cpa[0], parameters={'name': 'chl',
                                                    'long_name': 'Chlorophyl-a',
                                                    'units': 'mg m-3'})
        custom_n.add_band(array=cpa[1], parameters={'name': 'tsm',
                                                    'long_name': 'Total suspended matter',
                                                    'units': 'g m-3'})
        custom_n.add_band(array=cpa[2], parameters={'name': 'doc',
                                                    'long_name': 'Dissolved organic carbon',
                                                    'units': 'gC m-3'})
        custom_n.add_band(array=cpa[3], parameters={'name': 'mse',
                                                    'long_name': 'Root Mean Square Error',
                                                    'units': 'sr-1'})
        custom_n.add_band(array=cpa[4], parameters={'name': 'mask',
                                                    'long_name': 'L2 Boreali mask',
                                                    'units': '1'})

        #   custom_n.export(final_path + obj.split('/')[-1] + 'cpa_OSW.nc')
        return custom_n

    def get_r(self, coords, wavelengths, r_type='Rrs_'):
        y, x = coords
        band_numbers = list(map(lambda x: self.ifile._get_band_number(r_type + str(x)), wavelengths))
        r_list = [self.ifile[band][y][x] for band in band_numbers]
        return r_list

    def boreali_lm(self, coords, wavelengths_set='1x1km_bands', bottom_type=0, show='on', title=None, hydro_optic='michigan'):
        """
        :param wavelengths: list.
        :param coords: tuple. i, j 
        :param bottom_type: int.
        0 - sand; 1 - sargassum; 2 - silt; 3 - boodlea; 4 - limestone; 
        5 - enteromorpha; 6 - cladophora; 7 - chara; 8 - sand2; 9 - sand3; 10 - charasand;
        11 - cladforaosand; 12 - sandmichi; 13 - cladomichi; 
        """
        wavelengths = self.wavelengths['modis'][wavelengths_set]
        y, x = coords
        b = Boreali(hydro_optic, wavelengths)
        model = b.get_homodel()
        theta = 0
        albedoType = bottom_type
        depth = float(self.get_bottom()[y, x])
        # depth = 7
        albedo = b.get_albedo([albedoType])[0]

        cpa_limits = [0.01, 3,
                      0.01, 1,
                      0.01, 1, 10]

        r = self.get_r((y, x), 'Rrsw_', wavelengths)
        c_deep = lm.get_c_deep(cpa_limits, model, [r], [theta], 4)[1]
        c_osw = lm.get_c_shal(cpa_limits, model, [r], [theta], [depth], [albedo], 4)[1]

        # backward case
        rrsw_deep = lm.get_rrsw_deep(model, c_deep[0:3], theta, len(wavelengths))[1]
        rrsw_osw = lm.get_rrsw_shal(model, c_osw[0:3], theta, depth, albedo, len(wavelengths))[1]
        # rrsw_osw_mod = lm.get_rrsw_shal(model, np.array([0.5, 0, 0]), theta, depth, albedo, len(wavelens))[1]

        if show == 'on':
            print title
            print 'BOREALI DEEP: chl=%5.2f, tsm=%5.2f, doc=%5.2f, rmse=%5.2f' % tuple(c_deep)
            print 'BOREALI OSW: chl=%5.2f, tsm=%5.2f, doc=%5.2f, rmse=%5.2f' % tuple(c_osw)
            #
            # num_bands = range(0, len(r))
            #
            # plt.figure(figsize=(12, 5))
            # plt.title('%s / %s / Point: (x: %d, y:%d) / Depth: %5.2f' % (title, 'Rrsw', x, y, depth))
            # plt.tick_params(labelsize=14)
            # plt.plot(num_bands, r, label='rrsw')
            # plt.plot(num_bands, rrsw_osw, label='rrsw_osw')
            # plt.plot(num_bands, rrsw_osw_mod, label='rrsw_osw_mod')
            # plt.plot(num_bands, rrsw_deep, label='rrsw_deep')
            #
            # plt.xticks(num_bands, wavelens)
            # plt.xlabel('wavelength, nm', fontsize=14)
            # plt.ylabel('%s, sr^-1' % ('Rrsw_'), fontsize=14)
            # plt.grid(color='black')
            # plt.legend()
            # plt.ylim([0, 0.016])
            # # plot_r(obj, wavelens=wavelens, coords=(y, x), size=(10, 4), title='title')
        else:
            return c_osw, c_deep, depth

    def bottom_classification(self, clusters, wavelengths_set='1x1km_bands'):

        k_train_arr = []
        for band in self.wavelengths['modis'][wavelengths_set]:
            k_train_arr.append(self.ifile['Rrs_%s' % (band)].flatten())

        k_train_arr = np.array(k_train_arr).T

        kmeans = KMeans(n_clusters=clusters, random_state=0).fit(k_train_arr)
        arr = kmeans.labels_
        arr = arr.reshape(self.domain.shape())
        return arr
