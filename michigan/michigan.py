from fusion import Fusion
import numpy as np
from nansat import Nansat
from boreali import Boreali


class MichiganProcessing(Fusion):

    BATHYMETRY_PATH = './requirements/michigan_lld.grd'

    def __init__(self, m_file, s_file=None, fuse=False, wavelengths='full', domain=False, reproject=False):
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

        if fuse:
            try:
                Fusion.__init__(self, m_file, s_file)
                self.ifile = self.fusion()
            except TypeError:
                print 'Requested Sentinel-2 file was not found'

        else:
            self.ifile = Nansat(self.m_file)
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

    def boreali_processing(self, h, wavelenghts='full', bottom_type=0, osw_mod='on', bathymetry_path=BATHYMETRY_PATH):

        wavelenghts = self.wavelengths['modis'][wavelenghts]

        cpa_limits = [0.01, 3,
                      0.01, 1,
                      0.01, 1, 10]

        b = Boreali('michigan', wavelenghts)
        theta = np.zeros_like(self.ifile[2])
        custom_n = Nansat(domain=self.ifile)
        band_rrs_numbers = list(map(lambda x: self.ifile._get_band_number('Rrs_' + str(x)), wavelenghts))

        for index in range(0, len(wavelenghts)):
            rrsw = self.ifile[band_rrs_numbers[index]] / (0.52 + 1.7 * self.ifile[band_rrs_numbers[index]])
            custom_n.add_band(rrsw, parameters={'name': 'Rrsw_' + str(wavelenghts[index]),
                                                'units': 'sr-1',
                                                'wavelength': wavelenghts[index]})

            # If we want to use OSW mod, we will need to add Rrs data in custom_n obj
            if osw_mod == 'on':
                custom_n.add_band(self.ifile[band_rrs_numbers[index]], parameters={'name': 'Rrs_' + str(wavelenghts[index]),
                                                                              'units': 'sr-1',
                                                                              'wavelength': wavelenghts[index]})
        # Creating of the mask
        # All pixels marked as -0.015534 in img will marked as 0.0 in the mask
        mask = np.where(self.ifile[2] != np.float(-0.015534), np.array(64.0), np.array(0.0))
        # Validation of mask according to bathymetry data.
        # If in the bathymetry pixel was marked as np.nan, in mask he will marked as 0.0
        # else nothing
        mask = np.where(np.isnan(h), np.array(0.), mask)
        # Adding of mask into custom_n obj
        custom_n.add_band(mask, parameters={'name': 'mask'})

        # h is trigger for OWS processing mod.
        # If want to star OSW we will need to add h to Boreali.process
        # else marked it as None
        if osw_mod == 'on':
            depth = self.get_bottom(bathymetry_path=bathymetry_path)
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
