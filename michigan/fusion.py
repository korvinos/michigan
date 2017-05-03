from nansat import Nansat
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from ovl_plugins.fusion.fusion import fuse
from dataprep import Data


class Fusion(Data):
    BATHYMETRY_PATH = './requirements/michigan_lld.grd'

    # Constants for Sentinel-2 image correction
    bMax = 2300  # OK
    # bMax = 200 # MAYBE A BIT TOO LOW
    # for noise corrected bands
    bMin = -30
    # bMax = 20 ???????

    cutsize = 2000

    def __init__(self, m_file, s_file, domain=None, smooth=False, skip=True, log=False, mask=True, cut=True,
                 negative_px=True, h_mask=9999, prepare_m=False, prepare_s=False):
        """
        :param s_file: str, path to Sentinel-2 file of <hiresfile>
        :param m_file: str, path to MODISa file or <loresfile>
        :param smooth: bool 
        :param skip: bool
        :param log: bool
        :param mask: bool
        :param cut: bool
        """

        if not domain:
            self.domain = self.sbd_dom

        self.m_file = m_file

        if prepare_m:
            Data.__init__(self, m_file)
            self.loresfile = self.modis_geo_location()
        else:
            self.loresfile = Nansat(m_file)

        # Open inputted files by Nansat
        # Load low resolution file - MODISa
        self.loresfile.reproject(self.domain)
        if negative_px:
            self.negpix = self.loresfile[2] < 0

        self.index = self.loresfile['index']

        # load hi-res
        # Load file with sentinel data
        self.s_file = s_file

        if prepare_s:
            Data.__init__(self, s_file)
            hiresfile = self.s2_downscale()
        else:
            hiresfile = Nansat(s_file)

        hiresfile.reproject(self.domain)

        # Get numbers of of each band
        band_rrs_numbers = [hiresfile._get_band_number('Rrs_%s' % wavelength)
                            for wavelength in sorted(self.wavelengths['sentinel2'].values())]
        # Create array form
        hires_np_array = np.array([hiresfile[i] for i in band_rrs_numbers])
        # remove out-of-swath
        hires_np_array[:, hires_np_array[0] == 0] = np.nan

        self.cut = cut

        if mask:
            hires_np_array = self.mask(hires_np_array, h_mask)

        if smooth and negative_px:
            hires_np_array = self.smooth(hires_np_array)

        if skip:
            hires_np_array = hires_np_array[0:5]

        if log:
            for hrn in range(len(hires_np_array)):
                hires_np_array[hrn] = np.log10(hires_np_array[hrn] + 1)

        if self.cut:
            hires_np_array, self.negpix, self.index = self.crop(hires_np_array)

        self.hires = hires_np_array

    def get_bottom(self, bathymetry_path=BATHYMETRY_PATH):
        bathymetry = Nansat(bathymetry_path)
        bathymetry.reproject(self.domain)
        # preparing of bottom field
        h = bathymetry[1]
        # all points there h >= 0 will marked as np.nan
        h = np.where(h >= 0, np.nan, np.float32(h) * -1)
        return h

    def get_land_mask(self, bathymetry_path=BATHYMETRY_PATH):
        h = self.get_bottom(bathymetry_path=bathymetry_path)
        # the mask of land
        land_mask = np.where(np.isfinite(h), np.nan, np.array(1))
        return land_mask

    def get_h_mask(self, h_max, h_min=None, bathymetry_path=BATHYMETRY_PATH, mask_val=np.nan):
        h_mask = self.get_bottom(bathymetry_path=bathymetry_path)
        h_mask[h_mask > h_max] = mask_val

        if h_min is not None:
            h_mask[h_mask < h_min] = mask_val

        return h_mask

    def smooth(self, hires_arr):
        ws = 1
        hires_arr[:, self.negpix] = np.nan
        hires_arr = gaussian_filter(hires_arr, (0, ws, ws))
        return hires_arr

    def mask(self, hires_arr, h_mask):
        watter_mask = self.loresfile.watermask()[1]
        watter_mask_filtered = gaussian_filter(watter_mask.astype(np.float32), 1)
        hires_arr[:, watter_mask_filtered > 1] = np.nan
        if h_mask is not 9999:
            hm = self.get_h_mask(h_mask)
            # return hm
            hires_arr[:, np.isnan(hm) == True] = np.nan

        # mask clouds
        hires_arr[:, hires_arr[7] > self.bMax] = np.nan
        hires_arr[:, hires_arr[0] < self.bMin] = np.nan
        return hires_arr

    def crop(self, hires_arr):
        hires_arr = hires_arr[:, :self.cutsize, :self.cutsize]
        negpix = self.negpix[:self.cutsize, :self.cutsize]
        index = self.index[:self.cutsize, :self.cutsize]
        return hires_arr, negpix, index

    def fusion(self, m_wavelengths='full'):
        bands = ['Rrs_%s' % wavelength for wavelength in self.wavelengths['modis'][m_wavelengths]]
        n_hires = Nansat(domain=self.domain)
        n_lores = Nansat(domain=self.domain)

        for band in bands:
            lores = self.loresfile[band]

            if self.cut:
                lores = lores[:self.cutsize, :self.cutsize]

            lores[self.negpix] = np.nan

            n_lores.add_band(lores, parameters={'name': band})
            # hires_fused = fuse(hires, lores, network_name=rgb_band,
            # iterations=100, threads=7, nn_structure=[5, 10, 7, 3])
            # TODO: We should use less number of iterations: 15 - 16
            hires_fused = fuse(self.hires, lores, network_name=band, iterations=20, threads=7, index=self.index)
            n_hires.add_band(hires_fused, parameters={'name': band})

        return n_lores, n_hires
