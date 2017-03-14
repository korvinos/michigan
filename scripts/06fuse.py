import os
from nansat import *
from ovl_plugins.fusion.fusion import fuse
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import LogNorm

hiresfile = os.path.split(sfile)[1] + '_cor.npz'
loresfile = os.path.split(mfile)[1] + '_pro.nc'
pixel_size = 60

# TEST 1
bands = ['Rrs_667', 'Rrs_555', 'Rrs_488']
pref = 'test01_'
smooth = False
skip = True
log = False
mask = True
cut = True


# load lo-res
nnn = Nansat(loresfile)
negpix = nnn[2] < 0
index = nnn['index']

# load hi-res
hires = np.load(hiresfile)['arr_0'].astype(np.float64)

# remove out-of-swath
hires[:, hires[0]==0] = np.nan

bMax = 230 # OK
#bMax = 200 # MAYBE A BIT TOO LOW

# for noise corrected bands
bMin = -30
bMax = 20

cutsize = 2000

if mask:
    # mask land and coastal zone
    wm = nnn.watermask()[1]
    wmgf = gaussian_filter(wm.astype(np.float32), 1)
    hires[:, wmgf > 1] = np.nan

    # mask clouds
    hires[:, hires[7] > bMax] = np.nan
    hires[:, hires[0] < bMin] = np.nan

# smooth hires band for training
if smooth:
    ws = 1
    hires[:, negpix] = np.nan
    hires = gaussian_filter(hires, (0, ws, ws))

### show histograms
#for hr in hires:
#    plt.hist(hr[np.isfinite(hr)], 100);plt.show()

# remove SWIR bands
if skip:
    hires = hires[0:5]

# logscale
if log:
    for hrn in range(len(hires)):
        hires[hrn] = np.log10(hires[hrn] + 1)

if cut:
    hires = hires[:, :cutsize, :cutsize]
    negpix = negpix[:cutsize, :cutsize]
    index = index[:cutsize, :cutsize]

for rgb_band in bands:
    lores = nnn[rgb_band]
    if cut:
        lores = lores[:cutsize, :cutsize]
    lores[negpix] = np.nan
    np.savez_compressed('%s_hires_%s.npz' % (os.path.split(sfile)[1],
                                             rgb_band),
                        hires=hires,
                        lores=lores,
                        index=index
                    )

hires_rgb = []
lores_rgb = []
for rgb_band in bands:
    lores = nnn[rgb_band]
    if cut:
        lores = lores[:cutsize, :cutsize]
    lores[negpix] = np.nan
    lores_rgb.append(lores)
    #hires_fused = fuse(hires, lores, network_name=rgb_band, iterations=100, threads=7, nn_structure=[5, 10, 7, 3])
    hires_fused = fuse(hires, lores, network_name=rgb_band, iterations=20, threads=7, index=index)
    hires_rgb.append(hires_fused)

r0, r1 = 0, None
c0, c1 = 0, None


## VIS HI-RES RGB
wm = nnn.watermask()[1]
wm[np.isfinite(hires[0])] = 64

if cut:
    wm = wm[:cutsize, :cutsize]

wm_crp = wm[r0:r1, c0:c1]
hires_rgb_crp = np.array(hires_rgb)[:, r0:r1, c0:c1]

f = Figure(hires_rgb_crp,
            mask_array=wm_crp,
            mask_lut={1:[220,220,220], 2:[128,128,128]},
            ratio=0.99)
clim = f.clim_from_histogram()
f.process(cmin=clim[0], cmax=clim[1])
f.save(pref + os.path.split(mfile)[1] + '_rgb_hires.png')

## VIS LO-RES RGB
wm = nnn.watermask()[1]
wm[~negpix] = 64

if cut:
    wm = wm[:cutsize, :cutsize]

wm_crp = wm[r0:r1, c0:c1]
lores_rgb_crp = np.array(lores_rgb)[:, r0:r1, c0:c1]

f = Figure(lores_rgb_crp,
            mask_array=wm_crp,
            mask_lut={1:[220,220,220], 2:[128,128,128]},
            ratio=0.99)
f.process(cmin=clim[0], cmax=clim[1] )
f.save(pref + os.path.split(mfile)[1] + '_rgb_lores.png')

## VIS scatter plot comparison

vmax = [0.0007, 0.0035, 0.0045]

for rgbi in range(3):
    gpi = np.isfinite(lores_rgb[rgbi]) * np.isfinite(hires_rgb[rgbi])
    r = np.corrcoef([lores_rgb[rgbi][gpi], hires_rgb[rgbi][gpi]])[1,0]
    plt.hist2d(lores_rgb[rgbi][gpi], hires_rgb[rgbi][gpi], 100,
               range=[[0, vmax[rgbi]], [0, vmax[rgbi]]], norm=LogNorm())
    plt.colorbar(shrink=0.8)
    plt.title('R=%f' % r)
    plt.xlabel('$R_{rsMODIS}, sr^{-1}$')
    plt.ylabel('$R_{rsFUSED}, sr^{-1}$')
    plt.savefig(pref + 'hist2d_%02d.png' % rgbi, bbox_inches='tight', pad_inches=0)
    plt.close()

gpi = np.isfinite(lores_rgb[0]) * np.isfinite(hires_rgb[0])
plt.hist2d(lores_rgb[0][gpi], hires_rgb[0][gpi], 100);plt.show()

plt.hist2d(hires[0][gpi], hires_rgb[2][gpi], 100);plt.show()


#plt.hist2d(lores_rgb[0][gpi], hires_rgb_tmp[0][gpi], 1000, range=[[0, 0.02], [0, 0.02]]);plt.show()

#plt.hist2d(lores_rgb[1][gpi], hires_rgb[1][gpi], 1000, range=[[0, 0.02], [0, 0.02]]);plt.show()
#plt.hist2d(lores_rgb[2][gpi], hires_rgb[2][gpi], 1000, range=[[0, 0.02], [0, 0.02]]);plt.show()
#plt.hist2d(lores_rgb[2][gpi], hires[0][gpi], 1000);plt.show()

#hrgf = [gaussian_filter(hrrgb, 2) for hrrgb in hires_rgb]
#f = Figure(hrgf)
#f.process(cmin=clim[0], cmax=clim[1])
#f.save(pref + mfile + '_rgb_hires_gf.png')
