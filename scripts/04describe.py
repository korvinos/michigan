# Find angle that gives max difference between bands
    
hiresfile = os.path.split(ifile)[1] + '.npz'
corfile = hiresfile.replace('.npz', '_cor.npz')

mask = np.load(hiresfile)['arr_0'][5] > 200

bands = range(0, 11)

corrected = []
angle = None
for band in bands:
    print band
    # load hi-res
    hires = np.load(hiresfile)['arr_0'][band].astype(np.float64)
    hires[mask] = np.nan

    # find angle
    if angle is None:
        angles = np.arange(23, 24, 0.1)
        maxdif = []
        for angle in angles:
            print angle
            hr_rot = rotate(hires, angle, order=0, cval=np.nan)
            plt.imshow(hr_rot[::5, ::5], vmin=950, vmax=1250);plt.show()
            # For what we make filtration?
            hr_mean = gaussian_filter(np.nanmean(hr_rot, axis=0), 10)
            maxdif.append(np.nanmax(np.abs(np.diff(hr_mean))))
            plt.plot(np.abs(np.diff(hr_mean)));plt.show()

        plt.plot(angles, maxdif, '.-');plt.show()
        angle = angles[np.argmax(maxdif)]

    hr_rot = rotate(hires, angle, order=0, cval=np.nan)
    hr_mean = np.nanmean(hr_rot, axis=0)
    #plt.imshow(hr_rot, vmin=950, vmax=1250);plt.show()
    plt.plot(hr_mean);plt.show()

    hr_mean_grd = np.repeat(hr_mean[None], hr_rot.shape[0], axis=0)

    hr_mean_grd_rot = rotate(hr_mean_grd, -angle, order=0, cval=np.nan)

    #plt.imshow(hr_rot - hr_mean_grd, vmin=-20, vmax=20);plt.colorbar();plt.show()


    h0, w0 = hires.shape
    aa = np.radians(angle)
    h1 = np.round(w0*np.sin(aa) + h0*np.cos(aa))
    w1 = np.round(h0*np.sin(aa) + w0*np.cos(aa))

    h2 = np.round(w1*np.sin(aa) + h1*np.cos(aa))
    w2 = np.round(h1*np.sin(aa) + w1*np.cos(aa))
    dc = (w2 - w0) / 2.
    dr = (h2 - h0) / 2.

    # Exception of rows of noise from main image
    hires_cln = hires - hr_mean_grd_rot[dr:dr+h0, dc:dc+w0]
    plt.imshow(hires_cln, vmin=-20, vmax=20);plt.colorbar();plt.show()

    corrected.append(hires_cln)

np.savez_compressed(corfile, np.array(corrected))


# make RGB from corrected and non corrected data
#%cpaste
b5 = np.load(hiresfile)['arr_0'][5]
mask = ((b5 > 200) + np.isnan(b5)).astype(np.uint8)

f0 = Figure(np.load(hiresfile)['arr_0'][:3],
            mask_array=mask,
            mask_lut={1:[255,255,255]},
            ratio=0.99)
cmin, cmax = f0.clim_from_histogram()
f0.process(cmin=cmin, cmax=cmax)
f0.save(hiresfile + '_rgb.jpg')


f1 = Figure(np.load(corfile)['arr_0'][:3],
            mask_array=mask,
            mask_lut={1:[255,255,255]},
            ratio=0.97)
cmin, cmax = f1.clim_from_histogram()
f1.process(cmin=cmin, cmax=cmax)
f1.save(corfile + '_rgb.jpg')
