import numpy as np
import matplotlib.pyplot as plt


def show_it(obj, name='no name', color='jet', limits=None, z10='off', size=(9, 6), coord_grid='off', land_mask=None):

    plt.figure(figsize=size)
    plt.imshow(obj, clim=limits, cmap=color)
    plt.colorbar()
    plt.title(name, y=1.08, fontsize=16)
    plt.tick_params(labelsize=14)

    if land_mask is not None:
        plt.imshow(land_mask, cmap='gray')

    if coord_grid == 'on':
        y = np.arange(0.0, 1248.0, 1248.0 / 7)
        y_labels = np.arange(45.3, 44.5, -0.1)
        x = np.arange(0.0, 1952.0, 1952.0 / 11)
        x_labels = np.arange(86.3, 85.2, -0.1)
        plt.yticks(y, y_labels)
        plt.xticks(x, x_labels)
        plt.xlabel('west')
        plt.ylabel('north')

    # black field is area which depth less than 10 meters
    # if z10 == 'on': plt.imshow(h_10m, cmap='gray')
    # plt.imshow(h_mask, cmap='Dark2')

    plt.grid(color='black')
    plt.show()


def plot_r(obj, r_set, wavelens, coords, r_type='Rrs_', size=(9, 6), title=None):
    y, x = coords
    num_bands = range(0, len(r_set))

    plt.figure(figsize=size)
    plt.title('%s / %s / Point: (x: %d, y:%d) / Depth: %5.2f' % (title, r_type, x, y, h[y, x]))
    plt.tick_params(labelsize=14)
    plt.plot(num_bands, r_set)
    plt.xticks(num_bands, wavelens)
    plt.xlabel('wavelength, nm', fontsize=14)
    plt.ylabel('%s, sr^-1' % (r_type), fontsize=14)
    plt.grid(color='black')
    plt.ylim([0, 0.008])
