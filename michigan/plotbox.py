import numpy as np
import matplotlib.pyplot as plt


def show_it(obj, name, color='jet', limits=None, z10='off', size=(9, 6), coord_grid='off'):
    plt.figure(figsize=size)
    plt.imshow(obj, clim=limits, cmap=color)
    plt.colorbar()
    plt.title(name, y=1.08, fontsize=16)
    plt.tick_params(labelsize=14)

    if coord_grid == 'on':
        y = np.arange(0.0, 79.0, 79.0 / 7)
        y_labels = np.arange(45.3, 44.5, -0.1)
        x = np.arange(0.0, 123.0, 123.0 / 11)
        x_labels = np.arange(86.3, 85.2, -0.1)
        plt.yticks(y, y_labels)
        plt.xticks(x, x_labels)

    # black field is area which depth less than 10 meters
    # if z10 == 'on': plt.imshow(h_10m, cmap='gray')
    # plt.imshow(h_mask, cmap='Dark2')
    # plt.grid(color='black')

    plt.show()