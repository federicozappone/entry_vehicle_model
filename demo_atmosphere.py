import cv2
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mars_atmosphere import Mars_Atmosphere


atmosphere = Mars_Atmosphere()

image_rows = 1000
image_cols = 2000
altitude = 125000

heatmaps = {}
keys = ["T", "p", "rho", "mu"]

for key in keys:
    heatmaps[key] = atmosphere.get_altitude_heatmaps(altitude, image_rows, image_cols, key)

fig, axs = plt.subplots(len(keys), figsize=(15, 35))

x_ticks = np.linspace(0, image_cols, 13)
x_ticklabels = ["180°", "150°W", "120°W", "90°W", "60°W", "30°W", "0°", "30°E", "60°E", "90°E", "120°E", "150°E", "180°"]

y_ticks = np.linspace(0, image_rows, 13)
y_ticklabels = ["90°N", "75°N", "60°N", "45°N", "30°N", "15°N", "0°", "15°S", "30°S", "45°S", "60°S", "75°S", "90°S"]

for ax, key in zip(axs, keys):
    ax.set_title(key)
    im = ax.imshow(heatmaps[key], cmap="jet")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    fig.colorbar(im, cax=cax, orientation="vertical")

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticklabels)
    ax.grid()

plt.show()
