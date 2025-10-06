import numpy as np
import matplotlib.colors as mcolors

def execute(data, data_min, data_max, data_threshold):
    data = np.clip(data, data_min, data_max) # clip values
    hue = (data - data_min) / (data_max - data_min) * (240.0 / 360.0)

    # assign color using HSV colormap
    hsv = np.zeros((hue.size, 3))
    hsv[:, 0] = hue
    hsv[:, 1] = 1.0
    hsv[:, 2] = 1.0
    map_color = mcolors.hsv_to_rgb(hsv)

    # assign non-active to gray
    non_active_id = (data <= data_threshold) | np.isnan(data)
    map_color[non_active_id, :] = 0.5

    return map_color
