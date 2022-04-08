import pandas as pd
import numpy as np
import math

from scipy.interpolate import RegularGridInterpolator


class Mars_Atmosphere:

    def __init__(self, metadata_path="mars_atmosphere/dataframes"):
        self.interpolators = {}

        metadata = pd.read_csv(f"{metadata_path}/metadata.csv", index_col=0).sort_index()

        for key in ["T", "p", "rho", "mu"]:
            data = []

            for altitude, row in metadata.iterrows():
                current_altitude_data = pd.read_csv(row[key], index_col=0)

                data.append(current_altitude_data.to_numpy())

                longitudes = np.array(current_altitude_data.index.values)
                latitudes = np.array(current_altitude_data.columns.values).astype(np.float32)

            altitudes = np.array(metadata.index.values)

            self.interpolators[key] = RegularGridInterpolator((altitudes, longitudes, latitudes), np.array(data))

    def __call__(self, t, longitude, latitude, altitude):
        if altitude < 0.0:
            return np.array([0.0, 0.0, 0.0])

        pt = np.array([altitude, longitude, latitude])

        T = self.interpolators["T"](pt) # temperature (K)
        p = self.interpolators["p"](pt) # pressure (Pa)
        rho = self.interpolators["rho"](pt) # density (kg/m3)
        mu = self.interpolators["mu"](pt) # viscosity coefficient (kg/(m*s))
        cs = math.sqrt(1.29 * 191.8 * T) # speed of sound (m/s)

        return np.array([p, rho, cs, mu], dtype=np.float32)

    def get_T(self, altitude, longitude, latitude):
        pt = np.array([altitude, longitude, latitude])
        return self.interpolators["T"](pt)

    def get_p(self, altitude, longitude, latitude):
        pt = np.array([altitude, longitude, latitude])
        return self.interpolators["p"](pt)

    def get_rho(self, altitude, longitude, latitude):
        pt = np.array([altitude, longitude, latitude])
        return self.interpolators["rho"](pt)

    def get_mu(self, altitude, longitude, latitude):
        pt = np.array([altitude, longitude, latitude])
        return self.interpolators["mu"](pt)

    def get_altitude_heatmaps(self, altitude, rows, cols, key):
        longitudes = np.linspace(-180, 180, cols)
        latitudes = np.linspace(90, -90, rows)

        xyz_grid = np.meshgrid(altitude, longitudes, latitudes, indexing="ij")
        xyz_list = np.reshape(xyz_grid, (3, -1), order="C").T

        result = self.interpolators[key](xyz_list)
        heatmap = result.reshape([cols, rows])
        heatmap = heatmap.T

        return heatmap


