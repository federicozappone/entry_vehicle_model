import pandas as pd
import numpy as np

from scipy.interpolate import RegularGridInterpolator


class Mars_Atmosphere:

    def __init__(self, metadata_path="mars_atmosphere/dataframes"):
        self.interpolators = {}

        metadata = pd.read_csv(f"{metadata_path}/metadata.csv", index_col=0).sort_index()

        for key in ["T", "p", "rho"]:
            data = []

            for altitude, row in metadata.iterrows():
                current_altitude_data = pd.read_csv(row[key], index_col=0)

                data.append(current_altitude_data.to_numpy())

                longitudes = np.array(current_altitude_data.index.values)
                latitudes = np.array(current_altitude_data.columns.values).astype(np.float32)

            altitudes = np.array(metadata.index.values)

            data = np.array(data)

            self.interpolators[key] = RegularGridInterpolator((altitudes, longitudes, latitudes), data)

    def get(self, altitude, longitude, latitude):
        pt = np.array([altitude, longitude, latitude])

        T = self.interpolators["T"](pt)
        p = self.interpolators["p"](pt)
        rho = self.interpolators["rho"](pt)

        return T, p, rho
