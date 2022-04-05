import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator

cfd_tables_path = "cfd_tables"

cfd_table_dataframe = pd.read_csv(f"{cfd_tables_path}/msl.csv")

print(cfd_table_dataframe.head())
cfd_table = cfd_table_dataframe.to_numpy()

alpha = cfd_table[:, 0]
Ma = cfd_table[:, 1]
V = cfd_table[:, 2]
rho = cfd_table[:, 3]
T = cfd_table[:, 4]

data = cfd_table[:, 5:7]

print(cfd_table[:, 0:2])
print(data)

interpolator = RegularGridInterpolator((alpha, Ma), np.array(data))

print(interpolator(12.0, 20))