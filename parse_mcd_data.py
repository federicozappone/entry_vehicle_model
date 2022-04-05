import numpy as np
import pandas as pd
import os


base_path = "mars_atmosphere"
files = []

metadata = {}

for (dirpath, dirnames, filenames) in os.walk(base_path):
    files.extend(filenames)
    break

files.sort()

for filename in files:
    file = open(f"{base_path}/{filename}", "r")
    lines = file.readlines()

    altitude = int(filename.replace(".txt", ""))
    print(f"processing file {filename} for altitude {altitude}")

    keys = ["T", "p", "rho", "mu"]

    tables_dict = {}

    tables_dict["T"] = lines[11:77]
    tables_dict["p"] = lines[88:154]
    tables_dict["rho"] = lines[165:231]
    tables_dict["mu"] = lines[242:308]

    # remove empty line
    tables_dict["T"].pop(1)
    tables_dict["p"].pop(1)
    tables_dict["rho"].pop(1)
    tables_dict["mu"].pop(1)

    metadata_row = []

    for key in keys:
        data_dict = {}

        columns = tables_dict[key][0].replace("---- ||   ", "").split("   ")
        columns = [float(numeric_string) for numeric_string in columns]

        for line in tables_dict[key][1:]:
            data_with_row_index = line.strip().split(" ||    ")
            row_index = float(data_with_row_index[0])

            data = data_with_row_index[1].split("    ")
            data = [float(numeric_string) for numeric_string in data]

            data_dict[row_index] = data

        dataframe = pd.DataFrame.from_dict(data_dict, orient="index", columns=columns)

        csv_filename = f"{base_path}/dataframes/{altitude}_{key}.csv"
        dataframe.to_csv(csv_filename)

        metadata_row.append(csv_filename)

    metadata[altitude] = metadata_row


metadata_dataframe = pd.DataFrame.from_dict(metadata, orient="index", columns=["T", "p", "rho", "mu"]).sort_index()
metadata_dataframe.to_csv(f"{base_path}/dataframes/metadata.csv")

print(metadata_dataframe.head())
