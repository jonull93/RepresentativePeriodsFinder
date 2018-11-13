import os
import matplotlib.pyplot as plt
import gams_addon as ga
import numpy as np
import pandas as pd

csv_file = os.path.join(os.getcwd(), "..", "..", "..", "data", "Merged[Load,Wind,PV]_2015_CLEANED.csv")

df = pd.read_csv(csv_file, sep=';')
print df.head()
x = 1.0 * np.arange(len(df)) / len(df) * 100.
y = df.sort_values(by="ActualTotalLoad_BE", ascending=False)["ActualTotalLoad_BE"]
print x
plt.plot(x, y, "r")

csv_file = os.path.join(os.getcwd(), "..", "..", "..", "results", "test_plotting", "resulting_profiles.csv")
df2 = pd.read_csv(csv_file, sep=';')
print df2.head()

csv_file = os.path.join(os.getcwd(), "..", "..", "..", "results", "test_plotting", "decision_variables_short.csv")
df3 = pd.read_csv(csv_file, sep=';')
print df3.head()
df2["WEIGHT"] = list(np.repeat(df3["weights"], 24))
y2 = df2.sort_values(by="ActualTotalLoad_BE", ascending=False)["ActualTotalLoad_BE"]
x2 = df2.sort_values(by="ActualTotalLoad_BE", ascending=False)["WEIGHT"].cumsum() / 8760. * 100.

plt.plot(x2, y2, "b:")
plt.show()
