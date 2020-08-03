import matplotlib.pyplot as plt
import pandas as pd
import mops

case = mops.case("stagnation1")
print(case._outputFileName)
print(case.times())
print("End time = {}s".format(case.tfinal()))
tfinal = case.tfinal()
print(case.fileDict.keys())
psl_final = case["psl("+str(tfinal)+"s)"]
print(psl_final.dtypes)
psl_final["Equiv. Sphere Diameter (nm)"].plot.hist(bins=20, alpha=0.5)
plt.show()
