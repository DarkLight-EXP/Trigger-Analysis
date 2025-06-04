import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd  # Make sure to import pandas
from scipy.optimize import curve_fit, minimize
from scipy.signal import fftconvolve



with uproot.open("/Users/gabbygelinas/Desktop/Masters/RootFiles/output00049.root") as f:
    fileID = "t14ns_1458"  # change to this "t23_2367" and run the code for a different dataset to export  this means we are looking at the time of 23 with concidence with 67
    timeDiff = f["dltdc"][fileID].to_numpy()[0] #y data
    #t14 = f["dltdc"]["t14ns_1458"].to_numpy()[0]  #"t14ns_1458"
    #t58 = f["dltdc"]["t58ns_1458"].to_numpy()[0]
    #t23 = f["dltdc"]["t23ns_2367"].to_numpy()[0]
    #t67 = f["dltdc"]["t67ns_2367"].to_numpy()[0]
    x_root = f["dltdc"][fileID].to_numpy()[1][0:-1] #x data
    
#need timeDiff(# of counts per bin) and x_root(bins)

# Create a DataFrame
df = pd.DataFrame({
    'timeDiff': timeDiff,
    'x_root': x_root
})

# Export to CSV
csv_path = "/Users/gabbygelinas/Desktop/Masters/Convolved_fitting"
df.to_csv(csv_path, index=False)