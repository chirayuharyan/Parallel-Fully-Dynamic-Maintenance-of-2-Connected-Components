import os
import sys
import os.path
from os import path
import subprocess

dataset=str(sys.argv[1])
if(path.exists("RawDataset/"+dataset) == 0):
	print("Dataset doesn't exists. Please check the name.")
	sys.exit()

datasetExtension = dataset.split('.')[-1]



subprocess.call("make process -C Data-Processing_Codes/",shell = True)
subprocess.call("make processYifanhu -C Data-Processing_Codes/",shell = True)

if(datasetExtension == "mtx"):
	print("\n\n ******* Running Preprocessing Program (processYifanhu)*************")
	subprocess.call("./Data-Processing_Codes/processYifanhu RawDataset/"+dataset+" Dataset/"+dataset[:-4]+".txt",shell = True)
else:
	print("\n\n ******* Running Preprocessing Program (process)*************")
	subprocess.call("./Data-Processing_Codes/process RawDataset/"+dataset+" Dataset/"+dataset,shell = True)