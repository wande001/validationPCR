import ConfigParser
import io
import sys
import os
import numpy as np
from scipy.stats import cumfreq
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import netCDF4 as nc


configFile = sys.argv[1]
riverID = sys.argv[2:]

def readConfigFile(configFileName):
  global config
  with open(configFileName) as f:
    sample_config = f.read()
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.readfp(io.BytesIO(sample_config))
  return config

def plotCDF(forecast, validation, title):
  ax1 = plt.figure(figsize=(7,5))
  vals, x1, x2, x3 = cumfreq(forecast['modelled'], len(forecast['modelled']))
  ax1 = plt.plot(np.linspace(np.min(forecast['modelled']), np.max(forecast['modelled']), len(forecast['modelled'])), vals/len(forecast['modelled']), "r", label=str(config.get('Main options', 'RunName')))
  vals, x1, x2, x3 = cumfreq(validation['modelled'], len(validation['modelled']))
  ax2 = plt.plot(np.linspace(np.min(validation['modelled']), np.max(validation['modelled']), len(validation['modelled'])), vals/len(validation['modelled']), "b", label=str(config.get('Reference options', 'RunName')))
  vals, x1, x2, x3 = cumfreq(validation['observations'], len(validation['observations']))
  ax3 = plt.plot(np.linspace(np.min(validation['observations']), np.max(validation['observations']), len(validation['observations'])), vals/len(validation['observations']), "black", label="Observations")
  ax3 = plt.legend(prop={'size': 10}, loc=2)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Discharge (m3/s)")
  ax1 = plt.ylabel("ECDF")
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotScatter(forecast, title):
  ax1 = plt.figure(figsize=(7,5))
  ax1 = plt.plot(forecast['modelled'], forecast['observations'], "bo", markersize=8)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Modelled")
  ax1 = plt.ylabel("Observed")
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotTimeLines(forecast, reference, title):
  ax1 = plt.figure(figsize=(12,5))
  ax1 = plt.plot(forecast['times'], forecast['modelled'], "r", markersize=8, label=str(config.get('Main options', 'RunName')))
  ax1 = plt.plot(reference['times'], reference['modelled'], "b", markersize=8, label=str(config.get('Reference options', 'RunName')))
  ax1 = plt.plot(forecast['times'], forecast['observations'], "black", markersize=8, label="Observations", lw=2)
  print forecast['times']
  ax1 = plt.legend(prop={'size': 10}, loc=2)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("")
  ax1 = plt.ylabel("Discharge (m3/s)")
  #ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def getOptions(config, option = "general"):
  includeRef = config.get('Output options', 'includeReference') == str("True")
  return includeRef

config = readConfigFile(configFile)

run1 = str(config.get('Main options', 'RunName'))
run2 = str(config.get('Reference options', 'RunName'))
includeRef = getOptions(config, "general")

output, output2, fullOutput, fullOutput2, waterBalOutput, waterBalOutput2 = pickle.load(open('validationResultsPool_%s_%s.obj' %(run1, run2), 'rb') )

IDs = [i for i,x in enumerate(fullOutput["ID"]) if x == riverID[0]]

pdf = PdfPages('riverResults_%s_%s.pdf' %(run1, run2))
for ID in IDs:
  
  plotTimeLines(fullOutput["data"][ID], fullOutput2["data"][ID], "River")
  plotScatter(fullOutput["data"][ID], "Simulations %s, R= %.2f, AC = %.2f, KGE = %.2f" %(str(config.get('Main options', 'RunName')), output[IDs[0],3], output[IDs[0],4], output[IDs[0],5]))
  plotScatter(fullOutput2["data"][ID], "Simulations %s, R= %.2f, AC = %.2f, KGE = %.2f" %(str(config.get('Reference options', 'RunName')), output2[IDs[0],3], output2[IDs[0],4], output2[IDs[0],5]))
  plotCDF(fullOutput["data"][ID], fullOutput2["data"][ID], "River")

pdf.close()


