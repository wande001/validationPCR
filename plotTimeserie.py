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
import datetime


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
  ax1 = plt.plot(forecast['observations'], forecast['modelled'], "bo", markersize=8)
  xmin, xmax = plt.xlim()
  ymin, ymax = plt.ylim()
  ax1 = plt.plot([0,10**10],[0,10**10], 'k-')
  ax1 = plt.xlim(xmin, xmax)
  ax1 = plt.ylim(ymin, ymax)
  ax1 = plt.title(title)
  ax1 = plt.ylabel("Modelled")
  ax1 = plt.xlabel("Observed")
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotTimeLines(forecast, reference, title):
  ax1 = plt.figure(figsize=(12,5))
  plotTimes = forecast['times']
  if forecast['times'][0].month != forecast['times'][1].month:
    plotTimes = monthSeries(forecast['times'])
  ax1 = plt.plot(plotTimes, forecast['modelled'], "r", markersize=8, label=str(config.get('Main options', 'RunName')))
  ax1 = plt.plot(plotTimes, forecast['observations'], "black", markersize=8, lw=2)
  if reference['times'][0].month != reference['times'][1].month:
    plotTimes = monthSeries(reference['times'])
  ax1 = plt.plot(plotTimes, reference['modelled'], "b", markersize=8, label=str(config.get('Reference options', 'RunName')))
  ax1 = plt.plot(plotTimes, reference['observations'], "black", markersize=8, label="Observations", lw=2)
  ax1 = plt.legend(prop={'size': 10}, loc=2)
  ax1 = plt.title(title)
  xmin, xmax = plt.xlim()
  ax1 = plt.xlim(xmin*0.9995, xmax*1.0005)
  ax1 = plt.xlabel("")
  ax1 = plt.ylabel("Discharge (m3/s)")
  #ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def monthSeries(times):
  outputTimes = []
  for time in times:
    outputTimes.append(datetime.datetime(time.year, time.month, 15))
  return outputTimes

def getOptions(config, option = "general"):
  includeRef = config.get('Output options', 'includeReference') == str("True")
  return includeRef

config = readConfigFile(configFile)

run1 = str(config.get('Main options', 'RunName'))
run2 = str(config.get('Reference options', 'RunName'))
includeRef = getOptions(config, "general")

output, output2, fullOutput, fullOutput2, waterBalOutput, waterBalOutput2 = pickle.load(open('validationResultsPool_%s_%s.obj' %(run1, run2), 'rb') )

IDs = []
nameS = []
titleS = []
for ID in riverID:
  find = False
  for i,x in enumerate(fullOutput["ID"]):
    if x == ID:
      if find == False:
        find = True
        IDs.append(i)
        nameS.append(ID)
        if output[i,-1] == 1:
          titleS.append("Daily")
        else:
          titleS.append("Monthly")
      elif find and output[IDs[-1],-1] != output[i,-1]:
        IDs.append(i)
        nameS.append(ID)
        if output[i,-1] == 1:
          titleS.append("Daily")
        else:
          titleS.append("Monthly")

IDs = np.array(IDs).flatten()
print IDs

pdf = PdfPages('riverResults_%s_%s.pdf' %(run1, run2))
for ID, name, title in zip(IDs, nameS, titleS):
  plotTimeLines(fullOutput["data"][ID], fullOutput2["data"][ID], "GRDC River %s %s" %(name, title))
  plotScatter(fullOutput["data"][ID], "%s %s, R= %.2f, AC = %.2f, KGE = %.2f" %(title, str(config.get('Main options', 'RunName')), output[ID,3], output[ID,4], output[ID,5]))
  plotScatter(fullOutput2["data"][ID], "%s %s, R= %.2f, AC = %.2f, KGE = %.2f" %(title, str(config.get('Reference options', 'RunName')), output2[ID,3], output2[ID,4], output2[ID,5]))
  plotCDF(fullOutput["data"][ID], fullOutput2["data"][ID], "GRDC River %s %s" %(name, title))

pdf.close()


