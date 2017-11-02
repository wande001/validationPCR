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

configFile = sys.argv[1]

def readConfigFile(configFileName):
  global config
  with open(configFileName) as f:
    sample_config = f.read()
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.readfp(io.BytesIO(sample_config))
  return config

def stackedPlotHistogram(metric, catchmentSize, title, legendLoc = 2, ymax=3500):
  plotData = []
  for lim in [10**4,25000,50000,10**5,25*10**4,25*10**10]:
    sel = catchmentSize/10**6 < lim
    plotData.append(metric[sel])
  ax1 = plt.hist(plotData, bins=np.arange(-1,1.01,0.1), width = 0.1, stacked=True, color=plt.get_cmap("Blues")(np.linspace(0, 1, 6)), label = ["$<10*10^3$","$<25*10^3$","$<50*10^3$","$<100*10^3$","$<250*10^3$","$\geq250*10^3$"], edgecolor = "none")
  ax1 = plt.legend(prop={'size': 10}, title="Catchment size ($km^2$)", loc = legendLoc)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(-1, 1)
  ax1 = plt.ylim(0, ymax)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotHistogram(metric, title):
  ax1 = plt.hist(metric, bins=np.arange(-1,1.01,0.1))
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(-1, 1)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()


def plotCDF(forecast, validation, title, xlims = [-1,1]):
  forecast[forecast < -1.01] = -1.01
  vals, x1, x2, x3 = cumfreq(forecast, len(forecast))
  ax1 = plt.plot(np.linspace(np.min(forecast), np.max(forecast), len(forecast)), vals/len(forecast), label=str(config.get('Main options', 'RunName')))
  validation[validation < -1.01] = -1.01
  vals, x1, x2, x3 = cumfreq(validation, len(validation))
  ax2 = plt.plot(np.linspace(np.min(validation), np.max(validation), len(validation)), vals/len(validation), label=str(config.get('Reference options', 'RunName')))
  ax2 = plt.legend(prop={'size': 10}, loc=2)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("ECDF")
  ax1 = plt.xlim(xlims[0], xlims[1])
  ax1 = plt.ylim(0, 1)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()


def plotScatter(forecast, validation, title):
  ax1 = plt.plot(validation, forecast, "ro", markersize=8)
  ax1 = plt.plot([-100,100], [-100,100])
  ax1 = plt.title(title)
  ax1 = plt.xlabel(str(config.get('Reference options', 'RunName')))
  ax1 = plt.ylabel(str(config.get('Main options', 'RunName')))
  ax1 = plt.xlim(-1, 1)
  ax1 = plt.ylim(-1, 1)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotHexBin(forecast, validation, title):
  forecast[forecast < -1.1] = -1.1
  validation[validation < -1.1] = -1.1
  ax1 = plt.hexbin(validation, forecast, gridsize=20, vmin=1, vmax=20, cmap="OrRd")
  ax1 = plt.plot([-100,100], [-100,100])
  ax1 = plt.title(title)
  ax1 = plt.xlabel(str(config.get('Reference options', 'RunName')))
  ax1 = plt.ylabel(str(config.get('Main options', 'RunName')))
  ax1 = plt.xlim(-1, 1)
  ax1 = plt.ylim(-1, 1)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotWorldMap(data, lons, lats, title, vmin = -1., vmax = 1., s=5):
  plt.figure(figsize=(8, 4))
  m = Basemap(projection='mill',lon_0=0, llcrnrlon=-180., llcrnrlat=-59.,
    urcrnrlon=180., urcrnrlat=90.)
  x,y = m(lons, lats)

  m.drawcountries(zorder=0, color="white")
  #m.drawcoastlines(zorder=0, color="black")
  m.fillcontinents(color = 'black',zorder=-1)

  m.scatter(x,y, c=data, cmap='RdBu', vmin=vmin, vmax=vmax, s=s, edgecolors='none')
  m.colorbar()
  
  plt.title(title)
  plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()
  plt.figure(figsize=(8, 6))


config = readConfigFile(configFile)

output, output2 = pickle.load(open('validationResultsPool_20170905.obj', 'rb') )
sel1 = (np.isnan(output[:,3]+output[:,2]+output[:,4]+output2[:,2]+output2[:,3]+output2[:,4]) == False)
sel2 = np.sum(output[:,3:], axis=1) != 0.0
sel3 = np.sum(output2[:,3:], axis=1) != 0.0
sel = [x and y and z for x, y, z in zip(sel1, sel2, sel3)]
sel5Min = sel

pdf = PdfPages(str(config.get('Output options', 'outputFile')))
matplotlib.rcParams.update({'font.size': 12})


plotWorldMap(output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Correlation with observations (%s)' %(str(config.get('Main options', 'RunName'))))
print len(sel)
print np.percentile(output2[sel,0], [25,50,75])
print np.percentile(output2[sel,1], [25,50,75])
print np.percentile(output2[sel,3], [25,50,75])
plotWorldMap(output2[sel,3], output2[sel,0], output2[sel,1], 'Correlation with observations (%s)' %(str(config.get('Reference options', 'RunName'))))
plotWorldMap(output[sel,3]-output2[sel,3], output[sel,0], output[sel,1], 'Correlation difference 5min - 30min', vmin=-0.5, vmax=0.5)

plotWorldMap(output[sel5Min,4], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation (%s)' %(str(config.get('Main options', 'RunName'))))
plotWorldMap(output2[sel,4], output2[sel,0], output2[sel,1], 'Anomaly Correlation (%s)' %(str(config.get('Reference options', 'RunName'))))
plotWorldMap(output[sel,4]-output2[sel,4], output[sel,0], output[sel,1], 'Anomaly Correlation difference', vmin=-0.5, vmax=0.5)

plotWorldMap(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation - Correlation (%s)' %(str(config.get('Main options', 'RunName'))))

stackedPlotHistogram(output[sel5Min,3], output[sel5Min,2], "Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))), ymax=3600)
stackedPlotHistogram(output2[sel,3], output2[sel,2], "Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))), ymax=3600)

stackedPlotHistogram(output[sel5Min,4], output[sel5Min,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))), ymax=3600)
stackedPlotHistogram(output2[sel,4], output2[sel,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))), ymax=3600)

stackedPlotHistogram(output[sel5Min,5], output[sel5Min,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Main options', 'RunName'))), ymax=2000)
stackedPlotHistogram(output2[sel,5], output2[sel,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Reference options', 'RunName'))), ymax=2000)

stackedPlotHistogram(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,2], "AC - R (%s)" %(str(config.get('Main options', 'RunName'))))

plotCDF(output[sel,3], output2[sel,3], "R")
plotCDF(output[sel,4], output2[sel,4], "AC")
plotCDF(output[sel,5], output2[sel,5], "KGE")

#plotScatter(output[sel,3], output2[sel,3], "R")
#plotScatter(output[sel,4], output2[sel,4], "AC")
#plotScatter(output[sel,5], output2[sel,5], "KGE")

plotHexBin(output[sel,3], output2[sel,3], "R")
plotHexBin(output[sel,4], output2[sel,4], "AC")
plotHexBin(output[sel,5], output2[sel,5], "KGE")

pdf.close()
