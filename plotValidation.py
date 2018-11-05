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

def readConfigFile(configFileName):
  global config
  with open(configFileName) as f:
    sample_config = f.read()
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.readfp(io.BytesIO(sample_config))
  return config

def stackedPlotHistogram(metric, catchmentSize, title, legendLoc = 2):
  metric[np.isfinite(metric) == False] = -10000
  plotData = []
  lims = [0,10**4,25000,50000,10**5,25*10**4,25*10**10]
  for lim in range(1,len(lims)):
    sel1 = catchmentSize/10**6 < lims[lim]
    sel2 = catchmentSize/10**6 > lims[lim-1]
    sel = [x and y for x, y in zip(sel1, sel2)]
    plotData.append(metric[sel])
  ax1 = plt.hist(plotData, bins=np.arange(-1,1.01,0.1), width = 0.1, stacked=True, color=plt.get_cmap("Blues")(np.linspace(0, 1, 6)), label = ["$<10*10^3$","$<25*10^3$","$<50*10^3$","$<100*10^3$","$<250*10^3$","$\geq250*10^3$"], edgecolor = "none")
  ymax = np.max(ax1[0])*1.02
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
  forecast[np.isfinite(forecast) == False] = -1.01
  forecast[forecast < -1.01] = -1.01
  vals, x1, x2, x3 = cumfreq(forecast, len(forecast))
  ax1 = plt.plot(np.linspace(np.min(forecast), np.max(forecast), len(forecast)), vals/len(forecast), label=str(config.get('Main options', 'RunName')))
  validation[np.isfinite(validation) == False] = -1.01
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
  
def getNCData(ncFile, varName, xSel, ySel):    
  # Get netCDF file and variable name:
  data = np.zeros((len(xSel)))
  f = nc.Dataset(ncFile)
  for x,y,count in zip(xSel, ySel, range(len(xSel))):
    lon = int(np.minimum(np.floor(x * 2.0+360),719))
    lat = int(np.minimum(np.floor(y * -2.0+180), 359))
    data[count] = f.variables[varName][lat, lon]
  f.close()
  return data

def plotClimateHistogram(climateSel, title, output, output2, sel5Min, sel, ymax=125):
  selClim = [x and y and z for x, y, z in zip(sel5Min, climateSel, sel)]
  if len(output[selClim,2]) != 0:
    stackedPlotHistogram(output[selClim,5]-output2[selClim,5], output[selClim,2], "Kling-Gupta Efficiency %s (%s - %s)" %(title, run1, run2))

def plotClimateCDF(climateSel, title, output, output2, sel5Min, sel):
  selClim = [x and y and z for x, y, z in zip(sel5Min, climateSel, sel)]
  if len(output[selClim,5]) > 10:
    plotCDF(output[selClim,5], output2[selClim,5], "Kling-Gupta Efficiency %s" %(title))


config = readConfigFile(configFile)

run1 = str(config.get('Main options', 'RunName'))
run2 = str(config.get('Reference options', 'RunName'))

output, output2 = pickle.load(open('validationResultsPool_%s_%s.obj' %(run1, run2), 'rb') )

koeppenMask = getNCData("/Users/niko/Scripts/Misc/Koeppen-Geiger-Classification-Reclassfied.nc", varName = 'Koeppen_Classification', xSel = output[:,0], ySel = output[:,1])
demMask = getNCData("/Users/niko/Scripts/PCR-GLOBWB/input30min/routing/elev.nc", varName = 'Band1', xSel = output[:,0], ySel = -output[:,1])

for step in [1,30]:

  sel1 = (np.isnan(output[:,3]+output[:,2]+output[:,4]+output[:,5]+output2[:,2]+output2[:,3]+output2[:,4]+output2[:,5]) == False)
  sel2 = np.sum(output[:,3:-1], axis=1) != 0.0
  sel3 = np.sum(output2[:,3:5], axis=1) != 0.0
  sel4 = (np.isfinite(output[:,3]+output[:,2]+output[:,4]+output[:,5]+output2[:,2]+output2[:,3]+output2[:,4]+output2[:,5]) == True)
  if step != 1:
    sel5 = output[:,-1] > 1
  else:
    sel5 = output[:,-1] == 1
  sel = [x and y and z and w and v for x, y, z, w, v in zip(sel1, sel2, sel3, sel4, sel5)]
  sel5Min = sel # [x and y and z and w and v for x, y, z, w in zip(sel1, sel2, sel4, sel5)]

  print np.sum(sel5Min)

  if step != 1:
    pdf = PdfPages('plotResults_%s_%s_monthly.pdf' %(run1, run2))
  else:
    pdf = PdfPages('plotResults_%s_%s_daily.pdf' %(run1, run2))
  matplotlib.rcParams.update({'font.size': 12})

  plotWorldMap(output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Correlation with observations (%s)' %(str(config.get('Main options', 'RunName'))))
  plotWorldMap(output2[sel,3], output2[sel,0], output2[sel,1], 'Correlation with observations (%s)' %(str(config.get('Reference options', 'RunName'))))
  plotWorldMap(output[sel,3]-output2[sel,3], output[sel,0], output[sel,1], 'Correlation difference 5min - 30min', vmin=-0.5, vmax=0.5)

  plotWorldMap(output[sel5Min,4], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation (%s)' %(str(config.get('Main options', 'RunName'))))
  plotWorldMap(output2[sel,4], output2[sel,0], output2[sel,1], 'Anomaly Correlation (%s)' %(str(config.get('Reference options', 'RunName'))))
  plotWorldMap(output[sel,4]-output2[sel,4], output[sel,0], output[sel,1], 'Anomaly Correlation difference', vmin=-0.5, vmax=0.5)

  plotWorldMap(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation - Correlation (%s)' %(str(config.get('Main options', 'RunName'))))
  plotWorldMap(output2[sel,4]-output2[sel,3], output2[sel,0], output2[sel,1], 'Anomaly Correlation - Correlation (%s)' %(str(config.get('Reference options', 'RunName'))))

  stackedPlotHistogram(output[sel5Min,3], output[sel5Min,2], "Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))))
  stackedPlotHistogram(output2[sel,3], output2[sel,2], "Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))))
  stackedPlotHistogram(output[sel5Min,3]-output2[sel5Min,3], output2[sel5Min,2], "Correlation difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

  stackedPlotHistogram(output[sel5Min,4], output[sel5Min,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))))
  stackedPlotHistogram(output2[sel,4], output2[sel,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))))
  stackedPlotHistogram(output[sel5Min,4]-output2[sel5Min,4], output2[sel5Min,2], "Anomaly Correlation difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

  stackedPlotHistogram(output[sel5Min,5], output[sel5Min,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Main options', 'RunName'))))
  stackedPlotHistogram(output2[sel,5], output2[sel,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Reference options', 'RunName'))))
  stackedPlotHistogram(output[sel5Min,5]-output2[sel5Min,5], output2[sel5Min,2], "Kling-Gupta Efficiency difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

  stackedPlotHistogram(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,2], "AC - R (%s)" %(str(config.get('Main options', 'RunName'))))
  stackedPlotHistogram(output2[sel,4]-output2[sel,3], output2[sel,2], "AC - R (%s)" %(str(config.get('Reference options', 'RunName'))))

  plotCDF(output[sel,3], output2[sel,3], "R")
  plotCDF(output[sel,4], output2[sel,4], "AC")
  plotCDF(output[sel,5], output2[sel,5], "KGE")

  selA = koeppenMask <= 4
  if np.sum(selA) > 0:
    plotClimateHistogram(selA, "Tropical Climate", output, output2, sel5Min, sel)
    plotClimateCDF(selA, "Tropical Climate", output, output2, sel5Min, sel)

  selB = [x and y for x, y in zip(koeppenMask >= 5, koeppenMask <= 8)]
  if np.sum(selB) > 0:
    plotClimateHistogram(selB, "Desert Climate", output, output2, sel5Min, sel)
    plotClimateCDF(selB, "Desert Climate", output, output2, sel5Min, sel)

  selC = [x and y for x, y in zip(koeppenMask >= 9, koeppenMask <= 17)]
  if np.sum(selC) > 0:
    plotClimateHistogram(selC, "Temperate Climate", output, output2, sel5Min, sel)
    plotClimateCDF(selC, "Temperate Climate", output, output2, sel5Min, sel)

  selD = [x and y for x, y in zip(koeppenMask >= 18, koeppenMask <= 28)]
  if np.sum(selD) > 0:
    plotClimateHistogram(selD, "Continental Climate", output, output2, sel5Min, sel)
    plotClimateCDF(selD, "Continental Climate", output, output2, sel5Min, sel)

  selLow = demMask <= 1000
  if np.sum(selLow) > 0:
    plotClimateHistogram(selLow, "Below 1000m elevation", output, output2, sel5Min, sel)
    plotClimateCDF(selLow, "Below 1000m elevation", output, output2, sel5Min, sel)

  selHigh = demMask > 1000
  if np.sum(selHigh) > 0:
    plotClimateHistogram(selHigh, "Above 1000m elevation", output, output2, sel5Min, sel)
    plotClimateCDF(selHigh, "Above 1000m elevation", output, output2, sel5Min, sel)

  pdf.close()
