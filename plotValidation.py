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

def stackedPlotHistogram(metric, catchmentSize, title, legendLoc = 2, yMax=False, xMin=-1):
  metric[np.isfinite(metric) == False] = -10000
  plotData = []
  xWidth = 0.1
  if xMin == 0.0: xWidth = 0.05
  lims = [0,10**4,25000,50000,10**5,25*10**4,25*10**10]
  for lim in range(1,len(lims)):
    sel1 = catchmentSize/10**6 < lims[lim]
    sel2 = catchmentSize/10**6 > lims[lim-1]
    sel = [x and y for x, y in zip(sel1, sel2)]
    plotData.append(metric[sel])
  ax1 = plt.hist(plotData, bins=np.arange(-1,1.01,xWidth), width = xWidth, stacked=True, color=plt.get_cmap("Blues")(np.linspace(0, 1, 6)), label = ["$<10*10^3$","$<25*10^3$","$<50*10^3$","$<100*10^3$","$<250*10^3$","$\geq250*10^3$"], edgecolor = "none")
  ymax = np.max(ax1[0])*1.02
  #ax1 = plt.legend(prop={'size': 8}, fontsize=8, title="Catchment size ($km^2$)", loc = legendLoc)
  ax1 = plt.legend(prop={'size': 8}, fontsize=8, loc = legendLoc, frameon=False)
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(xMin, 1)
  ax1 = plt.ylim(0, ymax)
  if yMax != False:
    ax1 = plt.ylim(0,yMax*1.02)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def findPlotMax(forecast, validation, xMin = -1):
  xWidth = 0.1
  if xMin == 0: xWidth = 0.05
  ax1 = plt.hist(forecast, bins=np.arange(-1,1.01,xWidth))[0]
  ax2 = plt.hist(validation, bins=np.arange(-1,1.01,xWidth))[0]
  binMax = np.maximum(ax1, ax2)
  plt.clf()
  return(np.max(binMax))


def plotHeatmap(forecast, validation, title, lab1, lab2, nbin = 40, xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5):
  forecast[forecast < xmin] = xmin
  forecast[forecast > xmax] = xmax
  validation[validation < ymin] = ymin
  validation[validation > ymax] = ymax
  ax1 = plt.hist2d(forecast, validation, bins=(nbin,nbin), range=[[xmin,xmax],[ymin,ymax]], norm=colors.PowerNorm(0.2), cmap="Blues")
  heatmap, xedges, yedges = np.histogram2d(forecast, validation, bins=(nbin, nbin), range=[[-1.,1.],[-1.,1.]])
  mid = nbin/2
  ax1 = plt.text(xmin-xmin*0.05,ymin-ymin*0.1,"%.3f" %(np.sum(heatmap[0:mid,0:mid])/np.sum(heatmap)), horizontalalignment='left')
  ax1 = plt.text(xmax-xmax*0.05,ymin-ymin*0.1,"%.3f" %(np.sum(heatmap[mid:nbin,0:mid])/np.sum(heatmap)), horizontalalignment='right')
  ax1 = plt.text(xmin-xmin*0.05,ymax-ymax*0.1,"%.3f" %(np.sum(heatmap[0:mid,mid:nbin])/np.sum(heatmap)), horizontalalignment='left')
  ax1 = plt.text(xmax-xmax*0.05,ymax-ymax*0.1,"%.3f" %(np.sum(heatmap[mid:nbin,mid:nbin])/np.sum(heatmap)), horizontalalignment='right')
  ax1 = plt.title(title)
  ax1 = plt.plot([0.,0.],[ymin,ymax], "black")
  ax1 = plt.plot([xmin,xmax],[0.,0.], "black")
  ax1 = plt.xlim(xmin, xmax)
  ax1 = plt.ylim(ymin, ymax)
  ax1 = plt.xlabel(lab2)
  ax1 = plt.ylabel(lab1)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def plotHistogram(metric, title, yMax = False):
  ax1 = plt.hist(metric, bins=np.arange(-1,1.01,0.1))
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(-1, 1)
  if yMax != False:
    ax1 = plt.ylim(0,yMax*1.02)
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()


def plotCDF(forecast, validation, title, xlims = [-1,1]):
  forecast[np.isfinite(forecast) == False] = -1.01
  forecast[forecast < -1.01] = -1.01
  vals, x1, x2, x3 = cumfreq(forecast, len(forecast))
  ax1 = plt.plot(np.linspace(np.min(forecast), np.max(forecast), len(forecast)), vals/len(forecast), label=str(config.get('Main options', 'RunName')))
  ax1 = plt.plot([np.median(forecast),np.median(forecast)],[0.0,0.5], c="blue",linestyle=':')
  if includeRef:
    validation[np.isfinite(validation) == False] = -1.01
    validation[validation < -1.01] = -1.01
    vals, x1, x2, x3 = cumfreq(validation, len(validation))
    ax2 = plt.plot(np.linspace(np.min(validation), np.max(validation), len(validation)), vals/len(validation), label=str(config.get('Reference options', 'RunName')))
    ax2 = plt.plot([np.median(validation),np.median(validation)],[0.0,0.5], c="red", linestyle=':')
    ax2 = plt.text(np.max([np.median(validation),np.median(forecast)])+0.05, 0.44, "Med=%.3f" %(np.median(validation)), rotation=90, color="red", size=10)
  ax1 = plt.text(np.max([np.median(validation),np.median(forecast)])+0.05, 0.2,"Med=%.3f" %(np.median(forecast)), rotation=90, color="blue", size=10)
  ax1 = plt.legend(prop={'size': 10}, loc=2)
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
  plt.figure(figsize=(12, 6))
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
  plt.figure(figsize=(6, 5))

def plotImprovementMap(data, data2, lons, lats, title, s=5, lab1 = "Shape", lab2="Bias"):
  plt.figure(figsize=(12, 6))
  m = Basemap(projection='mill',lon_0=0, llcrnrlon=-180., llcrnrlat=-59.,
urcrnrlon=180., urcrnrlat=90.)
  m.drawcountries(zorder=0, color="white")
  #m.drawcoastlines(zorder=0, color="black")
  m.fillcontinents(color = 'black',zorder=-1)
  for i in np.arange(3,-1,-1):
    if i == 0 or i == 1:
      sel1 = data > 0.0
    elif i == 2 or i == 3:
      sel1 = data < 0.0
    if i == 0 or i == 2:
      sel2 = data2 > 0.0
    elif i == 1 or i == 3:
      sel2 = data2 < 0.0
    sel = [x and y for x, y in zip(sel1, sel2)]
    x,y = m(lons[sel], lats[sel])

    m.scatter(x,y, c=['#1b9e77','#ffff99','#7570b3','#d95f02'][i], s=s, edgecolors='none')

  leg = plt.legend(("Better %s, %s" %(lab1, lab2),\
  "Better %s, Worse %s" %(lab1, lab2),\
  "Worse %s, Better %s" %(lab1, lab2),\
  "Worse %s, %s" %(lab1, lab2)), loc = 9, bbox_to_anchor=(0.5, 0.0), ncol=4, prop={'size': 8},borderaxespad=0, frameon=False)
  leg.legendHandles[0].set_color('#1b9e77')
  leg.legendHandles[1].set_color('#ffff99')
  leg.legendHandles[2].set_color('#7570b3')
  leg.legendHandles[3].set_color('#d95f02')
  plt.title(title)
  #plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()
  plt.figure(figsize=(6, 5))

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

def plotWaterBalance(data, title):
  totalLen = len(data["year"])
  numMasks = len(data["precipitation"])/len(data["year"])
  if numMasks == 1:
    ax1 = plt.plot(data["year"], data["precipitation"][-totalLen:], label="Precipitation= %.2f" %(np.mean(data["precipitation"][-totalLen:])/1000.))
    ax2 = plt.plot(data["year"], data["actualET"][-totalLen:], label="Evaporation= %.2f" %(np.mean(data["actualET"][-totalLen:])/1000.))
    ax3 = plt.plot(data["year"], data["runoff"][-totalLen:], label="Discharge= %.2f" %(np.mean(data["runoff"][-totalLen:])/1000.))
    ax4 = plt.plot(data["year"], data["totalPotentialGrossDemand"][-totalLen:], label="Demand= %.2f" %(np.mean(data["totalPotentialGrossDemand"][-totalLen:])/1000.))
    ax5 = plt.plot(data["year"], data["storage"][-totalLen:], label="Storage Change= %.2f" %(np.mean(data["storage"][-totalLen:])/1000.))
    ax5 = plt.legend(prop={'size': 10}, loc=2)
    ax1 = plt.title("Global Water Balance %s" %(title))
    ax1 = plt.xlabel("Year")
    ax1 = plt.ylabel("km3")
    ax1 = plt.ylim(np.min(data["storage"]), np.max(data["precipitation"]))
  if numMasks > 1:
    toPlot1 = []
    toPlot2 = []
    toPlot3 = []
    toPlot4 = []
    toPlot5 = []
    for y in range(totalLen):
      selID = np.arange(y,totalLen*numMasks,numMasks)
      toPlot1.append(np.sum(np.array(data["precipitation"])[selID]))
      toPlot2.append(np.sum(np.array(data["actualET"])[selID]))
      toPlot3.append(np.sum(np.array(data["runoff"])[selID]))
      toPlot4.append(np.sum(np.array(data["totalPotentialGrossDemand"])[selID]))
      toPlot5.append(np.sum(np.array(data["storage"])[selID]))
    ax1 = plt.plot(data["year"], toPlot1, label="Precipitation = %.2f" %(np.mean(toPlot1)/1000.))
    ax2 = plt.plot(data["year"], toPlot2, label="Evaporation= %.2f" %(np.mean(toPlot2)/1000.))
    ax3 = plt.plot(data["year"], toPlot3, label="Discharge= %.2f" %(np.mean(toPlot3)/1000.))
    ax4 = plt.plot(data["year"], toPlot4, label="Demand= %.2f" %(np.mean(toPlot4)/1000.))
    ax5 = plt.plot(data["year"], toPlot5, label="Storage Change= %.2f" %(np.mean(toPlot5)/1000.))
    ax5 = plt.legend(prop={'size': 10}, loc=2)
    ax1 = plt.title("Global Water Balance %s" %(title))
    ax1 = plt.xlabel("Year")
    ax1 = plt.ylabel("km3")
    ax1 = plt.ylim(np.min(toPlot5), np.max(toPlot1))
  ax1 = plt.xlim(np.min(data["year"]), np.max(data["year"]))
  ax1 = plt.gcf().set_tight_layout(True)
  pdf.savefig()
  plt.clf()

def getOptions(config, option = "general"):
  dem = str(config.get(option, 'dem')) == str(True)
  demFile = str(config.get(option, 'demFile'))
  demVarName = str(config.get(option, 'demVarName'))
  koeppen = str(config.get(option, 'koeppen')) == str(True)
  koeppenFile = str(config.get(option, 'koeppenFile'))
  koeppenVarName = str(config.get(option, 'koeppenVarName'))
  reportWaterBalance = str(config.get(option, 'reportWaterBalance')) == str(True)
  worldMaps = str(config.get(option, 'worldMaps')) == str(True)
  plotHistogram = str(config.get(option, 'plotHistogram')) == str(True)
  includeRef = config.get('Output options', 'includeReference') == str("True")
  return dem, demFile, demVarName, koeppen, koeppenFile, koeppenVarName, reportWaterBalance, worldMaps, plotHistogram, includeRef

print configFile
config = readConfigFile(configFile)

run1 = str(config.get('Main options', 'RunName'))
run2 = str(config.get('Reference options', 'RunName'))
dem, demFile, demVarName, koeppen, koeppenFile, koeppenVarName, reportWaterBalance, worldMaps, plotHistogram, includeRef = getOptions(config, "general")

output, output2, fullOutput, fullOutput2, waterBalOutput, waterBalOutput2 = pickle.load(open('validationResultsPool_%s_%s.obj' %(run1, run2), 'rb') )

for step in [1,30]:

  uniques = np.zeros((output.shape[0]), dtype=bool)
  lons = []
  lats = []
  for i in range(output.shape[0]):
    unique = True
    if output[i,-1] != 1 and step != 1:
      for x, y in zip(lons == output[i,0], lats == output[i,1]):
        if x and y and output[i,-1] != 1 and step != 1:
          unique = False
    if output[i,-1] == 1 and step == 1:
      for x, y in zip(lons == output[i,0], lats == output[i,1]):
        if x and y and output[i,-1] != 1 and step != 1:
          unique = False        
    if unique:
      lons.append(output[i,0])
      lats.append(output[i,1])
      uniques[i] = True

  sel1 = (np.isnan(output[:,3]+output[:,2]+output[:,4]+output[:,5]+output2[:,2]+output2[:,3]+output2[:,4]+output2[:,5]) == False)
  sel2 = np.sum(output[:,3:-1], axis=1) != 0.0
  if includeRef:
    sel3 = np.sum(output2[:,3:5], axis=1) != 0.0
  else:
    sel3 = sel2
  sel4 = (np.isfinite(output[:,3]+output[:,2]+output[:,4]+output[:,5]+output2[:,2]+output2[:,3]+output2[:,4]+output2[:,5]) == True)
  if step != 1:
    sel5 = output[:,-1] > 1
  else:
    sel5 = output[:,-1] == 1
  sel = [x and y and z and w and v and u for x, y, z, w, v, u in zip(sel1, sel2, sel3, sel4, sel5, uniques)]
  sel5Min = sel # [x and y and z and w and v for x, y, z, w in zip(sel1, sel2, sel4, sel5)]
  print np.sum(sel5Min)
  print step
  if np.sum(sel5Min) > 0:
    if step != 1:
      pdf = PdfPages('plotResults_%s_%s_monthly.pdf' %(run1, run2))
    else:
      pdf = PdfPages('plotResults_%s_%s_daily.pdf' %(run1, run2))
    matplotlib.rcParams.update({'font.size': 12})
    if worldMaps:
      plotWorldMap(output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Correlation with observations (%s)' %(str(config.get('Main options', 'RunName'))))
      if includeRef: plotWorldMap(output2[sel,3], output2[sel,0], output2[sel,1], 'Correlation with observations (%s)' %(str(config.get('Reference options', 'RunName'))))
      if includeRef: plotWorldMap(output[sel,3]-output2[sel,3], output[sel,0], output[sel,1], 'Correlation difference', vmin=-0.5, vmax=0.5)

      plotWorldMap(output[sel5Min,4], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation (%s)' %(str(config.get('Main options', 'RunName'))))
      if includeRef: plotWorldMap(output2[sel,4], output2[sel,0], output2[sel,1], 'Anomaly Correlation (%s)' %(str(config.get('Reference options', 'RunName'))))
      if includeRef: plotWorldMap(output[sel,4]-output2[sel,4], output[sel,0], output[sel,1], 'Anomaly Correlation difference', vmin=-0.5, vmax=0.5)

      plotWorldMap(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,0], output[sel5Min,1], 'Anomaly Correlation - Correlation (%s)' %(str(config.get('Main options', 'RunName'))))
      if includeRef: plotWorldMap(output2[sel,4]-output2[sel,3], output2[sel,0], output2[sel,1], 'Anomaly Correlation - Correlation (%s)' %(str(config.get('Reference options', 'RunName'))))

    if plotHistogram:
      if includeRef: yMax = findPlotMax(output[sel5Min,3], output2[sel,3])
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,3], output[sel5Min,2], "Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output2[sel,3], output2[sel,2], "Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output[sel5Min,3]-output2[sel5Min,3], output2[sel5Min,2], "Correlation difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,4], output2[sel,4])
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,4], output[sel5Min,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output2[sel,4], output2[sel,2], "Anomaly Correlation with observations (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output[sel5Min,4]-output2[sel5Min,4], output2[sel5Min,2], "Anomaly Correlation difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,5], output2[sel,5])
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,5], output[sel5Min,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output2[sel,5], output2[sel,2], "Kling-Gupta Efficiency (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output[sel5Min,5]-output2[sel5Min,5], output2[sel5Min,2], "Kling-Gupta Efficiency difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,6], output2[sel,6])
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,6], output[sel5Min,2], "Kling-Gupta Efficiency Correlation (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output2[sel,6], output2[sel,2], "Kling-Gupta Correlation (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output[sel5Min,6]-output2[sel5Min,6], output2[sel5Min,2], "Kling-Gupta Efficiency Correlation difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,7], output2[sel,7], xMin = 0)
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,7], output[sel5Min,2], "Kling-Gupta Efficiency Alpha (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax, xMin = 0)
      if includeRef: stackedPlotHistogram(output2[sel,7], output2[sel,2], "Kling-Gupta Efficiency Alpha (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax, xMin = 0)
      if includeRef: stackedPlotHistogram(output[sel5Min,7]-output2[sel5Min,7], output2[sel5Min,2], "Kling-Gupta Efficiency Alpha difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,8], output2[sel,8], xMin = 0)
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,8], output[sel5Min,2], "Kling-Gupta Efficiency Beta (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax, xMin = 0)
      if includeRef: stackedPlotHistogram(output2[sel,8], output2[sel,2], "Kling-Gupta Efficiency Beta (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax, xMin = 0)
      if includeRef: stackedPlotHistogram(output[sel5Min,8]-output2[sel5Min,8], output2[sel5Min,2], "Kling-Gupta Efficiency Beta difference %s - %s" %(str(config.get('Main options', 'RunName')), str(config.get('Reference options', 'RunName'))))

      if includeRef: yMax = findPlotMax(output[sel5Min,4]-output[sel5Min,3], output2[sel,4]-output2[sel,3])
      else: yMax = False
      stackedPlotHistogram(output[sel5Min,4]-output[sel5Min,3], output[sel5Min,2], "AC - R (%s)" %(str(config.get('Main options', 'RunName'))), yMax= yMax)
      if includeRef: stackedPlotHistogram(output2[sel,4]-output2[sel,3], output2[sel,2], "AC - R (%s)" %(str(config.get('Reference options', 'RunName'))), yMax= yMax)

    plotCDF(output[sel5Min,3], output2[sel,3], "R")
    plotCDF(output[sel5Min,4], output2[sel,4], "AC")
    plotCDF(output[sel5Min,5], output2[sel,5], "KGE")

    if worldMaps and includeRef:
      plotImprovementMap(output[sel5Min,7]-output2[sel,7], output[sel5Min,8] - output2[sel,8], output2[sel,0], output2[sel,1], "Improv. in KGE alpha and KGE beta" , s=5, lab1 = "Shape", lab2="Bias")

    if includeRef:
      plotHeatmap(output[sel5Min,3]-output2[sel,3], output[sel5Min,4] - output2[sel,4], "Improv. in R and AC", "R", "AC")
      plotHeatmap(output[sel5Min,3]-output2[sel,3], output[sel5Min,5] - output2[sel,5], "Improv. in R and KGE", "R", "KGE")
      plotHeatmap(output[sel5Min,4]-output2[sel,4], output[sel5Min,5] - output2[sel,5], "Improv. in AC and KGE", "AC", "KGE")
      plotHeatmap(output[sel5Min,7]-output2[sel,7], output[sel5Min,8] - output2[sel,8], "Improv. in KGE alpha and KGE beta", "Alpha", "Beta")

    if koeppen:
      koeppenMask = getNCData(koeppenFile, varName = koeppenVarName, xSel = output[:,0], ySel = output[:,1])

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

    if dem:
      demMask = getNCData(demFile, varName = demVarName, xSel = output[:,0], ySel = -output[:,1])

      selLow = demMask <= 1000
      if np.sum(selLow) > 0:
        plotClimateHistogram(selLow, "Below 1000m elevation", output, output2, sel5Min, sel)
        plotClimateCDF(selLow, "Below 1000m elevation", output, output2, sel5Min, sel)

      selHigh = demMask > 1000
      if np.sum(selHigh) > 0:
        plotClimateHistogram(selHigh, "Above 1000m elevation", output, output2, sel5Min, sel)
        plotClimateCDF(selHigh, "Above 1000m elevation", output, output2, sel5Min, sel)

    if reportWaterBalance:
      plotWaterBalance(waterBalOutput, "05min")
      if includeRef: plotWaterBalance(waterBalOutput2, "30min")

    pdf.close()
