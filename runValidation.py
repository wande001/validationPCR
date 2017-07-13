import ConfigParser
import io
import sys
import os
#import pcraster as pcr
from osgeo import gdal
import netCDF4 as nc
import datetime
import numpy as np
from scipy.stats.stats import spearmanr
from scipy.stats import cumfreq
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp


configFile = sys.argv[1]

def readConfigFile(configFileName):
  global config
  with open(configFileName) as f:
    sample_config = f.read()
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.readfp(io.BytesIO(sample_config))
  return config

def getArguments(configFile, reference):
  print reference
  global inputDir, dischargeFileName, summary, full, dischargeDir	
  config = readConfigFile(configFile)
  if reference:
    inputDir = str(config.get('Reference options', 'inputRefDir'))
    dischargeFileName = str(config.get('Reference options', 'dischargeRefFileName'))
  else:
    inputDir = str(config.get('Main options', 'inputDir'))
    dischargeFileName = str(config.get('Main options', 'dischargeFileName'))
  summary = config.get('Main options', 'summaryAnalysis')
  full = config.get('Main options', 'FullAnalysis')
  dischargeDir = str(config.get('GRDC data', 'dischargeDir'))
  numCores = str(config.get('other', 'numCores'))
  return inputDir, dischargeFileName, summary, full, dischargeDir
  
def getCatchmentMap(config):
  catchmentAreaMap = str(config.get('other', 'catchmentAreaMap'))
  cellAreaMap = str(config.get('other', 'cellAreaMap'))
  cellAreaConstant = str(config.get('other', 'cellAreaConstant'))
  routingMap = str(config.get('other', 'routingMap'))
  if os.path.exists(catchmentAreaMap):
    f = gdal.Open(catchmentAreaMap)
    catchmentArea = f.GetRasterBand(1).ReadAsArray()
    return catchmentArea
  elif os.path_exists(cellAreaMap) and os.path_exists(routingMap):
    setclone(routingMap)
    catchmentAcc = pcr.catchmenttotal(cellAreaMap, routingMap)
    catchmentArea = pcr.pcr2numpy(catchmentAcc, -9999.)
    return catchmentArea
  elif os.path_exists(cellAreaMap) and os.path_exists(routingMap):
    setclone(routingMap)
    catchmentAcc = pcr.catchmenttotal(pcr.scalar(cellAreaConstant), routingMap)
    catchmentArea = pcr.pcr2numpy(catchmentAcc, -9999.)
    return catchmentArea
  else:
	print "No valid catchment size and or routing map provided"
	return False
	
def readModelProps(ncFile):    
  # Get netCDF file and variable name:
  f = nc.Dataset(ncFile)      
  try:
    lat = f.variables['lat'][:]
    lon = f.variables['lon'][:]
  except:
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
  nctime = f.variables['time'][:]
  nctimeUnit = f.variables['time'].units
  nctimeCalendar = f.variables['time'].calendar
  timeVar = nc.num2date(nctime,units = nctimeUnit, calendar = nctimeCalendar)
  return lon, lat, timeVar[0], timeVar[-1]

def getWindowSize(config):
  windowSize = int(config.get('other', 'windowSizeForMatching'))
  return windowSize

def getAreaMisMatch(config):
  misMatch = float(config.get('other', 'areaMisMatchTolerance'))
  return misMatch

def getCores(config):
  numCores = int(config.get('other', 'numCores'))
  return numCores
	
def readObservationsProps(fileName):
  f = open(fileName)
  lines=f.readlines()
  obsLon = float(lines[13][25:-1])
  obsLat = float(lines[12][25:-1])
  obsCatchArea = float(lines[14][25:-1])*10**6
  obsStart = datetime.datetime(int(lines[41][0:4]), int(lines[41][5:7]), int(lines[41][8:10]))
  obsEnd = datetime.datetime(int(lines[-1][0:4]), int(lines[-1][5:7]), int(lines[-1][8:10]))
  return obsLon, obsLat, obsCatchArea, obsStart, obsEnd

def findMin(obs, mod):
  cellSize = np.abs(mod[0] - mod[1])
  dif = np.abs(obs - mod)
  if np.min(dif) <= cellSize:
    return np.argmin(dif)
  else:
	return -999.

def getModArea(xSel, ySel, windowSize, modCatchArea):
  dy, dx = np.shape(modCatchArea)
  xLower = xSel - windowSize
  xUpper = xSel + windowSize + 1
  yLower = max(ySel - windowSize,0)
  yUpper = min(ySel + windowSize + 1,dy)
  if xLower >= 0 and xUpper <= dx:
    areaSelection = modCatchArea[yLower:yUpper,xLower:xUpper].flatten()
    xSelection = np.tile(range(xLower, xUpper), len(range(yLower, yUpper)))
  elif xLower < 0:
    areaSelection = np.concatenate([modCatchArea[range(yLower,yUpper),xLower:],modCatchArea[range(yLower,yUpper),:xUpper]]).flatten()
    xSelection = np.tile(np.concatenate([range(xLower+dx, dx), range(xUpper)]),len(range(yLower, yUpper)))
  elif xUpper > dx:
    areaSelection = np.concatenate([modCatchArea[range(yLower,yUpper),xLower:], modCatchArea[range(yLower,yUpper),:(xUpper - dx)]]).flatten()
    xSelection = np.tile(np.concatenate([range(xLower,dx), range(xUpper - dx)]),len(range(yLower, yUpper)))
  ySelection = np.repeat(range(yLower, yUpper), len(range(xLower, xUpper)))
  return areaSelection, xSelection, ySelection
    
def matchLocation(obsLon, obsLat, obsCatchArea, modLon, modLat, modCatchArea, windowSize, misMatch):
  misMatch = getAreaMisMatch(config)
  xSel = findMin(obsLon,modLon)
  ySel = findMin(obsLat,modLat)
  if xSel != -999. and ySel != -999.:
    areaSelection, xSelection, ySelection = getModArea(xSel, ySel, windowSize, modCatchArea)
    if np.min(np.abs(obsCatchArea - modCatchArea)) != modCatchArea[ySel, xSel] and np.min(np.abs(obsCatchArea - areaSelection))/obsCatchArea < misMatch:
      locSel = np.argmin(np.abs(obsCatchArea - areaSelection))
      xSel = xSelection[locSel]
      ySel = ySelection[locSel]
      return xSel, ySel
    else:
      return -999., -999.
  else:
    return -999., -999.

def getModelData(ncFile, xSel, ySel, varName = 'discharge'):    
  # Get netCDF file and variable name:
  f = nc.Dataset(ncFile)
  data = f.variables[varName][:,ySel,xSel]
  f.close()
  return data
  
def getObservationData(fileName):
  f = open(fileName)
  lines=f.readlines()
  numLines = range(41, len(lines))
  obsData = np.zeros(len(numLines)) - 999.
  lineCount = 0
  for line in numLines:
    obsData[lineCount] = lines[line].split(";")[3]
    lineCount += 1
  return obsData

def matchSeries(obs, mod, obsStart, obsEnd, modStart, modEnd):
  obsStartShift = 0
  modStartShift = 0
  obsEndShift = 0
  modEndShift = 0
  if obsStart < modStart:
	obsStartShift += (modStart- obsStart).days
  elif obsStart > modStart:
	modStartShift += (obsStart - modStart).days
  if obsEnd > modEnd:
	obsEndShift += (modEnd - obsEnd).days
  elif obsEnd < modEnd:
	modEndShift += (obsEnd - modEnd).days
  if obsStart + datetime.timedelta(days =obsStartShift) < obsEnd and obsEnd + datetime.timedelta(days =obsEndShift) > obsStart and modStart + datetime.timedelta(days =modStartShift) < modEnd and modEnd + datetime.timedelta(days =modEndShift) > modStart:
    if obsEndShift == 0:
      obsOut = obs[obsStartShift:]
    else:
      obsOut = obs[obsStartShift:obsEndShift]
    if modEndShift == 0:
      modOut = mod[modStartShift:]
    else:
      modOut = mod[modStartShift:modEndShift]
    obsOut[obsOut< 0.0] = np.nan
    modOut[modOut< 0.0] = np.nan
    return obsOut, modOut
  else:
	return None
	
def nashSutcliffe(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  MSE = np.sum((obs[sel]-mod[sel])**2)
  MSEclim = np.sum((obs[sel] - np.mean(obs[sel]))**2)
  NSE = 1-MSE/MSEclim
  return np.maximum(NSE,-100.)

def rmse(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  out = np.mean((obs[sel]-mod[sel])**2)**0.5
  return out

def bias(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  out = np.mean(obs[sel]) -np.mean(mod[sel])
  return out, len(sel)
  
def kge(obs, mod):
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  cc = np.corrcoef(obs[sel],mod[sel])[0,1]
  alpha = np.std(obs[sel])/np.std(mod[sel])
  beta  = np.sum(obs[sel])/np.sum(mod[sel])
  kge   = 1- np.sqrt( (cc-1)**2 + (alpha-1)**2 + (beta-1)**2 )
  return kge

def calculateMetrics(obs, mod, obsStart, obsEnd, modStart, modEnd):
  obs, mod = matchSeries(obs, mod, obsStart, obsEnd, modStart, modEnd)
  if len(obs) > 0 and len(mod) > 0:
    R = spearmanr(obs, mod)[0]
    NS = nashSutcliffe(obs, mod)
    RMSE = rmse(obs, mod)
    Bias, numPoints = bias(obs, mod)
    KGE = kge(obs, mod)
    return R, KGE, NS, RMSE, Bias, numPoints
  else:
    return None

def stackedPlotHistogram(metric, catchmentSize, title):
  plotData = []
  for lim in [10**4,25000,50000,10**5,25*10**4,25*10**10]:
    sel = catchmentSize/10**6 < lim
    plotData.append(metric[sel])
  ax1 = plt.hist(plotData, bins=np.arange(-1,1,0.1), stacked=True, color=plt.get_cmap("Blues")(np.linspace(0, 1, 6)), label = ["$<10*10^3$","$<25*10^3$","$<50*10^3$","$<100*10^3$","$<250*10^3$","$\geq250*10^3$"])
  ax1 = plt.legend(prop={'size': 10}, title="Catchment size ($km^2$)")
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(-1, 1)
  pdf.savefig()
  plt.clf()

def plotHistogram(metric, title):
  ax1 = plt.hist(metric, bins=np.arange(-1,1,0.1))
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("Frequency")
  ax1 = plt.xlim(-1, 1)
  pdf.savefig()
  plt.clf()


def plotCDF(forecast, validation, title, xlims = [-1,1]):
  vals, x1, x2, x3 = cumfreq(forecast, len(forecast))
  ax1 = plt.plot(np.linspace(np.min(forecast), np.max(forecast), len(forecast)), vals/len(forecast), label='Simulation')
  vals, x1, x2, x3 = cumfreq(validation, len(validation))
  ax2 = plt.plot(np.linspace(np.min(validation), np.max(validation), len(validation)), vals/len(validation), label='Reference')
  ax2 = plt.legend(prop={'size': 10}) 
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Value")
  ax1 = plt.ylabel("ECDF")
  ax1 = plt.xlim(xlims[0], xlims[1])
  ax1 = plt.ylim(0, 1)
  pdf.savefig()
  plt.clf()


def plotScatter(forecast, validation, title):
  ax1 = plt.plot(validation, forecast, "ro", markersize=8)
  ax1 = plt.plot([-100,100], [-100,100])
  ax1 = plt.title(title)
  ax1 = plt.xlabel("Reference")
  ax1 = plt.ylabel("Simulation")
  ax1 = plt.xlim(-1, 1)
  ax1 = plt.ylim(-1, 1)
  pdf.savefig()
  plt.clf()


def getGlobalProperties(configFile, reference):
  global modLon, modLat, modStart, modEnd, inputDir, modCatchArea, windowSize, misMatch, locations, numCores
  getArguments(configFile, reference)
  locations = os.listdir(dischargeDir)
  modLon, modLat, modStart, modEnd = readModelProps("%s/%s" %(inputDir, dischargeFileName))
  modCatchArea = getCatchmentMap(config)
  windowSize = getWindowSize(config)
  misMatch = getAreaMisMatch(config)
  numCores = getCores(config)
  return modLon, modLat, modStart, modEnd, inputDir, modCatchArea, windowSize, misMatch, locations, numCores
  
def extractLocation(location,inputDir, dischargeFileName):
  output = np.zeros((9))
  if location[-4:] == ".day":
    obsLon, obsLat, obsCatchArea, obsStart, obsEnd = readObservationsProps("%s/%s" %(dischargeDir, location))
    output[0:2] = obsLon, obsLat
    output[2] = obsCatchArea
    xSel, ySel = matchLocation(obsLon, obsLat, obsCatchArea, modLon, modLat, modCatchArea, windowSize, misMatch)
    if xSel != -999. and ySel != -999.:
      modValues = getModelData("%s/%s" %(inputDir, dischargeFileName), xSel, ySel)
      obsValues = getObservationData("%s/%s" %(dischargeDir, location))
      output[3:] = calculateMetrics(obsValues, modValues, obsStart, obsEnd, modStart, modEnd)
  return output

def f(location):
  print location


getGlobalProperties(configFile, reference=False)
print inputDir, dischargeFileName
pool = mp.Pool(processes=numCores)

results = [pool.apply_async(extractLocation,args=(location,inputDir, dischargeFileName)) for location in locations]
outputList = [p.get() for p in results]
output = np.array(outputList)

getGlobalProperties(configFile, reference=True)
print inputDir, dischargeFileName

results2 = [pool.apply_async(extractLocation,args=(location,inputDir, dischargeFileName)) for location in locations]
outputList2 = [p.get() for p in results2]
output2 = np.array(outputList2)

pdf = PdfPages(str(config.get('Output options', 'outputFile')))
matplotlib.rcParams.update({'font.size': 12})
 
m = Basemap(projection='mill',lon_0=0, llcrnrlon=-180., llcrnrlat=-59.,
    urcrnrlon=180., urcrnrlat=90.)

lons = output[:,0]
lats = output[:,1]
x,y = m(lons, lats)

m.drawcountries(zorder=0)
m.fillcontinents(color = 'coral',zorder=-1)
   
m.scatter(x,y, c=output[:,3])
m.colorbar()

pdf.savefig()
plt.clf()

stackedPlotHistogram(output[:,3], output[:,2], "R")
stackedPlotHistogram(output[:,4], output[:,2], "KGE")

plotCDF(output[:,3], output2[:,3], "R")
plotCDF(output[:,4], output2[:,4], "KGE")

plotScatter(output[:,3], output2[:,3], "R")
plotScatter(output[:,4], output2[:,4], "KGE")

pdf.close()

