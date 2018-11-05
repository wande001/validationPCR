import ConfigParser
import io
import sys
import os
import pickle
#import pcraster as pcr
from osgeo import gdal
import netCDF4 as nc
import datetime
from calendar import monthrange
import numpy as np
from scipy.stats.stats import spearmanr
from scipy.stats import cumfreq
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp
import pcraster as pcr

configFile = sys.argv[1]

def readConfigFile(configFileName):
  global config
  with open(configFileName) as f:
    sample_config = f.read()
  config = ConfigParser.RawConfigParser(allow_no_value=True)
  config.readfp(io.BytesIO(sample_config))
  return config

def getArguments(configFile, reference):
  global inputDir, dischargeFileName, summary, full, dischargeDir, runName, refName
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
  numCores = str(config.get('general', 'numCores'))
  runName = str(config.get('Main options', 'RunName'))
  refName = str(config.get('Reference options', 'RunName'))
  return inputDir, dischargeFileName, summary, full, dischargeDir,runName,refName
  
def getCatchmentMap(config, modLon, modLat, option = "otherReference"):
  catchmentAreaMap = str(config.get(option, 'catchmentAreaMap'))
  cellAreaMap = str(config.get(option, 'cellAreaMap'))
  cellAreaConstant = str(config.get(option, 'cellAreaConstant'))
  routingMap = str(config.get(option, 'routingMap'))
  dxMod = modLon[1] - modLon[0]
  dyMod = modLat[1] - modLat[0]
  if os.path.exists(catchmentAreaMap):
    f = gdal.Open(catchmentAreaMap)
    dx = f.GetGeoTransform()[1]
    dy = f.GetGeoTransform()[5]
    catchmentArea = f.GetRasterBand(1).ReadAsArray()
    if (dxMod > 0 and dx < 0) or (dxMod < 0 and dx > 0):
      catchmentArea = catchmentArea[:,::-1]
    if (dyMod > 0 and dy < 0) or (dyMod < 0 and dy > 0):
      catchmentArea = catchmentArea[::-1,:]    
    return catchmentArea
  elif os.path.exists(cellAreaMap) and os.path.exists(routingMap):
    pcr.setclone(routingMap)
    catchmentAcc = pcr.catchmenttotal(cellAreaMap, routingMap)
    catchmentArea = pcr.pcr2numpy(catchmentAcc, -9999.)
    return catchmentArea
  elif os.path.exists(cellAreaMap) and os.path.exists(routingMap):
    pcr.setclone(routingMap)
    catchmentAcc = pcr.catchmenttotal(pcr.scalar(cellAreaConstant), routingMap)
    catchmentArea = pcr.pcr2numpy(catchmentAcc, -9999.)
    return catchmentArea
  else:
	print "No valid catchment size and or routing map provided"
	return False
	
def readModelProps(ncFile):    
  # Get netCDF file and variable name:
  print ncFile
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
  modStep = (timeVar[1] - timeVar[0]).days
  return lon, lat, timeVar[0], timeVar[-1], modStep

def getWindowSize(config, option = "general"):
  windowSize = int(config.get(option, 'windowSizeForMatching'))
  return windowSize

def getAreaMisMatch(config, option = "general"):
  misMatch = float(config.get(option, 'areaMisMatchTolerance'))
  return misMatch

def getCores(config, option = "general"):
  numCores = int(config.get(option, 'numCores'))
  return numCores
	
def readObservationsProps(fileName):
  f = open(fileName)
  lines=f.readlines()
  obsLon = float(lines[13][25:-1])
  obsLat = float(lines[12][25:-1])
  obsCatchArea = float(lines[14][25:-1])*10**6
  if lines[41][8:10] == "00" and lines[41][5:7] != "12":
    obsStart = datetime.datetime(int(lines[41][0:4]), int(lines[41][5:7])+1, 1) - datetime.timedelta(days=1)
  elif lines[41][8:10] == "00" and lines[41][5:7] == "12":
    obsStart = datetime.datetime(int(lines[41][0:4])+1, 1, 1) - datetime.timedelta(days=1)
  else:
    obsStart = datetime.datetime(int(lines[41][0:4]), int(lines[41][5:7]), int(lines[41][8:10]))
  if lines[42][8:10] == "00" and lines[42][5:7] != "12":
    obsSecond = datetime.datetime(int(lines[42][0:4]), int(lines[42][5:7])+1, 1) - datetime.timedelta(days=1)
  elif lines[42][8:10] == "00" and lines[42][5:7] == "12":
    obsSecond = datetime.datetime(int(lines[42][0:4])+1, 1, 1) - datetime.timedelta(days=1)
  else:
    obsSecond = datetime.datetime(int(lines[42][0:4]), int(lines[42][5:7]), int(lines[42][8:10]))
  if lines[-1][8:10] == "00" and lines[-1][5:7] != "12":
    obsEnd = datetime.datetime(int(lines[-1][0:4]), int(lines[-1][5:7])+1, int(1)) - datetime.timedelta(days=1)
  elif lines[-1][8:10] == "00" and lines[-1][5:7] == "12":
    obsEnd = datetime.datetime(int(lines[-1][0:4])+1, 1, 1) - datetime.timedelta(days=1)
  else:
    obsEnd = datetime.datetime(int(lines[-1][0:4]), int(lines[-1][5:7]), int(lines[-1][8:10]))
  obsStep = (obsSecond - obsStart).days
  return obsLon, obsLat, obsCatchArea, obsStart, obsEnd, obsStep

def readObservationsPropsLatin(fileName):
  f = open(fileName)
  lines=f.readlines()
  obsLon = float(lines[1][6:-1])
  obsLat = float(lines[2][6:-1])
  obsCatchArea = float(lines[3][7:-1])*10**6
  month, day, year = lines[5].split("\t")[0].split("/")
  year = int(year) + 1900
  if year < datetime.datetime.now().year-100: year += 100
  obsStart = datetime.datetime(int(year), int(month), int(day))
  month, day, year = lines[6].split("\t")[0].split("/")
  year = int(year) + 1900
  if year < datetime.datetime.now().year-100: year += 100
  obsSecond = datetime.datetime(int(year), int(month), int(day))
  month, day, year = lines[-1].split("\t")[0].split("/")
  year = int(year) + 1900
  if year < datetime.datetime.now().year-100: year += 100
  obsEnd = datetime.datetime(int(year), int(month), int(day))
  obsStep = (obsSecond - obsStart).days
  return obsLon, obsLat, obsCatchArea, obsStart, obsEnd, obsStep

def aggregateToMonth(data, startDate, endDate):
  deltaDays = (endDate - startDate).days
  date_list = [startDate + datetime.timedelta(days=x) for x in range(0, deltaDays+1)]
  tempMonth = startDate.month
  dayCount = 0
  tempStart = 0
  monthlyData = []
  while dayCount <= deltaDays:
    if date_list[dayCount].month != tempMonth:
      monthlyData.append(np.mean(data[tempStart:dayCount]))
      tempStart = dayCount
      tempMonth = date_list[dayCount].month
    dayCount += 1
  monthlyData.append(np.mean(data[tempStart:dayCount]))
  startDate = datetime.datetime(startDate.year, startDate.month, 1)
  endDate = datetime.datetime(endDate.year, endDate.month, monthrange(endDate.year, endDate.month)[1])
  return np.array(monthlyData), startDate, endDate

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
    areaSelection = np.concatenate([modCatchArea[range(yLower,yUpper),(xLower-1):].flatten(),modCatchArea[range(yLower,yUpper),:xUpper].flatten()]).flatten()
    xSelection = np.tile(np.concatenate([range((xLower-1)+dx, dx), range(xUpper)]),len(range(yLower, yUpper)))
  elif xUpper > dx:
    areaSelection = np.concatenate([modCatchArea[range(yLower,yUpper),(xLower-1):].flatten(), modCatchArea[range(yLower,yUpper),:(xUpper - dx + 1)].flatten()]).flatten()
    xSelection = np.tile(np.concatenate([range(xLower,dx), range(xUpper - dx + 1)]),len(range(yLower, yUpper)))
  ySelection = np.repeat(range(yLower, yUpper), len(range(xLower, xUpper)))
  return areaSelection, xSelection, ySelection
    
def matchLocation(obsLon, obsLat, obsCatchArea, modLon, modLat, modCatchArea, windowSize, misMatch):
  misMatch = getAreaMisMatch(config)
  xSel = findMin(obsLon,modLon)
  ySel = findMin(obsLat,modLat)
  if xSel != -999. and ySel != -999.:
    areaSelection, xSelection, ySelection = getModArea(xSel, ySel, windowSize, modCatchArea)
    if obsCatchArea > 0.0:
      if np.min(np.abs(obsCatchArea - modCatchArea)) != modCatchArea[ySel, xSel] and np.min(np.abs(obsCatchArea - areaSelection))/obsCatchArea < misMatch:
        locSel = np.argmin(np.abs(obsCatchArea - areaSelection))
        xSel = xSelection[locSel]
        ySel = ySelection[locSel]
        return xSel, ySel
      else:
        return -999., -999.
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

def getObservationDataLatin(fileName):
  f = open(fileName)
  lines=f.readlines()
  numLines = range(5, len(lines))
  obsData = np.zeros(len(numLines)) - 999.
  lineCount = 0
  for line in numLines:
    obsData[lineCount] = lines[line].split("\t")[1]
    lineCount += 1
  return obsData

def monthDelta(d1, d2):
  delta = 0
  run = True
  d1 = datetime.datetime(d1.year, d1.month, 1)
  while run:
    mdays = monthrange(d1.year, d1.month)[1]
    d1 += datetime.timedelta(days=mdays)
    if d1 <= d2:
      delta += 1
    else:
      run = False
  return delta

def addMonth(date, months):
  deltaMonth = date.month + months
  year = date.year + (deltaMonth-1)/12
  month = (deltaMonth - ((deltaMonth-1)/12)*12)
  return datetime.datetime(year, month, date.day)

def matchSeriesMonth(obs, mod, obsStart, obsEnd, modStart, modEnd):
  obsStartShift = 0
  modStartShift = 0
  obsEndShift = 0
  modEndShift = 0
  if obsStart < modStart:
        obsStartShift += monthDelta(obsStart, modStart)
  elif obsStart > modStart:
        modStartShift += monthDelta(modStart, obsStart)
  if obsEnd > modEnd:
        obsEndShift -= monthDelta(modEnd, obsEnd)
  elif obsEnd < modEnd:
        modEndShift -= monthDelta(obsEnd, modEnd)
  if addMonth(obsStart, obsStartShift) < obsEnd and addMonth(obsEnd, obsEndShift) > obsStart and addMonth(modStart, modStartShift) < modEnd and addMonth(modEnd, modEndShift) > modStart:
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
        return [], []

def matchSeriesDay(obs, mod, obsStart, obsEnd, modStart, modEnd):
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
	return [], []
	
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

def anomalyCorrelation(obs, mod, timeScale = "month"):
  normObs = normalizeMonth(obs)
  normMod = normalizeMonth(mod)
  return spearmanr(normObs, normMod)[0]

def normalizeMonth(data):
  seasonCycle = np.zeros(len(data))
  for m in range(12):
    monthSelection = np.arange(m,len(data), 12)
    seasonCycle[monthSelection] = np.mean(data[monthSelection])
  return data - seasonCycle

def calculateMetrics(obs, mod, obsStart, obsEnd, modStart, modEnd, obsStep, modStep):
  if obsStep > 1 or modStep > 1: obs, mod = matchSeriesMonth(obs, mod, obsStart, obsEnd, modStart, modEnd)
  if obsStep <= 1 and modStep <= 1: obs, mod = matchSeriesDay(obs, mod, obsStart, obsEnd, modStart, modEnd)
  #print len(obs), len(mod)
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  if len(obs) > 12 and len(mod) > 12:
    obsSel = np.isnan(obs) == False
    modSel = np.isnan(mod) == False
    sel = obsSel & modSel
    if len(obs[sel]) > 12 and len(mod[sel]) > 12:
      R = spearmanr(obs, mod)[0]
      NS = nashSutcliffe(obs, mod)
      RMSE = rmse(obs, mod)
      Bias, numPoints = bias(obs, mod)
      KGE = kge(obs, mod)
      AC = anomalyCorrelation(obs, mod)
      print R, AC, KGE, NS, RMSE, Bias, numPoints
      return R, AC, KGE, NS, RMSE, Bias, numPoints
    else:
      return np.zeros((7))
  else:
    return np.zeros((7))

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
  global modLon, modLat, modStart, modEnd, modStep, inputDir, modCatchArea, windowSize, misMatch, locations, numCores
  getArguments(configFile, reference)
  locations = os.listdir(dischargeDir)
  modLon, modLat, modStart, modEnd, modStep = readModelProps("%s/%s" %(inputDir, dischargeFileName))
  if reference:
    option = "otherReference"
  else:
    option = "other"
  modCatchArea = getCatchmentMap(config, modLon, modLat, option)
  windowSize = getWindowSize(config, "general")
  misMatch = getAreaMisMatch(config, "general")
  numCores = getCores(config, "general")
  return modLon, modLat, modStart, modEnd, modStep, inputDir, modCatchArea, windowSize, misMatch, locations, numCores
  
def extractLocation(location,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep):
  #print location/float(len(locations)), locations[location]
  location = locations[location]
  output = np.zeros((11))
  #if location[-4:] == ".mon":
  if location[-4:] == ".mon" or location[-4:] == ".day":
    obsLon, obsLat, obsCatchArea, obsStart, obsEnd, obsStep = readObservationsProps("%s/%s" %(dischargeDir, location))
    output[0:2] = obsLon, obsLat
    output[2] = obsCatchArea
    output[10] = float(obsStep)
    xSel, ySel = matchLocation(obsLon, obsLat, obsCatchArea, modLon, modLat, modCatchArea, windowSize, misMatch)
    print location, obsStep, xSel, ySel
    if xSel != -999. and ySel != -999.:
      modValues = getModelData("%s/%s" %(inputDir, dischargeFileName), xSel, ySel)
      obsValues = getObservationData("%s/%s" %(dischargeDir, location))
      if obsStep == 1 and modStep != 1: obsValues, obsStart, obsEnd = aggregateToMonth(obsValues, obsStart, obsEnd)
      if modStep == 1 and obsStep != 1: modValues, modStart, modEnd = aggregateToMonth(modValues, modStart, modEnd)
      if obsStep == 1 and modStep == 1: print "Daily data"
      output[3:-1] = calculateMetrics(obsValues, modValues, obsStart, obsEnd, modStart, modEnd, obsStep, modStep)
  return np.array(output)

def f(location):
  print location


getGlobalProperties(configFile, reference=False)
print inputDir, dischargeFileName
pool = mp.Pool(processes=numCores)

#output = np.zeros((len(locations), 11))
#for location in range(len(locations)):
#  print location/float(len(locations)), locations[location]
#  output[location,:] = extractLocation(location,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep)

results = [pool.apply_async(extractLocation,args=(loc,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep)) for loc in range(len(locations))]
outputList = [p.get() for p in results]
output = np.array(outputList)

getGlobalProperties(configFile, reference=True)
print inputDir, dischargeFileName, modCatchArea.shape

#print len(locations)
#output2 = np.zeros((len(locations), 11))
#for location in range(len(locations)):
#  print location/float(len(locations)), locations[location]
#  output2[location,:] = extractLocation(locations[location],inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea)

results2 = [pool.apply_async(extractLocation,args=(loc,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep)) for loc in range(len(locations))]
outputList2 = [p.get() for p in results2]
output2 = np.array(outputList2)

with open('validationResultsPool_%s_%s.obj' %(runName, refName), 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump([output, output2], f)


