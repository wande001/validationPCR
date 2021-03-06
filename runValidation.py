import ConfigParser
import io
import sys
import os
import pickle
import glob
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
  global inputDir, dischargeFileName, summary, full, dischargeDir, runName, refName, logFile, includeRef
  config = readConfigFile(configFile)
  if reference:
    inputDir = str(config.get('Reference options', 'inputDir'))
    dischargeFileName = str(config.get('Reference options', 'dischargeFileName'))
    logFile = str(config.get('Reference options', 'logFile'))
  else:
    inputDir = str(config.get('Main options', 'inputDir'))
    dischargeFileName = str(config.get('Main options', 'dischargeFileName'))
    logFile = str(config.get('Main options', 'logFile'))
  full = config.get('Output options', 'FullAnalysis') == str("True")
  summary = config.get('Output options', 'summaryAnalysis') == str("True")
  includeRef = config.get('Output options', 'includeReference') == str("True")
  dischargeDir = str(config.get('GRDC data', 'dischargeDir'))
  numCores = str(config.get('general', 'numCores'))
  runName = str(config.get('Main options', 'RunName'))
  refName = str(config.get('Reference options', 'RunName'))
  return inputDir, dischargeFileName, summary, full, dischargeDir,runName,refName, logFile, includeRef

def matchMaps(f, modLon, modLat):
  dxMod = modLon[1] - modLon[0]
  dyMod = modLat[1] - modLat[0]
  dx = f.GetGeoTransform()[1]
  dy = f.GetGeoTransform()[5]
  startx = f.GetGeoTransform()[0] + dx/2.
  starty = f.GetGeoTransform()[3] + dy/2.
  catchmentArea = f.GetRasterBand(1).ReadAsArray()
  leny, lenx = catchmentArea.shape
  endx = startx + dx * lenx
  endy = starty + dy * leny
  lonCatch = np.arange(startx, endx, dx)
  latCatch = np.arange(starty, endy, dy)
  if (dxMod < 0 and dx > 0):
    lonCatch = lonCatch[::-1]
    catchmentArea = catchmentArea[:,::-1]
  elif (dxMod > 0 and dx < 0):
    lonCatch = lonCatch[::-1]
    catchmentArea = catchmentArea[:,::-1]
  xSel1 = np.argmax(lonCatch == modLon[0])
  xSel2 = np.argmax(lonCatch == modLon[-1])+1
  lonSel = np.minimum(xSel1,xSel2)
  if (dxMod > 0 and dx < 0):
    catchmentArea = catchmentArea[::-1,:]
  elif (dxMod < 0 and dx > 0):
    catchmentArea = catchmentArea[::-1,:]
  ySel1 = np.argmax(latCatch == modLat[0])
  ySel2 = np.argmax(latCatch == modLat[-1])+1
  #latSel = np.minimum(ySel1,ySel2)
  catchmentArea = catchmentArea[ySel1:ySel2, xSel1:xSel2]
  print catchmentArea.shape, len(modLon), len(modLat)
  return(catchmentArea)

def getCatchmentMap(config, modLon, modLat, option = "Reference options"):
  catchmentAreaMap = str(config.get(option, 'catchmentAreaMap'))
  cellAreaMap = str(config.get(option, 'cellAreaMap'))
  cellAreaConstant = str(config.get(option, 'cellAreaConstant'))
  routingMap = str(config.get(option, 'routingMap'))
  if os.path.exists(catchmentAreaMap):
    f = gdal.Open(catchmentAreaMap)
    catchmentArea = matchMaps(f, modLon, modLat)
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
  return lon, lat, timeVar[0], timeVar[-1], modStep, timeVar

def getWindowSize(config, option = "general"):
  windowSize = int(config.get(option, 'windowSizeForMatching'))
  return windowSize

def getTimeSize(config, option = "general"):
  timeSize = int(config.get(option, 'minimumLengthForComparison'))
  return timeSize

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
  return np.array(monthlyData), startDate, endDate, 28

def aggregateToMonthTime(data, startDate, endDate):
  deltaDays = (endDate - startDate).days
  date_list = [startDate + datetime.timedelta(days=x) for x in range(0, deltaDays+1)]
  tempMonth = startDate.month
  dayCount = 0
  tempStart = 0
  monthlyData = []
  while dayCount <= deltaDays:
    if date_list[dayCount].month != tempMonth:
      monthlyData.append(data[tempStart])
      tempStart = dayCount
      tempMonth = date_list[dayCount].month
    dayCount += 1
  monthlyData.append(data[tempStart])
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
  day = monthrange(year, month)[1]
  return datetime.datetime(year, month, day)

def matchSeriesMonth(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes):
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
      plotTimes = modTimes[modStartShift:]
    else:
      modOut = mod[modStartShift:modEndShift]
      plotTimes = modTimes[modStartShift:modEndShift]
    obsOut[obsOut< 0.0] = np.nan
    modOut[modOut< 0.0] = np.nan
    return obsOut, modOut, plotTimes
  else:
        return [], [], []

def matchSeriesDay(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes):
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
      plotTimes = modTimes[modStartShift:]
    else:
      modOut = mod[modStartShift:modEndShift]
      plotTimes = modTimes[modStartShift:modEndShift]
    obsOut[obsOut< 0.0] = np.nan
    modOut[modOut< 0.0] = np.nan
    return obsOut, modOut, plotTimes
  else:
	return [], [], []
	
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
  return kge, cc, alpha, beta

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

def calculateMetrics(obs, mod, obsStart, obsEnd, modStart, modEnd, obsStep, modStep, modTimes):
  if obsStep > 1 or modStep > 1: obs, mod, plotTimes = matchSeriesMonth(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes)
  if obsStep <= 1 and modStep <= 1: obs, mod, plotTimes = matchSeriesDay(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes)
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  if len(obs) > timeSize and len(mod) > timeSize:
    obsSel = np.isnan(obs) == False
    modSel = np.isnan(mod) == False
    sel = obsSel & modSel
    if len(obs[sel]) > timeSize and len(mod[sel]) > timeSize:
      R = spearmanr(obs[sel], mod[sel])[0]
      NS = nashSutcliffe(obs[sel], mod[sel])
      RMSE = rmse(obs[sel], mod[sel])
      Bias, numPoints = bias(obs[sel], mod[sel])
      KGE, CC, Alpha, Beta = kge(obs[sel], mod[sel])
      AC = anomalyCorrelation(obs[sel], mod[sel])
      return R, AC, KGE, CC, Alpha, Beta, NS, RMSE, Bias, numPoints
    else:
      return np.zeros((10))
  else:
    return np.zeros((10))

def calculateFullMetrics(obs, mod, obsStart, obsEnd, modStart, modEnd, obsStep, modStep, modTimes):
  if obsStep > 1 or modStep > 1: obs, mod, plotTimes = matchSeriesMonth(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes)
  if obsStep <= 1 and modStep <= 1: obs, mod, plotTimes = matchSeriesDay(obs, mod, obsStart, obsEnd, modStart, modEnd, modTimes)
  obsSel = np.isnan(obs) == False
  modSel = np.isnan(mod) == False
  sel = obsSel & modSel
  if len(obs) > timeSize and len(mod) > timeSize:
    obsSel = np.isnan(obs) == False
    modSel = np.isnan(mod) == False
    sel = obsSel & modSel
    if len(obs[sel]) > timeSize and len(mod[sel]) > timeSize:
      outData = {
		"observations" : obs[sel],
		"modelled": mod[sel],
		"times": plotTimes[sel],}
      return outData
    else:
      outData = {
		"observations" : [],
		"modelled": [],
		"times": [],}
      return outData
  else:
    outData = {
		"observations" : [],
		"modelled": [],
		"times": [],}
    return outData


def getWaterBalance(fileNames):
  varData = {
		"year" : [],
		"precipitation" : [],
		"actualET": [],
		"runoff": [],
		"totalPotentialGrossDemand": [],
		"baseflow": [],
		"storage": [],}
  
  varNames = varData.keys()
  for fileName in glob.glob(fileNames):
    print fileName
    try:
      f = open(fileName, "r")
      lines = f.readlines()
      f.close()
    except:
          lines = ""
    for line in lines:
      varFields = line.split(" ")
      if len(varFields) > 2:
        if varFields[2] == "pcrglobwb" and len(varFields) == 9:
          try:
            year = int(varFields[-1][:4])
          except:
            if varFields[2] == "pcrglobwb" and len(varFields) == 16: year = int(varFields[11][:4])
          if year not in varData["year"]: varData["year"].append(year)
        if len(varFields) > 15:
          if varFields[7] == "days" and varFields[8] == "1" and varFields[9] == "to":
            for var in varNames:
              if varFields[6] == var:
                varData[var].append(float(varFields[14]))
        if len(varFields) > 14:
          if varFields[6] == "days" and varFields[7] == "1" and varFields[8] == "to":
            for var in varNames:
              if varFields[5] == var:
                varData[var].append(float(varFields[13]))
  return(varData)

def getGlobalProperties(configFile, reference):
  global modLon, modLat, modStart, modEnd, modStep, modTimes, inputDir, modCatchArea, windowSize, timeSize, misMatch, locations, numCores
  getArguments(configFile, reference)
  locations = os.listdir(dischargeDir)
  modLon, modLat, modStart, modEnd, modStep, modTimes = readModelProps("%s/%s" %(inputDir, dischargeFileName))
  if reference:
    option = "Reference options"
  else:
    option = "Main options"
  modCatchArea = getCatchmentMap(config, modLon, modLat, option)
  windowSize = getWindowSize(config, "general")
  timeSize = getWindowSize(config, "general")
  misMatch = getAreaMisMatch(config, "general")
  numCores = getCores(config, "general")
  return modLon, modLat, modStart, modEnd, modStep, modTimes, inputDir, modCatchArea, windowSize, timeSize, misMatch, locations, numCores
  
def extractLocation(location,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep, modTimes):
  print location/float(len(locations)), locations[location]
  location = locations[location]
  output = np.zeros((14))
  series = []
  #if location[-4:] == ".mon":
  if location[-4:] == ".mon" or location[-4:] == ".day":
    obsLon, obsLat, obsCatchArea, obsStart, obsEnd, obsStep = readObservationsProps("%s/%s" %(dischargeDir, location))
    output[0:2] = obsLon, obsLat
    output[2] = obsCatchArea
    xSel, ySel = matchLocation(obsLon, obsLat, obsCatchArea, modLon, modLat, modCatchArea, windowSize, misMatch)
    if xSel != -999. and ySel != -999.:
      modValues = getModelData("%s/%s" %(inputDir, dischargeFileName), xSel, ySel)
      obsValues = getObservationData("%s/%s" %(dischargeDir, location))
      if obsStep == 1 and modStep != 1: obsValues, obsStart, obsEnd, obsStep = aggregateToMonth(obsValues, obsStart, obsEnd)
      elif modStep == 1 and obsStep != 1:
        modValues, modStart, modEnd, modStep = aggregateToMonth(modValues, modStart, modEnd)
        modTimes, modStart, modEnd = aggregateToMonthTime(modTimes, modStart, modEnd)
      elif obsStep == 1 and modStep == 1: print "Daily data"
      if summary:
        output[3:-1] = calculateMetrics(obsValues, modValues, obsStart, obsEnd, modStart, modEnd, obsStep, modStep, modTimes)
      if full:
        series = calculateFullMetrics(obsValues, modValues, obsStart, obsEnd, modStart, modEnd, obsStep, modStep, modTimes)
    output[13] = float(obsStep)
  return np.array(output), series

def f(location):
  print location

getGlobalProperties(configFile, reference=False)
pool = mp.Pool(processes=numCores)

#output = np.zeros((len(locations), 14))
#for location in range(len(locations)):
#  print location/float(len(locations)), locations[location]
#  output[location,:] = extractLocation(location,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep, modTimes)[0]

results = [pool.apply_async(extractLocation,args=(loc,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep, modTimes)) for loc in range(len(locations))]
outputList = [p.get() for p in results]

output = np.zeros((len(locations), 14))
for loc in range(len(locations)):
  output[loc,:] = outputList[loc][0]
	
if full:
  fullOutput = {"ID" : [],"data": [],}
  for loc in range(len(locations)):
    fullOutput["ID"].append(locations[loc][:-4])
    fullOutput["data"].append(outputList[loc][1])
else:
  fullOutput = []

waterBalOutput = getWaterBalance(logFile)

if includeRef:
  print "IncludeRef"

  getGlobalProperties(configFile, reference=True)

  #output2 = np.zeros((len(locations), 11))
  #for location in range(len(locations)):
  #  print location/float(len(locations)), locations[location]
  #  output2[location,:] = extractLocation(location,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep, modTimes)[0]

  results2 = [pool.apply_async(extractLocation,args=(loc,inputDir, dischargeFileName, modStart, modEnd, modLon, modLat, modCatchArea, modStep, modTimes)) for loc in range(len(locations))]
  outputList2 = [p.get() for p in results2]

  output2 = np.zeros((len(locations), 14))
  for loc in range(len(locations)):
    output2[loc,:] = outputList2[loc][0]

  if full:
    fullOutput2 = {"ID" : [],"data": [],}
    for loc in range(len(locations)):
      fullOutput2["ID"].append(locations[loc][:-4])
      fullOutput2["data"].append(outputList2[loc][1])
  else:
    fullOutput2 = []

  waterBalOutput2 = getWaterBalance(logFile)

if includeRef == False:
  output2 = output
  fullOutput2 = fullOutput
  waterBalOutput2 = waterBalOutput

with open('validationResultsPool_%s_%s.obj' %(runName, refName), 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump([output, output2, fullOutput, fullOutput2, waterBalOutput, waterBalOutput2], f)


