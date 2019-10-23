"""

"""

import arcpy
from arcpy.sa import *
import numpy
import os
import sys
import csv
import gc

def getScore(csvPath, ec, pl):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)
		plIndex = 0
		for cell in table[1]:
			if cell == pl:
				break
			plIndex += 1

		score = None
		for row in table:
			if row[1] == ec:
				score = row[plIndex]
				break
	return score

def getMean(csvPath, pl):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)
		plIndex = 0
		for cell in table[1]:
			if cell == pl:
				break
			plIndex += 1

		m = 0
		sum = 0
		count = 0
		for row in table:
			if row[1].startswith("EC_"):
				count += 1
				sum += float(row[plIndex])
		m = sum / count
	return m

def printError(message, e):
	arcpy.AddMessage("\n---- ERROR ----")
	arcpy.AddMessage("Information: " + message)
	if e is not None:
		arcpy.AddMessage("Message: " + str(e))
	arcpy.AddMessage("---------------\n")


if __name__ == '__main__':

	# ------------- Set up --------------

	sensitivityScoresCsv = arcpy.GetParameterAsText(0) # Sensitivity scores CSV file
	plRasterFolder = arcpy.GetParameterAsText(1) # Pressures raster folder
	method = arcpy.GetParameterAsText(2) # SUM, Weighted SUM
	plLabels = arcpy.GetParameterAsText(3) # Pressure rasters labels from the tool interface
	outputFolder = arcpy.GetParameterAsText(4) # Folder to save final raster and statistics matrix
	resultName = arcpy.GetParameterAsText(5) # Result raster name (tif)

	# Removing ' from input rasters tools labels
	plLabels = plLabels.replace("'","")

	# Tool labels arrays
	plLabelsArray = plLabels.split(";")

	# Input raster names and labels
	plLayerNames = []
	plLayerLabels = []

	for label in plLabelsArray:
		plLayerNames.append(label.split(":")[0])
		plLayerLabels.append(label.split(":")[1])

	layersDict = {} # Raster matrices storage

	# ------------- Set up --------------





	# ------------- Rasters transformation to matrices--------------

	arcpy.AddMessage("----- Transforming rasters to matrices. ")
	# transform all Pressure rasters to 2D matrices
	for pl in plLayerNames:
		plPath = plRasterFolder + "\\" + pl + ".tif"
		plNumArray = arcpy.RasterToNumPyArray(plPath, nodata_to_value = -100)
		plNumArray = plNumArray.astype('float')
		plNumArray[plNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
		layersDict[pl] = plNumArray

	# ------------- Rasters transformation to matrices--------------

	# Calculating BSPI
	if ((len(plLayerNames) > 1) or ((len(plLayerNames) == 1) and (plLayerNames[0] != ''))):
		arcpy.AddMessage("----- Calculating BSPI matrix. ")

		# final 2D matrix that will will be transformed to final raster
		finalArray = None
		# Describe raster object
		dsc = None
		# raster index for processing
		plIndex = 0

		for pl in plLayerNames:
			try:
				if plIndex == 0:
					# Get Describe object from the first processed raster
					dsc = arcpy.Describe(plRasterFolder + "\\" + pl + ".tif")

					if method == "SUM":
						# Init final 2D matrix
						finalArray = layersDict.get(pl)
					elif method == "Weighted SUM":
						s = getMean(sensitivityScoresCsv, pl)
						finalArray = layersDict.get(pl) * float(s)
				else:
					if method == "SUM":
						# Add matrix to final matrix
						finalArray += layersDict.get(pl)
					elif method == "Weighted SUM":
						s = getMean(sensitivityScoresCsv, pl)
						tmpArray = layersDict.get(pl) * float(s)
						finalArray += tmpArray


			except Exception as e:
				printError("Error calculating BSPI matrix.", e)

			plIndex += 1

		arcpy.AddMessage("----- Creating BSPI raster.")
		# transform final matrix to raster
		finalRaster = arcpy.NumPyArrayToRaster(finalArray, arcpy.Point(dsc.Extent.XMin, dsc.Extent.YMin), dsc.meanCellWidth, dsc.meanCellHeight)
		outFile = outputFolder + "\\" + resultName

		# save final raster to output folder
		finalRaster.save(outFile)

		arcpy.AddMessage("----- Defining projection of BSPI raster.")
		# definal final raster projection
		arcpy.DefineProjection_management(outFile, dsc.SpatialReference)

	else:
		printError("Calculation can't be performed. Select at least 1 Pressures layer. ", None)
