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

def printError(message, e):
	arcpy.AddMessage("\n---- ERROR ----")
	arcpy.AddMessage("Information: " + message)
	if e is not None:
		arcpy.AddMessage("Message: " + str(e))
	arcpy.AddMessage("---------------\n")


if __name__ == '__main__':

	# ------------- Set up --------------

	#method = arcpy.GetParameterAsText(0) # SUM or Halpern MEAN (Halpern MEAN not implemented yet)
	sensitivityScoresCsv = arcpy.GetParameterAsText(0) # Sensitivity scores CSV file
	ecRasterFolder = arcpy.GetParameterAsText(1) # Ecosystem component raster folder
	plRasterFolder = arcpy.GetParameterAsText(2) # Pressures raster folder
	ecLabels = arcpy.GetParameterAsText(3) # Ecosystem component rasters labels from the tool interface
	plLabels = arcpy.GetParameterAsText(4) # Pressure rasters labels from the tool interface
	calculateStatsMatrix = arcpy.GetParameterAsText(5) # True or False to calculate statistics matrix
	outputFolder = arcpy.GetParameterAsText(6) # Folder to save final raster and statistics matrix
	resultName = arcpy.GetParameterAsText(7) # Result raster name (tif)
	matrixFileName = arcpy.GetParameterAsText(8) # Statistics matrix name (csv)

	# Removing ' from input rasters tools labels
	ecLabels = ecLabels.replace("'","")
	plLabels = plLabels.replace("'","")

	# Tool labels arrays
	ecLabelsArray = ecLabels.split(";")
	plLabelsArray = plLabels.split(";")

	# Input raster names and labels
	ecLayerNames = []
	ecLayerLabels = []
	plLayerNames = []
	plLayerLabels = []

	for label in ecLabelsArray:
		ecLayerNames.append(label.split(":")[0])
		ecLayerLabels.append(label.split(":")[1])

	for label in plLabelsArray:
		plLayerNames.append(label.split(":")[0])
		plLayerLabels.append(label.split(":")[1])

	layersDict = {} # Raster matrices storage

	martrixArray = [] # Statistics result matrix
	plTmpSumDict = {} # Used for statistics matrix (Last row sums)

	# ------------- Set up --------------





	# ------------- Rasters transformation to matrices--------------

	arcpy.AddMessage("----- Transforming rasters to matrices. ")
	# transform all Ecosystem component rasters to 2D matrices
	for ec in ecLayerNames:
		ecPath = ecRasterFolder + "\\" + ec + ".tif"
		ecNumArray = arcpy.RasterToNumPyArray(ecPath, nodata_to_value = -100)
		ecNumArray = ecNumArray.astype('float')
		ecNumArray[ecNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
		layersDict[ec] = ecNumArray

	# transform all Pressure rasters to 2D matrices
	for pl in plLayerNames:
		plPath = plRasterFolder + "\\" + pl + ".tif"
		plNumArray = arcpy.RasterToNumPyArray(plPath, nodata_to_value = -100)
		plNumArray = plNumArray.astype('float')
		plNumArray[plNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
		layersDict[pl] = plNumArray
		if calculateStatsMatrix == "true":
			plTmpSumDict[pl] = 0 # Used for statistics matrix

	# ------------- Rasters transformation to matrices--------------

	# arcpy.env.cellSize = "MINOF"
	#
	# nodata = None
	# if ignoreNoData == "true":
	# 	nodata = "DATA"
	# elif ignoreNoData == "false":
	# 	nodata = "NODATA"

	# Calculating BSII
	if (((len(ecLayerNames) > 1) or ((len(ecLayerNames) == 1) and (ecLayerNames[0] != '')))) and (((len(plLayerNames) > 1) or ((len(plLayerNames) == 1) and (plLayerNames[0] != '')))):
		arcpy.AddMessage("----- Calculating BSII matrix. ")

		# Two first rows of CSV statistics matrix
		if calculateStatsMatrix == "true":
			matrixplLayerLabels = ["", ""] + plLayerLabels
			matrixplLayerLabels.append("")
			matrixplLayerLabels.append("SUM")
			matrixplLayerNames = ["", ""] + plLayerNames
			martrixArray.append(matrixplLayerLabels)
			martrixArray.append(matrixplLayerNames)

		# needed?
		#ecTmpRasters = []

		# final 2D matrix that will will be transformed to final raster
		finalArray = None
		# Describe raster object
		dsc = None
		# raster index for processing
		ecIndex = 0

		for ec in ecLayerNames:

			if calculateStatsMatrix == "true":
				# Create statistics matrix row with label ec layer and name at the beginning
				matrixEcRow = [ecLayerLabels[ecIndex], ec]

			# Sum for the statistics matrix last column of each row
			ecTmpSum = 0

			# Needed for Halpern method
			#if method == "Halpern MEAN":
			#	ecTmpRasters.append(Raster(ecPath))

			for pl in plLayerNames:

				try:
					# get sensitivity score for each layer combination
					s = getScore(sensitivityScoresCsv, ec, pl)

					# Temp matrix - Multiply 2D matrices and sensitivity score
					tmpArray = layersDict.get(ec) * layersDict.get(pl) * float(s)

					if ecIndex == 0:
						# Get Describe object from the first processed raster
						dsc = arcpy.Describe(ecRasterFolder + "\\" + ec + ".tif")
						# Init final 2D matrix
						finalArray = tmpArray

					else:
						# Add temp matrix to final matrix
						finalArray += tmpArray


					if calculateStatsMatrix == "true":
						# Sum of temp matrix values, treating NaNs as zeros
						tmpMatrixSum = numpy.nansum(tmpArray)
						# Add temp matrix sum to the statistics matrix
						matrixEcRow.append(tmpMatrixSum) # matrix
						# Add temp matrix sum to be written to last statistics matrix row
						plTmpSumDict[pl] += tmpMatrixSum
						# Add temp matrix sum to be written to last statistics matrix column
						ecTmpSum += tmpMatrixSum

				except Exception as e:
					printError("Error calculating BSII matrix or statistics matrix.", e)

			if calculateStatsMatrix == "true":
				matrixEcRow.append("")
				# Add temp matrices sum to last statistics matrix column
				matrixEcRow.append(ecTmpSum)
				# Statistics matrix row finalized, add to statistics matrix
				martrixArray.append(matrixEcRow) # matrix

			ecIndex += 1

		arcpy.AddMessage("----- Creating BSII raster.")
		# transform final matrix to raster
		finalRaster = arcpy.NumPyArrayToRaster(finalArray, arcpy.Point(dsc.Extent.XMin, dsc.Extent.YMin), dsc.meanCellWidth, dsc.meanCellHeight)
		outFile = outputFolder + "\\" + resultName

		# save final raster to output folder
		finalRaster.save(outFile)

		arcpy.AddMessage("----- Defining projection of BSII raster.")
		# definal final raster projection
		arcpy.DefineProjection_management(outFile, dsc.SpatialReference)

		if calculateStatsMatrix == "true":
			arcpy.AddMessage("----- Writing statistics matrix.")
			martrixArray.append("")
			plTmsSumsRow = ["SUM", ""]
			# Add last row to statistics matrix
			for pl in plLayerNames:
				plTmsSumsRow.append(plTmpSumDict.get(pl))
			plTmsSumsRow.append("")
			# Add sum of all raster values to last cell of statistics matrix
			plTmsSumsRow.append(numpy.nansum(finalArray))
			martrixArray.append(plTmsSumsRow)

			# Save statistics matrix in CSV file
			if (sys.version_info > (3, 0)): # Python 3
				with open(outputFolder + "//" + matrixFileName, mode='w', newline='') as matrixFile:
					matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
					for matrixRow in martrixArray:
						matrixWriter.writerow(matrixRow)
			else: # Python 2
				with open(outputFolder + "//" + matrixFileName, mode='wb') as matrixFile:
					matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
					for matrixRow in martrixArray:
						matrixWriter.writerow(matrixRow)

	else:
		printError("Calculation can't be performed. Select at least 1 Ecosystem component and 1 Pressures layer. ", None)
