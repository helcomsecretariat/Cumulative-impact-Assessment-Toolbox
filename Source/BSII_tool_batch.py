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

def getEcNames(csvPath):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	ecNames = []
	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)

		i = 0
		for cell in table[1]:
			if cell == "EC CODE":
				ecIndex = i
				break
			i += 1

		for row in table:
			if row[i].startswith("EC_"):
				ecNames.append(row[i])

	return ecNames

def getPlNames(csvPath):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	plNames = []
	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)

		for cell in table[1]:
			if cell.startswith("PL_"):
				plNames.append(cell)

	return plNames

def saveMatrixAsRaster(keyMatrix, dsc, path):
	# transform final matrix to raster
	keyRaster = arcpy.NumPyArrayToRaster(keyMatrix, arcpy.Point(dsc.Extent.XMin, dsc.Extent.YMin), dsc.meanCellWidth, dsc.meanCellHeight)

	# save final raster to output folder
	keyRaster.save(path)

	# definal final raster projection
	arcpy.DefineProjection_management(path, dsc.SpatialReference)

def printError(message, e):
	arcpy.AddMessage("\n---- ERROR ----")
	arcpy.AddMessage("Information: " + message)
	if e is not None:
		arcpy.AddMessage("Message: " + str(e))
	arcpy.AddMessage("---------------\n")


if __name__ == '__main__':

	method = arcpy.GetParameterAsText(0) # EV or ES method
	sensitivityScoresFolder = arcpy.GetParameterAsText(1) # Sensitivity scores folder with CSV files
	ecRasterFolder = arcpy.GetParameterAsText(2) # Ecosystem component raster folder
	plRasterFolder = arcpy.GetParameterAsText(3) # Pressures raster folder
	calculateStatsMatrix = arcpy.GetParameterAsText(4) # True or False to calculate statistics matrix
	outputFolder = arcpy.GetParameterAsText(5) # Folder to save final rasters and statistics matrices

	layersDict = {} # Raster matrices storage
	aggregatedLayers = {}

	# Describe raster object
	dsc = None

	for file in os.listdir(sensitivityScoresFolder):
		# Input raster names
		ecLayerNames = []
		plLayerNames = []

		if ((file.endswith(".csv")) and (file.startswith(method+"_"))):

			resultName = file.split(".")[0]
			arcpy.AddMessage("----- Processing " + resultName)
			ecLayerNames = getEcNames(sensitivityScoresFolder+"//"+file)
			plLayerNames = getPlNames(sensitivityScoresFolder+"//"+file)
			# arcpy.AddMessage(resultName)
			# arcpy.AddMessage("EC Layers " + str(ecLayerNames))
			# arcpy.AddMessage("PL Layers " + str(plLayerNames))

			martrixArray = [] # Statistics result matrix
			plTmpSumDict = {} # Used for statistics matrix (Last row sums)

			# ------------- Rasters transformation to matrices--------------
			#arcpy.AddMessage("----- Transforming rasters to matrices. ")
			# transform all Ecosystem component rasters to 2D matrices
			for ec in ecLayerNames:
				if not ec in layersDict:
					ecPath = ecRasterFolder + "\\" + ec + ".tif"
					ecNumArray = arcpy.RasterToNumPyArray(ecPath, nodata_to_value = -100)
					ecNumArray = ecNumArray.astype('float')
					ecNumArray[ecNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
					layersDict[ec] = ecNumArray

			# transform all Pressure rasters to 2D matrices
			for pl in plLayerNames:
				if not pl in layersDict:
					plPath = plRasterFolder + "\\" + pl + ".tif"
					plNumArray = arcpy.RasterToNumPyArray(plPath, nodata_to_value = -100)
					plNumArray = plNumArray.astype('float')
					plNumArray[plNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
					layersDict[pl] = plNumArray

				if calculateStatsMatrix == "true":
					plTmpSumDict[pl] = 0 # Used for statistics matrix

			# ------------- Rasters transformation to matrices--------------

			# Calculating BSII
			if (((len(ecLayerNames) > 1) or ((len(ecLayerNames) == 1) and (ecLayerNames[0] != '')))) and (((len(plLayerNames) > 1) or ((len(plLayerNames) == 1) and (plLayerNames[0] != '')))):
				#arcpy.AddMessage("----- Calculating BSII matrix. ")

				# Two first rows of CSV statistics matrix
				if calculateStatsMatrix == "true":
					matrixplLayerLabels = ["", ""]
					matrixplLayerNames = ["", ""] + plLayerNames
					matrixplLayerNames.append("")
					matrixplLayerNames.append("SUM")
					martrixArray.append(matrixplLayerLabels)
					martrixArray.append(matrixplLayerNames)


				# final 2D matrix that will will be transformed to final raster
				finalArray = None

				# raster index for processing
				ecIndex = 0

				# total BSII sum
				bsiiSum = 0

				for ec in ecLayerNames:

					if calculateStatsMatrix == "true":
						# Create statistics matrix row with label ec layer and name at the beginning
						matrixEcRow = ["", ec]

					# Sum for the statistics matrix last column of each row
					ecTmpSum = 0

					for pl in plLayerNames:

						try:
							# get sensitivity score for each layer combination
							s = getScore(sensitivityScoresFolder+"//"+file, ec, pl)

							# Temp matrix - Multiply 2D matrices and sensitivity score
							tmpArray = layersDict.get(ec) * layersDict.get(pl) * float(s)

							if ecIndex == 0:
								# Get Describe object from the first processed raster
								if dsc is None:
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
								# Add temp matrix sum to be written to last statistics matrix cell
								bsiiSum += tmpMatrixSum

						except Exception as e:
							printError("Error calculating BSII matrix or statistics matrix.", e)

					if calculateStatsMatrix == "true":
						matrixEcRow.append("")
						# Add temp matrices sum to last statistics matrix column
						matrixEcRow.append(ecTmpSum)
						# Statistics matrix row finalized, add to statistics matrix
						martrixArray.append(matrixEcRow) # matrix

					ecIndex += 1

				# Calculate EV aggregated matrices
				# Banthic habitats
				if resultName.endswith("_ECG_BH"):
					if "Aggregated_ECG_BH" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECG_BH"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECG_BH"] += finalArray
				# Fish
				elif resultName.endswith("_ECG_FH"):
					if "Aggregated_ECG_FH" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECG_FH"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECG_FH"] += finalArray
				# Bird
				elif resultName.endswith("_ECG_BD"):
					if "Aggregated_ECG_BD" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECG_BD"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECG_BD"] += finalArray
				# Mammal
				elif resultName.endswith("_ECG_MM"):
					if "Aggregated_ECG_MM" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECG_MM"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECG_MM"] += finalArray

				# Calculate aggregated ECSG matrices (Sum all same ES aspect matrices of each ECSG)
				# Benthic landscapes
				if resultName.endswith("_ECSG_BL"):
					if "Aggregated_ECSG_BL" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_BL"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_BL"] += finalArray
				# Aggregated Habitatforming Species
				elif resultName.endswith("_ECSG_HS"):
					if "Aggregated_ECSG_HS" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_HS"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_HS"] += finalArray
				# N2000
				elif resultName.endswith("_ECSG_N2"):
					if "Aggregated_ECSG_N2" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_N2"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_N2"] += finalArray
				# Birds
				elif resultName.endswith("_ECSG_BD"):
					if "Aggregated_ECSG_BD" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_BD"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_BD"] += finalArray
				# Mammals
				elif resultName.endswith("_ECSG_MM"):
					if "Aggregated_ECSG_MM" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_MM"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_MM"] += finalArray
				# Fish habitats
				elif resultName.endswith("_ECSG_EFH"):
					if "Aggregated_ECSG_EFH" not in aggregatedLayers:
						aggregatedLayers["Aggregated_ECSG_EFH"] = finalArray
					else:
						aggregatedLayers["Aggregated_ECSG_EFH"] += finalArray


				#arcpy.AddMessage("----- Creating BSII raster.")
				saveMatrixAsRaster(finalArray, dsc, outputFolder + "\\BSII_" + resultName + ".tif")

				if calculateStatsMatrix == "true":
					#arcpy.AddMessage("----- Writing statistics matrix.")
					martrixArray.append("")
					plTmsSumsRow = ["SUM", ""]
					# Add last row to statistics matrix
					for pl in plLayerNames:
						plTmsSumsRow.append(plTmpSumDict.get(pl))
					plTmsSumsRow.append("")
					# Add sum of all raster values to last cell of statistics matrix
					#plTmsSumsRow.append(numpy.nansum(finalArray))
					plTmsSumsRow.append(bsiiSum)
					martrixArray.append(plTmsSumsRow)

					# Save statistics matrix in CSV file
					if (sys.version_info > (3, 0)): # Python 3
						with open(outputFolder + "\\Matrix_" + resultName + ".csv", mode='w', newline='') as matrixFile:
							matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
							for matrixRow in martrixArray:
								matrixWriter.writerow(matrixRow)
					else: # Python 2
						with open(outputFolder + "\\Matrix_" + resultName + ".csv", mode='wb') as matrixFile:
							matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
							for matrixRow in martrixArray:
								matrixWriter.writerow(matrixRow)

			else:
				printError("Calculation can't be performed.", None)

	#Aggregated matrices
	finalMatrix = None
	finalIndex = 0
	matricesForMax = [] # used for ES only

	#arcpy.AddMessage("----- Calculating Aggregated matrices.")
	for key, matrix in aggregatedLayers.items():
		# Normalize aggregated ECSG matrices
		minv = numpy.nanmin(matrix)
		maxv = numpy.nanmax(matrix)
		matrix = (matrix - minv) / (maxv - minv)

		# Save normalized aggregated matrices as raster layers
		saveMatrixAsRaster(matrix, dsc, outputFolder + "\\BSII_" + key + ".tif")

		if method == "EV":
			# Sum aggregated ECG matrices
			if finalIndex == 0:
				finalMatrix = matrix
			else:
				finalMatrix += matrix
			finalIndex += 1

		elif method == "ES":
			# prepare matrices for MAX aggregation (ES only)
			if key == "Aggregated_ECSG_BL":
				matricesForMax.append(matrix)
			if key == "Aggregated_ECSG_HS":
				matricesForMax.append(matrix)
			if key == "Aggregated_ECSG_N2":
				matricesForMax.append(matrix)
			if key == "Aggregated_ECSG_EFH":
				matricesForMax.append(matrix)

			# add 2 matrices (if selected) to the final Aggregated_Ecosystem_Service_matrix.tif (ES only)
			if key == "Aggregated_ECSG_BD":
				if finalIndex == 0:
					finalMatrix = matrix
				else:
					finalMatrix += matrix
				finalIndex += 1
			if key == "Aggregated_ECSG_MM":
				if finalIndex == 0:
					finalMatrix = matrix
				else:
					finalMatrix += matrix
				finalIndex += 1

	if method == "ES":
		# calculate MAX pixel-wise (combine 4 (or less if not all selected) matrices into one MAX matrix)
		aggregatedMaxMatrix = None
		if len(matricesForMax) > 0:
			aggregatedMaxMatrix = matricesForMax[0]
		if len(matricesForMax) > 1:
			for i in range(1, len(matricesForMax)):
				aggregatedMaxMatrix = numpy.fmax(aggregatedMaxMatrix, matricesForMax[i])

		if aggregatedMaxMatrix is not None:
			saveMatrixAsRaster(aggregatedMaxMatrix, dsc, outputFolder + "\\BSII_Aggregated_MAX_matrix.tif")
			# add aggregated MAX matrix to final matrix
			if finalIndex == 0:
				finalMatrix = aggregatedMaxMatrix
			else:
				finalMatrix += aggregatedMaxMatrix
			finalIndex += 1

	# create final matrix
	if finalMatrix is not None:
		finalMatrix = finalMatrix / finalIndex
		minv = numpy.nanmin(finalMatrix)
		maxv = numpy.nanmax(finalMatrix)
		finalMatrix = (finalMatrix - minv) / (maxv - minv)
		if method == "EV":
			saveMatrixAsRaster(finalMatrix, dsc, outputFolder + "\\BSII_Aggregated_Ecological_Value_matrix.tif")
		elif method == "ES":
			saveMatrixAsRaster(finalMatrix, dsc, outputFolder + "\\BSII_Aggregated_Ecosystem_Service_matrix.tif")

	# Clean memory if possible
	gc.collect()
