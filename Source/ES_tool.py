"""

"""

import arcpy
from arcpy.sa import *
import numpy
import os
import sys
import csv
import gc

def getData(csvPath, es, ecg):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	dataArray = []
	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)

		# finding indices of EC layers (ec) and EC groups (ecg) in second row of csv file
		ecgIndex = None
		ecIndex = None
		i = 0
		for cell in table[1]:
			if cell == "EC CODE":
				ecIndex = i
			if cell == "GROUP BY":
				ecgIndex = i
			i += 1

		# finding indeces of ES aspects in second row of csv file
		giIndex = 0
		for cell in table[1]:
			if cell == es:
				break
			giIndex += 1

		# adding tuple ("EC layer code", "coeficient") to the array. Array has as many entries as there are EC layers in ECSG.
		for row in table:
			if row[ecgIndex] == ecg:
				dataArray.append((row[ecIndex], row[giIndex]))

	return dataArray

def getEcMatrices(csvPath, ecRasterFolder, ecgCodes):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	dict = {}

	with open(csvPath, mode) as csvFile:
		table = csv.reader(csvFile, delimiter=';')
		table = list(table)

		# finding indices of EC layers (ec) and EC groups (ecg) in second row of csv file
		ecgIndex = None
		ecIndex = None
		i = 0
		for cell in table[1]:
			if cell == "EC CODE":
				ecIndex = i
			if cell == "GROUP BY":
				ecgIndex = i
			i += 1

		# finding rasters in the folder and transforming them to numpy arrays. Using only those rasters where ECG is present in csv file
		for row in table:
			ecg = row[ecgIndex]
			ec = row[ecIndex]
			if ((ecg in ecgCodes) and (ec.startswith("EC_")) and (ec not in dict)):
				ecPath = ecRasterFolder + "\\" + ec + ".tif"
				try:
					ecNumArray = arcpy.RasterToNumPyArray(ecPath, nodata_to_value = -100)
					ecNumArray = ecNumArray.astype('float')
					ecNumArray[ecNumArray == -100] = numpy.nan # NoData values stored as Python nan values in the matrices
					dict[ec] = ecNumArray
				except Exception as e:
					arcpy.AddError("Can't read " + ec + " raster. Check if raster exists. --- " + str(e))
					sys.exit(0)
					#arcpy.AddWarning("Can't read " + ec + " raster. Check if raster exists. " + ec + " not included in calculations.")

	return dict

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

	# ------------- Set up --------------

	esCoeficientsCsv = arcpy.GetParameterAsText(0) # ES coeficients CSV file
	ecRasterFolder = arcpy.GetParameterAsText(1) # Ecosystem component raster folder
	esLabels = arcpy.GetParameterAsText(2) # ES aspects labels from the tool interface
	ecgLabels = arcpy.GetParameterAsText(3) # Ecosystem component subgroups labels from the tool interface
	outputFolder = arcpy.GetParameterAsText(4) # Location to save rasters folder
	rasterFolderName = arcpy.GetParameterAsText(5) # Name of the folder to save result rasters

	# Removing ' from input input data tools labels
	ecgLabels = ecgLabels.replace("'","")
	esLabels = esLabels.replace("'","")

	# Tool labels arrays
	ecgLabelsArray = ecgLabels.split(";")
	esLabelsArray = esLabels.split(";")

	# Input data codes and labels
	ecgCodes = []
	#ecgDesriptions = [] # not in use
	esCodes = []
	#giDescriptions = [] # not in use

	for label in ecgLabelsArray:
		ecgCodes.append(label.split(":")[0])
		#ecgDesriptions.append(label.split(":")[1])

	for label in esLabelsArray:
		esCodes.append(label.split(":")[0])
		#giDescriptions.append(label.split(":")[1])

	# Input raster storage folders
	#ecRasterFolder = os.path.dirname(os.path.realpath(__file__ + "/../")) + "/rasters/EC_for_GI/"

	# Create folder for output rasters
	outputFolderPath = outputFolder + "\\" + rasterFolderName + "\\"
	try:
		arcpy.CreateFolder_management(outputFolder, rasterFolderName)
	except Exception as e:
		arcpy.AddError("Can't create output folder. --- " + str(e))
		sys.exit(0)

	inputDataDict = {} # Input data example: ( inputDataDict["ECSG_BL_ES_FINU"] = [("EC_02", 0), ("EC_03", 1), ...] )
	arcpy.AddMessage("----- Transforming input rasters to matrices.")
	try:
		ecLayersDict = getEcMatrices(esCoeficientsCsv, ecRasterFolder, ecgCodes) #EC rasters transformed to matrices
	except Exception as e:
		arcpy.AddError("Can't read rasters. Check EC CODE and GROUP BY columns in csv file and raster names in the folder. --- " + str(e))
		sys.exit(0)

	aggregatedLayers = {}

	# ------------- Set up --------------



	# Read input data from ES input matrix and store in dictionary
	# Each dictionary entry is an array of tuples. Entry key is Ecosystem component subgroup code _ ES code. Array has as many tuples as there are EC layers in the subgroup.
	# Each tuple stores EC layer name and coeficient to multiply this layer by.
	arcpy.AddMessage("----- Reading input data from ES matrix.")
	try:
		for es in esCodes:
			for ecg in ecgCodes:
				inputDataDict[es+"_"+ecg] = getData(esCoeficientsCsv, es, ecg)
	except Exception as e:
		arcpy.AddError("Can't read coeficient values . Check EC CODE and GROUP BY columns in csv file. --- " + str(e))
		sys.exit(0)


	# Calculate final matrices and save in Rasters

	# dsc - Raster properties
	dsc = None
	countTotal = 0
	# Create one ES aspect matrix for each EC group. Total 60 matrices according to current csv file (10 ES aspects x 6 ECSG)
	arcpy.AddMessage("----- Calculating ES x ECSG matrices.")
	for key, value in inputDataDict.items():
		countKey = 0
		keyMatrix = None

		# Add all EC layers to the final EV x ECG matrix
		for item in value:
			# EC layer name
			ec = item[0]
			# coeficient
			s = float(item[1])

			if ec in ecLayersDict:
				if countTotal == 0:
					# Get Describe object from the first processed raster
					dsc = arcpy.Describe(ecRasterFolder + "\\" + ec + ".tif")

				# Multiply matrix by coeficient and add to final GI x ECSG matrix
				if countKey == 0:
					keyMatrix = ecLayersDict[ec] * s
				else:
					keyMatrix += ecLayersDict[ec] * s

				countTotal += 1
				countKey += 1

		# Minimum and maximum matrices values for normalization
		minv = numpy.nanmin(keyMatrix)
		maxv = numpy.nanmax(keyMatrix)
		# If matrices are not empty
		if (not ((minv == 0) and (maxv == 0))):
			# Normalize GI x ECG matrix
			keyMatrix = (keyMatrix - minv) / (maxv - minv)

			# Calculate aggregated ECSG matrices (Sum all same ES aspect matrices of each ECSG)
			# Benthic landscapes
			if key.endswith("_ECSG_BL"):
				if "Aggregated_ECSG_BL" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_BL"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_BL"] += keyMatrix
			# Aggregated Habitatforming Species
			elif key.endswith("_ECSG_HS"):
				if "Aggregated_ECSG_HS" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_HS"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_HS"] += keyMatrix
			# N2000
			elif key.endswith("_ECSG_N2"):
				if "Aggregated_ECSG_N2" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_N2"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_N2"] += keyMatrix
			# Birds
			elif key.endswith("_ECSG_BD"):
				if "Aggregated_ECSG_BD" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_BD"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_BD"] += keyMatrix
			# Mammals
			elif key.endswith("_ECSG_MM"):
				if "Aggregated_ECSG_MM" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_MM"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_MM"] += keyMatrix
			# Fish habitats
			elif key.endswith("_ECSG_EFH"):
				if "Aggregated_ECSG_EFH" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECSG_EFH"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECSG_EFH"] += keyMatrix
			#saveMatrixAsRaster(keyMatrix, dsc, outputFolderPath + key + ".tif")
		# Save normalized ES x ECSG matrix as raster layer
		saveMatrixAsRaster(keyMatrix, dsc, outputFolderPath + key + ".tif")

	#Final "Aggregated ecological value" matrix is an average of 4 aggregated matrices
	finalMatrix = None
	finalIndex = 0
	matricesForMax = []
	arcpy.AddMessage("----- Calculating Aggregated ECSG matrices.")
	for key, matrix in aggregatedLayers.items():
		# Normalize aggregated ECSG matrices
		minv = numpy.nanmin(matrix)
		maxv = numpy.nanmax(matrix)
		matrix = (matrix - minv) / (maxv - minv)

		# Save normalized aggregated ECSG matrices as raster layers
		saveMatrixAsRaster(matrix, dsc, outputFolderPath + key + ".tif")

		# prepare 4 matrices (or less if selected less) for MAX aggregation
		if key == "Aggregated_ECSG_BL":
			matricesForMax.append(matrix)
		if key == "Aggregated_ECSG_HS":
			matricesForMax.append(matrix)
		if key == "Aggregated_ECSG_N2":
			matricesForMax.append(matrix)
		if key == "Aggregated_ECSG_EFH":
			matricesForMax.append(matrix)

		# add 2 matrices (if selected) to the final Aggregated_Ecosystem_Service_matrix.tif
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

	# calculate MAX pixel-wise (combine 4 (or less if not all selected) matrices into one MAX matrix)
	aggregatedMaxMatrix = None
	if len(matricesForMax) > 0:
		aggregatedMaxMatrix = matricesForMax[0]
	if len(matricesForMax) > 1:
		for i in range(1, len(matricesForMax)):
			aggregatedMaxMatrix = numpy.fmax(aggregatedMaxMatrix, matricesForMax[i])

	if aggregatedMaxMatrix is not None:
		saveMatrixAsRaster(aggregatedMaxMatrix, dsc, outputFolderPath +"Aggregated_MAX_matrix.tif")
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
		saveMatrixAsRaster(finalMatrix, dsc, outputFolderPath + "Aggregated_Ecosystem_Service_matrix.tif")

	# Clean memory if possible
	gc.collect()
