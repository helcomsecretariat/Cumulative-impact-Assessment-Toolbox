"""

"""

import arcpy
from arcpy.sa import *
import numpy
import os
import sys
import csv
import gc

def getData(csvPath, ev, ecg):
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

		# finding indeces of EV aspects in second row of csv file
		giIndex = 0
		for cell in table[1]:
			if cell == ev:
				break
			giIndex += 1

		# adding tuple ("EC layer code", "coeficient") to the array. Array has as many entries as there are EC layers in ECG.
		for row in table:
			if row[ecgIndex] == ecg:
				dataArray.append((row[ecIndex], row[giIndex]))

	return dataArray

def getEcMatrices(csvPath, ecRasterFolder, ecgCodes):
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	arcpy.AddMessage("-----. " + str(ecgCodes))
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

	evCoeficientsCsv = arcpy.GetParameterAsText(0) # EV coeficients CSV file
	ecRasterFolder = arcpy.GetParameterAsText(1) # Ecosystem component raster folder
	evLabels = arcpy.GetParameterAsText(2) # EV aspects labels from the tool interface
	ecgLabels = arcpy.GetParameterAsText(3) # Ecosystem component groups labels from the tool interface
	outputFolder = arcpy.GetParameterAsText(4) # Location to save rasters folder
	rasterFolderName = arcpy.GetParameterAsText(5) # Name of the folder to save result rasters

	# Removing ' from input input data tools labels
	ecgLabels = ecgLabels.replace("'","")
	evLabels = evLabels.replace("'","")

	# Tool labels arrays
	ecgLabelsArray = ecgLabels.split(";")
	evLabelsArray = evLabels.split(";")

	# Input data codes and labels
	ecgCodes = []
	#ecgDesriptions = [] # not in use
	evCodes = []
	#giDescriptions = [] # not in use

	for label in ecgLabelsArray:
		ecgCodes.append(label.split(":")[0])
		#ecgDesriptions.append(label.split(":")[1])

	for label in evLabelsArray:
		evCodes.append(label.split(":")[0])
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

	inputDataDict = {} # Input data example: ( inputDataDict["ECG_BH_EV_BD"] = [("EC_02", 0), ("EC_03", 0), ...] )
	arcpy.AddMessage("----- Transforming input rasters to matrices.")
	try:
		ecLayersDict = getEcMatrices(evCoeficientsCsv, ecRasterFolder, ecgCodes) #EC rasters transformed to matrices
	except Exception as e:
		arcpy.AddError("Can't read rasters. Check EC CODE and GROUP BY columns in csv file and raster names in the folder. --- " + str(e))
		sys.exit(0)

	aggregatedLayers = {}

	# ------------- Set up --------------



	# Read input data from EV input matrix and store in dictionary
	# Each dictionary entry is an array of tuples. Entry key is Ecosystem component group code _ EV code. Array has as many tuples as there are EC layers in the group.
	# Each tuple stores EC layer name and coeficient to multiply this layer by.
	arcpy.AddMessage("----- Reading input data from EV matrix.")
	try:
		for ev in evCodes:
			for ecg in ecgCodes:
				inputDataDict[ev+"_"+ecg] = getData(evCoeficientsCsv, ev, ecg)
	except Exception as e:
		arcpy.AddError("Can't read coeficient values . Check EC CODE and GROUP BY columns in csv file. --- " + str(e))
		sys.exit(0)


	# Calculate final matrices and save in Rasters

	# dsc - Raster properties
	dsc = None
	countTotal = 0
	# Create one EV aspect matrix for each EC group. Total 24 matrices according to current csv file (6 EV aspects x 4 ECG)
	arcpy.AddMessage("----- Calculating EV x ECG matrices.")
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

				# Multiply matrix by coeficient and add to final GI x ECG matrix
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

			# Calculate 4 aggregated matrices (Sum all same EV aspect matrices of each ECG)
			# Banthic habitats
			if key.endswith("_ECG_BH"):
				if "Aggregated_ECG_BH" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECG_BH"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECG_BH"] += keyMatrix
			# Fish
			elif key.endswith("_ECG_FH"):
				if "Aggregated_ECG_FH" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECG_FH"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECG_FH"] += keyMatrix
			# Bird
			elif key.endswith("_ECG_BD"):
				if "Aggregated_ECG_BD" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECG_BD"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECG_BD"] += keyMatrix
			# Mammal
			elif key.endswith("_ECG_MM"):
				if "Aggregated_ECG_MM" not in aggregatedLayers:
					aggregatedLayers["Aggregated_ECG_MM"] = keyMatrix
				else:
					aggregatedLayers["Aggregated_ECG_MM"] += keyMatrix

		# Save normalized GI x ECG matrix as raster layer
		saveMatrixAsRaster(keyMatrix, dsc, outputFolderPath + key + ".tif")

	# Final "Aggregated ecological value" matrix is an average of all aggregated matrices
	finalMatrix = None
	finalIndex = 0
	arcpy.AddMessage("----- Calculating Aggregated ECG matrices.")
	for key, matrix in aggregatedLayers.items():
		#arcpy.AddMessage("In agg loop: " + str(finalIndex) + " " + str(key))
		# Normalize aggregated ECG matrices
		#arcpy.AddMessage("Before " + key + " min = " + str(numpy.nanmin(matrix)) + ", max = " + str(numpy.nanmax(matrix)))
		minv = numpy.nanmin(matrix)
		maxv = numpy.nanmax(matrix)
		matrix = (matrix - minv) / (maxv - minv)
		#arcpy.AddMessage("After " + key + " min = " + str(numpy.nanmin(matrix)) + ", max = " + str(numpy.nanmax(matrix)))

		# Save normalized aggregated ECG matrices as raster layers
		saveMatrixAsRaster(matrix, dsc, outputFolderPath + key + ".tif")

		# Sum aggregated ECG matrices
		if finalIndex == 0:
			finalMatrix = matrix
		else:
			finalMatrix += matrix
		finalIndex += 1

	#arcpy.AddMessage("After agg loop: " + str(finalIndex))
	# Calculate average (dividing sum by amount of aggregated matrices) and save "Aggregated ecological value" raster layer
	arcpy.AddMessage("----- Calculating Aggregated ecological value matrix.")
	#finalMatrix = finalMatrix / 4
	finalMatrix = finalMatrix / finalIndex
	minv = numpy.nanmin(finalMatrix)
	maxv = numpy.nanmax(finalMatrix)
	finalMatrix = (finalMatrix - minv) / (maxv - minv)
	saveMatrixAsRaster(finalMatrix, dsc, outputFolderPath +"Aggregated_Ecological_Value_matrix.tif")

	# Clean memory if possible
	gc.collect()
