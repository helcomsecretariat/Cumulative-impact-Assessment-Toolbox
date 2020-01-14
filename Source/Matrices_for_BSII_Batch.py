"""

"""

import arcpy
from arcpy.sa import *
import numpy
import os
import sys
import csv
import gc

if __name__ == '__main__':

	# ------------- Set up --------------

	sensitivityScoresCsv = arcpy.GetParameterAsText(0) # Sensitivity scores CSV file
	coefMatrixCsv = arcpy.GetParameterAsText(1) # Coeficients matrix CSV file
	outputFolder = arcpy.GetParameterAsText(2) # Folder to save final matrices

	ssTable = None
	ssEcIndex = None
	if (sys.version_info > (3, 0)): # Python 3
		mode = 'r'
	else: # Python 2
		mode = 'rb'

	# Read sensitivity scores table
	with open(sensitivityScoresCsv, mode) as csvFileSS:
		ssTable = csv.reader(csvFileSS, delimiter=';')
		ssTable = list(ssTable)

		# finding index of EC layers (ec) in second row of sensitivity scores csv file
		i = 0
		for cell in ssTable[1]:
			if cell == "EC CODE":
				ssEcIndex = i
				break
			i += 1

		if ssEcIndex is None:
			arcpy.AddError("EC CODE header is not present in second row of sensitivity scores file. --- ")
			sys.exit(0)

	# Read coeficients table
	with open(coefMatrixCsv, mode) as csvFileCoef:
		coefTable = csv.reader(csvFileCoef, delimiter=';')
		coefTable = list(coefTable)

		# finding indices of EC layers (ec) and EC groups (ecg) in second row of csv file
		ecgIndex = None
		ecIndex = None
		i = 0
		for cell in coefTable[1]:
			if cell == "EC CODE":
				ecIndex = i
			if cell == "GROUP BY":
				ecgIndex = i
			i += 1

		if ((ecIndex is None) or (ecgIndex is None)):
			arcpy.AddError("EC CODE or GROUP BY header is not present in second row of coeficients file. --- ")
			sys.exit(0)

		ecgCodeList = []
		for row in coefTable:
			ecgCode = row[ecgIndex]
			if (((ecgCode.startswith("ECG_")) or (ecgCode.startswith("ECSG_"))) and (ecgCode not in ecgCodeList)):
				ecgCodeList.append(ecgCode)


		for ecg in ecgCodeList:
			outputFileName = None
			index = 0
			for cell in coefTable[1]:
				# coeficients are stored in columns with headers srarting with EV_ or ES_
				# headers will be a csv file names
				if ((cell.startswith("EV_")) or (cell.startswith("ES_"))):
					outputFileName = cell + "_" + ecg + ".csv"

					outputMatrix = []
					outputMatrix.append(ssTable[0])
					outputMatrix.append(ssTable[1])

					for row in coefTable:
						if ((row[ecIndex].startswith("EC_")) and (row[ecgIndex] == ecg)):
							ec = row[ecIndex]
							score = row[index]

							for ssRow in ssTable:
								if ssRow[ssEcIndex] == ec:
									newRow = []
									ssI = 0
									for ssCell in ssRow:
										if ssI < 2:
											newRow.append(ssCell)
										else:
											newRow.append(float(ssCell) * float(score))
										ssI += 1
									outputMatrix.append(newRow)
									break

					# Save new sensitivity scores in CSV file
					if (sys.version_info > (3, 0)): # Python 3
						with open(outputFolder + "//" + outputFileName, mode='w', newline='') as matrixFile:
							matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
							for outputMatrixRow in outputMatrix:
								matrixWriter.writerow(outputMatrixRow)
					else: # Python 2
						with open(outputFolder + "//" + outputFileName, mode='wb') as matrixFile:
							matrixWriter = csv.writer(matrixFile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
							for outputMatrixRow in outputMatrix:
								matrixWriter.writerow(outputMatrixRow)
				index += 1
