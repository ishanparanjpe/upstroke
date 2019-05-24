import operator
from pandas import read_csv
from pathlib import Path
import operator

import numpy as np
import pandas as pd

import sys
import getopt
from os import path
import os
import PeakFinder
import csv
from optparse import OptionParser
import os

#########
# Parameters
# peak_start_percentage - what percentage of the maximum slope to be considered 
#						start of upstroke
# peak_frame_cutoff - peaks detected before this frame number will be considered
#						noise and discarded
#######
peak_start_percentage = .8
peak_frame_cutoff = 10
upstroke_frame_cutoff = 10


"""
data= pd.read_csv('data_sep20/170503_6.csv')

cell1= PeakFinder.DataSeries('data_sep20', '170503_3','Cell10', data['Cell33'])
cell1.smooth()
cell1.findPeaks()
cell1.removeEarlyLatePeaks()

cell1.findUpstroke()
cell1.findMaxPeak()
cell1.findAmplitude()
cell1.findDuration()

#cell1.removeShortPeaks()
cell1.plotPeaks_transformed(display=True)
"""

def pipeline(current_cell):
	smooth_data= current_cell.smooth()

	current_cell.findPeaks()
	current_cell.removeEarlyLatePeaks()
	upstroke_frame, upstroke_height=current_cell.findUpstroke()
	max_peak = current_cell.findMaxPeak()
	amplitude, peak_height= current_cell.findAmplitude()
	duration=current_cell.findDuration()
	area = current_cell.integrate()
	return (smooth_data,upstroke_frame,upstroke_height,max_peak,amplitude,peak_height,duration, area)


def single_cell(path,cell, save):
	data = pd.read_csv(path)
	data = data[cell]
	
	folder_name = path.split('/')[0]
	filename= path.split('/')[1]
	new_path = folder_name+'/'+cell+'_'+filename
	data.to_csv(new_path, index=False, header=True)
	print(new_path)
	batchmode(new_path, save, display=True, single_cell=True)

	"""
	print('filename:', cell+filename)
	failed_cells=[]
	output= folder_name+str('/output_'+cell+'_'+filename)

	with open(output,'w') as f1:
		writer= csv.writer(f1, delimiter=',',lineterminator = '\n')
		writer.writerow(['Peak Start Percentage',peak_start_percentage])
		writer.writerow(['Peak frame cutoff',peak_frame_cutoff])
		writer.writerow(['Upstroke frame cutoff',upstroke_frame_cutoff])

		writer.writerow(['Cell #', 'Upstroke Frame','Upstroke Value', 'Peak Frame Number', 'Peak Value' ,'Amplitude', 'Duration', 'Area'])
		current_cell = PeakFinder.DataSeries(folder_name, filename, cell, data[cell],peak_start_percentage,peak_frame_cutoff, upstroke_frame_cutoff)

		try:
			
			
			smooth_data,upstroke_frame,upstroke_height,max_peak,amplitude,peak_height,duration, area = pipeline(current_cell)
			writer.writerow([cell, upstroke_frame,upstroke_height, max_peak, peak_height, amplitude,duration,area])

		except:
			print('error')
			failed_cells.append(cell)
			
		writer.writerow(['Failed Cells:']+ failed_cells)
		f1.close()
		current_cell.plotPeaks_transformed(display=True, save=save)
			
		"""

			

def batchmode(path, save, display=False, upstroke_path= None, single_cell=False):
	data = pd.read_csv(path)
	dir_name = path.split('/')[0]+'/'+path.split('/')[1].split('.')[0]
	print('dir:', dir_name)
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	
	filename= path.split('/')[1]
	cells= list(data)
	output= dir_name+str('/output_'+filename)
	failed_cells=[]
	transformed_data =[]
	if(single_cell):
		os.remove(path)
	upstrokes= pd.DataFrame({'0':[0], '1':[0]})
	if(upstroke_path!=None):
		upstrokes = pd.read_csv(upstroke_path, header=None)
		print(upstrokes)
		print(upstrokes.ix[:,0].tolist())

	with open(output,'w') as f1:
		writer= csv.writer(f1, delimiter=',',lineterminator = '\n')
		writer.writerow(['Peak Start Percentage',peak_start_percentage])
		writer.writerow(['Peak frame cutoff',peak_frame_cutoff])
		writer.writerow(['Upstroke frame cutoff',upstroke_frame_cutoff])

		writer.writerow(['Cell #', 'Upstroke Frame','Upstroke Value', 'Peak Frame Number', 'Peak Value' ,'Amplitude', 'Duration', 'Area'])
		plot=None
		for cell in cells:
			
			try:
			
				current_cell = PeakFinder.DataSeries(dir_name, filename, cell, data[cell],peak_start_percentage,peak_frame_cutoff, upstroke_frame_cutoff)
				smooth_data= current_cell.smooth()
				transformed_data.append(np.concatenate([[cell],smooth_data]))

				current_cell.findPeaks()
				current_cell.removeEarlyLatePeaks()
				max_peak = current_cell.findMaxPeak()
				#upstroke_frame =0
				
				if(cell in upstrokes.ix[:,0].tolist()):
					current_cell.upstroke = upstrokes.loc[upstrokes.iloc[:,0]==cell,1].values[0]
					current_cell.start_positions=[current_cell.upstroke]
					upstroke_frame= current_cell.upstroke
				else:
					upstroke_frame=current_cell.findUpstroke()
				upstroke_height= current_cell.findUpstrokeHeight()
				amplitude, peak_height= current_cell.findAmplitude()
				duration=current_cell.findDuration()
				area = current_cell.integrate()
				#smooth_data,upstroke_frame,upstroke_height,max_peak,amplitude,peak_height,duration, area= pipeline(current_cell)
			
			except:	
				print('error')
				failed_cells.append(cell)
				continue
			plot= current_cell.plotPeaks_transformed( save=save)
		
			writer.writerow([cell, upstroke_frame,upstroke_height, max_peak, peak_height, amplitude,duration,area])
			print('written')
		
			
		writer.writerow(['Failed Cells:']+ failed_cells)

	output_transformed_name = dir_name +str('/transformed_data_'+filename)
	transformed_data = pd.DataFrame(transformed_data).transpose()
	transformed_data.to_csv(output_transformed_name,index=False, header=False)
	
	if(display and plot is not None):
		plot.show()
			



def main(argv):
	num_arguments = len(argv)
	print (num_arguments)
	if(num_arguments==1):
		batchmode(argv[0])
	

if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option('-c','--cell', dest= 'cell_num', help='Invidual cell to process. Do not use this flag if you want to process entire file')
	parser.add_option('-s','--savefig',type=None, action= 'store_true',dest ='save' ,help='Type this if flag if you want to save the figure.')
	#parser.add_option('-u', 'upstrokeframe', type=None,dest = 'upstroke_frame', help ='Compute parameters from designated upstroke frame.')
	(options,args)=parser.parse_args()
	num_args = len(args)
	print(args)
	if(len(args)>2):
		print('Error: Incorrect number of arguments. filename should be the last parameter')
	else:
		if(options.cell_num==None):
			if(num_args ==2):
				batchmode(args[0],options.save, upstroke_path = args[1])

			else:
				batchmode(args[0],options.save)

		else:
			print(options.cell_num)
			
			single_cell(args[0],options.cell_num,options.save)
			


