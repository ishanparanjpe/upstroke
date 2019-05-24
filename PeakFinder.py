import numpy as np
from scipy.signal import butter, filtfilt,argrelextrema
from residual_analysis import residual_analysis
import peakutils
from matplotlib import pyplot 
from scipy.integrate import quad
import matplotlib.patches as mpatches


#peak_start_percentage = .8
#frame_start_position_cutoff = 10
min_peak_height=0.01


class DataSeries:
	def __init__(self, foldername, filename, cell_name, data, peak_start_percentage, frame_start_position_cutoff, upstroke_frame_cutoff):
		self.peak_start_percentage= peak_start_percentage
		self.frame_start_position_cutoff= frame_start_position_cutoff
		self.min_peak_height=min_peak_height
		self.upstroke_frame_cutoff = upstroke_frame_cutoff
		self.foldername = foldername #name of folder that data file is in
		self.filename=filename #filename 
		self.cell_name=cell_name #cell name (should be first row in csv file)
		self.data=data #single column with all values for a single cell
		self.smooth_data=[] #transformed data
		self.peaks=[] #list of peaks found in smooth data 
		self.start_positions= [] #upstroke start positions (one for each peak)
		self.peakHeights=[] 
		self.max_peak=0 #frame number of peak with largest value 
		self.upstroke=0 #earliest start position
		self.duration=0 #frame number where curve comes back to upstroke - upstroke frame number
		self.area = 0 # integrated area under peak
		self.amplitude=0
	def smooth(self):
		num_frames = len(self.data)
		trimmed_data = np.array(self.data.dropna())	
		offset = num_frames -len(trimmed_data)
		# print(trimmed_data)
		# fc_opt = residual_analysis(trimmed_data,show=False)
		# print('fc_opt',fc_opt)	
		# C = 0.802  # for dual pass; C = (2**(1/npasses) - 1)**0.25
		# freq = 1

		order=2 # order of the filter
		btype='lowpass' # type
		F_sampling=5   # my sampling freq
		Nyqvist_freq=F_sampling/2 # Nyqvist frequency in function of my sampling frequency
		cut_off_Hz=.3 # my cutoff frequency in Hz
		cutoff_frequency=cut_off_Hz/Nyqvist_freq    # normalized cut_off frequency
		analog=False #digital filter
		b, a= butter(order,cutoff_frequency,btype, analog)
		

		# b, a = butter(2, (fc_opt/C) / (freq / 2),'low',analog=True)
		yf = filtfilt(b, a, trimmed_data)
		diff = np.diff(yf, n=1)
		
		self.smooth_data= yf
		print('done smoothing')
		return self.smooth_data
	def findPeaks(self):
		self.peaks = peakutils.indexes(self.smooth_data, thres=.3)
		print('peaks:',self.peaks)
		#print(self.peaks)

	def removeEarlyLatePeaks(self):
		# remove all peaks that are before the frame cutoff (peaks due to noise)
		self.peaks= [peak for peak in self.peaks if peak>self.frame_start_position_cutoff and peak<len(self.data)-1]

	def findUpstroke(self):
		diff= np.diff(self.smooth_data) #1st derivatives
		for peak in self.peaks:
			i = peak-2
			while diff[i]>0 and i>=0:
			
				i = i-1
			base_frame= i+1 #this is the point when derivative = 0
			print('base frame', base_frame)

			max_slope = max(diff[base_frame:peak]) #maximum slope obtained between base frame and peak
			print('max slope', max_slope)
			
			i=i+1
			while(diff[i]<max_slope*self.peak_start_percentage):
				
				i=i+1
				
		#	if(i>frame_start_position_cutoff-10):	
			self.start_positions.append(i)
		print('start positions',self.start_positions, self.start_positions==[0])
		#self.start_positions= argrelextrema(self.smooth_data, np.less)[0]  #argrelextreme returns a tuple so this flattens it
		if len(self.start_positions)==1:
			self.upstroke=self.start_positions[0]
		else:
			self.upstroke= min([x for x in self.start_positions if x>self.upstroke_frame_cutoff and x<self.max_peak])
		return self.upstroke
	def findUpstrokeHeight(self):
		return self.smooth_data[self.upstroke]
	def findMaxPeak(self):
		peaks = [self.smooth_data[peak] for peak in self.peaks]
		self.max_peak = self.peaks[max( (v, i) for i, v in enumerate(peaks) )[1]] # index is the maximum frame number in peaks array
		print('max peak:',self.max_peak)
		return self.max_peak
	def findAmplitude(self):
	
		#self.peakHeights= list(map(lambda x,y: self.smooth_data[x]-self.smooth_data[y],self.peaks, self.start_positions))
		self.amplitude = self.smooth_data[self.max_peak]- self.smooth_data[self.upstroke]
		return self.amplitude, self.smooth_data[self.max_peak]
	def findDuration(self):
		i = self.max_peak
		while(i<len(self.data) and self.smooth_data[i]>self.smooth_data[self.upstroke]):
			i=i+1
		
		self.duration = i-self.upstroke
		return self.duration
	# remove all peaks which are smaller than min_peak_height
	def removeShortPeaks(self):
		temp= [(x,y,z) for (x,y,z) in zip(self.peaks, self.start_positions,self.peakHeights) if z>self.min_peak_height ]

		if(len(temp)>0):
			self.peaks,self.start_positions,self.peakHeights= zip(*temp)
			print(self.peaks,self.start_positions,self.peakHeights)
	def integrate(self):
		self.area = sum((self.smooth_data[x] - self.smooth_data[self.upstroke]) for x in list(range(self.upstroke,self.upstroke+self.duration+1)))
		return self.area
	def plotPeaks_transformed(self, save=True):
		fig=pyplot.figure()
		ax=pyplot.subplot(2,1,1)
		#plot legend
		peak_patch = mpatches.Patch(color='blue', label='Peak')
		upstroke_patch = mpatches.Patch(color='red', label='Main Upstroke(used in calculations)')
		other_upstroke_patch = mpatches.Patch(color='green', label='Other Upstrokes')

		pyplot.legend(handles= [peak_patch,upstroke_patch,other_upstroke_patch],prop={'size': 6},loc=1,ncol=1,shadow=True, fancybox=True)

		pyplot.title('Transformed ({}) {}'. format(self.filename, self.cell_name), fontsize= 15)
		print('Upstroke',self.upstroke)
		# plot all upstrokes with frame number labels
		for i in self.start_positions:
			pyplot.scatter(i, self.smooth_data[i],color='green',s=20)
			ax.annotate('%s' % i,xy=(i,self.smooth_data[i]),textcoords='data')
		pyplot.scatter(self.max_peak, self.smooth_data[self.max_peak], color= 'blue', s = 20)
		pyplot.scatter(self.upstroke, self.smooth_data[self.upstroke], color= 'red', s = 20)
		
		#pyplot.xticks(np.arange(min(x), max(x)+1, 5.0), fontsize=5)
		x = range(len(self.data))
		pyplot.plot(x,self.smooth_data , color=[1, 0, 0, .5],linewidth=2, label='Opt. filtered')
		#pyplot.scatter(peaks[0], [yf[i] for i in sorted_start_positions][0], color= 'blue', s = 200)
		ax=pyplot.subplot(2,1,2)
		pyplot.title('Original Data ({}) {}'. format(self.filename, self.cell_name), fontsize= 15)
		# plot all upstrokes with frame number labels
		for i in self.start_positions:
			pyplot.scatter(i, self.data[i],color='green',s=20)
			ax.annotate('%s' % i,xy=(i,self.data[i]),textcoords='data')
		pyplot.plot(x, self.data)
		pyplot.scatter(self.max_peak, self.data[self.max_peak], color= 'blue', s = 20)
		pyplot.scatter(self.upstroke, self.data[self.upstroke], color= 'red', s = 20)

		#plotting table at the bottom
		
		colnames= ('cell', 'upstroke_frame','upstroke_height', 'max_peak', 'peak_height', 'amplitude','duration','area')	
		table_data = [self.cell_name, self.upstroke, round(self.smooth_data[self.upstroke],3), self.max_peak,round(self.smooth_data[self.max_peak],3),round(self.amplitude,3),self.duration,round(self.area,3)]
		#table_data = [[i] for i in table_data]
		
		table= pyplot.table(cellText= [table_data],colLabels= colnames, cellLoc='bottom',bbox=[0, -0.27, 1, 0.15])
		table.auto_set_font_size(False)
		table.set_fontsize(5)
		table.scale(2,2)

		if save:	
			fig.savefig(self.foldername+'/'+self.foldername+'_'+self.cell_name,dpi=500)		
		
		return pyplot
