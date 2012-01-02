#! /usr/bin/python
from scipy import signal
import numpy as np

# An adaptor class for use in Lattice Wave Digital Filters
class Adaptor():
	type = 0
	dtype = np.float
	alpha = np.zeros((1,),dtype=np.float)
	coeff_scale = 0
	def __init__(self, gamma=0, coeff_scale=1, dtype=np.float):
		self.set_params(gamma, coeff_scale, dtype)
		
	def set_params(self, gamma=0, coeff_scale=1, dtype=np.float):
		# Define the adaptor
		self.dtype = dtype
		self.alpha = np.zeros((1,),self.dtype)
		self.coeff_scale = coeff_scale
		# Choose the adaptor type
		if (gamma > 0.5) & (gamma < 1):
			self.type = 0
			temp = (1<<self.coeff_scale)*(1 - gamma)
		elif (gamma > 0) & (gamma <= 0.5):
			self.type = 1
			temp = (1<<self.coeff_scale)*gamma
		elif (gamma >= -0.5) & (gamma < 0):
			self.type = 2
			temp = (1<<self.coeff_scale)*abs(gamma)
		elif (gamma > -1) & (gamma < -0.5):
			self.type = 3
			temp = (1<<self.coeff_scale)*(1 + gamma)
		else:
			print('Invalid gamma.')
		if self.dtype != np.float:
			temp = temp.round()
		self.alpha[0] = temp
		print gamma, self.alpha, self.type
		
	def process(self, wave_in):
		_diff = np.zeros((1,),dtype=self.dtype)
		_waveout = np.zeros((2,),self.dtype)
		_wavein = np.array((wave_in),self.dtype)
		if self.type == 0:
			_diff = _wavein[0] - _wavein[1]
			_waveout[1] = self.multiply_by_alpha(_diff) + _wavein[1]
			_waveout[0] = _waveout[1] - _diff
		elif self.type == 1:
			_diff = self.multiply_by_alpha(_wavein[1] - _wavein[0])
			_waveout[1] = _diff + _wavein[0]
			_waveout[0] = _diff + _wavein[1]
		elif self.type == 2:
			_diff = self.multiply_by_alpha(_wavein[0] - _wavein[1])
			_waveout[1] = _diff - _wavein[0]
			_waveout[0] = _diff - _wavein[1]
		elif self.type == 3:
			_diff = _wavein[1] - _wavein[0]
			_waveout[1] = self.multiply_by_alpha(_diff) - _wavein[1]
			_waveout[0] = _waveout[1] - _diff
		else:
			print('Invalid type.')
		return _waveout
		
	def multiply_by_alpha(self, input):
		if self.dtype == np.float:
			return (self.alpha[0]*input)/(1<<self.coeff_scale)
		else:
			return ((np.cast[np.int64](self.alpha[0])*input)+(1<<(self.coeff_scale-1)))>>self.coeff_scale
		
# Class for constructing and building a Lattice Wave Digital Filter
class Filter():
	order = 0
	adaptors = []
	registers = []
	output = [0,0]
	dtype = np.float
	coeff_scale = 0
	def __init__(self, order=3, wn=0.01, coeff_scale=1, rp=1, rs=60, ftype='butter', dtype=np.float):
		if order%2 == 0:
			print 'No even order filters allowed.'
			return
		# Filter settings
		self.order = order
		self.wn = wn
		self.rp = rp
		self.rs = rs
		self.dtype = dtype
		self.ftype = ftype
		self.coeff_scale = coeff_scale
		gammas = []
		psi = np.tan(np.pi*self.wn/2.0)
		psi2 = psi*psi
		if self.ftype == 'butter':
			(z,p,k) = signal.iirfilter(self.order, psi, btype='lowpass', analog=1, ftype='butter', output='zpk')
			filter_roots = np.sort(p)
		elif self.ftype == 'bessel':
			print 'Please avoid using Bessel filters as they don\'t translate well to LWDFs.'
			(z,p,k) = signal.iirfilter(self.order, psi, btype='lowpass', analog=1, ftype='bessel', output='zpk')
			filter_roots = np.sort(p)
		elif self.ftype == 'cheby1':
			(z,p,k) = signal.iirfilter(self.order, psi, rp=1, btype='lowpass', analog=1, ftype='cheby1', output='zpk')
			filter_roots = np.sort(p)
		elif self.ftype == 'cheby2':
			(z,p,k) = signal.iirfilter(self.order, psi, rs=self.rs, btype='lowpass', analog=1, ftype='cheby2', output='zpk')
			filter_roots = np.sort(p)
		else:
			print 'Invalid filter type.'
			return
		# Separate the real pole from the complex poles
		real_index = 0;
		for i in range(self.order):
			if abs(filter_roots[i].imag) <= 1e-16:
				real_index = i
				break
		complex_roots = np.concatenate((filter_roots[0:real_index],filter_roots[real_index+1:]))
		# Put in the real pole's gamma value
		h_B = -1.0*filter_roots[real_index].real
		gammas.append((1.0 - h_B) / (1.0 + h_B))
		# Calculate coefficients of the individual Hurwitz polynomials
		for i in (range((order-1)/2)):
			h_A = -2.0*(complex_roots[2*i].real)
			h_B = abs(complex_roots[2*i])**2
			gammas.append((h_A - h_B - 1.0)/(h_A + h_B + 1.0))
			gammas.append((1.0 - h_B) / (1.0 + h_B))
		# Construct filter
		for i in range(self.order):
			self.adaptors.append(Adaptor(gammas[i],self.coeff_scale,self.dtype))
			self.registers.append(0)
			
	def push(self, value):
		# Insert the new value into both branches
		if self.dtype != np.float:
			output = [value<<self.coeff_scale,value<<self.coeff_scale]
		else:
			output = [value*(1<<self.coeff_scale),value*(1<<self.coeff_scale)]
		
		# Process real pole in the top path
		wave_in = [output[0], self.registers[0]]
		wave_out = self.adaptors[0].process(wave_in)
		self.registers[0] = wave_out[1]
		output[0] = wave_out[0]
		
		# Process top branch
		for i in range(int(self.order/4)):
			x = [4*i+3, 4*i+4]
			
			wave_in = [self.registers[x[0]], self.registers[x[1]]]
			wave_out = self.adaptors[x[1]].process(wave_in)
			self.registers[x[1]] = wave_out[1]
			
			wave_in = [output[0], wave_out[0]]
			wave_out = self.adaptors[x[0]].process(wave_in)
			self.registers[x[0]] = wave_out[1]
			
			output[0] = wave_out[0]
		
		# Process bottom branch
		for i in range(int((self.order+1)/4)):
			x = [4*i+1, 4*i+2]
			
			wave_in = [self.registers[x[0]], self.registers[x[1]]]
			wave_out = self.adaptors[x[1]].process(wave_in)
			self.registers[x[1]] = wave_out[1]
			
			wave_in = [output[1], wave_out[0]]
			wave_out = self.adaptors[x[0]].process(wave_in)
			self.registers[x[0]] = wave_out[1]
			
			output[1] = wave_out[0]
		
		# Save the output and return it
		self.output = output
		return output

if __name__ == "__main__":
	from scipy import fftpack, signal
	import matplotlib.pyplot as plt
	# Fetch example data
	[ir,red,amb] = np.load('test_data.npy')
	
	# Make signal zero mean and scale it up to avoid rounding errors
	# ir = ir-(sum(ir)/len(ir))
	# ir = ir<<10
	
	# Define filter (cutt-off in pi radians per seconds; limit is 0.5pi)
	order = 5
	wn = 50.0/500.0
	coeff_scale = 16
	dtype = np.int32
	rs = 70
	dtype=np.int32
	ftype = 'cheby2'
	# ftype = 'butter'
	lwdfilter = Filter(order,wn,coeff_scale,rs=70,ftype=ftype,dtype=dtype)
	# lwdfilter = Filter(order,wn,coeff_scale,'butter',dtype=np.int32)
	# lwdfilter = Filter(order,wn,coeff_scale,'butter',dtype=dtype)
	
	# Perform sample-by-sample filtering
	ir_out = []
	for x in ir:
		temp = lwdfilter.push(x)
		if dtype==np.float:
			ir_out.append((temp[0]+temp[1])/(2.0*(1<<coeff_scale)))
		else:
			ir_out.append((temp[0]+temp[1]+(1<<coeff_scale))>>(coeff_scale+1))
		
	# Compare with SciPy equivalent
	if ftype=='cheby1':
		(B,A) = signal.iirfilter(order, wn, rp=1,btype='lowpass', analog=0, ftype='cheby1', output='ba')
	elif ftype=='cheby2':
		(B,A) = signal.iirfilter(order, wn, rs=rs,btype='lowpass', analog=0, ftype='cheby2', output='ba')
	else:
		(B,A) = signal.iirfilter(order, wn, btype='lowpass', analog=0, ftype=ftype, output='ba')
	zf = signal.lfiltic(B,A,np.zeros(len(A)-1), np.zeros(len(B)-1))
	ir_comp, zf = signal.lfilter(B,A,ir,zi=zf)
	
	# Plot the data
	plt.plot(ir, 'k-', linewidth=2.0)
	plt.plot(ir_comp, 'r-', linewidth=2.0)
	plt.plot(ir_out, 'b-', linewidth=2.0)
	plt.grid()
	plt.legend(('Raw data', 'SciPy filtered - float32', 'LWDF - int32'), loc='upper left')
	plt.show()
