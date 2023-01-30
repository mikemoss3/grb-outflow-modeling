import numpy as np
import matplotlib.pyplot as plt

c_cm = 2.99792*np.power(10.,10.); # speed of light, in cm/s

# class SHELL(object):

# 	def __init__(self,Gamma):
# 		self.Gamma = Gamma


class PLOT(object):

	def __init__(self, kappa=0.2, E_dot_iso=1e52, sigma=0.3, delte=1,Gamma_min=100,Gamma_max=400,num=100):
		self.ax = plt.figure().gca()

		self.gamma_range = np.linspace(Gamma_min,Gamma_max,num)

		self.kappa = kappa
		self.E_dot_iso = E_dot_iso
		self.sigma = sigma
		self.delte = delte

	def plot_aesthetics(self,ax,fontsize=14,fontweight='bold'):
		"""
		This function is used to make bold and increase the font size of all plot tick markers
		"""

		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

		for tick in ax.xaxis.get_minor_ticks():
			tick.label1.set_fontweight(fontweight)
			tick.label2.set_fontweight(fontweight)


		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontweight(fontweight)
			tick.label2.set_fontweight(fontweight)

	def plot_radii_comp(self, save_pref=None, fontsize=14,fontweight='bold',logscale=True):
		"""
		Make Plot 
		"""
		if logscale is True:
			self.ax.set_yscale('log')
			# self.ax.set_xscale('log')

		# Make plots look good
		self.plot_aesthetics(self.ax,fontsize=fontsize,fontweight=fontweight)

		self.ax.set_xlabel('Gamma',fontsize=fontsize,fontweight=fontweight)
		self.ax.set_ylabel('Radius (cm)',fontsize=fontsize,fontweight=fontweight)
		
		# ax.grid(True)
		self.ax.margins(x=0,y=0)
		plt.tight_layout()

		if save_pref is not None:
			plt.savefig('{}-is-vs-phot-rad.png'.format(save_pref)) 

	def phot_rad(self, Gamma, E_dot_iso = None, kappa=None, sigma=None):
		"""
		Calculate photospheric radius of a shell with a Lorentz factor Gamma 
		"""
		if sigma is None:
			sigma = self.sigma
		if kappa is None:
			kappa = self.kappa
		if E_dot_iso is None:
			E_dot_iso = self.E_dot_iso

		r_phot = (2.9 * 10**(15)) * (kappa/0.2) * (E_dot_iso / 1e52) / (1+sigma) / (Gamma/10)**3

		return r_phot

	def add_lines_phot_rad(self, linestyle="solid", alpha=1,color="r"):
		self.ax.plot(self.gamma_range, self.phot_rad(Gamma = self.gamma_range), alpha=alpha,linestyle=linestyle,color=color)

	def is_coll_rad(self,Gamma_1,Gamma_2,delte):
		"""
		Calculate the radius of collision of two shells where Gamma_2 > Gamma_1 with a seperation length dell
		"""

		r_coll = c_cm * delte * (2*Gamma_1**2 - 1)*(2*Gamma_2**2 -1) / 2 / (Gamma_2**2 - Gamma_1**2)
		return r_coll

	def add_lines_is_coll_rad(self, linestyle="solid", alpha=1):

		vals = self.is_coll_rad(Gamma_1=self.gamma_range[0], Gamma_2=self.gamma_range[1:], delte=self.delte)

		self.ax.plot(self.gamma_range[1:], vals , alpha=alpha,linestyle=linestyle)



if __name__ == "__main__":

	# Make radius comparison plot
	plot = PLOT(Gamma_min=50,Gamma_max=500,num=100)

	plot.add_lines_phot_rad(linestyle="dashed")

	plot.delte=1
	plot.add_lines_is_coll_rad()
	plot.delte=0.1
	plot.add_lines_is_coll_rad()
	plot.delte=2
	plot.add_lines_is_coll_rad()

	plot.delte=1
	plot.E_dot_iso = 1e54
	plot.add_lines_phot_rad(color="C4",linestyle="dashed")

	plot.plot_radii_comp(save_pref="2022-01-03")

	plt.show()