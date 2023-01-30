import numpy as np
import matplotlib.pyplot as plt 


def plot_aesthetics(ax,fontsize=14,fontweight='bold',xax=True,yax=True):
	"""
	This function is used to make bold and increase the font size of all plot tick markers
	"""

	if xax is True:
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

	if yax is True:
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize=fontsize)
			tick.label1.set_fontweight(fontweight)

			tick.label2.set_fontsize(fontsize=fontsize)
			tick.label2.set_fontweight(fontweight)

c = 3*10**10 # cm / s
mp = 1.67e-24 # grams 


E_dot_iso = 10**(51) # erg/s
tw = 30 # seconds
E_iso = E_dot_iso * tw # erg
# E_iso = 10**(53) # erg
Gamma_bar = 100

R_coll = 10**(15) # cm

gamma_rs = 5
mass_rs = 10**(30) # g
gamma_ej =  15 # g
mass_ej = 10**(31)

rho_est = E_dot_iso / (4*np.pi*R_coll**2 * Gamma_bar**2 * c**3)
rho_obs_frame_est = rho_est / gamma_ej**2
n_est = rho_est / mp
n_obs_frame_est = n_est / gamma_ej**2

m_ej = 10**(30) # g
l_ej_est = m_ej / rho_est / (np.pi * R_coll**2)
l_ej_obs_frame_est = l_ej_est / gamma_ej


print("rho_est (obs frame) = {:.3e}".format(rho_obs_frame_est))
print("n_est (obs frame) = {:.3e}".format(n_obs_frame_est))
print("l_ej_est (obs frame) = {:.3e}".format(l_ej_obs_frame_est))


def l_ej(n,mass,rad):
	return mass/(np.pi * rad**2 * n * mp)

def t_shock_cross(width,gamma_1,gamma_2):
	return width * 2 * gamma_1**2 * gamma_2**2 / c / (gamma_2**2 - gamma_1**2)

def gamma_m(gamma_1,mass_1,gamma_2,mass_2):
	return np.sqrt(gamma_1*gamma_2)*np.sqrt(( (mass_1*gamma_1) + (mass_2*gamma_2) ) / ( (mass_1*gamma_2) + (mass_2*gamma_1) ))

n_arr = np.logspace(np.log10(n_est)-2,np.log10(n_est)+2)
# print("n_arr =", n_arr)
n_obs_frame_arr = n_arr/ (gamma_ej**2)

l_ej_arr = l_ej(n_arr, m_ej, R_coll)
# print("l_ej_arr =", l_ej_arr)
l_ej_obs_frame_arr = l_ej_arr / gamma_ej

t_shock_cross_arr = t_shock_cross(l_ej_arr,gamma_rs,gamma_ej)
# print("t_shock_cross_arr =", t_shock_cross_arr)
t_shock_cross_obs_frame_arr = t_shock_cross_arr / (2*gamma_m(gamma_rs, mass_rs, gamma_ej, mass_ej)**2)

ax = plt.figure().gca()
at = ax.twinx()

ax.plot(l_ej_obs_frame_arr, t_shock_cross_obs_frame_arr,color="C0",alpha=0.7)
at.plot(l_ej_obs_frame_arr, n_obs_frame_arr,color="C1",alpha=0.7)

ax.hlines(y=1e5,xmin=l_ej_obs_frame_arr[0],xmax=l_ej_obs_frame_arr[-1],color="C0",alpha=0.7,linestyle="dashed")
at.hlines(y=n_obs_frame_est,xmin=l_ej_obs_frame_arr[0],xmax=l_ej_obs_frame_arr[-1],color="C1",alpha=0.6,linestyle="dotted")
# ax.vlines(x=5.5e13,ymin=t_shock_cross_obs_frame_arr[0],ymax=t_shock_cross_obs_frame_arr[-1],color="grey",alpha=0.6,linestyle="dashed")

ax.margins(x=0,y=0)
ax.set_yscale('log')
ax.set_xscale('log')
at.set_yscale('log')
at.set_xscale('log')

fontsize = 14
fontweight = "bold"
ax.set_xlabel('Shell Width \n (observer frame, cm)',fontsize=fontsize,fontweight=fontweight)
ax.set_ylabel('Shock Crossing Time \n (observer frame, sec)',fontsize=fontsize,fontweight=fontweight,color="C0")
at.set_ylabel("Shell Density\n" r"(observer frame, cm$^{-3}$)",fontsize=fontsize,fontweight=fontweight,color="C1")
at.yaxis.set_label_position("right")

# ax.yaxis.offsetText.set_fontsize(fontsize)
# ax.xaxis.offsetText.set_fontsize(fontsize)
# ax.ticklabel_format(style='sci', scilimits=(0,0))

plot_aesthetics(ax,fontsize=fontsize,fontweight=fontweight)
plot_aesthetics(at,fontsize=fontsize,fontweight=fontweight)
plt.tight_layout()

# plt.savefig("shock-crossing-time.png")

plt.show()
