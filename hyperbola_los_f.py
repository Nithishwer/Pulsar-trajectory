import matplotlib.pyplot as plt
import numpy as np

M0 = 1.989e30								#Solar mass
Mp = 1.4*M0 								#Mass of pulsar
Mc = 0.3*M0								#Mass of companion
i = 60*np.pi/180							#Angle of inclination
P0 = 0.5*60*60*24							#Orbital Period	
G = 6.67408e-11								#Gravitaional constant
v_inf = 5000.0

#Calculation of projected periastron distance
mu = G * (Mp + Mc)
#a = mu/np.power(v_inf, 2)
ap = 336431967.08858365 #a * np.sin(i)
a = ap/np.sin(i)
nu = np.power(mu/np.power(a,3), 0.5)

#omega is the longitude of periastron
#e is the eccentricity
def hyperbola_f(omega, e):
	eterm = np.power(e,2) - 1

	print('ap: ', ap)
	print('nu: ', nu)

	fmax = np.arccos(-1/e)
	f = np.arange(-fmax+fmax/100.0,fmax-fmax/100.0, fmax/100.0)
	fdeg = f*180.0/np.pi

	#calculate line of sight radius
	rl =  ap * eterm * np.sin(f + omega) / (1 + e*np.cos(f))

	#calculate the line of sight velocity
	vl = nu * ap / np.sqrt(eterm) * (np.cos(f + omega) + e*np.cos(omega))

	#calculate line of sight acceleration
	al = - np.power(nu, 2) * ap * np.sin(f + omega) * np.power((1 + e*np.cos(f)), 2) * np.power(eterm, -2)

	#calculate line of sight jerk
	jl = - np.power(nu, 3) * ap * np.power(eterm, -7.0/2.0) * np.power(1 + e*np.cos(f), 3) * (np.cos(f + omega) + e*np.cos(omega) - 3.0*e*np.sin(f + omega)*np.sin(f))



	# Plot the points using matplotlib
	plt.figure(1)
	plt.plot(fdeg, rl)

	plt.figure(2)
	plt.plot(fdeg, vl)

	plt.figure(3)
	plt.plot(fdeg, al)

	plt.figure(4)
	plt.plot(fdeg, jl)
#End of function ellipse_f

#Call function
hyperbola_f(0.0, 1.2)
hyperbola_f(0.0, 1.5)

#Set plot properties for each figure
plt.figure(1)
plt.legend(['e = 1.2', 'e = 1.5'])
plt.title("Variation of line-of-sight-radius with true anomaly")
plt.xlabel('f (deg)')
plt.ylabel('r$_l$ (m)')
plt.figure(2)
plt.legend(['e = 1.2', 'e = 1.5'])
plt.title("Variation of line-of-sight-velocity with true anomaly")
plt.xlabel('f (deg)')
plt.ylabel('v$_l$ (m s$-2$)')
plt.figure(3)
plt.legend(['e = 1.2', 'e = 1.5'])
plt.title("Variation of line-of-sight-acceleration with true anomaly")
plt.xlabel('f (deg)')
plt.ylabel('a$_l$ (m s$^-2$)')
plt.figure(4)
plt.legend(['e = 1.2', 'e = 1.5'])
plt.title("Variation of line-of-sight-jerk with true anomaly")
plt.xlabel('f (deg)')
plt.ylabel('j$_l$ (m s$^-3$)')

#Display plot
plt.show()
