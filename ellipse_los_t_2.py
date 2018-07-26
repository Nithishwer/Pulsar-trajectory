import matplotlib.pyplot as plt
import numpy as np

M0 = 1.989e30								#Solar mass
Mp = 1.4*M0 								#Mass of pulsar
Mc = 0.3*M0								#Mass of companion
i = 60*np.pi/180							#Angle of inclination
P0 = 0.5*60*60*24							#Orbital Period	
G = 6.67408e-11								#Gravitaional constant
tau = 0.0

#omega is the longitude of periastron
#e is the eccentricity
def ellipse_t(omega, e):
	#Calculation of projected periastron distance
	ap = np.power(np.power(P0/(2*np.pi), 2) * G * (Mp + Mc), 1.0/3.0) * Mc * np.sin(i) / (Mp + Mc)
	print(ap);
	eterm = 1 - np.power(e,2)

	f = np.arange(0, 2*np.pi, 0.05)
	fdot = (2*np.pi/P0)*np.power(eterm,-1.5)*np.power(1.0+e*np.cos(f),2)
	fdotdot = (2*np.pi/P0) * np.power(eterm, -1.5) * 2.0 * (1 + e*np.cos(f)) * (-e*np.sin(f))
	fdotdotdot = -2*e*(2*np.pi/P0) * np.power(eterm, -1.5) * (e*np.cos(2*f) + np.cos(f))
	expr = np.power((1-e)/(1+e), 0.5) * np.tan(f/2.0)
	E = 2.0 * np.arctan2(np.absolute(expr), np.sign(expr))
	M = E - e*np.sin(E)
	t = M/(2*np.pi/P0) + tau;

	#calculate line of sight radius
	rl = ap * eterm * np.sin(f + omega) / (1 + e*np.cos(f))
	r = ap * eterm  / (1 + e*np.cos(f))

	#calculate the line of sight velocity
	#vrl = (2*np.pi/P0) * ap / np.sqrt(eterm) * (np.cos(f + omega) + e*np.cos(omega))
	vr = (2*np.pi/P0)*ap*e*np.sin(f) / np.sqrt(eterm)
	vt = r*fdot
	vl = np.sqrt(vr*vr + vt*vt) * np.sin(f + omega)

	#calculate line of sight acceleration
	#arl = - np.power((2.0*np.pi/P0), 2) * ap * np.sin(f + omega) * np.power((1 + e*np.cos(f)), 2) * np.power(eterm, -2)
	ar = (2*np.pi/P0)*ap*e*np.cos(f) * fdot / np.sqrt(eterm)
	at = vr*fdot + r*fdotdot
	al = np.sqrt(ar*ar + at*at) * np.sin(f + omega)

	#calculate line of sight jerk
	#jrl = - np.power(2*np.pi/P0, 3) * ap * np.power(eterm, -7.0/2.0) * np.power(1 + e*np.cos(f), 3) * (np.cos(f + omega) + e*np.cos(omega) - 3.0*e*np.sin(f + omega)*np.sin(f))
	jr = - np.power((2*np.pi/P0), 3) * ap * e * np.sin(f) * (1 - e*np.cos(f)) * np.power(1 + e*np.cos(f), 3) / np.power(eterm, 7.0/2.0)
	jt = ar*fdot + 2.0*vr*fdotdot + r*fdotdotdot
	jl = np.sqrt(jr*jr + jt*jt) * np.sin(f + omega)



	# Plot the points using matplotlib
	plt.figure(1)
	plt.plot(t, rl)

	plt.figure(2)
	plt.plot(t, vl)

	plt.figure(3)
	plt.plot(t, al)

	plt.figure(4)
	plt.plot(t, jl)
#End of function ellipse_t


#Call function
ellipse_t(0.0, 0.5)
ellipse_t(0.0, 0.65)

#Set plot properties for each figure
plt.figure(1)
plt.legend(['e = 0.5', 'e = 0.65'])
plt.xlabel('t (s)')
plt.ylabel('r$_l$ (m)')
plt.figure(2)
plt.legend(['e = 0.5', 'e = 0.65'])
plt.xlabel('t (s)')
plt.ylabel('v$_l$ (m s${-2}$)')
plt.figure(3)
plt.legend(['e = 0.5', 'e = 0.65'])
plt.xlabel('t (s)')
plt.ylabel('a$_l$ (m s$^{-2}$)')
plt.figure(4)
plt.legend(['e = 0.5', 'e = 0.65'])
plt.xlabel('t (s)')
plt.ylabel('j$_l$ (m s$^{-3}$)')

#Display plot
plt.show()
