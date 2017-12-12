import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math

tList   = [] # Time
xList   = [] # x
yList   = [] # y
vList   = [] # v
aList   = [] # a
hList   = [] # Height above Earth surface
TList   = [] # Temperature
pList   = [] # (Static) Pressure
qList   = [] # Dynamic Pressure
taList  = [] # True Anomaly
rhoList = [] # Density

# Init
G    = 6.67e-11 	# Gravitational constant

# Sun
MSun = 1.989e30
RSun = 696000000

# Earth
MEarth = 5.976e24
REarth = 6378000
DEarth = 149.6e9

# Mars
MMars = 6.4171e23
RMars = 3396200
DMars = 227.92e9

# Central object
M = MEarth
R = REarth
D = DEarth

t    = 0
y    = 0
ta   = 0
vx   = 0
done = False
i    = 0
SOI  = D*(M/MSun)**(2/5) # Gravitational Sphere of Influence - for escape determination

lastTa = ta
numberOfOrbits = 0

m     = 4200        # Mass of satellite, Dragon spacecraft m=4200

A     = math.pi * (3.7/2)**2 # Cross section area of spacecraft, Dragon spacecraft Ø=3.7 m
Cd    = 0.65 # Aerodynamic drag coefficient
h     = 400000		# Initial orbit height
v     = 7600		# Initial satellite orbit velocity
r     = R + h
x     = r
vy    = v
g0    = math.sqrt(G*M/R**2)
dt    = 1 # Time step in seconds

apoapsis  = False
periapsis = False
r_old = r

print("Start distance ",r/1000," km")

while not done:

	if not apoapsis and not periapsis:
		if r<r_old:
			apoapsis=True
		if r>r_old:
			periapsis=True

	if periapsis and not apoapsis and r<r_old:
		print("Orbit " , numberOfOrbits , "apoapsis " , ("%.0f" % (r/1000)) , " km (h =",("%.0f" % ((r-R)/1000))," km)")
		apoapsis  = True
		periapsis = False

	if apoapsis and not periapsis and r>r_old:
		print("Orbit " , numberOfOrbits , "periapsis" , ("%.0f" % (r/1000)) , " km (h =",("%.0f" % ((r-R)/1000))," km), arg. of periapsis=",("%.2f" % (math.degrees(ta)))," deg")
		apoapsis  = False
		periapsis = True
		numberOfOrbits = numberOfOrbits + 1

	if numberOfOrbits > 100:
		done = True

	if r < R:
		print("Impact speed:   " , ("%.2f" % (v)) , " m/s")
		done = True

	if r > SOI:
		print("Escape speed:   " , ("%.2f" % (v)) , " m/s")
		done = True

	t = i*dt
	
	# Gravity
	g = -G*M/r**2

	# Modified gravity from general relativity - from http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node116.html
	#g = -G*M/r**2 - 3 * G * M * ( 3.3e23*r**2*() )**2 / (299792458**2*r**4)

	# Atmosphere density

	# Earth
	# https://physics.stackexchange.com/questions/121809/why-e-in-the-formula-for-air-density
	rho = 1.1225 * math.e**(h*-0.0001)

	# Mars
	# https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
	#rho = 0.014517 * math.e**(h*-0.000073385)

	# Velocity
	v = math.sqrt(vx**2+vy**2)

	# Dynamic pressure
	q = 0.5 * rho * v**2

	# Drag
	D = -q * Cd * A
	Dx = D * vx/v
	Dy = D * vy/v

	# Gravity
	gx = g * math.cos(ta)
	gy = g * math.sin(ta)

	# Acceleration
	ax = gx + Dx/m
	ay = gy + Dy/m

	# Velocity
	vx = vx + ax*dt
	vy = vy + ay*dt

	# Position
	x = x + vx*dt
	y = y + vy*dt

	r_old = r
	r = (x**2+y**2)**0.5
	h = r - R
	
	# True Anomaly
	lastTa = ta
	ta = math.atan2(y,x)
	if ta < 0:
		ta = ta + 2*math.pi

	# Flight Path Angle
	fpa = math.atan2(vy,vx) - ta - math.pi/2

	a = math.sqrt(ax**2+ay**2)

	# Append to lists
	tList.append(t)
	xList.append(x)
	yList.append(y)
	hList.append(h)
	aList.append(a/g0)
	vList.append(v)
	qList.append(q)
	taList.append(math.degrees(ta))

	i = i + 1

#print("Orbital period:    " , ("%.2f" % (t/3600)) , " hours")
#print("Perigee height:    " , ("%.2f" % ( ( min( [ max(xList) , -min(xList) ] ) - R ) / 1000 ) ) , " km")
#print("Apogee height:     " , ("%.2f" % ( ( max( [ max(xList) , -min(xList) ] ) - R ) / 1000 ) ) , " km")
print("Max. acceleration: " , ("%.2f" % ( ( max( aList )) ) ) , " G")

# Planet
f, ax = plt.subplots(figsize=(6, 6))
circle = Circle((0,0),R,color='#0000FF')
ax.add_patch(circle)

# Orbit
plt.figure(1)
plt.plot(xList,yList)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(False)
plt.axis('scaled')

# Graphs
plt.figure(2)

plt.subplot(411)
plt.plot(tList,hList)
plt.ylabel('Altitude / [m]')
plt.grid(True)

plt.subplot(412)
plt.plot(tList,vList)
plt.ylabel('Speed / [m/s]')
plt.grid(True)

plt.subplot(413)
plt.plot(tList,aList)
plt.ylabel('Acceleration / [G]')
plt.grid(True)

plt.subplot(414)
plt.plot(tList,qList)
plt.xlabel('Time / [s]')
plt.ylabel('Dynamic Pressure / [Pa]')
plt.grid(True)

plt.show()
