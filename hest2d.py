#import modules
import numpy as np
import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d  import Axes3D  # for 3D plots
import math


#Analytical solution 
def trueu(z, x, t):
  u = np.sin(2*np.pi*z/Lz)*np.sin(4*np.pi*x/Lx)*np.cos(t*omega)
  return u

#seperate points into different groups
def getInside(zleftwall,zrightwall,xupwall,xdownwall,zvals,xvals,dz,dx):

  inside = []
  edge = []
  corner = []
  outside = []
  gaml = np.ones(len(zleftwall))
  gamr = np.ones(len(zrightwall))
  gamu = np.ones(len(xupwall))
  gamd = np.ones(len(xdownwall))

  
  
  for i in range(zvals.shape[0]):
    for j in range(zvals.shape[1]):

      count = 0

	#if the point is inside the boundary
      if zvals[i][j] >= zleftwall[i] and zvals[i][j] <= zrightwall[i] and    xvals[i][j] >= xdownwall[j] and xvals[i][j] <= xupwall[j]:

	inside.append((i,j))

	#if the point is the closest to a boundary
	if zvals[i][j] - zleftwall[i] < dz/2.001:
          edge.append((i,j,'l'))
	  gaml[i] = (zvals[i][j] - zleftwall[i])/dz
	  count+=1

        if zrightwall[i] - zvals[i][j] < dz/2.001: 
	  edge.append((i,j,'r'))
	  gamr[i] =  (zrightwall[i] - zvals[i][j])/dz
	  count+=1

	if xvals[i][j] - xdownwall[j] < dx/2.001:
	  edge.append((i,j,'d'))
	  gamd[j] = (xvals[i][j] - xdownwall[j])/dx
	  count+=1

	if xupwall[j] - xvals[i][j] < dx/2.001:
	  edge.append((i,j,'u'))
	  gamu[j] = (xupwall[j] - xvals[i][j])/dx
	  count+=1


	#if the point was along an edge, remove it from interior         
	if count >= 1:
	  inside.pop()

	#if the point was along two edges, it is in the corner.
	if count >= 2:
	  v = edge.pop()[2]
	  h = edge.pop()[2]
	  corner.append((i,j,v,h))
  

      else:
        outside.append((i, j))



  
  return inside, edge, corner, gaml, gamr, gamu, gamd, outside


def updateBz(oldE, oldBz):
	newBz = oldBz.copy()
	l, w = np.shape(oldBz)
	#first update everything with yee
	newBz[0:l, 0:w] = oldBz[0:l, 0:w] - dt/dx*(oldE[1:l + 1, 0:w] - oldE[0:l, 0:w])

	#now overwrite the edge and corner values with hesthaven scheme values

	#first the B values on the edge themselves
	for i, j, k in Bzedge:
	  if k == 'u':
	    #Bz on the top
	    newBz[i][j] = oldBz[i][j] + dt/dx*(2*oldE[i][j]/(2*Bzgamu[j] + 1))
	    if abs(Bzgamu[j]) < .000001:
	      newBz[i][j] = 0
	  elif k == 'd':
	    #Bz on the bottom
	    newBz[i][j] = oldBz[i][j] + dt/dx*(-2*oldE[i + 1][j]/(2*Bzgamd[j] + 1))
	    if abs(Bzgamd[j]) < .000001:
	      newBz[i][j] = 0

	  elif k == 'l' and abs(Bzgaml[i]) < .000001:
	    newBz[i][j] = 0
	  elif k == 'r' and abs(Bzgamr[i]) < .000001:
	    newBz[i][j] = 0

	#Now overwrite corners
	for i, j, v, h in Bzcorner:
	  if v == 'u':
	    #Bz on the top
	    newBz[i][j] = oldBz[i][j] + dt/dx*(2*oldE[i][j]/(2*Bzgamu[j] + 1))
	    if h == 'l' and (abs(Bzgaml[i]) < .000001 or abs(Bzgamu[j]) < .000001):
	      newBz[i][j] = 0
	    if h == 'r' and (abs(Bzgamr[i]) < .000001 or abs(Bzgamu[j]) < .000001):
	      newBz[i][j] = 0
	  elif v == 'd':
	    #Bz on the bottom
	    newBz[i][j] = oldBz[i][j] + dt/dx*(-2*oldE[i + 1][j]/(2*Bzgamd[j] + 1))
	    if h == 'l' and (abs(Bzgaml[i]) < .000001 or abs(Bzgamd[j]) < .000001):
	      newBz[i][j] = 0
	    if h == 'r' and (abs(Bzgamr[i]) < .000001 or abs(Bzgamd[j]) < .000001):
	      newBz[i][j] = 0
	

	#second, the B values next to E values on the edge

	for i, j, k in Eedge:
	  if k == 'u':
	    #E on top
	    newBz[i - 1][j] = oldBz[i - 1][j] + dt/dx*(oldE[i - 1][j]/(1 + Egamu[j]))
	  elif k == 'd':
	    #E on bottom
	    newBz[i][j] = oldBz[i][j] + dt/dx*(-oldE[i + 1][j]/(1 + Egamd[j]))


	#Now to the same for the corner values
	for i, j, v, h in Ecorner:
	  if v == 'u':
	    #E on top
	    newBz[i - 1][j] = oldBz[i - 1][j] + dt/dx*(oldE[i - 1][j]/(1 + Egamu[j]))
	  elif v == 'd':
	    #E on bottom
	    newBz[i][j] = oldBz[i][j] + dt/dx*(-oldE[i + 1][j]/(1 + Egamd[j]))


	return newBz


def updateBx(oldE, oldBx):
	newBx = oldBx.copy()
	l, w = np.shape(oldBx)

	#update all points using yee
	newBx[0:l, 0:w] = oldBx[0:l, 0:w] + dt/dz*(oldE[0:l, 1:w + 1] - oldE[0:l, 0:w])

	#  for i in range(1, l - 1):
	#    for j in range(w):  
	#      newBx[i][j] = oldBx[i][j] - dt/dz*(oldE[i][j + 1] - oldE[i][j])

	#now overwrite the edge and corner values with hesthaven scheme values

	#first the B values on the edge themselves
	for i, j, k in Bxedge:
	  if k == 'l':
	    #Bx on the left
	    newBx[i][j] = oldBx[i][j] + dt/dz*(2*oldE[i][j + 1]/(2*Bxgaml[i] + 1))
	    if abs(Bxgaml[i]) < .000001:#change to check within certain tolerance 10^-14, anywhere checking for == 0
	      newBx[i][j] = 0
	  elif k == 'r':
	    #Bx on the right
	    newBx[i][j] = oldBx[i][j] + dt/dz*(-2*oldE[i][j]/(2*Bxgamr[i] + 1))
	    if abs(Bxgamr[i]) < .000001:
	      newBx[i][j] = 0

	  elif k == 'u' and abs(Bzgamu[j]) < .000001:
	    newBx[i][j] = 0
	  elif k == 'd' and abs(Bzgamd[j]) < .000001:
	    newBx[i][j] = 0



	for i, j, v, h in Bxcorner:
	  if h == 'l':
	    #Bx on the left
	    newBx[i][j] = oldBx[i][j] + dt/dz*(2*oldE[i][j + 1]/(2*Bxgaml[i] + 1))
	    if v == 'u' and (abs(Bxgamu[j]) < .000001 or abs(Bxgaml[i]) < .000001):
	      newBx[i][j] = 0
	    if v == 'd' and (abs(Bxgamd[j]) < .000001 or abs(Bxgaml[i]) < .000001):
	      newBx[i][j] = 0

	  elif h == 'r':
	    #Bx on the right
	    newBx[i][j] = oldBx[i][j] + dt/dz*(-2*oldE[i][j]/(2*Bxgamr[i] + 1))
	    if v == 'u' and (abs(Bxgamu[j]) < .000001 or abs(Bxgamr[i]) < .000001):
	      newBx[i][j] = 0
	    if v == 'd' and (abs(Bxgamd[j]) < .000001 or abs(Bxgamr[i]) < .000001):
	      newBx[i][j] = 0


  	#second, the B values next to E values on the edge
	for i, j, k in Eedge:
	  if k == 'l':
	    #E on the left
	    newBx[i][j] = oldBx[i][j] + dt/dz*(oldE[i][j + 1]/(1 + Egaml[i]))

	  elif k == 'r':
	    #E on the right
	    newBx[i][j - 1] = oldBx[i][j - 1] + dt/dz*(-oldE[i][j - 1]/(1 + Egamr[i]))


	#Now do the same for the corners
	for i, j, v, h in Ecorner:
	  if h == 'l':
	    #E on the left
	    newBx[i][j] = oldBx[i][j] + dt/dz*(oldE[i][j + 1]/(1 + Egaml[i]))

	  elif h == 'r':
	    #E on the right
	    newBx[i][j - 1] = oldBx[i][j - 1] + dt/dz*(-oldE[i][j - 1]/(1 + Egamr[i]))


	return newBx


def updateE(oldE, oldBz, oldBx):
	newE = oldE.copy()
	l, w = np.shape(oldE)
	
	#First, update all with yee scheme
	newE[1:l - 1, 1:w - 1] = oldE[1:l - 1, 1:w - 1] + dt*c**2*(-(oldBz[1:l - 1, 1:w - 1] - oldBz[:l - 2, 1:w - 1])/dx + (oldBx[1:l - 1, 1:w - 1] - oldBx[1:l - 1, :w - 2])/dz)

	#Now update edges
	for i, j, k in Eedge:
	  
	  if k == 'l':
	    #E on the left
	    newE[i][j] = oldE[i][j] + dt*c**2*(Egaml[i]/((Egaml[i] + 1)*dz)*(oldBx[i][j + 1] - oldBx[i][j]) - (oldBz[i][j] - oldBz[i - 1][j])/dx)
	    if abs(Egaml[i]) < .000001:
	      newE[i][j] = 0

	  elif k == 'r':
	    #E on the right
	    newE[i][j] = oldE[i][j] + dt*c**2*(-Egamr[i]/((Egamr[i] + 1)*dz)*(oldBx[i][j - 2] - oldBx[i][j - 1]) - (oldBz[i][j] - oldBz[i - 1][j])/dx)
	    if abs(Egamr[i]) < .000001:
	      newE[i][j] = 0


	  elif k == 'u':
	    #E on the top
	    newE[i][j] = oldE[i][j] + dt*c**2*((oldBx[i][j] - oldBx[i][j - 1])/dz + Egamu[j]/((Egamu[j] + 1)*dx)*(oldBz[i - 2][j] - oldBz[i - 1][j]))
	    if abs(Egamu[j]) < .000001:
	      newE[i][j] = 0

	    
	  elif k == 'd':
	    #E on the bottom
	    newE[i][j] = oldE[i][j] + dt*c**2*((oldBx[i][j] - oldBx[i][j - 1])/dz - Egamd[j]/((Egamd[j] + 1)*dx)*(oldBz[i + 1][j] - oldBz[i][j]))
	    if abs(Egamd[j]) < .000001:
	      newE[i][j] = 0


	#now update corners
	for i, j, v, h in Ecorner:

	  if v == 'u' and h == 'l':
	    #E in upper left corner
	    newE[i][j] = oldE[i][j] + dt*c**2*(Egaml[i]/((Egaml[i] + 1)*dz)*(oldBx[i][j + 1] - oldBx[i][j]) + Egamu[j]/((Egamu[j] + 1)*dx)*(oldBz[i - 2][j] - oldBz[i - 1][j]))
	    if abs(Egaml[i]) < .000001 or abs(Egamu[j]) < .000001:
	      newE[i][j] = 0


	  elif v == 'd' and h == 'l':
	    #E in lower left corner
	    newE[i][j] = oldE[i][j] + dt*c**2*(Egaml[i]/((Egaml[i] + 1)*dz)*(oldBx[i][j + 1] - oldBx[i][j]) - Egamd[j]/((Egamd[j] + 1)*dx)*(oldBz[i + 1][j] - oldBz[i][j]))
	    if abs(Egaml[i]) < .000001 or abs(Egamd[j]) < .000001:
	      newE[i][j] = 0

	  elif v == 'u' and h == 'r':
	    #E in upper right corner
	    newE[i][j] = oldE[i][j] + dt*c**2*(-Egamr[i]/((Egamr[i] + 1)*dz)*(oldBx[i][j - 2] - oldBx[i][j - 1]) + Egamu[j]/((Egamu[j] + 1)*dx)*(oldBz[i - 2][j] - oldBz[i - 1][j]))
	    if abs(Egamr[i]) < .000001 or abs(Egamu[j]) < .000001:
	      newE[i][j] = 0

	  elif v == 'd' and h == 'r':
	    #E in lower right corner
	    newE[i][j] = oldE[i][j] + dt*c**2*(-Egamr[i]/((Egamr[i] + 1)*dz)*(oldBx[i][j - 2] - oldBx[i][j - 1]) - Egamd[j]/((Egamd[j] + 1)*dx)*(oldBz[i + 1][j] - oldBz[i][j]))
	    if abs(Egamr[i]) < .000001 or abs(Egamd[j]) < .000001:
	      newE[i][j] = 0

	for i, j in Eoutside:
	  newE[i][j] = 0

	return newE
	

#constants
startoffset = 1/(3*np.e)
endoffset = 1/(3*np.e)
zstart = -2#0 - startoffset #makes e on edge but not corner
zend = 15#4*np.pi + endoffset
xstart = -2#0 - startoffset
xend = 15#4*np.pi + endoffset
numEzstart = 30
numExstart = 30
Lz = zend - zstart
Lx = xend - xstart
T = 15.0#1.*2*np.pi/np.sqrt(5)#4*np.pi/np.sqrt(5)#1 full period
t = 0.0
c = 1.0
omega = c*np.sqrt((2*np.pi/(4*np.pi))**2 + (4*np.pi/(4*np.pi))**2)



zE, dz = np.linspace(zstart, zend, numEzstart + 1, retstep = True)
xE, dx = np.linspace(xstart, xend, numExstart + 1, retstep = True)



Lz = 4*np.pi
Lx = 4*np.pi


#dtstart = 1/(c*np.sqrt(1/dz**2 + 1/dx**2))/2
dtstart = dx*4/13
numstepsstart = math.ceil(T/dtstart)
errors = []
dts = []


#####Start trial loop here

for trialno in range(2, 3):

	t = 0.0
	#deal with changing variables
	numEz = numEzstart*2**trialno + 1
	numEx = numExstart*2**trialno + 1
	dt = dtstart/2**trialno
	numsteps = int(numstepsstart*2**trialno)

	numBz = numEz - 1
	numBx = numEx - 1

	#make E values
	zE, dz = np.linspace(zstart, zend, numEz, retstep = True)
	xE, dx = np.linspace(xstart, xend, numEx, retstep = True)
	#x coordinates of Bz s
	xBz = np.zeros(numBx)
	#z coordinates of Bx s
	zBx = np.zeros(numBz)


	for i in range(numBx): 
	  xBz[i] = (xE[i] + xE[i + 1])/2

	for i in range(numBz):
	  zBx[i] = (zE[i] + zE[i + 1])/2


	Eleftwall = np.ones(len(xE))
	Erightwall = np.ones(len(xE))
	Edownwall = np.ones(len(zE))
	Eupwall = np.ones(len(zE))

	Eleftwall *= 0
	Erightwall *= 4*np.pi
	Edownwall *= 0
	Eupwall *= 4*np.pi


	Bxleftwall = np.ones(len(xE))
	Bxrightwall = np.ones(len(xE))
	Bxdownwall = np.ones(len(zE) - 1)
	Bxupwall = np.ones(len(zE) - 1)

	Bxleftwall *= 0
	Bxrightwall *= 4*np.pi
	Bxdownwall *= 0
	Bxupwall *= 4*np.pi

	Bzleftwall = np.ones(len(xE) - 1)
	Bzrightwall = np.ones(len(xE) - 1)
	Bzdownwall = np.ones(len(zE))
	Bzupwall = np.ones(len(zE))

	Bzleftwall *= 0
	Bzrightwall *= 4*np.pi
	Bzdownwall *= 0
	Bzupwall *= 4*np.pi



	#make grid of E values
	zzE, xxE = np.meshgrid(zE, xE)
	#make grid of Bz values, shares z values of E
	zzBz, xxBz = np.meshgrid(zE, xBz)
	#make grid of Bx values, shares x values of E
	zzBx, xxBx = np.meshgrid(zBx, xE)

	#seperate points

	Einterior, Eedge, Ecorner, Egaml, Egamr, Egamu, Egamd, Eoutside = getInside(Eleftwall, Erightwall, Eupwall, Edownwall, zzE, xxE, dz, dx)

	Bzinterior, Bzedge, Bzcorner, Bzgaml, Bzgamr, Bzgamu, Bzgamd, Bzoutside = getInside(Bzleftwall, Bzrightwall, Bzupwall, Bzdownwall, zzBz, xxBz, dz, dx)
	
	Bxinterior, Bxedge, Bxcorner, Bxgaml, Bxgamr, Bxgamu, Bxgamd, Bxoutside = getInside(Bxleftwall, Bxrightwall, Bxupwall, Bxdownwall, zzBx, xxBx, dz, dx)

#	Egamu = Egamu*0
#	Bxgamu = Egamu*0




	#Now, to deal with base cases of E0, and B12

	
	#E0 = np.sin(2*np.pi*zzE/Lz)*np.sin(4*np.pi*xxE/Lx)
	E0 = np.e**(-80*(zzE/(4*np.pi) - .4)**2)*np.e**(-80*(xxE/(4*np.pi) - .4)**2)
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1, projection = '3d')

	ax.plot_surface(zzE, xxE, E0, vmin = -1.1, vmax = 1.1, cmap=plt.cm.jet)
	ax.set_zlim3d(-1.1, 1.1)

	fig.show()
	  
	##Now deal with B0 and B1/2
	Bz0 = np.zeros(np.shape(zzBz))
	Bx0 = np.zeros(np.shape(zzBx))


	Bz12 = Bz0.copy()

	for i,j in Bzinterior:
	  Bz12[i][j] = Bz0[i][j] - dt/(dx*2) * (E0[i+1][j] - E0[i][j])


	#hestavhen for edges
	#for Bz points, use hest if on up or down edge, otherwise use yee
	for i,j,k in Bzedge:

	  if k == 'd':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * ( -2*E0[i+1][j] / ( (2*Bzgamd[j] +1)) )
	  elif k == 'u':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * ( 2*E0[i][j] / ( (2*Bzgamu[j] + 1)) )
	  else:
            Bz12[i][j] = Bz0[i][j] - dt/(dx*2) * (E0[i+1][j] - E0[i][j])

	
	#update corner points
	for i,j,v,h in Bzcorner:
	  if v == 'd':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * ( -2*E0[i+1][j] / ( (2*Bzgamd[j] +1)) )
	  elif v == 'u':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * ( 2*E0[i][j] / ( (2*Bzgamu[j] + 1)) )


	#now use Hest on B points next to E points on edges
  	for i,j,k in Eedge:
	  if k == 'd':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * (-E0[i+1][j]/(1+Egamd[j]))
          elif k == 'u':
            Bz12[i-1][j] = Bz0[i-1][j] + dt/(2*dx) * (E0[i-1][j]/(1+Egamu[j]))

	#do the same for E corner points
	for i,j,v,h in Ecorner:
	  if v == 'd':
            Bz12[i][j] = Bz0[i][j] + dt/(2*dx) * (-E0[i+1][j]/(1+Egamd[j]))
          elif v == 'u':
            Bz12[i-1][j] = Bz0[i-1][j] + dt/(2*dx) * (E0[i-1][j]/(1+Egamu[j]))


	##################################################Now Bx12

	Bx12 = Bx0.copy()

	for i,j in Bxinterior:
	  Bx12[i][j] = Bx0[i][j] + dt/(dz*2) *(E0[i][j+1] - E0[i][j])

	for i,j,k in Bxedge:
	  if k == 'l':
            Bx12[i][j] = Bx0[i][j] + dt/(2*dz) * ( 2*E0[i][j+1] / ( (2*Bxgaml[i] + 1)) )
	  elif k == 'r':
            Bx12[i][j] = Bx0[i][j] + dt/(2*dz) * ( -2*E0[i][j] / ( (2*Bxgamr[i] + 1)) )
	  else:
            Bx12[i][j] = Bx0[i][j] + dt/(dz*2) *(E0[i][j+1] - E0[i][j])


	for i,j,v,h in Bxcorner:
	  if h == 'l':
            Bx12[i][j] = Bx0[i][j] + dt/(2*dz) * ( 2*E0[i][j+1] / ( (2*Bxgaml[i] + 1)) )
	  elif h == 'r':
            Bx12[i][j] = Bx0[i][j] + dt/(2*dz) * ( -2*E0[i][j] / ( (2*Bxgamr[i] + 1)) )

	
	for i,j,k in Eedge:
	  if k == 'l':
            Bx12[i][j] = Bx0[i][j] +  dt/(2*dz) *(E0[i][j+1]/(1 + Egaml[i]))
	  elif k == 'r':
            Bx12[i][j-1] = Bx0[i][j-1] + dt/(2*dz) * (-E0[i][j-1]/(1+Egamr[i]))

	for i,j,v,h in Ecorner:
	  if h == 'l':
            Bx12[i][j] = Bx0[i][j] +  dt/(2*dz) *(E0[i][j+1]/(1 + Egaml[i]))
	  elif h == 'r':
            Bx12[i][j-1] = Bx0[i][j-1] + dt/(2*dz) * (-E0[i][j-1]/(1+Egamr[i]))


	
	


	newE = E0
	oldBz = Bz12.copy()
	oldBx = Bx12.copy()



	#main animation loop
	for step in range(numsteps):
	  
	  newE = updateE(newE, oldBz, oldBx) 
	  newBz = updateBz(newE, oldBz)
	  newBx = updateBx(newE, oldBx)
	  oldBz = newBz
	  oldBx = newBx


	  if(step%5 == 0):
	    ax.collections = []  
	    ax.plot_surface(zzE, xxE, newE, vmin = -1.1, vmax = 1.1, cmap=plt.cm.jet)
	 
	    ax.set_zlim3d(-1.1, 1.1)
	    ax.set_xlabel('"z" axis')
	    ax.set_ylabel('"x" axis')
	    plt.draw()
	    fig.show()
	    

	  t += dt




	for i, j in Eoutside:
	  newE[i][j] = 0

	ax.collections = []  
	ax.plot_surface(zzE, xxE, newE, vmin = -1.1, vmax = 1.1, cmap=plt.cm.jet)
	ax.set_zlim3d(-1.1, 1.1)
	ax.set_xlabel('"z" axis')
	ax.set_ylabel('"x" axis')
	plt.draw()
	fig.show()

	#error analysis
#	exactu = trueu(zzE, xxE, t)
#	
#	e = np.max(abs(exactu - newE))

	#only test for error in the interior and border 
	e = 0
	maxi = -100
	maxj = -100
	for l in range(len(Einterior)):
	  i, j = Einterior[l]
	  temp = abs(trueu(zzE[i][j], xxE[i][j], t) - newE[i][j])
	  if temp > e: 
	    e = temp
	    maxi = i
	    maxj = j

	for l in range(len(Eedge)):
	  i, j, k = Eedge[l]
	  temp = abs(trueu(zzE[i][j], xxE[i][j], t) - newE[i][j])
	  if temp > e: 
	    e = temp

	for l in range(len(Ecorner)):
	  i, j, v, h = Ecorner[l]
	  temp = abs(trueu(zzE[i][j], xxE[i][j], t) - newE[i][j])
	  if temp > e: 
	    e = temp


#	#This is L1 error as opposed to L infinity error
#	for l in range(len(Einterior)):
#	  i, j = Einterior[l]
#	  e += dx*dz*abs(trueu(zzE[i][j], xxE[i][j], t) - newE[i][j])


	#This is first order relative order
#	top = 0
#	bottom = 0
#	for i, j in Einterior:
#	  top += abs(trueu(zzE[i][j], xxE[i][j], t) - newE[i][j])
#	  bottom += abs(trueu(zzE[i][j], xxE[i][j], t))
#	e = top/bottom


	errors.append(e)
	dts.append(dt)
	print 't is ', t
	print 'dt is ', dt
	print 'e is ', e
	print dt/dx
	
#	fig5 = plt.figure()
#	xx = fig5.add_subplot(1, 1, 1, projection = '3d')

#	xx.set_zlim3d(-1.1, 1.1)
#	xx.set_xlabel('"z" axis')
#	xx.set_ylabel('"x" axis')

#	truesol = trueu(zzE, xxE, t)

#	xx.plot_surface(zzE, xxE, truesol, vmin = -1.1, vmax = 1.1, cmap=cm.jet)
#	xx.set_zlim3d(-1.1, 1.1)
#	plt.draw()
#	fig5.show()
	print 'number of elements in Eedge is: ', len(Eedge)

	

	


if len(errors) > 1:

  fig = plt.figure()
  ax = fig.add_subplot(1, 1, 1)

  ax.plot(np.log2(dts), np.log2(errors))
  ax.set_xlabel('log(dt)')
  ax.set_ylabel('log(error)')
  ax.set_xlim(-6, 6)
  ax.set_ylim(-12, 0)
  

  plt.show()

  print 'log fraction of errors is ', np.log2(errors[-2]/errors[-1])


#fig2 = plt.figure()
#bx = fig2.add_subplot(1,1,1,projection = '3d')
#bx.plot_surface(zzBz,xxBz,Bz12, cmap=cm.jet)
#bx.set_xlabel('"z" axis')
#bx.set_ylabel('"x" axis')
#plt.show()


#fig3 = plt.figure()
#cx = fig3.add_subplot(1,1,1)
#cx.plot(zzBx[1],newBx[1])

#fig3.show()


