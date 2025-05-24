"""Brill Wave initial data via spectral methods. """
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#create grid
def Cheb_grid_and_D(xmin,xmax,Nx):
     x=np.zeros(Nx)
     cx=(xmax-xmin)/2
     dx=(xmax+xmin)/2
     for j in range(0,Nx):
         x[j]=-np.cos(j*np.pi/(Nx-1))*cx+dx
     Dx=np.zeros((Nx,Nx))
     Dx[0][0]=-(2*(Nx-1)**2+1)/6
     Dx[Nx-1][Nx-1]=(2*(Nx-1)**2+1)/6
     for i in range(0,Nx):
         if(i==0 or i==Nx-1):
             ci=2
         else:
             ci=1
         for j in range(0,Nx):
             if(j==0 or j==Nx-1):
                 cj=2
             else:
                 cj=1
             if(i==j and i!=0 and i!=Nx-1):
                Dx[i][i]=-((x[i]-dx)/cx)/(1-((x[i]-dx)/cx)**2)*0.5
             if(i!=j):
                Dx[i][j]=ci/cj*(-1)**(i+j)/(x[i]/cx-x[j]/cx)

     Dx=Dx/cx
     return x,Dx,np.dot(Dx,Dx)

