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
#Holz data
sigma_rho=1
sigma_z=1
Amp=10
rho0=0.0
z0=0.0

def fq(rho,z):

    return Amp*(rho)**2*np.exp(-(rho-rho0)**2/sigma_rho**2-(z-z0)**2/sigma_z**2)
#Initial values
Nz=41
Nrho=51

zmin=-3.5
zmax=3.5

rhomax=3.5
rhomin=0.0




rho,Drho,D2rho=Cheb_grid_and_D(rhomin,rhomax,Nrho)
z,Dz,D2z=Cheb_grid_and_D(zmin,zmax,Nz)

#create linear operator L
dz=np.kron(np.identity((Nrho)),Dz)
d2z=np.kron(np.identity((Nrho)),D2z)
drho=np.kron(Drho,np.identity((Nz)))
d2rho=np.kron(D2rho,np.identity((Nz)))

d2rhoq=np.dot(d2rho,q)
d2zq=np.dot(d2z,q)

fneu=rhov**2*(d2rhoq+d2zq)/4

c1rho=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')
c2rho=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')
c2z=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')
cf=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')

b1rho=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')
b2rho=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')
b2z=np.zeros([Nrho*Nz,Nrho*Nz], dtype='double')

for j in range(0,Nrho):
        for i in range(0,Nz):
            cf[i+j*Nz][i+j*Nz]=fneu[i+j*Nz]


            c1rho[i+j*Nz][i+j*Nz]=rhov[i+j*Nz]

            c2rho[i+j*Nz][i+j*Nz]=rhov[i+j*Nz]**2
            c2z[i+j*Nz][i+j*Nz]=rhov[i+j*Nz]**2


            b1rho[i+j*Nz][i+j*Nz]=(rho[j]**2+z[i]**2)
            b2rho[i+j*Nz][i+j*Nz]=rho[j]
            b2z[i+j*Nz][i+j*Nz]=z[i]

boundrho=np.matmul(b1rho,drho)+b2rho
boundz=np.matmul(b1rho,dz)+b2z

L=np.matmul(c1rho,drho)+np.matmul(c2rho,d2rho)+np.matmul(c2z,d2z)+cf

q=np.zeros(Nrho*Nz)


rhov=np.zeros(Nrho*Nz)

for j in range(0,Nrho):
          for i in range(0,Nz):

                  q[i+j*Nz]=fq(rho[j],z[i])
                  rhov[i+j*Nz]=rho[j]
