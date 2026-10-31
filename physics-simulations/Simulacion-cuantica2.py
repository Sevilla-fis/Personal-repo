import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from functools import partial
from scipy.linalg import expm
from mpl_toolkits.mplot3d import Axes3D

#Constantes

L=1
h=1.05*10**-34
m1=9.11e-31 #Electrón
m2=1.62e-25 #Boson Z0
E1=(np.pi)**2*h**2/(2*m1*L**2)
E2=(np.pi)**2*h**2/(2*m2*L**2)
FSt1prima=np.array([1/np.sqrt(2),np.sqrt(1/2)])
FSt2prima=np.array([1/np.sqrt(4),np.sqrt(3/4)])
BSt0prima=np.array([1/2,np.sqrt(1/2),1/2])

#Funciones de onda

def phi_x(x,n):
    phi_x=(2/L)**(1/2)*np.sin(n*np.pi*x/L)
    return phi_x
def PB(x1,x2,n1,n2):
    phisim=1/np.sqrt(2)*(phi_x(x1,n1)*phi_x(x2,n2)+phi_x(x1,n2)*phi_x(x2,n1))
    return abs(phisim)**2
def PF(x1,x2,n1,n2):
    phiantisim=1/np.sqrt(2)*(phi_x(x1,n1)*phi_x(x2,n2)-phi_x(x1,n2)*phi_x(x2,n1))
    return abs(phiantisim)**2

#Gráfica
x1,x2 = np.linspace(0,L,10),np.linspace(0,L,10)
MFF=np.eye(len(x1),len(x2))
MFB=np.eye(len(x1),len(x2))
for j in range(len(x1)):
    for k in range(len(x2)):
        MFF[j][k]=PF(x1[j],x2[k],1,2)
for j in range(len(x1)):
    for k in range(len(x2)):
        MFB[j][k]=PB(x1[j],x2[k],1,2)
plt.imshow(#MODIFICAR AQUI SI QUEREMOS MOSTRAR FERMION O BOSON (MFB-BOSON/MFF-FERMION)
    MFB
    ####################################
    , cmap='hot',origin='lower', interpolation='spline36')
plt.xlabel('x1'),plt.ylabel('x2')

plt.show()