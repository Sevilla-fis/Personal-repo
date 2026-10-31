import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from functools import partial
from scipy.linalg import expm
from mpl_toolkits.mplot3d import Axes3D

plt.ion()

#Constantes

tnorm1=1e2
tnorm2=1e8

#He normalizado el tiempo para que se pueda apreciar el movimiento de cada onda sin que sea o demsaido rapido o demasiado lento, ya que el Bosón al ser mucho más pesado se mueve mucho más lento

n=1
L=1
h=1.05*10**-34
m1=9.11e-31 #Electrón
m2=1.62e-25 #Boson Z0
E1=(np.pi)**2*h**2/(2*m1*L**2)
E2=(np.pi)**2*h**2/(2*m2*L**2)
#MODIFICAR ESTADO INICIAL AQUI



FSt0prima=np.array([1/np.sqrt(2),np.sqrt(1/2)])
BSt0prima=np.array([1/2,np.sqrt(1/2),1/2])



################################################
sigma_xF=np.array([[0,1],[1,0]])
sigma_yF=np.array([[0,-1j],[1j,0]])
sigma_zF=np.array([[1,0],[0,-1]])
sigma_xB=np.array([[0,1,0],[1,0,1],[0,1,0]])
sigma_yB=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
sigma_zB=np.array([[1,0,0],[0,0,0],[0,0,-1]])
#Funciones temporales, función espacial y esfera de Bloch

def phi_t1(t):
    phi_t1=np.exp(-E1*t*tnorm1*1j/h)
    return phi_t1
def phi_t2(t):
    phi_t2=np.exp(-E2*t*tnorm2*1j/h)
    return phi_t2
def phi_x(x):
    phi_x=(2/L)**(1/2)*np.sin(n*np.pi*x/L)
    return phi_x
def esfera(a,b):
    x=2*(a.conjugate()*b).real
    y=2*(a.conjugate()*b).imag
    z=abs(a)**2-abs(b)**2
    return x,y,z
(x0,y0,z0)=esfera(FSt0prima[0],FSt0prima[1])
#Gráfica


fig4 = plt.figure(figsize=(7,7))
ax4 = fig4.add_subplot(111, projection='3d')
u = np.linspace(0, 2*np.pi, 80)
v = np.linspace(0, np.pi, 40)
X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))
ax4.plot_surface(X, Y, Z, rstride=1, cstride=1, color='#e6f2ff', alpha=0.3)
#--- Ecuadores y meridianos ---
# Ecuador (z=0)
eq_u = np.linspace(0, 2*np.pi, 200)
ax4.plot(np.cos(eq_u), np.sin(eq_u), 0*eq_u, color='gray', linewidth=1)
# Meridiano en φ=0
mer_v = np.linspace(0, np.pi, 200)
ax4.plot(np.sin(mer_v), 0*mer_v, np.cos(mer_v), color='gray', linewidth=1)
# Meridiano en φ=π/2
ax4.plot(0*mer_v, np.sin(mer_v), np.cos(mer_v), color='gray', linewidth=1)

# --- Ejes X, Y, Z ---


((fig1,ax1),(fig2,ax2))=(plt.subplots(),plt.subplots())
ax1.set_ylim(-5,5)
ax2.set_ylim(-5,5)
x = np.arange(0, L, 0.001)

fig3=plt.figure(figsize=(6, 2))
axtheta=fig3.add_axes([0.1, 0.4, 0.8, 0.2])
axtheta.set_xticks([])
axtheta.set_yticks([])
for spine in axtheta.spines.values():
    spine.set_visible(False)

ax4.plot([-1,1],[0,0],[0,0], color='black', linewidth=1)
ax4.plot([0,0],[-1,1],[0,0], color='black', linewidth=1)
ax4.plot([0,0],[0,0],[-1,1], color='black', linewidth=1)
ax4.text(1.05,0,0,'X', fontsize=12)
ax4.text(0,1.05,0,'Y', fontsize=12)
ax4.text(0,0,1.05,'Z', fontsize=12)

#Animación


vector_line, = ax4.plot([0, x0], [0, y0], [0, z0], color='#1f77b4', lw=2)
punto = ax4.scatter([x0], [y0], [z0], color='#1f77b4', s=50)

angle = Slider(
    ax=axtheta,
    label='Angle [Rad]',
    valmin=0,
    valmax=2*np.pi,
    valinit=0,
    )
RFSt1, = ax1.plot(x, phi_x(x).real,'r',alpha=abs(FSt0prima[0])**2)
IFSt1, = ax1.plot(x,phi_x(x).imag,'r--',alpha=abs(FSt0prima[0])**2)
RFSt2, = ax1.plot(x, phi_x(x).real,'b',alpha=abs(FSt0prima[1])**2)
IFSt2, = ax1.plot(x,phi_x(x).imag,'b--',alpha=abs(FSt0prima[1])**2)
RBSt1, = ax2.plot(x, phi_x(x).real,'r',alpha=abs(BSt0prima[0])**2)
IBSt1, = ax2.plot(x,phi_x(x).imag,'r--',alpha=abs(BSt0prima[0])**2)
RBSt2, = ax2.plot(x, phi_x(x).real,'b',alpha=abs(BSt0prima[1])**2)
IBSt2, = ax2.plot(x,phi_x(x).imag,'b--',alpha=abs(BSt0prima[1])**2)
RBSt3, = ax2.plot(x, phi_x(x).real,'g',alpha=abs(BSt0prima[2])**2)
IBSt3, = ax2.plot(x,phi_x(x).imag,'g--',alpha=abs(BSt0prima[2])**2)
def animate1(i,anguloF,FSt0):
    valinit=0
    Rx=expm(-1j*anguloF/2*sigma_xF)
    Ry=expm(-1j*anguloF/2*sigma_yF)
    Rz=expm(-1j*anguloF/2*sigma_zF)
    if anguloF!=valinit:
        #MODIFICAR AQUI EJES DE ROTACION


        FSt0=Rz@FSt0


        ################################
        opacidadFSt1=abs(FSt0[0])**2
        opacidadFSt2=abs(FSt0[1])**2
        phi1=phi_x(x)*phi_t1(i)
        WFSt1=phi1*FSt0[0]
        WFSt2=phi1*FSt0[1]
        RFSt1.set_ydata(WFSt1.real)
        RFSt2.set_ydata(WFSt2.real)
        IFSt1.set_ydata(WFSt1.imag)
        IFSt2.set_ydata(WFSt2.imag)
        RFSt1.set_alpha(opacidadFSt1)
        RFSt2.set_alpha(opacidadFSt2)
        IFSt1.set_alpha(opacidadFSt1)
        IFSt2.set_alpha(opacidadFSt2)
        rx,ry,rz=esfera(FSt0[0],FSt0[1])
        
        # Actualiza Bloch
        vector_line.set_data_3d([0, rx], [0, ry], [0, rz])
        punto._offsets3d = (np.array([rx]), np.array([ry]), np.array([rz]))
        valinit=anguloF
        return RFSt1,IFSt1,RFSt2,IFSt2,vector_line,punto,
    else:
        FSt0=FSt0
        opacidadFSt1=abs(FSt0[0])**2
        opacidadFSt2=abs(FSt0[1])**2
        phi1=phi_x(x)*phi_t1(i)
        WFSt1=phi1*FSt0[0]
        WFSt2=phi1*FSt0[1]
        RFSt1.set_ydata(WFSt1.real)
        RFSt2.set_ydata(WFSt2.real)
        IFSt1.set_ydata(WFSt1.imag)
        IFSt2.set_ydata(WFSt2.imag)
        RFSt1.set_alpha(opacidadFSt1)
        RFSt2.set_alpha(opacidadFSt2)
        IFSt1.set_alpha(opacidadFSt1)
        IFSt2.set_alpha(opacidadFSt2)
        r=esfera(FSt0[0],FSt0[1])
        valinit=anguloF
        return RFSt1,IFSt1,RFSt2,IFSt2,vector_line,punto,
    
def animate2(i,anguloB,BSt0):
    valinit=0
    Rx=expm(-1j*anguloB/np.sqrt(2)*sigma_xB)
    Ry=expm(-1j*anguloB/np.sqrt(2)*sigma_yB)
    Rz=expm(-1j*anguloB*sigma_zB)
    if anguloB!=valinit:
        #MODIFICAR AQUI EJES DE ROTACION


        BSt0=Rz@BSt0

        
        ################################
        opacidadBSt1=abs(BSt0[0])**2
        opacidadBSt2=abs(BSt0[1])**2
        opacidadBSt3=abs(BSt0[2])**2
        phi1=phi_x(x)*phi_t2(i)
        WBSt1=phi1*BSt0[0]
        WBSt2=phi1*BSt0[1]
        WBSt3=phi1*BSt0[2]
        RBSt1.set_ydata(WBSt1.real)
        RBSt2.set_ydata(WBSt2.real)
        RBSt3.set_ydata(WBSt3.real)
        IBSt1.set_ydata(WBSt1.imag)
        IBSt2.set_ydata(WBSt2.imag)
        IBSt3.set_ydata(WBSt3.imag)
        RBSt1.set_alpha(opacidadBSt1)
        RBSt2.set_alpha(opacidadBSt2)
        RBSt3.set_alpha(opacidadBSt3)
        IBSt1.set_alpha(opacidadBSt1)
        IBSt2.set_alpha(opacidadBSt2)
        IBSt3.set_alpha(opacidadBSt3)
        valinit=anguloB
        return RBSt1,IBSt1,RBSt2,IBSt2,RBSt3,IBSt3
    else:
        BSt0=BSt0
        opacidadBSt1=abs(BSt0[0])**2
        opacidadBSt2=abs(BSt0[1])**2
        opacidadBSt3=abs(BSt0[2])**2
        phi1=phi_x(x)*phi_t2(i)
        WBSt1=phi1*BSt0[0]
        WBSt2=phi1*BSt0[1]
        WBSt3=phi1*BSt0[2]
        RBSt1.set_ydata(WBSt1.real)
        RBSt2.set_ydata(WBSt2.real)
        RBSt3.set_ydata(WBSt3.real)
        IBSt1.set_ydata(WBSt1.imag)
        IBSt2.set_ydata(WBSt2.imag)
        IBSt3.set_ydata(WBSt3.imag)
        RBSt1.set_alpha(opacidadBSt1)
        RBSt2.set_alpha(opacidadBSt2)
        RBSt3.set_alpha(opacidadBSt3)
        IBSt1.set_alpha(opacidadBSt1)
        IBSt2.set_alpha(opacidadBSt2)
        IBSt3.set_alpha(opacidadBSt3)
        valinit=anguloB
        return RBSt1,IBSt1,RBSt2,IBSt2,RBSt3,IBSt3

ani1 = animation.FuncAnimation(
    fig1, lambda i: animate1(i,angle.val,FSt0prima), interval=20, blit=False)
ani2 = animation.FuncAnimation(
    fig2, lambda i: animate2(i,angle.val,BSt0prima), interval=20, blit=False)
plt.show()


while plt.fignum_exists(fig1.number) or plt.fignum_exists(fig2.number):
    plt.pause(0.01)
