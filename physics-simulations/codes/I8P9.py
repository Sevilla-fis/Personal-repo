from sympy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from sympy.plotting import plot
from matplotlib.widgets import TextBox, Slider
import tabulate as tabulate
mpl.use('TkAgg')
from typing import Tuple
def ajuste_linear(x:NDArray[np.float64],y:NDArray[np.float64],errorx:NDArray[np.float64] | float=0,errory:NDArray[np.float64] | float=0,xlabel:str='',ylabel:str='',titulo:str='',color:str='orange',colorerr:str='black',colordot:str='red'):
    """
    Realiza un ajuste lineal de los datos x e y introducidos en la función, y devuelve la pendiente, la ordenada en origen y el coeficiente de correlación al cuadrado.
    """
    n=len(x)
    xm=np.mean(x)
    ym=np.mean(y)
    sx=np.sum(x)
    sy=np.sum(y)
    sxy=np.sum(x*y)
    sx2=np.sum(x**2)
    sy2=np.sum(y**2)
    m=(n*sxy-sx*sy)/(n*sx2-sx**2)
    b=(sx2*sy-sx*sxy)/(n*sx2-sx**2)
    xi=symbols('xi')
    f=b+m*xi
    fx=lambdify(xi,f)
    fi=fx(x)
    numerador=n*sxy-sx*sy
    raiz1=np.sqrt(n*sx2-sx**2)
    raiz2=np.sqrt(n*sy2-sy**2)
    r=numerador/(raiz1*raiz2)
    r2=r**2
    r2_porcentaje=np.around(r2*100,2)
    sigma=(sum((y-m*x-b)**2)/(len(x)-2))**(1/2)
    errorm=(((len(x)*sigma**2)/(len(x)*sum(x**2)-(sum(x))**2)))**(1/2)
    errorn=((sigma**2*sum(x**2))/(len(x)*sum(x**2)-(sum(x))**2))**(1/2)
    m1='m='+str(m)+'+/-'+str(errorm)
    b1='n='+str(b)+'+/-'+str(errorn)
    plt.plot(x,y,'o',label='Datos')
    plt.plot(x,fi,color=color,label='Ajuste')
    plt.errorbar(x,y,color=colorerr,yerr=errory,fmt='.')
    plt.errorbar(x,y,color=colorerr,xerr=errorx,fmt='.')
    for i in range(0,n,1):
        y0=np.min([y[i],fi[i]])
        y1=np.max([y[i],fi[i]])
        plt.vlines(x[i],y0,y1, color=colordot,
            linestyle ='dotted')
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titulo)
    plt.show
    f1='f='+str(f)
    r21='r^2='+str(r2)
    m1val=m
    return f1,m1,b1,r21,m1val

def properrores(variables:list[Symbol],f:Symbol,errores:list[float],valores1:NDArray[np.float64]):
    """
    Es necesario definir de manera simbólica antes de realizar la propagación de errores.
    En los valores se ha de utilizar un array bidimensional (una matriz).
    Realiza la propagación de errores de una función dadas unas variables y unos valores con sus respectivos errores, para devolver la funcion de error y los valores de estos.
    """
    ecuacion=0
    errorindirecto1=[]
    valores=valores1.T
    for i in range(len(variables)):
        ecuacion+=abs(diff(f,variables[i]))*errores[i]
    print(ecuacion)
    eqlatex=latex(ecuacion)
    print(eqlatex)
    for j in range(len(valores)):
        errorindirecto2=ecuacion.subs(dict(zip(variables,valores[j])))
        errorindirecto1.append(errorindirecto2)
    return np.array(errorindirecto1)
def csvtopython(csvfile,dtype=None,delimiter=';'):

    data=np.genfromtxt(csvfile,delimiter=delimiter,unpack=True,dtype=dtype) #We obtain data with np.genfromtxt
    #It's important that data is separated by . in excel (0.23), so type is float
    #Delimiter is set to ; (Change if necessary)

    datos=[]

    for i in range(len(data)): #Remove first element of array (name of variable) so we obtain only data (type:float)

        dat=data[i][1:]

        datos.append(dat)

    return datos
#Error Omega
m=1.1
g=9.8
L,k=symbols('L,k')
errk1=properrores([L,k],2*np.pi/(sqrt(g/L+2*k/m)),[0.001,0.13],np.array([[0.982],[12.84]]))
errk2=properrores([L,k],2*np.pi/(sqrt(g/L+2*k/m)),[0.001,1.11],np.array([[0.982],[3.07]]))
errk3=properrores([L,k],2*np.pi/(sqrt(g/L+2*k/m)),[0.001,0.85],np.array([[0.982],[20.28]]))
errk4=properrores([L,k],2*np.pi/(sqrt(g/L+2*k/m)),[0.001,0.09],np.array([[0.982],[2.83]]))
print(errk1)
print(errk2)
print(errk3)
print(errk4)
#Comrpobacion k
m_compro=4*np.pi**2*np.array([0.250,0.350,0.550,0.750])
T1=20**2/np.array([22**2,19**2,15**2,13**2])
T2=10**2/np.array([15**2,13**2,10**2,9**2])
m_compro3=4*np.pi**2*np.array([0.150,0.250,0.350,0.450])
T3=10**2/np.array([15**2,11**2,9**2,8**2])
m_compro4=4*np.pi**2*np.array([0.03,0.03,0.04,0.06])
T4=20**2/np.array([29**2,25**2,23**2,21**2])
print('k1=',ajuste_linear(T1,m_compro))
print('k2=',ajuste_linear(T2,m_compro))
print('k3=',ajuste_linear(T3,m_compro3))
print('k4=',ajuste_linear(T4,m_compro4))
plt.show()
 #Muelle 1              #x [cm], m [g]
Tabla1=np.array([[14.3,10],
                 [15.0,20],
                 [15.8,30],
                 [16.6,40],
                 [17.3,50],
                 [18.1,60]])
d1=12.3*10**-2
print(Tabla1.T[0])
a1=ajuste_linear(10**-2*Tabla1.T[0]-d1,9.8*10**-3*Tabla1.T[1],xlabel='Elongación [m]',ylabel='Masa [kg]',titulo='Cte elástica muelle 1')
m1=a1[4]
print(a1)
print('m1=',m1)
print(tabulate.tabulate(Tabla1,["Longitud del muelle [cm]","Masa newp[g]",],tablefmt="latex"))
g=9.8
m=symbols('m')
properrores([m],g/m,[0.013289773277664185],np.array([[m1]]))
plt.show()
#Muelle 2
d2=7.3*10**-2
Tabla2=np.array([[8.7,10],
                 [9.4,20],
                 [9.9,30],
                 [10.9,40],
                 [11.9,50],
                 [21.7,60]])
a2=ajuste_linear(10**-2*Tabla2.T[0]-d2,9.8*10**-3*Tabla2.T[1],xlabel='Elongación [m]',ylabel='Masa [kg]',titulo='Cte elástica muelle 2')
m2=a2[4]
print(a2)
print('m2=',m2)
print(tabulate.tabulate(Tabla2,["Longitud del muelle [cm]","Masa [g]",],tablefmt="latex"))
m=symbols('m')
properrores([m],g/m,[0.11296578621130098],np.array([[m2]]))
plt.show()
#Muelle 3
d3=18.5*10**-2
Tabla3=np.array([[19.2,10],
                 [19.7,20],
                 [20.2,30],
                 [20.6,40],
                 [21.0,50],
                 [21.7,60]])
a3=ajuste_linear(10**-2*Tabla3.T[0]-d3,9.8*10**-3*Tabla3.T[1],xlabel='Elongación [m]',ylabel='Masa [kg]',titulo='Cte elástica muelle 3')
m3=a3[4]
print(a3)
print('m3=',m3)
print(tabulate.tabulate(Tabla3,["Longitud del muelle [cm]","Masa [g]",],tablefmt="latex"))
m=symbols('m')
properrores([m],g/m,[0.08620689655172416],np.array([[m3]]))
plt.show()
#Muelle 4
d4=11.7*10**-2
Tabla4=np.array([[18.5,10],
                 [22.3,20],
                 [25.6,30],
                 [28.9,40],
                 [33.2,50],
                 [35.4,60]])
a4=ajuste_linear(10**-2*Tabla4.T[0]-d4,9.8*10**-3*Tabla4.T[1],xlabel='Elongación [m]',ylabel='Masa [kg]',titulo='Cte elástica muelle 4')
m4=a4[4]
print(a4)
print('m4=',m4)
print(tabulate.tabulate(Tabla4,["Longitud del muelle [cm]","Masa [g]",],tablefmt="latex"))
m=symbols('m')
properrores([m],g/m,[0.009398999979574149],np.array([[m4]]))
plt.show()
#Periodo Pendulos (40s de medidas, 27 promedio móvil, tasa de muestro 200Hz)
def calculo_periodos(Pos,Tiempo):
    t1=[]
    t0=[]
    cont0 = 0
    for i in range(len(Pos)):
        if Pos[i]==0 and cont0==0 and Pos[i+1]!=0:
            t0.append(Tiempo[i])
            cont0+=1
            continue
        if Pos[i]==0 and cont0==1 and Pos[i+1]!=0:
            cont0+=1
            continue
        if Pos[i]==0 and cont0==2 and Pos[i+1]!=0:
            t1.append(Tiempo[i])
            t0.append(Tiempo[i])
            cont0=1
            continue
    t0=t0[:-1]
    T=np.array(t1)-np.array(t0)
    T_medio=sum(T)/len(T)
    return t1,t0,T_medio
#Muelle 1
print('Muelle 1')
datos21=csvtopython('DatosParte21.csv')
t21=[]
for i in range(len(datos21[0])):
    datos21[0][i]=datos21[0][i].replace(',','.')
    t21.append(float(datos21[0][i]))
t21=np.array(t21)
print(t21)
U21=[]
for i in range(len(datos21[1])):
    datos21[1][i]=datos21[1][i].replace(',','.')
    U21.append(float(datos21[1][i]))
U21=np.array(U21)
print(U21)
t11,t01,T1=calculo_periodos(U21,t21)
print('t11=',t11)
print('t01=',t01)
print('T1_medio=',T1)
#Muelle 2
print('Muelle 2')
datos22=csvtopython('DatosParte22.csv')
t22=[]
for i in range(len(datos22[0])):
    datos22[0][i]=datos22[0][i].replace(',','.')
    t22.append(float(datos22[0][i]))
t22=np.array(t22)
print(t22)
U22=[]
for i in range(len(datos22[1])):
    datos22[1][i]=datos22[1][i].replace(',','.')
    U22.append(float(datos22[1][i]))
U22=np.array(U22)
print(U22)
t12,t02,T2=calculo_periodos(U22,t22)
print('t12=',t12)
print('t02=',t02)
print('T2_medio=',T2)
#Oscilacion en fase
#Muelle 1
def calculo_periodos_final3(Data):
    datos = csvtopython(Data)
    t = []
    for i in range(len(datos[0])):
        datos[0][i] = datos[0][i].replace(',', '.')
        t.append(float(datos[0][i]))
    t = np.array(t)
    print('t=',t)
    U1 = []
    for i in range(len(datos[1])):
        datos[1][i] = datos[1][i].replace(',', '.')
        U1.append(float(datos[1][i]))
    U1 = np.array(U1)
    print('U1=',U1)
    U2 = []
    for i in range(len(datos[2])):
        datos[2][i] = datos[2][i].replace(',', '.')
        U2.append(float(datos[2][i]))
    U2 = np.array(U2)
    print('U2=',U2)



    t1=[]
    t0=[]
    cont0 = 0
    for i in range(len(U1)):
        if U1[i]==0 and cont0==0 and U1[i+1]!=0:
            t0.append(t[i])
            cont0+=1
            continue
        if U1[i]==0 and cont0==1 and U1[i+1]!=0:
            cont0+=1
            continue
        if U1[i]==0 and cont0==2 and U1[i+1]!=0:
            t1.append(t[i])
            t0.append(t[i])
            cont0=1
            continue
    t0=t0[:-1]
    T=np.array(t1)-np.array(t0)
    T_medio=sum(T)/len(T)



    t12=[]
    t02=[]
    cont02 = 0
    for i in range(len(U2)):
        if U2[i]==0 and cont02==0 and U2[i+1]!=0:
            t02.append(t[i])
            cont02+=1
            continue
        if U2[i]==0 and cont02==1 and U2[i+1]!=0:
            cont02+=1
            continue
        if U2[i]==0 and cont02==2 and U2[i+1]!=0:
            t12.append(t[i])
            t02.append(t[i])
            cont02=1
            continue
    t02=t02[:-1]
    T2=np.array(t12)-np.array(t02)
    T_medio2=sum(T2)/len(T2)
    return T_medio, T_medio2

def calculo_periodos_final003(Data):
    datos = csvtopython(Data)
    t = []
    for i in range(len(datos[0])):
        datos[0][i] = datos[0][i].replace(',', '.')
        t.append(float(datos[0][i]))
    t = np.array(t)
    print('t=',t)
    U1 = []
    for i in range(len(datos[1])):
        datos[1][i] = datos[1][i].replace(',', '.')
        U1.append(float(datos[1][i]))
    U1 = np.array(U1)+0.05
    print('U1=',U1)
    U2 = []
    for i in range(len(datos[2])):
        datos[2][i] = datos[2][i].replace(',', '.')
        U2.append(float(datos[2][i]))
    U2 = np.array(U2)-0.05
    print('U2=',U2)



    t1=[]
    t0=[]
    cont0 = 0
    for i in range(len(U1)):
        if U1[i]==0 and cont0==0 and U1[i+1]!=0:
            t0.append(t[i])
            cont0+=1
            continue
        if U1[i]==0 and cont0==1 and U1[i+1]!=0:
            cont0+=1
            continue
        if U1[i]==0 and cont0==2 and U1[i+1]!=0:
            t1.append(t[i])
            t0.append(t[i])
            cont0=1
            continue
    t0=t0[:-1]
    T=np.array(t1)-np.array(t0)
    T_medio=sum(T)/len(T)



    t12=[]
    t02=[]
    cont02 = 0
    for i in range(len(U2)):
        if U2[i]==0 and cont02==0 and U2[i+1]!=0:
            t02.append(t[i])
            cont02+=1
            continue
        if U2[i]==0 and cont02==1 and U2[i+1]!=0:
            cont02+=1
            continue
        if U2[i]==0 and cont02==2 and U2[i+1]!=0:
            t12.append(t[i])
            t02.append(t[i])
            cont02=1
            continue
    t02=t02[:-1]
    T2=np.array(t12)-np.array(t02)
    T_medio2=sum(T2)/len(T2)
    return T_medio, T_medio2

def calculo_periodos_final008(Data):
    datos = csvtopython(Data)
    t = []
    for i in range(len(datos[0])):
        datos[0][i] = datos[0][i].replace(',', '.')
        t.append(float(datos[0][i]))
    t = np.array(t)
    print('t=',t)
    U1 = []
    for i in range(len(datos[1])):
        datos[1][i] = datos[1][i].replace(',', '.')
        U1.append(float(datos[1][i]))
    U1 = np.array(U1)+0.05
    print('U1=',U1)
    U2 = []
    for i in range(len(datos[2])):
        datos[2][i] = datos[2][i].replace(',', '.')
        U2.append(float(datos[2][i]))
    U2 = np.array(U2)-0.05
    print('U2=',U2)



    t1=[]
    t0=[]
    cont0 = 0
    for i in range(len(U1)):
        if U1[i]==0 and cont0==0 and U1[i+1]!=0:
            t0.append(t[i])
            cont0+=1
            continue
        if U1[i]==0 and cont0==1 and U1[i+1]!=0:
            cont0+=1
            continue
        if U1[i]==0 and cont0==2 and U1[i+1]!=0:
            t1.append(t[i])
            t0.append(t[i])
            cont0=1
            continue
    t0=t0[:-1]
    T=np.array(t1)-np.array(t0)
    T_medio=sum(T)/len(T)



    t12=[]
    t02=[]
    cont02 = 0
    for i in range(len(U2)):
        if U2[i]==0 and cont02==0 and U2[i+1]!=0:
            t02.append(t[i])
            cont02+=1
            continue
        if U2[i]==0 and cont02==1 and U2[i+1]!=0:
            cont02+=1
            continue
        if U2[i]==0 and cont02==2 and U2[i+1]!=0:
            t12.append(t[i])
            t02.append(t[i])
            cont02=1
            continue
    t02=t02[:-1]
    T2=np.array(t12)-np.array(t02)
    T_medio2=sum(T2)/len(T2)
    return T_medio, T_medio2

print('31',calculo_periodos_final3('DatosParte31.csv'))
print('32',calculo_periodos_final3('DatosParte32.csv'))
print('33',calculo_periodos_final3('DatosParte33.csv'))
print('34',calculo_periodos_final3('DatosParte34.csv'))
print('41',calculo_periodos_final003('DatosParte41.csv'))
print('42',calculo_periodos_final003('DatosParte42.csv'))
print('43',calculo_periodos_final008('DatosParte43.csv'))
print('44',calculo_periodos_final3('DatosParte44.csv'))
print('54',calculo_periodos_final3('DatosParte54.csv'))
#Antifase (tiempo de medida 15s)









#datos31=csvtopython('DatosParte31.csv')
#t31=[]
#for i in range(len(datos31[0])):
    #datos31[0][i]=datos31[0][i].replace(',','.')
    #t31.append(float(datos31[0][i]))
#t31=np.array(t31)
#print(t31)
#U311=[]
#for i in range(len(datos31[1])):
    #datos31[1][i]=datos31[1][i].replace(',','.')
    #U311.append(float(datos31[1][i]))
#U311=np.array(U311)
#print(U311)
#U312=[]
#for i in range(len(datos31[2])):
    #datos31[2][i]=datos31[2][i].replace(',','.')
    #U312.append(float(datos31[2][i]))
#U312=np.array(U312)
#print(U312)
#t131,t031,T31=calculo_periodos(U311,t31)
#('t12=',t131)
#print('t02=',t031)
#print('T2_medio1=',T31)
#t132,t032,T32=calculo_periodos(U312,t31)
#print('t12=',t12)
#print('t02=',t02)
#print('T2_medio2=',T32)
