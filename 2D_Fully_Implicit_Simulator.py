import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import legend
import matplotlib.dates as mdates
import seaborn as sns
import matplotlib.cm as cmap
plt.rcParams['text.usetex']=True

kx=200 #mD
ky=200 #mD
Bo=np.array([[400,1.0120],[800,1.0255],[1200,1.0380],[1600,1.0510],[2000,1.0630],[2400,1.0750],[2800,1.0870],
             [3200,1.0985],[3600,1.1100],[4000,1.1056],[4400,1.1012],[4800,1.0968],[5200,1.0925],[5600,1.0882]])
Bof=Bo[:,1]
Bop=Bo[:,0]
mu=np.array([[400,1.17],[800,1.14],[1200,1.11],[1600,1.08],[2000,1.06],[2400,1.03],[2800,1.00],
             [3200,0.98],[3600,0.95],[4000,0.97],[4400,0.99],[4800,1.01],[5200,1.03],[5600,1.05]])
muf=mu[:,1]
mup=mu[:,0]
deno=np.array([[400,46.497],[800,48.1000],[1200,49.372],[1600,50.726],[2000,52.072],[2400,53.318],[2800,54.399],
             [3200,55.424],[3600,56.203],[4000,56.421],[4400,56.646],[4800,56.871],[5200,57.096],[5600,57.320]])
denof=deno[:,1]
denop=deno[:,0]
cr= 4e-6 #psi-1
co= 1.26e-5 #psi-1
ct=cr+co #psi-1
Lx=4000 #ft
Ly=4000 #ft
h=50
Ix=21
Iy=21
dt=5 #dias
dtt=0
dx=Lx/Ix
dy=Ly/Iy
q=np.zeros((Ix,Iy))
q[10,10]=-2500
pi=5500
pi=np.ones((Ix,Iy))*pi
tsim=180 #días
Tx=np.zeros((Ix,Iy))
Ty=np.zeros((Ix,Iy))
TermA=np.zeros((Ix,Iy))
por=0.2
por=np.ones((Ix,Iy,tsim))*por
dp=1e-6
error=1
n=0
tol=1e-3
b=np.zeros((Ix,Iy))
bx=np.zeros((Ix,Iy))
by=np.zeros((Ix,Iy))
b1=np.zeros(Ix*Iy)
a=np.zeros((Ix,Iy))
a1=np.zeros(Ix*Iy)
c=np.zeros((Ix,Iy))
c1=np.zeros(Ix*Iy)
e=np.zeros((Ix,Iy-1))
e1=np.zeros(Ix*Iy)
f=np.zeros((Ix,Iy-1))
f1=np.zeros(Ix*Iy)
d=np.zeros((Ix,Iy))
Fx=np.zeros((Ix,Iy))
Fy=np.zeros((Ix,Iy))
d1=np.zeros(Ix*Iy)
Derpor=np.zeros((Ix,Iy))
DerTx=np.zeros((Ix,Iy))
Derdenx=np.zeros((Ix,Iy))
DerTy=np.zeros((Ix,Iy))
Derdeny=np.zeros((Ix,Iy))
f1=np.zeros((Ix,Iy))
f2=np.zeros((Ix,Iy))
f3=np.zeros((Ix,Iy))
deltap=np.zeros((Ix,Iy,tsim))
deltapp=np.zeros((Ix,Iy,tsim))
delta=([])
p=np.zeros((Ix,Iy,tsim))
psl=np.zeros((Ix,Iy,tsim))
dd=np.zeros((Ix,Iy,tsim))
BAL=np.zeros((Ix,Iy,tsim))
nodo=np.zeros((Ix,Iy))
znodo=np.zeros((Ix,Iy))
time=np.zeros(tsim)
zref=9035
ang=0
NIT=0

for j in range(Iy):
    for i in range(Ix):
        nodo[i,j]=i
        znodo[i,j]=zref-(nodo[i,j]*dx*np.tan(np.radians(ang)))
        p[i,j,0]=pi[i,j]+((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-zref))/144)
        psl[i,j,0]=pi[i,j]+((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-zref))/144)
        pi[i,j]=p[i,j,0]

while dtt<tsim:
    n=n+1
    while error>tol:
        NIT=NIT+1
        p[:,:,n]=pi[:,:]
        for j in range(Iy):
            for i in range(Ix):
                por[i,j,n]=por[i,j,n-1]*(1+(cr*(p[i,j,n]-p[i,j,n-1])))
                
        for j in range(Iy):
            for i in range(Ix):
                Tx[i,j]=0.001127*((h*dy)/dx)*(kx/(np.interp(p[i,j,n], mup, muf)*np.interp(p[i,j,n], Bop, Bof)))
                DerTx[i,j]=0.001127*((h*dy)/dx)*kx*((1/(np.interp(p[i,j,n]+dp, mup, muf)*np.interp(p[i,j,n]+dp, Bop, Bof)))-
                                                    (1/(np.interp(p[i,j,n]-dp, mup, muf)*np.interp(p[i,j,n]-dp, Bop, Bof))))/(2*dp)
                Derdenx[i,j]=(((np.interp(p[i,j,n]+dp, denop, denof)))-((np.interp(p[i,j,n]-dp, denop, denof))))/(2*dp)
                
        for i in range(Ix):
            for j in range(Iy):
                Ty[i,j]=0.001127*((h*dx)/dy)*(ky/(np.interp(p[i,j,n], mup, muf)*np.interp(p[i,j,n], Bop, Bof)))
                DerTy[i,j]=0.001127*((h*dx)/dy)*ky*((1/(np.interp(p[i,j,n]+dp, mup, muf)*np.interp(p[i,j,n]+dp, Bop, Bof)))-
                                                    (1/(np.interp(p[i,j,n]-dp, mup, muf)*np.interp(p[i,j,n]-dp, Bop, Bof))))/(2*dp)
                Derdeny[i,j]=(((np.interp(p[i,j,n]+dp, denop, denof)))-((np.interp(p[i,j,n]-dp, denop, denof))))/(2*dp)
        for j in range(Iy):
            for i in range(1,Ix):
                a[i,j]=Ty[i-1,j]+(((Ty[i-1,j]*(znodo[i,j]-znodo[i-1,j]))*Derdeny[i-1,j])/144)-((p[i,j,n]-p[i-1,j,n]-
                            ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144))*(DerTx[i-1,j]))
                
        a1=a        
        a1=np.transpose(a1)
        a1=np.reshape(a1, -1)

        for j in range(Iy):
             for i in range(0,Ix-1):
                c[i,j]=Ty[i+1,j]-(((Ty[i+1,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144)+((p[i+1,j,n]-p[i,j,n]-
                            ((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144))*(DerTx[i+1,j]))
                
        c1=c
        c1=np.transpose(c1)
        c1=np.reshape(c1, -1)
        for j in range(1,Iy):
            for i in range(Ix):
                e[i,j-1]=Tx[i,j-1]+(((Tx[i,j-1]*(znodo[i,j]-znodo[i,j-1]))*Derdenx[i,j-1])/144)-((p[i,j,n]-p[i,j-1,n]-
                            ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144))*(DerTy[i,j-1]))
                
        e1=e
        e1=np.transpose(e1)
        e1=np.reshape(e1, -1)
        for j in range(0,Iy-1):
            for i in range(Ix):
                f[i,j]=Tx[i,j+1]-(((Tx[i,j+1]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144)+((p[i,j+1,n]-p[i,j,n]-
                            ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144))*(DerTy[i,j+1]))
                
        f1=f
        f1=np.transpose(f1)
        f1=np.reshape(f1, -1)
        for i in range(Ix):
            for j in range(Iy):
            
                if (i==0)&(j==0):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])-((((Ty[i,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(DerTy[i,j]*((p[i+1,j,n]-p[i,j,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                    bx[i,j]=-(Tx[i,j])-((((Tx[i,j]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(DerTx[i,j]*((p[i,j+1,n]-p[i,j,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (i==(Ix-1))&(j==(Iy-1)):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])+((((Ty[i,j]*(znodo[i,j]-znodo[i-1,j]))*Derdenx[i-1,j])/144))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                    bx[i,j]=-(Tx[i,j])+((((Tx[i,j]*(znodo[i,j]-znodo[i,j-1]))*Derdeny[i,j-1])/144))-(DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (i==0)&(j==(Iy-1)):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])-((((Ty[i,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(DerTy[i,j]*((p[i+1,j,n]-p[i,j,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                    bx[i,j]=-(Tx[i,j])+((((Tx[i,j]*(znodo[i,j]-znodo[i,j-1]))*Derdeny[i,j-1])/144))-(DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (i==(Ix-1))&(j==0):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])+((((Ty[i,j]*(znodo[i,j]-znodo[i-1,j]))*Derdeny[i-1,j])/144))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                    bx[i,j]=-(Tx[i,j])-((((Tx[i,j]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(DerTx[i,j]*((p[i,j+1,n]-p[i,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (i==0)&(0<j<(Iy-1)):    
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])-((((Ty[i,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(DerTy[i,j]*((p[i+1,j,n]-p[i,j,n])-
                                         ((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                    bx[i,j]=-(Tx[i,j+1]+Tx[i,j])-((((Tx[i,j+1]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(((((Tx[i,j]*
                                        (znodo[i,j]-znodo[i,j-1]))*Derdenx[i,j]))/144))+((((DerTx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-
                                        ((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144))))-(DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (0<i<(Ix-1))&(j==0):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    bx[i,j]=-(Tx[i,j])-((((Tx[i,j]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(((DerTx[i,j]*((p[i,j+1,n]-p[i,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))))
                    by[i,j]=-(Ty[i+1,j]+Ty[i,j])-((((Ty[i+1,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(((((Ty[i,j]*
                                        (znodo[i,j]-znodo[i-1,j]))*Derdeny[i,j]))/144))+((((DerTy[i+1,j]*((p[i+1,j,n]-p[i,j,n])-
                                        ((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144))))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (i==(Ix-1))&(0<j<(Iy-1)):
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i,j])+(((((Ty[i,j]*(znodo[i,j]-znodo[i-1,j]))*Derdeny[i,j]))/144))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                    bx[i,j]=-(Tx[i,j+1]+Tx[i,j])-((((Tx[i,j+1]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(((((Tx[i,j]*
                                        (znodo[i,j]-znodo[i,j-1]))*Derdenx[i,j]))/144))+(((DerTx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-
                                        ((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144))))-(DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144))))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
                elif (0<i<(Ix-1))&(j==(Iy-1)): 
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i+1,j]+Ty[i,j])-((((Ty[i+1,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(((((Ty[i,j]*
                                        (znodo[i,j]-znodo[i-1,j]))*Derdeny[i,j]))/144))+(((DerTy[i+1,j]*((p[i+1,j,n]-p[i,j,n])-
                                        ((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144))))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                        ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144))))
                    bx[i,j]=-(Tx[i,j])+(((((Tx[i,j]*(znodo[i,j]-znodo[i,j-1]))*Derdenx[i,j]))/144))-((DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                          ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144))))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))

                else:
                    Derpor[i,j]=((((por[i,j,n]/np.interp(p[i,j,n]+dp, Bop, Bof))-(por[i,j,n]/np.interp(p[i,j,n-1]-dp, Bop, Bof)))/(2*dp)))
                    by[i,j]=-(Ty[i+1,j]+Ty[i-1,j])-((((Ty[i+1,j]*(znodo[i+1,j]-znodo[i,j]))*Derdeny[i+1,j])/144))+(((((Ty[i,j]*(znodo[i,j]-znodo[i-1,j]))
                                    *Derdeny[i,j]))/144))+((((DerTy[i+1,j]*((p[i+1,j,n]-p[i,j,n])-
                                    ((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144))))-(DerTy[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                    ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))))
                    bx[i,j]=-(Tx[i,j+1]+Tx[i,j-1])-((((Tx[i,j+1]*(znodo[i,j+1]-znodo[i,j]))*Derdenx[i,j+1])/144))+(((((Tx[i,j]*(znodo[i,j]-znodo[i,j-1]))
                                    *Derdenx[i,j]))/144))+((((DerTx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-
                                    ((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144))))-(DerTx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                    ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))))
                    b[i,j]=bx[i,j]+by[i,j]-(((h*dy*dx)/(5.615*dt))*(Derpor[i,j]))
        b1=b
        b1=np.transpose(b1)
        b1=np.reshape(b1,-1)  

        for i in range(Ix):
            for j in range(Iy):
                
                    if (i==0)&(j==0):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=(Ty[i,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                        Fx[i,j]=(Tx[i,j]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (i==(Ix-1))&(j==(Iy-1)):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=-(Ty[i,j]*((p[i,j,n]-p[i-1,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        Fx[i,j]=-(Tx[i,j]*((p[i,j,n]-p[i,j-1,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (i==0)&(j==(Iy-1)):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=(Ty[i,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                        Fx[i,j]=-(Tx[i,j]*((p[i,j,n]-p[i,j-1,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (i==(Ix-1))&(j==0):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=-(Ty[i,j]*((p[i,j,n]-p[i-1,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        Fx[i,j]=(Tx[i,j]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (i==0)&(0<j<(Iy-1)):  
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=(Ty[i,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))
                        Fx[i,j]=(Tx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))-(Tx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                                                  ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (0<i<(Ix-1))&(j==0):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))  
                        Fy[i,j]=(Ty[i+1,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))-(Ty[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                                                  ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        Fx[i,j]=(Tx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (i==(Ix-1))&(0<j<(Iy-1)):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=-(Ty[i,j]*((p[i,j,n]-p[i-1,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        Fx[i,j]=(Tx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))-(Tx[i,j]*((p[i,j,n]-p[i,j-1,n])-
                                                                  ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    elif (0<i<(Ix-1))&(j==(Iy-1)):
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fy[i,j]=(Ty[i+1,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))-(Ty[i,j]*((p[i,j,n]-p[i-1,j,n])-
                                                                  ((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        Fx[i,j]=-(Tx[i,j]*((p[i,j,n]-p[i,j-1,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                    else:
                        f3[i,j]=((h*dy*dx)/(5.615*dt))*por[i,j,n]*((((1/np.interp(p[i,j,n], Bop, Bof))-(1/np.interp(p[i,j,n-1], Bop, Bof)))))
                        Fx[i,j]=(Tx[i,j+1]*((p[i,j+1,n]-p[i,j,n])-((np.interp(p[i,j+1,n], denop, denof)*(znodo[i,j+1]-znodo[i,j]))/144)))-(Tx[i,j]*
                                           ((p[i,j,n]-p[i,j-1,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i,j-1]))/144)))
                        Fy[i,j]=(Ty[i+1,j]*((p[i+1,j,n]-p[i,j,n])-((np.interp(p[i+1,j,n], denop, denof)*(znodo[i+1,j]-znodo[i,j]))/144)))-(Ty[i,j]*
                                           ((p[i,j,n]-p[i-1,j,n])-((np.interp(p[i,j,n], denop, denof)*(znodo[i,j]-znodo[i-1,j]))/144)))
                        d[i,j]=-(Fx[i,j]+Fy[i,j]-f3[i,j]+q[i,j])
                

        d1=d
        d1=np.transpose(d1)
        d1=np.reshape(d1,-1)      
        A = np.diag(b1,0) + np.diag(a1[1:],-1) + np.diag(c1[:(Ix*Iy)-1],1) + np.diag(e1[:],-Ix) + np.diag(f1[:],Ix)
        
        delta=np.linalg.solve(A,d1)
        error=np.linalg.norm(delta)
        delta=np.transpose(delta)
        delta=np.reshape(delta,(Ix,Iy))
        pi=pi+delta
        

        if NIT==1:
            psl[:,:,n]=pi

        if NIT>=6:
            dt=dt/2
            break
    
    if 3<=NIT<=5:
        dtt=dtt+dt
        time[n]=dtt
        dt=dt
        if dt>30.5:
            dt=30.5
        else:
            dt=dt

    elif NIT<3:
        dtt=dtt+dt
        time[n]=dtt
        dt=2*dt
        if dt>30.5:
            dt=30.5
        else:
            dt=dt

    elif (NIT>=6)&(error<tol):
        dtt=dtt+dt
        time[n]=dtt
        dt=5
    elif (NIT>=6)&(error>tol):
        n=n-1
        pi=p[:,:,n]

    NIT=0
    error=1

pn=np.zeros((Ix,Iy,n+1))
pnsl=np.zeros((Ix,Iy,n+1))
t=np.zeros(n+1)

for k in range(n+1):
    for j in range(Iy):
        for i in range(Ix):
            pn[i,j,k]=p[i,j,k]
            pnsl[i,j,k]=psl[i,j,k]
    t[k]=time[k]


fig, ax = plt.subplots(figsize = (12, 9))
im = ax.imshow(pn[:,:,0], interpolation='gaussian', cmap='jet')
ax.invert_yaxis()
ax.set_xlabel('Nodo X')
ax.set_ylabel('Nodo Y')
plt.title(label=f'Distribución de presiones, t={t[0]} días', fontsize='20')
cbar = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cbar.set_label('Presión, psia')


fig, ax = plt.subplots(figsize = (12, 9))
im = ax.imshow(pn[:,:,7], interpolation='nearest', cmap='jet')
ax.invert_yaxis()
ax.set_xlabel('Nodo X')
ax.set_ylabel('Nodo Y')
plt.title(label=f'Distribución de presiones, t={np.round(t[n],0)} días', fontsize='20')
cbar = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cbar.set_label('Presión, psia')
plt.show()
#for j in range(0,Ix,10):
plt.plot(t,pn[10,10,:],linewidth=2.5,label=f'Totalmente Implicito, i={10},j={10}', color='blue')
plt.plot(t,pnsl[10,10,:],linewidth=2.5,label=f'SI Linealizado, i={10},j={10}', color='red', linestyle='--')
plt.xlabel('Tiempo [días]')
plt.ylabel('Presión [psia]')
plt.legend()
plt.grid()
plt.ylim(1000,6000)
plt.xlim(0,np.round(t[n]))
plt.show()

