import pickle
import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

print("Comecou!")

dps_value = 100

mp.dps = dps_value


K = 100
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(5)     
N = mpf(0.01)    
L = mpf(5000)

subdivisions = 100

pizao = mp.pi(dps=dps_value)

def H(y):
    return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * (y-L/mpf(2))/L)) )


def Bsfc(y):
    return mpf(0.1) * mp.sin(mpf(2) * pizao * y/L)


z = np.arange(0,2010,1) 
y = np.arange(-float(L),float(L)+5,5) 
Y,Z = np.meshgrid(y,z)
Y = Y * mpf(1)
Z = Z * mpf(1)



with open('/home/gome0004/kataprogramms/output/Ek.pickle', 'rb') as handle:
    Ek = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/Eki.pickle', 'rb') as handle:
    Eki = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/Ck.pickle', 'rb') as handle:
    Ck = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/Cki.pickle', 'rb') as handle:
    Cki = pickle.load(handle)   
    

with open('/home/gome0004/kataprogramms/output/Dk.pickle', 'rb') as handle:
    Dk = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/Dki.pickle', 'rb') as handle:
    Dki = pickle.load(handle)


with open('/home/gome0004/kataprogramms/output/B.pickle', 'rb') as handle:
    B = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/V.pickle', 'rb') as handle:
    V = pickle.load(handle)
    

# with open('/home/gome0004/kataprogramms/output/VlargeY.pickle', 'rb') as handle:
#     V = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/U.pickle', 'rb') as handle:
    U = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/Eq.pickle', 'rb') as handle:
    Eq = pickle.load(handle)
    
    
with open('/home/gome0004/kataprogramms/output/Eqi.pickle', 'rb') as handle:
    Eqi = pickle.load(handle)
 
    
with open('/home/gome0004/kataprogramms/output/W.pickle', 'rb') as handle:
    W = pickle.load(handle)


# with open('/home/gome0004/kataprogramms/output/WlargeY.pickle', 'rb') as handle:
#     W = pickle.load(handle)


with open('/home/gome0004/kataprogramms/output/P.pickle', 'rb') as handle:
    P = pickle.load(handle)
    

with open('/home/gome0004/kataprogramms/output/psi.pickle', 'rb') as handle:
    psi = pickle.load(handle)
    
print("Comecou!")

fig,ax1 = plt.subplots(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.close()

    
Yplot,Zplot = np.meshgrid(y,z)
Bplot = np.ones_like(B)*[mpf(0)]
for k in range(0,len(B)):
    for t in range(0,len(B[0])):
        Bplot[k][t] = float(B[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)
       

    
##Plotting the buoyancy
fig,ax1 = plt.subplots(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
maxB = np.nanmax(Bplot)
minB = np.nanmin(Bplot)
plt.title(r'B$_\max$ = $%(number).2f$ m $\rms^{-2}$        B$_\min$ = $%(number2).2f$ m $\rms^{-2}$' % {"number": maxB, "number2": minB},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Bplot,np.arange(-0.1,0.11,0.01),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
#plt.colorbar(label='[$ms^{-2}$]')
cbar = plt.colorbar(ticks=np.arange(-0.1,0.12,0.02))
cbar.ax.tick_params(length=14, width=1)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

nameoffigure = 'buoyancy.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()   






Yplot,Zplot = np.meshgrid(y,z)
Vplot = np.ones_like(V)*[mpf(0)]
for k in range(0,len(V)):
    for t in range(0,len(V[0])):
        Vplot[k][t] = float(V[k][t].real) 

##Plotting the V wind
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
maxV = np.nanmax(Vplot)
minV = np.nanmin(Vplot)
plt.title(r'V$_\max$ = $%(number).2f$ m $\rms^{-1}$        V$_\min$ = $%(number2).2f$ m $\rms^{-1}$' % {"number": maxV, "number2": minV},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Vplot,np.arange(-4,4.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
cbar = plt.colorbar()
cbar.ax.tick_params(length=14, width=1)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Vwind.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()  




Yplot,Zplot = np.meshgrid(y,z)
Uplot = np.ones_like(U)*[mpf(0)]
for k in range(0,len(U)):
    for t in range(0,len(U[0])):
        Uplot[k][t] = float(U[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)

print(U.shape)
#Plotting the U wind
fig,ax1 = plt.subplots(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
maxU = np.nanmax(Uplot)
minU = np.nanmin(Uplot)
plt.title(r'U$_\max$ = $%(number).2f$ m $\rms^{-1}$        U$_\min$ = $%(number2).2f$ m $\rms^{-1}$' % {"number": maxU, "number2": minU},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Uplot,np.arange(-2,2.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
cbar = plt.colorbar()
cbar.ax.tick_params(length=14, width=1)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])


ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Uwind.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()






Yplot,Zplot = np.meshgrid(y,z)
Wplot = np.ones_like(W)*[mpf(0)]
for k in range(0,len(W)):
    for t in range(0,len(W[0])):
        Wplot[k][t] = float(W[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)
           

##Plotting the W wind
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
maxW = np.nanmax(Wplot)
minW = np.nanmin(Wplot)
plt.title(r'W$_\max$ = $%(number).2f$ m $\rms^{-1}$        W$_\min$ = $%(number2).2f$ m $\rms^{-1}$' % {"number": maxW, "number2": minW},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Wplot,np.arange(-2,2.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
cbar = plt.colorbar()
cbar.ax.tick_params(length=14, width=1)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Wwind.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()




#Calculating the U* and plotting it  
Ustar = np.ones_like(y)*[mpf(0)]
ympf = y*[mpf(1)]
for k in range(-K,K+1):
    for t in range(0,len(y)):
        Ustar[t] = Ustar[t] +  1j*Eq[Eqi.index(k)] * mp.sin(mpf(2) * mpf(k) * pizao * ympf[t] / L)/mp.cos(alpha)


Ustarplot = np.ones_like(Ustar)*[mpf(0)]
for t in range(0,len(y)):
    Ustarplot[t] = float(Ustar[t].real)

        
##Plotting U* infinity
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


#plt.title('U star plot')
plt.plot(y,Ustarplot,linewidth=2.0)
plt.xlabel("Y [m]")
plt.ylabel('U$^{\u2605}_\infty$ [m s$^{-1}$]')
#plt.xlim([-5000,5000])
plt.xlim([float(-L),float(L)])
#plt.ylim([-6,2])
plt.ylim([-4,4])
plt.yticks(np.arange(-4,6,2))
plt.grid('True')

nameoffigure = 'Ustar.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()





#Plotting the theta angle
Ustar_real = U * mp.cos(alpha) + W * mp.sin(alpha)
Wstar_real = - U * mp.sin(alpha) + W * mp.cos(alpha)
theta = np.ones_like(U)*[mpf(0)]
for k in range(0,len(theta)):
    for t in range(0,len(theta[0])):
        theta[k][t] = mp.acos( abs(Ustar_real[k][t] / mp.sqrt(U[k][t]**2 + W[k][t]**2 )))
theta = theta * mpf(180) / pizao

Yplot,Zplot = np.meshgrid(y,z)
thetaplot = np.ones_like(theta)*[mpf(0)]
thetaplot_sign = np.ones_like(theta)*[mpf(0)]
for k in range(0,len(theta)):
    for t in range(0,len(theta[0])):
        thetaplot[k][t] = float(theta[k][t].real) #+ 1j*float(theta[k][t].imag)
        # thetaplot_sign[k][t] = float(theta[k][t].real * abs(Wstar_real[k][t].real)/Wstar_real[k][t].real)
        thetaplot_sign[k][t] = np.sign(float(Wstar_real[k][t].real)) * float(theta[k][t].real)
        
        
#Plotting the theta angle vs height at y=0 
thetaplotf = (thetaplot).T
thetaplotf2 = (thetaplot_sign).T 
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
ax1.plot([0,0],[0,1500],linewidth=2,color="k")
ax1.plot(thetaplotf2[50][1:],z[1:],linewidth=3,color="r") #Number might depends on how many points y has
plt.xlabel(r"$\rm\theta$ [$\rm^{o}$]")
plt.ylabel("Z [m]")
plt.xlim([-0.5,6])
ax1.set_ylim([0,1500])
ax1.set_xticks(np.arange(-0.5,6.5,0.5))
ax1.set_yticks(np.arange(0,1600,100))
ax1.tick_params('both', length=10, width=1, which='major')
plt.grid(True)

nameoffigure = 'theta.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()


#Plotting the Wstar_real 
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

interval = 12
boundsl = []
for k in range(-interval,interval+1):
    if k < 0:
        boundsl.append(0.3/interval*k)
    if k > 0:
        boundsl.append(3/interval*k)
    if k ==0:
        boundsl.append(k)
bounds = np.array(boundsl)


Yplot,Zplot = np.meshgrid(y,z)
Wstar_realplot = np.ones_like(Wstar_real)*[mpf(0)]
for k in range(0,len(Wstar_real)):
    for t in range(0,len(Wstar_real[0])):
        Wstar_realplot[k][t] = float(Wstar_real[k][t].real) 
 
# Yploy = Yplot
# for k in range(0,W.shape[0]):
#     for t in range(0,W.shape[1]):
#         if Wstarplot[k][t] == nan:
#             Yploy[k][t] = 0
#         else:
#             Yploy[k][t] = Wstarplot[k][t]


Wstarplot = Wstar_realplot
fig, ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
maxWstar = np.nanmax(Wstarplot)
minWstar = np.nanmin(Wstarplot)
plt.title(r'W$^{\bigstar}_\max$ = $%(number).2f$ m $\rms^{-1}$        W$^{\bigstar}_\min$ = $%(number2).2f$ m $\rms^{-1}$' % {"number": maxWstar, "number2": minWstar},x=0.5, y=1.02)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# pcm = plt.pcolormesh(Yplot,Zplot,Yploy,norm=norm,cmap='seismic')

plt.contourf(Yplot,Zplot,Wstarplot,np.arange(-2,2.2,0.2),cmap='seismic')

cbar = plt.colorbar()
cbar.ax.tick_params(length=14, width=1)

plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

nameoffigure = 'Wstar.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()


            

Yplot,Zplot = np.meshgrid(y,z)
Pplot = np.ones_like(P)*[mpf(0)]
for k in range(0,len(P)):
    for t in range(0,len(P[0])):
        Pplot[k][t] = float(P[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)
        
            

##Plotting the pressure
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
maxP = np.nanmax(Pplot)
minP = np.nanmin(Pplot)
plt.title(r'$\Pi_\max$ = $%(number).2f$ $\rmm^{2}$ $\rms^{-2}$        $\Pi_\min$ = $%(number2).2f$ $\rmm^{2}$ $\rms^{-2}$' % {"number": maxP, "number2": minP},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Pplot,np.arange(-3,3.3,0.3),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
cbar = plt.colorbar()
cbar.ax.tick_params(length=14, width=1)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
#plt.xlim([-10000,10000])
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'pressure.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=300)
plt.show()
plt.close()






Yplot,Zplot = np.meshgrid(y,z)
psiplot = np.ones_like(psi)*[mpf(0)]
for k in range(0,len(psi)):
    for t in range(0,len(psi[0])):
        psiplot[k][t] = float(psi[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)


##Plotting the streamlines
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.rcParams['contour.negative_linestyle'] = 'solid'
#plt.title('Streamfunction')
#fig.add_subplot(1,1,1)
#plt.contourf(Y,Z,psi,np.arange(-300,305,5),cmap='seismic')
CS = plt.contour(Yplot,Zplot,psiplot,15,colors='k',linewidths=1.2)
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Yplot,Zplot,Uplot,np.arange(-100000,110000,10000),cmap='seismic')
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
#plt.contourf(Y,Z,psi,cmap='seismic')
jk = 20
# q = plt.quiver(Yplot[::jk,::jk],Zplot[::jk,::jk],Vplot[::jk,::jk],Wplot[::jk,::jk],scale=50,angles="xy")
# plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
# plt.streamplot(Yplot,Zplot,Vplot,Wplot,density = 3,arrowstyle='->',arrowsize = 1.5)
# plt.xlabel("Y axis [m]")
# plt.ylabel("Height [m]")
#plt.xlim([-5000,5000])
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.subplots_adjust(bottom=0.08, top=0.99, hspace=0.1,right=0.97,left=0.13)
nameoffigure = 'streamlines.pdf'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/gome0004/kataprogramms/output/'+string_in_string,dpi=500)
#plt.show()
plt.close()

print("Calculando Vinfinity")
#Getting Vinfinity
Vinf=mpf(0)
for k in range(-K,K+1):
    
    R = mpf(2) * N**2 * mp.cos(alpha)**2 / (visc * diff) * (mpf(k) * pizao / L)**2
    
    Q = N**2 * mp.sin(alpha)**2 / (mpf(3) * visc * diff)
    
    S1 = abs(R + mp.sqrt(Q**3 + R**2) )**(1/3)
    S2 = - abs( mp.sqrt(Q**3 + R**2) -R )**(1/3)
    
    phi = mp.sqrt(S1**2 + S2**2 - S1*S2)
    Lk = mp.acos(- (S1 + S2)/ (2 * phi) )
    
    m1 = - mp.sqrt(S1 + S2)
    m2 = - mp.sqrt(phi) * mp.exp(1j * Lk/2)
    m3 = mp.conj(m2)
    
    
    def f1r(y):
        return (mp.exp(m1 * H(y)) * mp.cos(mpf(2) * mpf(- k) * pizao * y / L) )
    gamma1k = mpf(2)/L * mp.quad(f1r,[0,L/mpf(2)])
        
    def f2r(y):
        return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * mpf(- k) * pizao * y / L) )
    gamma2k = mpf(2)/L * mp.quad(f2r,[0,L/mpf(2)])
    
    if k != 0:
        Vinf = Vinf + ( mp.cos(alpha) / visc * mpf(2) * mpf(k) * pizao / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + mpf(2) * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) )

print("V infinity =", Vinf)

print(V[2000][0])
