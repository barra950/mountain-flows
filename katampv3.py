import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


#Solving the equations for the even case
dps_value = 100

mp.dps = dps_value


K = 70
alpha = mpf(0.1) 
visc = mpf(1)     
diff = mpf(1)     
N = mpf(0.01)    
L = mpf(1000)   

subdivisions = 100

pizao = mp.pi(dps=dps_value)

# def H(y):
#     return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )

def H(y):
    return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * (y-L/mpf(2))/L)) )

def Bsfc(y):
    return mpf(0.1)


final_system = []
b=[]
for q in range(-K,K+1):
    equation1 = []
    equation2 = []
    equation3 = []
    Aki = []
    Cki = []
    Dki = []
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
            return (mp.exp(m1 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma1 = 2/L * mp.quad(f1r,[mpf(0),mpf(L/2)]) 
        
        def f2r(y):
            return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma2 = 2/L * mp.quad(f2r,[mpf(0),mpf(L/2)])
        
        
        if k == 0:
            equation1.append(mpf(2) * gamma2.real)
            Cki.append(k)
            equation1.append(mpf(-2) * gamma2.imag)
            Dki.append(k)
        else:
            equation1.append(gamma1)
            Aki.append(k)
            equation1.append(mpf(2) * gamma2.real)
            Cki.append(k)
            equation1.append(mpf(-2) * gamma2.imag)
            Dki.append(k)
            
        if q != 0:
            
            if k == 0:
                equation2.append(mpf(0))
                equation2.append(mpf(0))
            else:
                equation2.append(mpf(k) * gamma1 / (m1**3) )
                equation2.append(mpf(2) * mpf(k) * (gamma2 / (m2**3) ).real)
                equation2.append(mpf(-2) * mpf(k) * (gamma2 / (m2**3) ).imag)
        
        if k == 0:
            equation3.append(mpf(2) * (m2**2 * gamma2).real)
            equation3.append(mpf(-2) * (m2**2 * gamma2).imag)
        else:
            equation3.append(m1**2 * gamma1)
            equation3.append(mpf(2) * (m2**2 * gamma2).real)
            equation3.append(mpf(-2) * (m2**2 * gamma2).imag)
            
            
    final_system.append(equation1)
    def f4r(y):
        return (Bsfc(y) * mp.cos(mpf(2) * mpf(q) * pizao * y / L) )
    b.append(2/L * mp.quad(f4r,[mpf(0),mpf(L/2)]))
    
    if q != 0:
        final_system.append(equation2)
        b.append(mpf(0))
    
    
    final_system.append(equation3)
    b.append(mpf(0))


final_system = matrix(final_system)
b=matrix(b)
           
            
#Qr solver 
X = qr_solve(final_system, b)

difference = b - (final_system * X[0])

difference2 = []
for k in difference:
    difference2.append(float(k))
    
difference2 = np.array(difference2)

if abs(np.amax(difference2)) < 0.00000000000000001 and abs(np.amin(difference2)) < 0.00000000000000001:
    print(True)
else:
    print(False)



#Getting the values for Ak, Ck and Dk


solution = []
for k in X:
    solution.append(k)

strings = []

for k in range(-K,K+1):
    if k != 0:
        strings.append('A')     
            
    strings.append('R')           
    strings.append('I')             
        


Ak = []
Rk = []
Ik = []
for k in range(0,len(solution[0])):
     if 'A' in strings[k]:
         Ak.append(solution[0][k])
         
     if 'R' in strings[k]:
         Rk.append(solution[0][k])
        
     if 'I' in strings[k]:
         Ik.append(solution[0][k])
         
Ck=[]
      
for k in range(0,len(Rk)):
    Ck.append(Rk[k] + Ik[k] * 1j)

Ck = np.array(Ck)

Dk = Ck.conjugate()

Ak = np.array(Ak)

# #Now, we start converting out of mpf to float type
# for k in range(0,len(Ak)):
#     Ak[k] = float(Ak[k].real) + float(Ak[k].imag)*1j 
            
# for k in range(0,len(Ck)):
#     Ck[k] = float(Ck[k].real) + float(Ck[k].imag)*1j
    
# for k in range(0,len(Dk)):
#     Dk[k] = float(Dk[k].real) + float(Dk[k].imag)*1j
            
  
# Xpost = []
# for k in range(0,len(solution[0])):
#     Xpost.append(solution[0][k])

# Xpost = np.array(Xpost)

# for k in range(0,len(Xpost)):
#     Xpost[k] = float(Xpost[k].real) + float(Xpost[k].imag)*1j
    
   
#Getting the Buoyancy value

z = np.arange(0,2010,10) 
y = np.arange(-float(L),float(L)+10,10) 
Y,Z = np.meshgrid(y,z)
Y = Y * mpf(1)
Z = Z * mpf(1)
B = np.ones_like(Y)*[mpf(0)]

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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                B[i][t] = B[i][t] + ( Ak[Aki.index(k)] * mp.exp(m1 * Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)  )
            B[i][t] = B[i][t] + ( ( Ck[Cki.index(k)] * mp.exp(m2 * Z[i][t]) + Dk[Dki.index(k)] * mp.exp(m3 * Z[i][t]) )  * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) )


for k in range(0,B.shape[0]):
    for t in range(0,B.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            B[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (B[k][t], "B value at the ground")
#         if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#             if B[k][t] > 0.101:
#                 print (B[k][t],'fudeu geral -------------------------------------------------')
# #            print (B[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )
    
#Bp = Bsfc(Y) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.cos(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )

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
plt.title(r'B$_\max$ = %(number).2f m $\rms^{-2}$        B$_\min$ = %(number2).2f m $\rms^{-2}$' % {"number": maxB, "number2": minB},x=0.5, y=1.02)
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

nameoffigure = 'buoyancy.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()   

#Buoyancy for the prandtl case 
Bp = np.ones_like(B)*[mpf(0)]
for i in range(0,len(Y)):
    for t in range(0,len(Y[0])):  
        Bp[i][t] = Bsfc(Y[i][t]) * mp.exp(-Z[i][t] * mp.sqrt(N * mp.sin(alpha) ) / (mpf(4)*visc*diff)**(1/4) ) * mp.cos(mp.sqrt(N*mp.sin(alpha)) /((mpf(4)*visc*diff)**(1/4))*Z[i][t] )

   


#Getting the value of the V wind
# z = np.arange(0,2020,20) 
# y = np.arange(-float(L),float(L)+10,10) 
# Y,Z = np.meshgrid(y,z)
# Y = Y * mpf(1)
# Z = Z * mpf(1)
V = np.ones_like(Y)*[mpf(0)] 


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
    
    if k != 0:
        for i in range(0,len(Y)):
            for t in range(0,len(Y[0])):
                
                V[i][t] = V[i][t] + mp.cos(alpha)/visc * 2j*mpf(k)*pizao/L *  ( Ak[Aki.index(k)]*mp.exp(m1*Z[i][t])/(m1**3) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**3)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**3) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
    

for k in range(0,V.shape[0]):
    for t in range(0,V.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            V[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (V[k][t], "V value at ground")
#         if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#             if V[k][t] > 0.1:
#                 print (V[k][t],'fudeu geral -------------------------------------------------')
# #            print (V[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


Yplot,Zplot = np.meshgrid(y,z)
Vplot = np.ones_like(V)*[mpf(0)]
for k in range(0,len(V)):
    for t in range(0,len(V[0])):
        Vplot[k][t] = float(V[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)


##Plotting the V wind
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
maxV = np.nanmax(Vplot)
minV = np.nanmin(Vplot)
plt.title(r'V$_\max$ = %(number).2f m $\rms^{-1}$        V$_\min$ = %(number2).2f m $\rms^{-1}$' % {"number": maxV, "number2": minV},x=0.5, y=1.02)
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
nameoffigure = 'Vwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()  


#Getting the value of the U wind
#We first need the value of Eq
Eq=[]
Eqi=[]
for q in range(-K,K+1):
    E = mpf(0)
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
            return (mp.exp(m1 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma1 = 2/L * mp.quad(f1r,[mpf(0),mpf(L/2)])
        
        def f2r(y):
            return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma2 = 2/L * mp.quad(f2r,[mpf(0),mpf(L/2)])
        
            
        if k != 0:
            E = E - mp.sin(alpha)/visc * Ak[Aki.index(k)]*gamma1/(m1**2) 
        E = E - mpf(2)*mp.sin(alpha)/visc * ( Ck[Cki.index(k)]*gamma2/(m2**2) ).real
    
    Eq.append(E)
    Eqi.append(q)

Eq= np.array(Eq)




#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
U = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                U[i][t] = U[i][t] + mp.sin(alpha)/visc * Ak[Aki.index(k)]/(m1**2) * mp.exp(m1*Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
            U[i][t] = U[i][t] + mp.sin(alpha)/visc * ( Ck[Cki.index(k)]/(m2**2)*mp.exp(m2*Z[i][t]) + Dk[Dki.index(k)]/(m3**2)*mp.exp(m3*Z[i][t]) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) + Eq[Eqi.index(k)] * mp.cos(mpf(2) * mpf(k) * pizao * Y[i][t] / L)


for k in range(0,U.shape[0]):
    for t in range(0,U.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            U[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (U[k][t], "U value at ground")
#         if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#             if U[k][t] > 0.1:
#                 print (U[k][t],'fudeu geral -------------------------------------------------')
# #            print (U[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


#U for prandtl case:
#Up = -Bsfc(Y)/N * np.sqrt(diff/visc) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.sin(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )


Yplot,Zplot = np.meshgrid(y,z)
Uplot = np.ones_like(U)*[mpf(0)]
for k in range(0,len(U)):
    for t in range(0,len(U[0])):
        Uplot[k][t] = float(U[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)

#Plotting the U wind
fig,ax1 = plt.subplots(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
maxU = np.nanmax(Uplot)
minU = np.nanmin(Uplot)
plt.title(r'U$_\max$ = %(number).2f m $\rms^{-1}$        U$_\min$ = %(number2).2f m $\rms^{-1}$' % {"number": maxU, "number2": minU},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Uplot,np.arange(-14,15,1),cmap='seismic')
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
nameoffigure = 'Uwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()


#Uwind for the prandtl case 
Up = np.ones_like(U)*[mpf(0)]
for i in range(0,len(Y)):
    for t in range(0,len(Y[0])):  
        Up[i][t] = -Bsfc(Y[i][t])/N * mp.sqrt(diff/visc) * mp.exp(-Z[i][t] * mp.sqrt(N * mp.sin(alpha) ) / (mpf(4)*visc*diff)**(1/4) ) * mp.sin(mp.sqrt(N*mp.sin(alpha)) /((mpf(4)*visc*diff)**(1/4))*Z[i][t] )

   

#PLotting the W wind
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
W = np.ones_like(Y)*[mpf(0)]


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
    
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                W[i][t] = W[i][t] + mp.cos(alpha)/visc * mpf(4)*mpf(k)**2 *pizao**2 / L**2 * (Ak[Aki.index(k)]*mp.exp(m1*Z[i][t])/(m1**4) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**4)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**4)   ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) +  mp.tan(alpha)*(Eq[Eqi.index(k)] * mp.cos(mpf(2) * mpf(k) * pizao * Y[i][t] / L) )
            else:
                W[i][t] = W[i][t] + mp.tan(alpha)*(Eq[Eqi.index(k)] * mp.cos(mpf(2) * mpf(k) * pizao * Y[i][t] / L) )
        


for k in range(0,W.shape[0]):
    for t in range(0,W.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            W[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (W[k][t], "W value at ground")
#         if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#             if W[k][t] > 0.1:
#                 print (W[k][t],'fudeu geral -------------------------------------------------')
# #            print (W[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )

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
plt.title(r'W$_\max$ = %(number).2f m $\rms^{-1}$        W$_\min$ = %(number2).2f m $\rms^{-1}$' % {"number": maxW, "number2": minW},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Wplot,np.arange(-6.2,6.3,0.2),cmap='seismic')
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
nameoffigure = 'Wwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()



#Now, we need to get the values of the pressure 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
P = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                P[i][t] = P[i][t] + mp.cos(alpha) * Ak[Aki.index(k)] / m1 * mp.exp(m1*Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
            P[i][t] = P[i][t] + mp.cos(alpha) * (Ck[Cki.index(k)] / m2 * mp.exp(m2*Z[i][t]) + Dk[Dki.index(k)] / m3 * mp.exp(m3*Z[i][t])) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)


for k in range(0,P.shape[0]):
    for t in range(0,P.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            P[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (P[k][t], "P value at ground")
            
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
plt.title(r'$\Pi_\max$ = %(number).2f $\rmm^{2}$ $\rms^{-2}$        $\Pi_\min$ = %(number2).2f $\rmm^{2}$ $\rms^{-2}$' % {"number": maxP, "number2": minP},x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Pplot,np.arange(-10,10.5,0.5),cmap='seismic')
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
nameoffigure = 'pressure.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
plt.show()
plt.close()


#Calculating the streamlines
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
psi = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                psi[i][t] = psi[i][t] - mp.tan(alpha) * 1j*L* Eq[Eqi.index(k)] / (mpf(2)*mpf(k)*pizao) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) + Eq[Eqi.index(0)] * Y[i][t] * mp.tan(alpha) - mp.cos(alpha)/visc * 2j*mpf(k)*pizao/L *  ( Ak[Aki.index(k)]*mp.exp(m1*Z[i][t])/(m1**4) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**4)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**4) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)


for k in range(0,psi.shape[0]):
    for t in range(0,psi.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            psi[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (psi[k][t], "psi value at ground")
            
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
CS = plt.contour(Yplot,Zplot,psiplot,50,colors='k',linewidths=0.7)
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
nameoffigure = 'streamlines.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi = 500)
#plt.show()
plt.close()



#Calculating the U* and plotting it  
Ustar = np.ones_like(y)*[mpf(0)]
ympf = y*[mpf(1)]
for k in range(-K,K+1):
    for t in range(0,len(y)):
        Ustar[t] = Ustar[t] +  Eq[Eqi.index(k)] * mp.cos(mpf(2) * mpf(k) * pizao * ympf[t] / L)/mp.cos(alpha)
        
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
plt.ylim([-14,4])
plt.yticks(np.arange(-14,6,2))
plt.grid('True')

nameoffigure = 'Ustar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()


#Testing whether V integral is constant with height
Vintegral = []
Hmax = 500
for k in range(0,len(V)):
    if Zplot[k][0] > Hmax:
        Vsum = 0
        for t in range(0,len(V[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Vsum = Vsum + Vplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Vintegral.append(Vsum)
        #print (Z[k][0])
print ("Max and min values of Vintegral", np.amax( Vintegral), np.amin( Vintegral))

jato=10
for k in range(0,len(Vintegral)):
    print(Vintegral[k].real,500+jato,"m" )
    jato = jato + 10
    

#Testing whether W integral is constant with height
Wintegral = []
for k in range(0,len(W)):
    if Zplot[k][0] > Hmax:
        Wsum = 0
        for t in range(0,len(Wplot[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Wsum = Wsum + Wplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Wintegral.append(Wsum)
        #print (Z[k][0])
print ("Max and min values of Wintegral", np.amax( Wintegral), np.amin( Wintegral))

jato=10
for k in range(0,len(Wintegral)):
    print(Wintegral[k].real,500+jato,"m" )
    jato = jato + 10
    
#Testing whether Ustar integral is constant with height
Ustsum = 0
for t in range(0,len(y)):
    if (y[t] >= float(-L/2) and y[t] <= float(L/2)):
        Ustsum = Ustsum + Ustarplot[t]*abs(y[0]-y[1])
print ("Ustar integral value", Ustsum)

#Getting the U integral (not constant with height)
Uintegral = []
integralZ = []
for k in range(0,len(U)):
    if Zplot[k][0] > Hmax:
        Usum = 0
        for t in range(0,len(U[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Usum = Usum + Uplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Uintegral.append(Usum)
        integralZ.append(z[k])
        #print (Z[k][0])
print ("Max and min values of Uintegral", np.amax( Uintegral), np.amin( Uintegral))


#Getting the B integral
Bintegral = []
for k in range(0,len(B)):
    if Zplot[k][0] > Hmax:
        Bsum = 0
        for t in range(0,len(B[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Bsum = Bsum + Bplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Bintegral.append(Bsum)
        #print (Z[k][0])
print ("Max and min values of Bintegral", np.amax( Bintegral), np.amin( Bintegral))


#Testing wether C is a constant
Uintegral = np.array(Uintegral)
integralZ = np.array(integralZ)
Bintegral = np.array(Bintegral)
delta = float((4*visc*diff)**(1/4)) / np.sqrt(float(N) * float(mp.sin(alpha)))
Cconst = np.sqrt( np.exp(2 * integralZ/delta) * (Bintegral**2 + float(visc/diff) * float(N)**2 * Uintegral**2 ) )


#Plotting the prandtl solution
Bpplot = np.ones_like(Bp)*[mpf(0)]
Upplot = np.ones_like(Up)*[mpf(0)]
for k in range(0,len(U)):
    for t in range(0,len(U[0])):
        Bpplot[k][t] = float(Bp[k][t].real)
        Upplot[k][t] = float(Up[k][t].real)

Bpplot2 = np.array(Bpplot).T
Upplot2 = np.array(Upplot).T
fig,ax1=plt.subplots(figsize=(10,10))
#plt.xticks(maxtime, time2plot, rotation='vertical')
plt.rcParams.update({'font.size':16})
fsize = 20

ax1.set_ylim([-0.02,0.1])
plt.xlabel('Z [m]',name='Arial')
plt.ylabel(r'B [m $\rms^{-2}$]',name='Arial')
plt.plot(z,Bpplot2[0][:],linewidth=3,color='b',label='B')
ax1.tick_params('both', length=10, width=1, which='major')
ax2=ax1.twinx()
ax2.set_ylim([-4,4])
plt.plot(z,Upplot2[0][:],linewidth=3,color='r',label='U')
plt.ylabel(r'U [m $\rms^{-1}$]',name='Arial')
ax1.set_xlim([0,1000])
ax1.set_xticks(np.arange(0,1100,100))
ax2.tick_params('both', length=10, width=1, which='major')
ax1.legend(bbox_to_anchor=(0.7, 0.95), loc='upper left', borderaxespad=0)
ax2.legend(bbox_to_anchor=(0.7, 0.88), loc='upper left', borderaxespad=0)

nameoffigure = 'prandtl.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
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
        thetaplot_sign[k][t] = float(theta[k][t].real * abs(Wstar_real[k][t].real)/Wstar_real[k][t].real)
        
#Plotting the theta angle vs height at y=0 
thetaplotf = (thetaplot).T
thetaplotf2 = (thetaplot_sign).T 
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
ax1.plot([0,0],[0,1500],linewidth=2,color="k")
ax1.plot(thetaplotf2[100][1:],z[1:],linewidth=3,color="r") #Number might depends on how many points y has
plt.xlabel(r"$\rm\theta$ [$\rm^{o}$]")
plt.ylabel("Z [m]")
plt.xlim([-0.5,6])
ax1.set_ylim([0,1500])
ax1.set_xticks(np.arange(-0.5,6.5,0.5))
ax1.set_yticks(np.arange(0,1600,100))
ax1.tick_params('both', length=10, width=1, which='major')
plt.grid(True)

nameoffigure = 'theta.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
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
        boundsl.append(0.6/interval*k)
    if k > 0:
        boundsl.append(6/interval*k)
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
plt.title(r'W$^{\bigstar}_\max$ = %(number).2f m $\rms^{-1}$        W$^{\bigstar}_\min$ = %(number2).2f m $\rms^{-1}$' % {"number": maxWstar, "number2": minWstar},x=0.5, y=1.02)

norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# pcm = plt.pcolormesh(Yplot,Zplot,Yploy,norm=norm,cmap='seismic')

pcm = plt.contourf(Yplot,Zplot,Wstarplot,norm=norm,cmap='seismic',levels=bounds)

cbar = fig.colorbar(pcm)
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

nameoffigure = 'Wstar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string,dpi=300)
#plt.show()
plt.close()






#%%

import pickle
import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


dps_value = 100

mp.dps = dps_value


K = 70
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(10)     
N = mpf(0.01)    
L = mpf(1000)

subdivisions = 100

pizao = mp.pi(dps=dps_value)

# def H(y):
#     return ( mpf(200) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )

def H(y):
    return ( mpf(200) * (mpf(1) + mp.cos(mpf(2) * pizao * (y-L/mpf(2))/L)) )

def Bsfc(y):
    return mpf(0.1)

z = np.arange(0,2010,10) 
y = np.arange(-float(L),float(L)+10,10) 
Y,Z = np.meshgrid(y,z)
Y = Y * mpf(1)
Z = Z * mpf(1)


with open('/home/owner/Documents/katabatic_flows/variables/Ak.pickle', 'wb') as handle:
    pickle.dump(Ak, handle, protocol=pickle.HIGHEST_PROTOCOL)


    


with open('/home/owner/Documents/katabatic_flows/variables/Aki.pickle', 'wb') as handle:
    pickle.dump(Aki, handle, protocol=pickle.HIGHEST_PROTOCOL)


    
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Ck.pickle', 'wb') as handle:
    pickle.dump(Ck, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Cki.pickle', 'wb') as handle:
    pickle.dump(Cki, handle, protocol=pickle.HIGHEST_PROTOCOL)





with open('/home/owner/Documents/katabatic_flows/variables/Dk.pickle', 'wb') as handle:
    pickle.dump(Dk, handle, protocol=pickle.HIGHEST_PROTOCOL)


    

with open('/home/owner/Documents/katabatic_flows/variables/Dki.pickle', 'wb') as handle:
    pickle.dump(Dki, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    
       

with open('/home/owner/Documents/katabatic_flows/variables/B.pickle', 'wb') as handle:
    pickle.dump(B, handle, protocol=pickle.HIGHEST_PROTOCOL)

 
    
 
with open('/home/owner/Documents/katabatic_flows/variables/V.pickle', 'wb') as handle:
    pickle.dump(V, handle, protocol=pickle.HIGHEST_PROTOCOL)




with open('/home/owner/Documents/katabatic_flows/variables/VlargeY.pickle', 'wb') as handle:
    pickle.dump(V, handle, protocol=pickle.HIGHEST_PROTOCOL)



with open('/home/owner/Documents/katabatic_flows/variables/U.pickle', 'wb') as handle:
    pickle.dump(U, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eq.pickle', 'wb') as handle:
    pickle.dump(Eq, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eqi.pickle', 'wb') as handle:
    pickle.dump(Eqi, handle, protocol=pickle.HIGHEST_PROTOCOL)

 
with open('/home/owner/Documents/katabatic_flows/variables/W.pickle', 'wb') as handle:
    pickle.dump(W, handle, protocol=pickle.HIGHEST_PROTOCOL)


    
with open('/home/owner/Documents/katabatic_flows/variables/WlargeY.pickle', 'wb') as handle:
    pickle.dump(W, handle, protocol=pickle.HIGHEST_PROTOCOL)



with open('/home/owner/Documents/katabatic_flows/variables/P.pickle', 'wb') as handle:
    pickle.dump(P, handle, protocol=pickle.HIGHEST_PROTOCOL)

    

with open('/home/owner/Documents/katabatic_flows/variables/psi.pickle', 'wb') as handle:
    pickle.dump(psi, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
####################################################################################################################################

with open('/home/owner/Documents/katabatic_flows/variables/Ak.pickle', 'rb') as handle:
    Ak = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Aki.pickle', 'rb') as handle:
    Aki = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Ck.pickle', 'rb') as handle:
    Ck = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Cki.pickle', 'rb') as handle:
    Cki = pickle.load(handle)   
    

with open('/home/owner/Documents/katabatic_flows/variables/Dk.pickle', 'rb') as handle:
    Dk = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Dki.pickle', 'rb') as handle:
    Dki = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/B.pickle', 'rb') as handle:
    B = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/V.pickle', 'rb') as handle:
    V = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/VlargeY.pickle', 'rb') as handle:
    V = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/U.pickle', 'rb') as handle:
    U = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Eq.pickle', 'rb') as handle:
    Eq = pickle.load(handle)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eqi.pickle', 'rb') as handle:
    Eqi = pickle.load(handle)
 
    
with open('/home/owner/Documents/katabatic_flows/variables/W.pickle', 'rb') as handle:
    W = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/WlargeY.pickle', 'rb') as handle:
    W = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/P.pickle', 'rb') as handle:
    P = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/psi.pickle', 'rb') as handle:
    psi = pickle.load(handle)




#%%
#Calculating some equations
zmp = z * mpf(1)
du2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):   
        du2dz2[k,t] = ( U[k+1,t] - mpf(2)*U[k,t] + U[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

term1 = -B[1:-1,:]*mp.sin(alpha)
term2 = visc*du2dz2[1:-1,:]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '666666666666666666666666666666666666666666')
    

########################################################################################################################################################################

zmp = z * mpf(1)
dv2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dv2dz2[k,t] = ( V[k+1,t] - mpf(2)*V[k,t] + V[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

ymp = y * mpf(1)
dpdy=np.ones_like(Y)*[mpf(0)]
for k in range(1,len(y)-1):
    for t in range(0,len(z)):
        dpdy[t,k] = ( P[t,k+1] - P[t,k-1] ) / abs(ymp[k-1]-ymp[k+1])
    
term1 = -dpdy[1:-1,1:-1]
term2 = visc*dv2dz2[1:-1,1:-1]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '77777777777777777777777777777777777777')   
             
########################################################################################################################################################################        

zmp = z * mpf(1)        
dpdz=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dpdz[k,t] = ( P[k+1,t] - P[k-1,t] ) / abs(zmp[k-1]-zmp[k+1])

term1 = -dpdz[1:-1,:] 
term2 = B[1:-1,:]*mp.cos(alpha)

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '99999999999999999999999999999')       
        
#############################################################################################################################################################        

zmp = z * mpf(1)         
db2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        db2dz2[k,t] = ( B[k+1,t] - mpf(2)*B[k,t] + B[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

term1 = N**2 * (U[1:-1,:]*mp.sin(alpha) - W[1:-1,:]*mp.cos(alpha))
term2 = diff*db2dz2[1:-1,:]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '8888888888888888888888888888888888888' )       


#############################################################################################################################################################       
        
ymp = y * mpf(1)       
dvdy=np.ones_like(Y)*[mpf(0)]
for k in range(1,len(y)-1):
    for t in range(0,len(z)):
        dvdy[t,k] = ( V[t,k+1] - V[t,k-1] ) / abs(ymp[k-1]-ymp[k+1])

zmp = z * mpf(1)    
dwdz=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dwdz[k,t] = ( W[k+1,t] - W[k-1,t] ) / abs(zmp[k-1]-zmp[k+1])


term1 = dvdy[1:-1,1:-1]
term2 = dwdz[1:-1,1:-1]  

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum) 

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '55555555555555555555555555555555555555555555555555') 

#%%

import pickle
import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

#Solving the equations for the odd case
dps_value = 100

mp.dps = dps_value


K = 70
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(5)     
N = mpf(0.01)    
L = mpf(1000)

subdivisions = 100

pizao = mp.pi(dps=dps_value)

def H(y):
    return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )


def Bsfc(y):
    return mpf(0.1) * mp.sin(mpf(2) * pizao * y/L)


final_system = []
b=[]
for q in range(-K,K+1):
    equation1 = []
    equation2 = []
    equation3 = []
    Eki = []
    Cki = []
    Dki = []
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
            return (mp.exp(m1 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma1 = 2/L * mp.quad(f1r,[mpf(0),mpf(L/2)]) 
        
        def f2r(y):
            return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma2 = 2/L * mp.quad(f2r,[mpf(0),mpf(L/2)])  
        
        
        if k == 0:
            equation1.append(mpf(2) * gamma2.imag)
            Cki.append(k)
            equation1.append(mpf(2) * gamma2.real)
            Dki.append(k)
        else:
            equation1.append(gamma1)
            Eki.append(k)
            equation1.append(mpf(2) * gamma2.imag)
            Cki.append(k)
            equation1.append(mpf(2) * gamma2.real)
            Dki.append(k)
            
        if q != 0:
            
            if k == 0:
                equation2.append(mpf(0))
                equation2.append(mpf(0))
            else:
                equation2.append(mpf(k) * gamma1 / (m1**3) )
                equation2.append(mpf(2) * mpf(k) * (gamma2 / (m2**3) ).imag)
                equation2.append(mpf(2) * mpf(k) * (gamma2 / (m2**3) ).real)
        
        if k == 0:
            equation3.append(mpf(2) * (m2**2 * gamma2).imag)
            equation3.append(mpf(2) * (m2**2 * gamma2).real)
        else:
            equation3.append(m1**2 * gamma1)
            equation3.append(mpf(2) * (m2**2 * gamma2).imag)
            equation3.append(mpf(2) * (m2**2 * gamma2).real)
    
    final_system.append(equation1)
    def f4r(y):
        return (Bsfc(y) * mp.sin(mpf(2) * mpf(q) * pizao * y / L) )
    b.append(-2/L * mp.quad(f4r,[mpf(0),mpf(L/2)]))
    
    if q != 0:
        final_system.append(equation2)
        b.append(mpf(0))
    
    
    final_system.append(equation3)
    b.append(mpf(0))


final_system = matrix(final_system)
b=matrix(b)


#Qr solver 
X = qr_solve(final_system, b)

difference = b - (final_system * X[0])

difference2 = []
for k in difference:
    difference2.append(float(k))
    
difference2 = np.array(difference2)

if abs(np.amax(difference2)) < 0.00000000000000001 and abs(np.amin(difference2)) < 0.00000000000000001:
    print(True)
else:
    print(False)
    
    
#Getting the values for Ek, Ck and Dk


solution = []
for k in X:
    solution.append(k)

strings = []

for k in range(-K,K+1):
    if k != 0:
        strings.append('E')     
            
    strings.append('R')           
    strings.append('I')             
        


Ek = []
Rk = []
Ik = []
for k in range(0,len(solution[0])):
     if 'E' in strings[k]:
         Ek.append(solution[0][k])
         
     if 'R' in strings[k]:
         Rk.append(solution[0][k])
        
     if 'I' in strings[k]:
         Ik.append(solution[0][k])
         
Ck=[]
      
for k in range(0,len(Rk)):
    Ck.append(Rk[k] + Ik[k] * 1j)

Ck = np.array(Ck)

Dk = - Ck.conjugate()

Ek = np.array(Ek)



#Getting the Buoyancy value

z = np.arange(0,2001,1) 
y = np.arange(-float(L),float(L)+10,10) 
Y,Z = np.meshgrid(y,z)
Y = Y * mpf(1)
Z = Z * mpf(1)
B = np.ones_like(Y)*[mpf(0)]

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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                B[i][t] = B[i][t] + ( 1j*Ek[Eki.index(k)] * mp.exp(m1 * Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)  )
            B[i][t] = B[i][t] + ( ( Ck[Cki.index(k)] * mp.exp(m2 * Z[i][t]) + Dk[Dki.index(k)] * mp.exp(m3 * Z[i][t]) )  * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) )


for k in range(0,B.shape[0]):
    for t in range(0,B.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            B[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (B[k][t], "B value at the ground")
#         if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#             if B[k][t] > 0.101:
#                 print (B[k][t],'fudeu geral -------------------------------------------------')
# #            print (B[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )
    
#Bp = Bsfc(Y) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.cos(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )

Yplot,Zplot = np.meshgrid(y,z)
Bplot = np.ones_like(B)*[mpf(0)]
for k in range(0,len(B)):
    for t in range(0,len(B[0])):
        Bplot[k][t] = float(B[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)
       

    
##Plotting the buoyancy
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title(r'B$_\max$ 5 $\rmms^{-2}$        B$_\min$ 5 $\rmms^{-2}$',x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Bplot,np.arange(-0.1,0.11,0.01),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
#plt.colorbar(label='[$ms^{-2}$]')
plt.colorbar()
plt.xlabel("Y axis [m]")
plt.ylabel("Height [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'buoyancy.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



#Getting the value of the V wind
# Y,Z = np.meshgrid(y,z)
# Y = Y * mpf(1)
# Z = Z * mpf(1)
V = np.ones_like(Y)*[mpf(0)]


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
    gamma1k = mpf(2)/L * mp.quad(f1r,[mpf(0),L/mpf(2)])
        
    def f2r(y):
        return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * mpf(- k) * pizao * y / L) )
    gamma2k = mpf(2)/L * mp.quad(f2r,[mpf(0),L/mpf(2)])
    
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                V[i][t] = V[i][t] + mp.cos(alpha)/visc * 2j*mpf(k)*pizao/L *  ( 1j * Ek[Eki.index(k)]*mp.exp(m1*Z[i][t])/(m1**3) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**3)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**3) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
                V[i][t] = V[i][t] + mp.cos(alpha) / visc * mpf(2) * mpf(k) * pizao / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + mpf(2) * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) 

for k in range(0,V.shape[0]):
    for t in range(0,V.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            V[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (V[k][t], "V value at ground")
#        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if V[k][t] > 0.1:
#                print (V[k][t],'fudeu geral -------------------------------------------------')
#            print (V[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


Yplot,Zplot = np.meshgrid(y,z)
Vplot = np.ones_like(V)*[mpf(0)]
for k in range(0,len(V)):
    for t in range(0,len(V[0])):
        Vplot[k][t] = float(V[k][t].real) 

##Plotting the V wind
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title(r'V$_\max$ = 5 m $\rms^{-2}$        V$_\min$ = 5 m $\rms^{-2}$',x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Vplot,np.arange(-7,7.03,0.03),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
plt.colorbar()
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Vwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()  


#Getting the value of the U wind
#We first need the value of Eq
Eq=[]
Eqi=[]
for q in range(-K,K+1):
    E = mpf(0)
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
            return (mp.exp(m1 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma1 = mpf(2)/L * mp.quad(f1r,[mpf(0),L/mpf(2)])
        
        def f2r(y):
            return (mp.exp(m2 * H(y)) * mp.cos(mpf(2) * (mpf(q) - mpf(k)) * pizao * y / L) )
        gamma2 = mpf(2)/L * mp.quad(f2r,[0,L/2])
        
#        def f3r(y):
#            return (np.exp(m3 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
#        def f3i(y):
#            return (np.exp(m3 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
#        gamma3 = 2/L * (quad(f3r,0,L/2)[0] + quad(f3i,0,L/2)[0]*1j) 
        
#        gamma1 = 0.0
#        gamma2 = 0.0
#        gamma3 = 0.0
#        for y in points:
#            gamma1 = 2/L * f1(y)*tick + gamma1
#            gamma2 = 2/L * f2(y)*tick + gamma2
#            #gamma3 = 2/L * f3(y)*tick + gamma3
            
        if k != 0:
            E = E - 1j * mp.sin(alpha)/visc * ( Ek[Eki.index(k)]*gamma1/(m1**2) + mpf(2) * ( Ck[Cki.index(k)] * gamma2 / (m2**2) ).imag )
        
    
    Eq.append(E)
    Eqi.append(q)

Eq= np.array(Eq)



# Y,Z = np.meshgrid(y,z)
# Y = Y * mpf(1)
# Z = Z * mpf(1)
U = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                U[i][t] = U[i][t] + mp.sin(alpha)/visc * 1j*Ek[Eki.index(k)]/(m1**2) * mp.exp(m1*Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
            U[i][t] = U[i][t] + mp.sin(alpha)/visc * ( Ck[Cki.index(k)]/(m2**2)*mp.exp(m2*Z[i][t]) + Dk[Dki.index(k)]/(m3**2)*mp.exp(m3*Z[i][t]) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) + 1j *Eq[Eqi.index(k)] * mp.sin(mpf(2) * mpf(k) * pizao * Y[i][t] / L)


for k in range(0,U.shape[0]):
    for t in range(0,U.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            U[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (U[k][t], "U value at ground")
#        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if U[k][t] > 0.1:
#                print (U[k][t],'fudeu geral -------------------------------------------------')
##            print (U[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )

Yplot,Zplot = np.meshgrid(y,z)
Uplot = np.ones_like(U)*[mpf(0)]
for k in range(0,len(U)):
    for t in range(0,len(U[0])):
        Uplot[k][t] = float(U[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)


#Plotting the U wind
fig,ax1 = plt.subplots(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title(r'U$_\max$ = 5 m $\rms^{-1}$        U$_\min$ = 5 m $\rms^{-1}$',x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Uplot,np.arange(-25,26,1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
plt.colorbar()
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Uwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



#Now, we need to get the values of the W wind 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
W = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                W[i][t] = W[i][t] + mp.cos(alpha)/visc * mpf(4)*mpf(k)**2 *pizao**2 / L**2 * (1j*Ek[Eki.index(k)]*mp.exp(m1*Z[i][t])/(m1**4) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**4)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**4)   ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L) +  mp.tan(alpha)*(1j*Eq[Eqi.index(k)] * mp.sin(mpf(2) * mpf(k) * pizao * Y[i][t] / L) )
            else:
                W[i][t] = W[i][t] + mp.tan(alpha)*(1j*Eq[Eqi.index(k)] * mp.sin(mpf(2) * mpf(k) * pizao * Y[i][t] / L) )
        


for k in range(0,W.shape[0]):
    for t in range(0,W.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            W[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (W[k][t], "W value at ground")
#        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if W[k][t] > 0.1:
#                print (W[k][t],'fudeu geral -------------------------------------------------')
##            print (W[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


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
plt.title(r'W$_\max$ = 5 m $\rms^{-1}$        W$_\min$ = 5 m $\rms^{-1}$',x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Wplot,np.arange(-6,6.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
#plt.colorbar(label='[$ms^{-1}$]')
plt.colorbar()
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])

ax1.tick_params('both', length=14, width=1, which='major')
ax1.tick_params('both', length=7, width=1, which='minor')

ax1.minorticks_on()
ax1.xaxis.set_tick_params(which='minor', bottom=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
nameoffigure = 'Wwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



#Now, we need to get the values of the pressure 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
P = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                P[i][t] = P[i][t] + mp.cos(alpha) * 1j *Ek[Eki.index(k)] / m1 * mp.exp(m1*Z[i][t]) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
            P[i][t] = P[i][t] + mp.cos(alpha) * (Ck[Cki.index(k)] / m2 * mp.exp(m2*Z[i][t]) + Dk[Dki.index(k)] / m3 * mp.exp(m3*Z[i][t])) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)


for k in range(0,P.shape[0]):
    for t in range(0,P.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            P[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (P[k][t], "P value at ground")
        
            

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
plt.title(r'$\Pi_\max$ = 5 $\rmm^{2}$ $\rms^{-2}$        $\Pi_\min$ = 5 $\rmm^{2}$ $\rms^{-2}$',x=0.5, y=1.02)
plt.contourf(Yplot,Zplot,Pplot,np.arange(-15,15.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar()
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
nameoffigure = 'pressure.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
plt.show()
plt.close()



#Calculating the streamlines 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
#Y,Z = np.meshgrid(y,z)
psi = np.ones_like(Y)*[mpf(0)]


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
    
    for i in range(0,len(Y)):
        for t in range(0,len(Y[0])):
            if k != 0:
                psi[i][t] = psi[i][t] - mp.tan(alpha) * 1j*L* Eq[Eqi.index(k)] / (mpf(2)*mpf(k)*pizao) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)  - mp.cos(alpha)/visc * 2j*mpf(k)*pizao/L *  ( 1j*Ek[Eki.index(k)]*mp.exp(m1*Z[i][t])/(m1**4) + Ck[Cki.index(k)]*mp.exp(m2*Z[i][t])/(m2**4)  + Dk[Dki.index(k)]*mp.exp(m3*Z[i][t])/(m3**4) ) * mp.exp(2j * mpf(k) * pizao * Y[i][t] / L)
                psi[i][t] = psi[i][t] - ( mp.cos(alpha) / visc * mpf(2) * mpf(k) * pizao / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + mpf(2) * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) ) * Z[i][t]
psi[i][t] = psi[i][t] + Eq[Eqi.index(0)] * Y[i][t] * mp.tan(alpha)

for k in range(0,psi.shape[0]):
    for t in range(0,psi.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            psi[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (psi[k][t], "psi value at ground")
#        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if psi[k][t] > 0.1:
#                print (psi[k][t],'fudeu geral -------------------------------------------------')
##            print (psi[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


Yplot,Zplot = np.meshgrid(y,z)
psiplot = np.ones_like(psi)*[mpf(0)]
for k in range(0,len(psi)):
    for t in range(0,len(psi[0])):
        psiplot[k][t] = float(psi[k][t].real) #+ 1j*float(B[k][t].imag)
        # if B[k][t].real < 0 and abs(B[k][t].real) > 0.1:
        #     print(B[k][t].real)


##Plotting the streamlines
fig,ax1 = plt.subplots(figsize=(20,20)) 
plt.rcParams.update({'font.size':16})
plt.rcParams['contour.negative_linestyle'] = 'solid'
#plt.title('Streamfunction')
#fig.add_subplot(1,1,1)
#plt.contourf(Y,Z,psi,np.arange(-300,305,5),cmap='seismic')
CS = plt.contour(Yplot,Zplot,psiplot,50,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Yplot,Zplot,Uplot,np.arange(-100000,110000,10000),cmap='seismic')
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
plt.subplots_adjust(bottom=0.07, top=0.99, hspace=0.1,right=0.98,left=0.05)
nameoffigure = 'streamlines.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
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

        
##Plotting U* and streamlines
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


fig.add_subplot(2,1,1)
#plt.title('U star plot')
plt.plot(y,Ustarplot)
plt.xlabel("Y [m]")
plt.ylabel('U$^{\u2605}_\infty$ [m s$^{-1}$]')
#plt.xlim([-5000,5000])
plt.xlim([float(-L),float(L)])
#plt.ylim([-6,2])
plt.ylim([-18,10])
plt.grid('True')

fig.add_subplot(2,1,2)
#CS = plt.contour(Y,Z,psi,50,colors='k')
CS = plt.contour(Yplot,Zplot,psiplot,30,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Yplot,Zplot,Uplot,np.arange(-100000,110000,10000),cmap='seismic')
#plt.contourf(Y,Z,psi,cmap='seismic')
#jk = 20
#q = plt.quiver(Y[::jk,::jk],Z[::jk,::jk],V[::jk,::jk],W[::jk,::jk],scale=50,angles="xy")
#plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
#plt.streamplot(Y,Z,V,W,density = 3,arrowstyle='->',arrowsize = 1.5)
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
#plt.xlim([-5000,5000])
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'Ustar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


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



#Testing whether V integral is constant with height
Vintegral = []
Hmax = 400
for k in range(0,len(V)):
    if Zplot[k][0] > Hmax:
        Vsum = 0
        for t in range(0,len(V[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Vsum = Vsum + Vplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Vintegral.append(Vsum)
        #print (Z[k][0])
print ("Max and min values of Vintegral", np.amax( Vintegral), np.amin( Vintegral))

jato=10
for k in range(0,len(Vintegral)):
    print(Vintegral[k].real,Hmax+jato,"m" )
    jato = jato + 10
        
        
#Testing whether W integral is constant with height
Wintegral = []
for k in range(0,len(W)):
    if Zplot[k][0] > Hmax:
        Wsum = 0
        for t in range(0,len(Wplot[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Wsum = Wsum + Wplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Wintegral.append(Wsum)
        #print (Z[k][0])
print ("Max and min values of Wintegral", np.amax( Wintegral), np.amin( Wintegral))

jato=10
for k in range(0,len(Wintegral)):
    print(Wintegral[k].real,Hmax+jato,"m" )
    jato = jato + 10      
        
        
        
#Testing whether Ustar integral is constant with height
Ustsum = 0
for t in range(0,len(y)):
    if (y[t] >= float(-L/2) and y[t] <= float(L/2)):
        Ustsum = Ustsum + Ustarplot[t]*abs(y[0]-y[1])
print ("Ustar integral value", Ustsum)

#Getting the U integral (not constant with height)
Uintegral = []
integralZ = []
for k in range(0,len(U)):
    if Zplot[k][0] > Hmax:
        Usum = 0
        for t in range(0,len(U[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Usum = Usum + Uplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Uintegral.append(Usum)
        integralZ.append(z[k])
        #print (Z[k][0])
print ("Max and min values of Uintegral", np.amax( Uintegral), np.amin( Uintegral))


#Getting the B integral
Bintegral = []
for k in range(0,len(B)):
    if Zplot[k][0] > Hmax:
        Bsum = 0
        for t in range(0,len(B[0])):
            if (Yplot[k][t] >= float(-L/2) and Yplot[k][t] <= float(L/2)):
                Bsum = Bsum + Bplot[k][t]*abs(Yplot[0][0]-Yplot[0][1])
        Bintegral.append(Bsum)
        #print (Z[k][0])
print ("Max and min values of Bintegral", np.amax( Bintegral), np.amin( Bintegral))


#Testing wether C is a constant
Uintegral = np.array(Uintegral)
integralZ = np.array(integralZ)
Bintegral = np.array(Bintegral)
delta = float((4*visc*diff)**(1/4)) / np.sqrt(float(N) * float(mp.sin(alpha)))
Cconst = np.sqrt( np.exp(2 * integralZ/delta) * (Bintegral**2 + float(visc/diff) * float(N)**2 * Uintegral**2 ) )


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
        thetaplot_sign[k][t] = float(theta[k][t].real * abs(Wstar_real[k][t].real)/Wstar_real[k][t].real)
        
#Plotting the theta angle vs height at y=0 
thetaplotf = (thetaplot).T
thetaplotf2 = (thetaplot_sign).T 
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
ax1.plot(thetaplotf2[50][1:],z[1:],linewidth=3,color="r") #Number might depends on how many points y has
plt.xlabel(r"$\rm\theta$ [$\rm^{o}$]")
plt.ylabel("Z [m]")
plt.xlim([-0.5,6])
ax1.set_ylim([0,1500])
ax1.set_xticks(np.arange(-0.5,6.5,0.5))
ax1.set_yticks(np.arange(0,1600,100))
ax1.tick_params('both', length=10, width=1, which='major')
plt.grid(True)


#Plotting the Wstar_real 

interval = 10
boundsl = []
for k in range(-interval,interval+1):
    if k < 0:
        boundsl.append(0.5/interval*k)
    if k > 0:
        boundsl.append(7/interval*k)
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
plt.title(r'W$^{\bigstar}_\max$ = 5 m $\rms^{-1}$        W$^{\bigstar}_\min$ = 5 m $\rms^{-1}$',x=0.5, y=1.02)
plt.rcParams.update({'font.size':16})

norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# pcm = plt.pcolormesh(Yplot,Zplot,Yploy,norm=norm,cmap='seismic')

pcm = plt.contourf(Yplot,Zplot,Wstarplot,norm=norm,cmap='seismic',levels=bounds)

cbar = fig.colorbar(pcm)
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

nameoffigure = 'Wstar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


#%%

import pickle
import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

dps_value = 100

mp.dps = dps_value


K = 70
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(5)     
N = mpf(0.01)    
L = mpf(1000)

subdivisions = 100

pizao = mp.pi(dps=dps_value)

def H(y):
    return ( mpf(250) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )


def Bsfc(y):
    return mpf(0.1) * mp.sin(mpf(2) * pizao * y/L)

z = np.arange(0,2010,10) 
y = np.arange(-float(L),float(L)+10,10) 
Y,Z = np.meshgrid(y,z)
Y = Y * mpf(1)
Z = Z * mpf(1)


with open('/home/owner/Documents/katabatic_flows/variables/Ek.pickle', 'wb') as handle:
    pickle.dump(Ek, handle, protocol=pickle.HIGHEST_PROTOCOL)


    


with open('/home/owner/Documents/katabatic_flows/variables/Eki.pickle', 'wb') as handle:
    pickle.dump(Eki, handle, protocol=pickle.HIGHEST_PROTOCOL)


    
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Ck.pickle', 'wb') as handle:
    pickle.dump(Ck, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Cki.pickle', 'wb') as handle:
    pickle.dump(Cki, handle, protocol=pickle.HIGHEST_PROTOCOL)





with open('/home/owner/Documents/katabatic_flows/variables/Dk.pickle', 'wb') as handle:
    pickle.dump(Dk, handle, protocol=pickle.HIGHEST_PROTOCOL)


    

with open('/home/owner/Documents/katabatic_flows/variables/Dki.pickle', 'wb') as handle:
    pickle.dump(Dki, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    
       

with open('/home/owner/Documents/katabatic_flows/variables/B.pickle', 'wb') as handle:
    pickle.dump(B, handle, protocol=pickle.HIGHEST_PROTOCOL)

 
    
 
with open('/home/owner/Documents/katabatic_flows/variables/V.pickle', 'wb') as handle:
    pickle.dump(V, handle, protocol=pickle.HIGHEST_PROTOCOL)




with open('/home/owner/Documents/katabatic_flows/variables/VlargeY.pickle', 'wb') as handle:
    pickle.dump(V, handle, protocol=pickle.HIGHEST_PROTOCOL)



with open('/home/owner/Documents/katabatic_flows/variables/U.pickle', 'wb') as handle:
    pickle.dump(U, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eq.pickle', 'wb') as handle:
    pickle.dump(Eq, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eqi.pickle', 'wb') as handle:
    pickle.dump(Eqi, handle, protocol=pickle.HIGHEST_PROTOCOL)

 
with open('/home/owner/Documents/katabatic_flows/variables/W.pickle', 'wb') as handle:
    pickle.dump(W, handle, protocol=pickle.HIGHEST_PROTOCOL)


    
with open('/home/owner/Documents/katabatic_flows/variables/WlargeY.pickle', 'wb') as handle:
    pickle.dump(W, handle, protocol=pickle.HIGHEST_PROTOCOL)



with open('/home/owner/Documents/katabatic_flows/variables/P.pickle', 'wb') as handle:
    pickle.dump(P, handle, protocol=pickle.HIGHEST_PROTOCOL)

    

with open('/home/owner/Documents/katabatic_flows/variables/psi.pickle', 'wb') as handle:
    pickle.dump(psi, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
####################################################################################################################################

with open('/home/owner/Documents/katabatic_flows/variables/Ek.pickle', 'rb') as handle:
    Ek = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Eki.pickle', 'rb') as handle:
    Eki = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Ck.pickle', 'rb') as handle:
    Ck = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Cki.pickle', 'rb') as handle:
    Cki = pickle.load(handle)   
    

with open('/home/owner/Documents/katabatic_flows/variables/Dk.pickle', 'rb') as handle:
    Dk = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Dki.pickle', 'rb') as handle:
    Dki = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/B.pickle', 'rb') as handle:
    B = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/V.pickle', 'rb') as handle:
    V = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/VlargeY.pickle', 'rb') as handle:
    V = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/U.pickle', 'rb') as handle:
    U = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/Eq.pickle', 'rb') as handle:
    Eq = pickle.load(handle)
    
    
with open('/home/owner/Documents/katabatic_flows/variables/Eqi.pickle', 'rb') as handle:
    Eqi = pickle.load(handle)
 
    
with open('/home/owner/Documents/katabatic_flows/variables/W.pickle', 'rb') as handle:
    W = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/WlargeY.pickle', 'rb') as handle:
    W = pickle.load(handle)


with open('/home/owner/Documents/katabatic_flows/variables/P.pickle', 'rb') as handle:
    P = pickle.load(handle)
    

with open('/home/owner/Documents/katabatic_flows/variables/psi.pickle', 'rb') as handle:
    psi = pickle.load(handle)


#%%
#Calculating some equations
zmp = z * mpf(1)
du2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):   
        du2dz2[k,t] = ( U[k+1,t] - mpf(2)*U[k,t] + U[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

term1 = -B[1:-1,:]*mp.sin(alpha)
term2 = visc*du2dz2[1:-1,:]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '666666666666666666666666666666666666666666')
    
########################################################################################################################################################################

zmp = z * mpf(1)
dv2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dv2dz2[k,t] = ( V[k+1,t] - mpf(2)*V[k,t] + V[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

ymp = y * mpf(1)
dpdy=np.ones_like(Y)*[mpf(0)]
for k in range(1,len(y)-1):
    for t in range(0,len(z)):
        dpdy[t,k] = ( P[t,k+1] - P[t,k-1] ) / abs(ymp[k-1]-ymp[k+1])
    
term1 = -dpdy[1:-1,1:-1]
term2 = visc*dv2dz2[1:-1,1:-1]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '77777777777777777777777777777777777777')   
             
########################################################################################################################################################################        

zmp = z * mpf(1)        
dpdz=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dpdz[k,t] = ( P[k+1,t] - P[k-1,t] ) / abs(zmp[k-1]-zmp[k+1])

term1 = -dpdz[1:-1,:] 
term2 = B[1:-1,:]*mp.cos(alpha)

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '99999999999999999999999999999')       
        
#############################################################################################################################################################        

zmp = z * mpf(1)         
db2dz2=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        db2dz2[k,t] = ( B[k+1,t] - mpf(2)*B[k,t] + B[k-1,t] ) / abs(zmp[k]-zmp[k+1])**2

term1 = N**2 * (U[1:-1,:]*mp.sin(alpha) - W[1:-1,:]*mp.cos(alpha))
term2 = diff*db2dz2[1:-1,:]

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum)

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '8888888888888888888888888888888888888' )       


#############################################################################################################################################################       
        
ymp = y * mpf(1)       
dvdy=np.ones_like(Y)*[mpf(0)]
for k in range(1,len(y)-1):
    for t in range(0,len(z)):
        dvdy[t,k] = ( V[t,k+1] - V[t,k-1] ) / abs(ymp[k-1]-ymp[k+1])

zmp = z * mpf(1)    
dwdz=np.ones_like(Z)*[mpf(0)]
for k in range(1,len(z)-1):
    for t in range(0,len(y)):
        dwdz[k,t] = ( W[k+1,t] - W[k-1,t] ) / abs(zmp[k-1]-zmp[k+1])


term1 = dvdy[1:-1,1:-1]
term2 = dwdz[1:-1,1:-1]  

sumterms = term1 + term2
final_sum = []
for k in range(0,len(sumterms)):
    for t in range(0,len(sumterms[0])):
        final_sum.append(float(sumterms[k][t].real) + 1j*float(sumterms[k][t].imag))
final_sum = np.array(final_sum) 

print (np.nanmax(final_sum) , np.nanmin(final_sum) , '55555555555555555555555555555555555555555555555555') 
        
        
