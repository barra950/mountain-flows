import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np

#Solving the equations for the Prandtl case
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
    
#%%    
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
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('Buoyancy')
plt.contourf(Yplot,Zplot,Bplot,np.arange(-0.1,0.11,0.01),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
plt.colorbar(label='1/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'buoyancy.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()      


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
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('V Wind')
plt.contourf(Yplot,Zplot,Vplot,np.arange(-7,7.03,0.03),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'Vwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()  

#%%
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
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('U Wind')
plt.contourf(Yplot,Zplot,Uplot,np.arange(-25,26,1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-float(L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'Uwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()

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
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('W Wind')
plt.contourf(Yplot,Zplot,Wplot,np.arange(-6,6.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])
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
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('Pressure')
plt.contourf(Yplot,Zplot,Pplot,np.arange(-15,15.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar(label='hPa')
plt.xlabel("Y axis")
plt.ylabel("Height")
#plt.xlim([-10000,10000])
plt.xlim([float(-L),float(L)])
plt.ylim([0,1500])
nameoffigure = 'pressure.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
plt.show()
plt.close()


#Calculating the U* and plotting it  
Ustar = np.ones_like(y)*[mpf(0)]
ympf = y*[mpf(1)]
for k in range(-K,K+1):
    for t in range(0,len(y)):
        Ustar[t] = Ustar[t] +  Eq[Eqi.index(k)] * mp.cos(mpf(2) * mpf(k) * pizao * ympf[t] / L)/mp.cos(alpha)
        
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
for k in range(0,len(theta)):
    for t in range(0,len(theta[0])):
        thetaplot[k][t] = float(theta[k][t].real) #+ 1j*float(theta[k][t].imag)
        
#Plotting the theta angle vs height at y=0 
thetaplotf = (thetaplot).T
fig,ax1 = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
ax1.plot(thetaplotf[50][1:],z[1:],linewidth=3,color="r")
plt.xlabel(r"$\rm\theta$ ($\rm^{o}$)")
plt.ylabel("Height (m)")
plt.xlim([0,6])
ax1.set_ylim([0,1500])
ax1.set_xticks(np.arange(0,6.5,0.5))
ax1.set_yticks(np.arange(0,1600,100))
ax1.tick_params('both', length=10, width=1, which='major')


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
        
        
        
        
        
        
        
        