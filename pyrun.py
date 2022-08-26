import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import numpy as np
from scipy.linalg import lu , lu_factor, lu_solve
from scipy.integrate import quad
import matplotlib.pyplot as plt

#Solving the equations for the Prandtl case

K = 100
alpha = 0.1 
visc = 5 
diff = 5   
N = 0.01     
L = 5000 

def H(y):
    return ( 300 * (1 + np.cos(2 * np.pi * y/L)) )

for q in range(-K,K+1):
    Aki = []
    Cki = []
    Dki = []
    for k in range(-K,K+1):
    
         
        if k == 0:
            Cki.append(k)
            Dki.append(k)
        else:
            Aki.append(k)
            Cki.append(k)
            Dki.append(k)
        
    


import datetime
import re

#copy the ihop sounding into the ihopsounding file
f = open('chupar')
r = f.readlines()


for k in range(0,len(r)):
    r[k] = re.sub('\s+', '', r[k])
    r[k] = r[k].split('oo')

b=[]
final_system = []
for k in r:
    if len(k) > 1:
        equations = []
        for t in k:
            if len(t) > 2:
                equations.append(float(t))
        final_system.append(equations)
    else:
        b.append(float(k[0]))
            
final_system = np.array(final_system)
b=np.array(b)

#LU solver 2
P, Ls, U = scipy.linalg.lu(final_system)

Bl = np.linalg.inv(P) @ b 

Z = np.linalg.solve(Ls,Bl)

X = np.linalg.solve(U,Z)

print (np.allclose(final_system @ X, b))


#Getting the values for Ak, Ck and Dk
strings = []

for k in range(-K,K+1):
    if k != 0:
        strings.append('A')     
            
    strings.append('R')           
    strings.append('I')             
        


Ak = []
Rk = []
Ik = []
for k in range(0,len(X)):
     if 'A' in strings[k]:
         Ak.append(X[k])
         
     if 'R' in strings[k]:
         Rk.append(X[k])
        
     if 'I' in strings[k]:
         Ik.append(X[k])
         
Ck=[]
      
for k in range(0,len(Rk)):
    Ck.append(Rk[k] + Ik[k] * 1j)

Ck = np.array(Ck)

Dk = Ck.conjugate()

Ak = np.array(Ak)


#Getting the Buoyancy value
z = np.arange(0,2002,2) 
y = np.arange(-5000,5010,10) 
Y,Z = np.meshgrid(y,z)
B = np.ones_like(Y)*[0]

for k in range(-K,K+1):
    
    R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

    Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
    S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
    S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
    phi = np.sqrt(S1**2 + S2**2 - S1*S2)
    Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
    m1 = - np.sqrt(S1 + S2)
    m2 = -np.sqrt(phi) * np.exp(1j * Lk/2)
    m3 = m2.conjugate()
    
    if k != 0:
        B = B + ( Ak[Aki.index(k)] * np.exp(m1 * Z) * np.exp(2j * (k) * np.pi * Y / L)  )
    B = B + ( ( Ck[Cki.index(k)] * np.exp(m2 * Z) + Dk[Dki.index(k)] * np.exp(m3 * Z) )  * np.exp(2j * (k) * np.pi * Y / L) )
    

for k in range(0,B.shape[0]):
    for t in range(0,B.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            B[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (B[k][t], "B value at the ground")
#        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if B[k][t] > 0.101:
#               print (B[k][t],'fudeu geral -------------------------------------------------')
#            print (B[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )
    


    
##Plotting the buoyancy
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('Buoyancy')
plt.contourf(Y,Z,B,np.arange(-0.2,0.201,0.001),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
plt.colorbar(label='1/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])  
#plt.show()

plt.savefig('/home/owner/rato.png')
plt.close()

#Getting the value of the V wind
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
V = np.ones_like(Y)*[0]


for k in range(-K,K+1):
    
    R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

    Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
    S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
    S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
    phi = np.sqrt(S1**2 + S2**2 - S1*S2)
    Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
    m1 = - np.sqrt(S1 + S2)
    m2 = -np.sqrt(phi) * np.exp(1j * Lk/2)
    m3 = m2.conjugate()
    
    if k != 0:
        V = V + np.cos(alpha)/visc * 2j*k*np.pi/L *  ( Ak[Aki.index(k)]*np.exp(m1*Z)/(m1**3) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**3)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**3) ) * np.exp(2j * (k) * np.pi * Y / L)
    

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


##Plotting the V wind
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('V Wind')
plt.contourf(Y,Z,V,np.arange(-7,7.05,0.05),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])
#plt.show()


plt.savefig('/home/owner/rato2.png')
plt.close()

print(np.nanmax(V.imag))
