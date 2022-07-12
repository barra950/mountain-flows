import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import numpy as np
from scipy.linalg import lu , lu_factor, lu_solve
from scipy.integrate import quad
import matplotlib.pyplot as plt

#Solving the equations for the Prandtl case

K = 20
alpha = 0.1  
visc = 5  
diff = 5   
N = 0.01     
L = 5000 

tick = 10
points = np.arange(0,L/2+tick,tick)

def H(y):
    return ( 200 * (1 + np.cos(2 * np.pi * y/L)) )
    #return 0
#    return 700 * 2 * abs(y) / L

def Bsfc(y):
    #return ( 0.1 * (1 + np.cos(2 * np.pi * y/L)) )
    return 0.1
    #return 0.2 * 2 * abs(y) / L

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
    
        R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

        Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
        S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
        S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
        phi = np.sqrt(S1**2 + S2**2 - S1*S2)
        Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
        m1 = - np.sqrt(S1 + S2)
        m2 = - np.sqrt(phi) * np.exp(1j * Lk/2)
        m3 = m2.conjugate()
        
        
        def f1r(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f1i(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma1 = 2/L * (quad(f1r,0,L/2)[0] + quad(f1i,0,L/2)[0]*1j)
        
        def f2r(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f2i(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma2 = 2/L * (quad(f2r,0,L/2)[0] + quad(f2i,0,L/2)[0]*1j)
        
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
        
         
        if k == 0:
            equation1.append(2 * gamma2.real)
            Cki.append(k)
            equation1.append(-2 * gamma2.imag)
            Dki.append(k)
        else:
            equation1.append(gamma1)
            Aki.append(k)
            equation1.append(2 * gamma2.real)
            Cki.append(k)
            equation1.append(-2 * gamma2.imag)
            Dki.append(k)
        
        if q != 0:
            
            if k == 0:
                equation2.append(0)
                equation2.append(0)
            else:
                equation2.append(k * gamma1 / (m1**3) )
                equation2.append(2 * k * (gamma2 / (m2**3) ).real)
                equation2.append(-2 * k * (gamma2 / (m2**3) ).imag)
        
        if k == 0:
            equation3.append(2 * (m2**2 * gamma2).real)
            equation3.append(-2 * (m2**2 * gamma2).imag)
        else:
            equation3.append(m1**2 * gamma1)
            equation3.append(2 * (m2**2 * gamma2).real)
            equation3.append(-2 * (m2**2 * gamma2).imag)
    
    final_system.append(equation1)
    def f4r(y):
        return (Bsfc(y) * np.cos(2 * q * np.pi * y / L) ).real
    def f4i(y):
        return (Bsfc(y) * np.cos(2 * q * np.pi * y / L) ).imag
    b.append(2/L * (quad(f4r,0,L/2)[0] + quad(f4i,0,L/2)[0]*1j))
    
    
    if q != 0:
        final_system.append(equation2)
        b.append(0)
    
    
    final_system.append(equation3)
    b.append(0)


final_system = np.array(final_system)
b=np.array(b)


#Normal solver 
#X = np.linalg.solve(final_system,b)
#    
#print (np.allclose(final_system @ X, b), '-----------------')

#LU solver 1
#LU, P = lu_factor(final_system)
#X = lu_solve((LU,P),b) 
#
#print (np.allclose(final_system @ X, b))

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
        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
            if B[k][t] > 0.101:
                print (B[k][t],'fudeu geral -------------------------------------------------')
#            print (B[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )
    
Bp = Bsfc(Y) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.cos(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )


    
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
        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
            if V[k][t] > 0.1:
                print (V[k][t],'fudeu geral -------------------------------------------------')
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
plt.show()




#Getting the value of the U wind
#We first need the value of Eq
Eq=[]
Eqi=[]
for q in range(-K,K+1):
    E = 0
    for k in range(-K,K+1):
    
        R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

        Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
        S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
        S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
        phi = np.sqrt(S1**2 + S2**2 - S1*S2)
        Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
        m1 = - np.sqrt(S1 + S2)
        m2 = - np.sqrt(phi) * np.exp(1j * Lk/2)
        m3 = m2.conjugate()
        
        
        def f1r(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f1i(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma1 = 2/L * (quad(f1r,0,L/2)[0] + quad(f1i,0,L/2)[0]*1j)
        
        def f2r(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f2i(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma2 = 2/L * (quad(f2r,0,L/2)[0] + quad(f2i,0,L/2)[0]*1j)
        
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
            E = E - np.sin(alpha)/visc * Ak[Aki.index(k)]*gamma1/(m1**2) 
        E = E - 2*np.sin(alpha)/visc * ( Ck[Cki.index(k)]*gamma2/(m2**2) ).real
    
    Eq.append(E)
    Eqi.append(q)

Eq= np.array(Eq)




#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
U = np.ones_like(Y)*[0]


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
        U = U + np.sin(alpha)/visc * Ak[Aki.index(k)]/(m1**2) * np.exp(m1*Z) * np.exp(2j * (k) * np.pi * Y / L)
    U = U + np.sin(alpha)/visc * ( Ck[Cki.index(k)]/(m2**2)*np.exp(m2*Z) + Dk[Dki.index(k)]/(m3**2)*np.exp(m3*Z) ) * np.exp(2j * (k) * np.pi * Y / L) + Eq[Eqi.index(k)] * np.cos(2 * (k) * np.pi * Y / L)


for k in range(0,U.shape[0]):
    for t in range(0,U.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            U[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (U[k][t], "U value at ground")
        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
            if U[k][t] > 0.1:
                print (U[k][t],'fudeu geral -------------------------------------------------')
#            print (U[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


#Plotting the U wind
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('U Wind')
plt.contourf(Y,Z,U,np.arange(-10,10.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])
plt.show()



#Now, we need to get the values of the W wind 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
W = np.ones_like(Y)*[0]


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
        W = W + np.cos(alpha)/visc * 4*(k)**2 *np.pi**2 / L**2 * (Ak[Aki.index(k)]*np.exp(m1*Z)/(m1**4) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**4)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**4)   ) * np.exp(2j * (k) * np.pi * Y / L) +  np.tan(alpha)*(Eq[Eqi.index(k)] * np.cos(2 * (k) * np.pi * Y / L) )
    else:
        W = W + np.tan(alpha)*(Eq[Eqi.index(k)] * np.cos(2 * (k) * np.pi * Y / L) )
        


for k in range(0,W.shape[0]):
    for t in range(0,W.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            W[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (W[k][t], "W value at ground")
        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
            if W[k][t] > 0.1:
                print (W[k][t],'fudeu geral -------------------------------------------------')
#            print (W[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )

           
##Plotting the W wind
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('W Wind')
plt.contourf(Y,Z,W,np.arange(-4,4.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])
plt.show()



#Now, we need to get the values of the pressure 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
P = np.ones_like(Y)*[0]


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
        P = P + np.cos(alpha) * Ak[Aki.index(k)] / m1 * np.exp(m1*Z) * np.exp(2j * (k) * np.pi * Y / L)
    P = P + np.cos(alpha) * (Ck[Cki.index(k)] / m2 * np.exp(m2*Z) + Dk[Dki.index(k)] / m3 * np.exp(m3*Z)) * np.exp(2j * (k) * np.pi * Y / L)


for k in range(0,P.shape[0]):
    for t in range(0,P.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            P[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (P[k][t], "P value at ground")
        
            

##Plotting the pressure
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('Pressure')
plt.contourf(Y,Z,P,np.arange(-4,4.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar(label='hPa')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
plt.ylim([0,1000])
plt.show()


#%%
#Calculating some equations

du2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    du2dz2[k,:] = ( U[k+1,:] - 2*U[k,:] + U[k-1,:] ) / abs(z[k]-z[k+1])**2

term1 = -B[1:-1,:]*np.sin(alpha)
term2 = visc*du2dz2[1:-1,:]

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '666666666666666666666666666666666666666666')
    
##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

########################################################################################################################################################################


dv2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dv2dz2[k,:] = ( V[k+1,:] - 2*V[k,:] + V[k-1,:] ) / abs(z[k]-z[k+1])**2

dpdy=np.ones_like(Y)*np.nan
for k in range(1,len(y)-1):
    dpdy[:,k] = ( P[:,k+1] - P[:,k-1] ) / abs(y[k-1]-y[k+1])
    
term1 = -dpdy[1:-1,1:-1]
term2 = visc*dv2dz2[1:-1,1:-1]

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '77777777777777777777777777777777777777')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


#############################################################################################################################################################

dpdz=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dpdz[k,:] = ( P[k+1,:] - P[k-1,:] ) / abs(z[k-1]-z[k+1])

term1 = -dpdz[1:-1,:] 
term2 = B[1:-1,:]*np.cos(alpha)

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '99999999999999999999999999999')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])



#############################################################################################################################################################

db2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    db2dz2[k,:] = ( B[k+1,:] - 2*B[k,:] + B[k-1,:] ) / abs(z[k]-z[k+1])**2

term1 = N**2 * (U[1:-1,:]*np.sin(alpha) - W[1:-1,:]*np.cos(alpha))
term2 = diff*db2dz2[1:-1,:]


print (np.nanmax(term1 + term2),np.nanmin(term1 + term2) , '8888888888888888888888888888888888888' )

    
##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

#############################################################################################################################################################

dvdy=np.ones_like(Y)*np.nan
for k in range(1,len(y)-1):
    dvdy[:,k] = ( V[:,k+1] - V[:,k-1] ) / abs(y[k-1]-y[k+1])
    
dwdz=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dwdz[k,:] = ( W[k+1,:] - W[k-1,:] ) / abs(z[k-1]-z[k+1])


term1 = dvdy[1:-1,1:-1]
term2 = dwdz[1:-1,1:-1]   

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '55555555555555555555555555555555555555555555555555')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

#############################################################################################################################################################

#%%
#Now solving for when Bsfc is an odd function of Y
import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import numpy as np
from scipy.linalg import lu , lu_factor, lu_solve
from scipy.integrate import quad
import matplotlib.pyplot as plt

#Solving the equations for the Prandtl case

K = 20
alpha = 0.1   
visc = 5  
diff = 5   
N = 0.01     
L = 5000 

tick = 10
points = np.arange(0,L/2+tick,tick)

def H(y):
    return ( 200 * (1 + np.cos(2 * np.pi * y/L)) )
    #return 0
#    return 700 * 2 * abs(y) / L

def Bsfc(y):
    return 0.1 * np.sin(2 * np.pi * y/L)
    #return 0.1
    #return 0.2 * 2 * abs(y) / L

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
    
        R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

        Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
        S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
        S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
        phi = np.sqrt(S1**2 + S2**2 - S1*S2)
        Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
        m1 = - np.sqrt(S1 + S2)
        m2 = - np.sqrt(phi) * np.exp(1j * Lk/2)
        m3 = m2.conjugate()
        
        
        def f1r(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f1i(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma1 = 2/L * (quad(f1r,0,L/2)[0] + quad(f1i,0,L/2)[0]*1j)
        
        def f2r(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f2i(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma2 = 2/L * (quad(f2r,0,L/2)[0] + quad(f2i,0,L/2)[0]*1j)
        
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
        
         
        if k == 0:
            equation1.append(2 * gamma2.imag)
            Cki.append(k)
            equation1.append(2 * gamma2.real)
            Dki.append(k)
        else:
            equation1.append(gamma1)
            Eki.append(k)
            equation1.append(2 * gamma2.imag)
            Cki.append(k)
            equation1.append(2 * gamma2.real)
            Dki.append(k)
        
        if q != 0:
            
            if k == 0:
                equation2.append(0)
                equation2.append(0)
            else:
                equation2.append(k * gamma1 / (m1**3) )
                equation2.append(2 * k * (gamma2 / (m2**3) ).imag)
                equation2.append(2 * k * (gamma2 / (m2**3) ).real)
        
        if k == 0:
            equation3.append(2 * (m2**2 * gamma2).imag)
            equation3.append(2 * (m2**2 * gamma2).real)
        else:
            equation3.append(m1**2 * gamma1)
            equation3.append(2 * (m2**2 * gamma2).imag)
            equation3.append(2 * (m2**2 * gamma2).real)
    
    final_system.append(equation1)
    def f4r(y):
        return (Bsfc(y) * np.sin(2 * q * np.pi * y / L) ).real
    def f4i(y):
        return (Bsfc(y) * np.sin(2 * q * np.pi * y / L) ).imag
    b.append(-2/L * (quad(f4r,0,L/2)[0] + quad(f4i,0,L/2)[0]*1j))
    
    
    if q != 0:
        final_system.append(equation2)
        b.append(0)
    
    
    final_system.append(equation3)
    b.append(0)


final_system = np.array(final_system)
b=np.array(b)


#Normal solver 
#X = np.linalg.solve(final_system,b)
#    
#print (np.allclose(final_system @ X, b), '-----------------')

#LU solver 1
#LU, P = lu_factor(final_system)
#X = lu_solve((LU,P),b) 
#
#print (np.allclose(final_system @ X, b))

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
        strings.append('E')     
            
    strings.append('R')           
    strings.append('I')             
        


Ek = []
Rk = []
Ik = []
for k in range(0,len(X)):
     if 'E' in strings[k]:
         Ek.append(X[k])
         
     if 'R' in strings[k]:
         Rk.append(X[k])
        
     if 'I' in strings[k]:
         Ik.append(X[k])
         
Ck=[]
      
for k in range(0,len(Rk)):
    Ck.append(Rk[k] + Ik[k] * 1j)

Ck = np.array(Ck)

Dk = - Ck.conjugate()

Ek = np.array(Ek)


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
        B = B + ( 1j*Ek[Eki.index(k)] * np.exp(m1 * Z) * np.exp(2j * (k) * np.pi * Y / L)  )
    B = B + ( ( Ck[Cki.index(k)] * np.exp(m2 * Z) + Dk[Dki.index(k)] * np.exp(m3 * Z) )  * np.exp(2j * (k) * np.pi * Y / L) )
    

for k in range(0,B.shape[0]):
    for t in range(0,B.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            B[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (B[k][t], "B value at the ground")
#       if abs(Z[k][t] - H(Y[k][t])) < 0.1:
#            if B[k][t] > 0.101:
#                print (B[k][t],'fudeu geral -------------------------------------------------')
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
    
    def f1r(y):
        return (np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.pi * y / L) ).real
    def f1i(y):
        return (np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.pi * y / L) ).imag
    gamma1k = 2/L * (quad(f1r,0,L/2)[0] + quad(f1i,0,L/2)[0]*1j)
        
    def f2r(y):
        return (np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.pi * y / L) ).real
    def f2i(y):
        return (np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.pi * y / L) ).imag
    gamma2k = 2/L * (quad(f2r,0,L/2)[0] + quad(f2i,0,L/2)[0]*1j)
    
    
    
    if k != 0:
        V = V + np.cos(alpha)/visc * 2j*k*np.pi/L *  ( 1j * Ek[Eki.index(k)]*np.exp(m1*Z)/(m1**3) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**3)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**3) ) * np.exp(2j * (k) * np.pi * Y / L)
        V = V + np.cos(alpha) / visc * 2 * k * np.pi / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + 2 * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) 

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
plt.show()  


#Getting the value of the U wind
#We first need the value of Eq
Eq=[]
Eqi=[]
for q in range(-K,K+1):
    E = 0
    for k in range(-K,K+1):
    
        R = 2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.pi / L)**2

        Q = N**2 * np.sin(alpha)**2 / (3 * visc * diff)
        
        S1 = abs(R + np.sqrt(Q**3 + R**2) )**(1/3)
        S2 = - abs( np.sqrt(Q**3 + R**2) -R )**(1/3)
        
        phi = np.sqrt(S1**2 + S2**2 - S1*S2)
        Lk = np.arccos(- (S1 + S2)/ (2 * phi) )
        
        m1 = - np.sqrt(S1 + S2)
        m2 = - np.sqrt(phi) * np.exp(1j * Lk/2)
        m3 = m2.conjugate()
        
        
        def f1r(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f1i(y):
            return (np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma1 = 2/L * (quad(f1r,0,L/2)[0] + quad(f1i,0,L/2)[0]*1j)
        
        def f2r(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real
        def f2i(y):
            return (np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag
        gamma2 = 2/L * (quad(f2r,0,L/2)[0] + quad(f2i,0,L/2)[0]*1j)
        
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
            E = E - 1j * np.sin(alpha)/visc * ( Ek[Eki.index(k)]*gamma1/(m1**2) + 2 * ( Ck[Cki.index(k)] * gamma2 / (m2**2) ).imag )
        
    
    Eq.append(E)
    Eqi.append(q)

Eq= np.array(Eq)





#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
U = np.ones_like(Y)*[0]


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
        U = U + np.sin(alpha)/visc * 1j*Ek[Eki.index(k)]/(m1**2) * np.exp(m1*Z) * np.exp(2j * (k) * np.pi * Y / L)
    U = U + np.sin(alpha)/visc * ( Ck[Cki.index(k)]/(m2**2)*np.exp(m2*Z) + Dk[Dki.index(k)]/(m3**2)*np.exp(m3*Z) ) * np.exp(2j * (k) * np.pi * Y / L) + 1j *Eq[Eqi.index(k)] * np.sin(2 * (k) * np.pi * Y / L)


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


#Plotting the U wind
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('U Wind')
plt.contourf(Y,Z,U,np.arange(-10,10.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])
plt.show()



#Now, we need to get the values of the W wind 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
W = np.ones_like(Y)*[0]


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
        W = W + np.cos(alpha)/visc * 4*(k)**2 *np.pi**2 / L**2 * (1j*Ek[Eki.index(k)]*np.exp(m1*Z)/(m1**4) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**4)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**4)   ) * np.exp(2j * (k) * np.pi * Y / L) +  np.tan(alpha)*(1j*Eq[Eqi.index(k)] * np.sin(2 * (k) * np.pi * Y / L) )
    else:
        W = W + np.tan(alpha)*(1j*Eq[Eqi.index(k)] * np.sin(2 * (k) * np.pi * Y / L) )
        


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

           
##Plotting the W wind
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('W Wind')
plt.contourf(Y,Z,W,np.arange(-4,4.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
#plt.ylim([1000,10000])
plt.show()





#Now, we need to get the values of the pressure 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
P = np.ones_like(Y)*[0]


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
        P = P + np.cos(alpha) * 1j *Ek[Eki.index(k)] / m1 * np.exp(m1*Z) * np.exp(2j * (k) * np.pi * Y / L)
    P = P + np.cos(alpha) * (Ck[Cki.index(k)] / m2 * np.exp(m2*Z) + Dk[Dki.index(k)] / m3 * np.exp(m3*Z)) * np.exp(2j * (k) * np.pi * Y / L)


for k in range(0,P.shape[0]):
    for t in range(0,P.shape[1]):
        if Z[k][t] < H(Y[k][t]):
            P[k][t] = np.nan
        if Z[k][t] == H(Y[k][t]):
            print (P[k][t], "P value at ground")
        
            

##Plotting the pressure
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('Pressure')
plt.contourf(Y,Z,P,np.arange(-4,4.1,0.1),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar(label='hPa')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-10000,10000])
plt.ylim([0,1000])
plt.show()


#%%
#Calculating some equations

du2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    du2dz2[k,:] = ( U[k+1,:] - 2*U[k,:] + U[k-1,:] ) / abs(z[k]-z[k+1])**2

term1 = -B[1:-1,:]*np.sin(alpha)
term2 = visc*du2dz2[1:-1,:]

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '666666666666666666666666666666666666666666')
    
##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

########################################################################################################################################################################


dv2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dv2dz2[k,:] = ( V[k+1,:] - 2*V[k,:] + V[k-1,:] ) / abs(z[k]-z[k+1])**2

dpdy=np.ones_like(Y)*np.nan
for k in range(1,len(y)-1):
    dpdy[:,k] = ( P[:,k+1] - P[:,k-1] ) / abs(y[k-1]-y[k+1])
    
term1 = -dpdy[1:-1,1:-1]
term2 = visc*dv2dz2[1:-1,1:-1]

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '77777777777777777777777777777777777777')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


#############################################################################################################################################################

dpdz=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dpdz[k,:] = ( P[k+1,:] - P[k-1,:] ) / abs(z[k-1]-z[k+1])

term1 = -dpdz[1:-1,:] 
term2 = B[1:-1,:]*np.cos(alpha)

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '99999999999999999999999999999')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

plt.subplot(2,1,2)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])



#############################################################################################################################################################

db2dz2=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    db2dz2[k,:] = ( B[k+1,:] - 2*B[k,:] + B[k-1,:] ) / abs(z[k]-z[k+1])**2

term1 = N**2 * (U[1:-1,:]*np.sin(alpha) - W[1:-1,:]*np.cos(alpha))
term2 = diff*db2dz2[1:-1,:]


print (np.nanmax(term1 + term2),np.nanmin(term1 + term2) , '8888888888888888888888888888888888888' )

    
##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,:],Z[1:-1,:],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.contourf(Y[1:-1,:],Z[1:-1,:],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

#############################################################################################################################################################

dvdy=np.ones_like(Y)*np.nan
for k in range(1,len(y)-1):
    dvdy[:,k] = ( V[:,k+1] - V[:,k-1] ) / abs(y[k-1]-y[k+1])
    
dwdz=np.ones_like(Z)*np.nan
for k in range(1,len(z)-1):
    dwdz[k,:] = ( W[k+1,:] - W[k-1,:] ) / abs(z[k-1]-z[k+1])


term1 = dvdy[1:-1,1:-1]
term2 = dwdz[1:-1,1:-1]   

print (np.nanmax(term1 + term2) , np.nanmin(term1 + term2) , '55555555555555555555555555555555555555555555555555')

##Plotting the individual terms vs sum
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})

plt.subplot(2,1,1)
plt.title('Terms vs sum')
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term1/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term1/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])


plt.subplot(2,1,2)
plt.contourf(Y[1:-1,1:-1],Z[1:-1,1:-1],term2/(abs(term1) + abs(term2)),np.arange(-1,1.05,0.05),cmap='seismic')
plt.colorbar(label='term2/(term1+term2)')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-5000,5000])
plt.ylim([0,1000])

#############################################################################################################################################################

    




