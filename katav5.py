import pprint
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import numpy as np
from scipy.linalg import lu , lu_factor, lu_solve
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def lu(A):
    
    #Get the number of rows
    n = A.shape[0]
    
    U = A.copy()
    L = np.eye(n, dtype=np.longdouble)
    
    #Loop over rows
    for i in range(n):
            
        #Eliminate entries below i with row operations 
        #on U and reverse the row operations to 
        #manipulate L
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]
        
    return L, U


def forward_substitution(L, b):
    
    #Get number of rows
    n = L.shape[0]
    
    #Allocating space for the solution vector
    y = np.zeros_like(b, dtype=np.longdouble);
    
    #Here we perform the forward-substitution.  
    #Initializing  with the first row.
    y[0] = b[0] / L[0, 0]
    
    #Looping over rows in reverse (from the bottom  up),
    #starting with the second to last row, because  the 
    #last row solve was completed in the last step.
    for i in range(1, n):
        y[i] = (b[i] - np.dot(L[i,:i], y[:i])) / L[i,i]
        
    return y



def back_substitution(U, y):
    
    #Number of rows
    n = U.shape[0]
    
    #Allocating space for the solution vector
    x = np.zeros_like(y, dtype=np.longdouble);

    #Here we perform the back-substitution.  
    #Initializing with the last row.
    x[-1] = y[-1] / U[-1, -1]
    
    #Looping over rows in reverse (from the bottom up), 
    #starting with the second to last row, because the 
    #last row solve was completed in the last step.
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]
        
    return x


def lu_solve(A, b):
    
    L, U = lu(A)
    
    y = forward_substitution(L, b)
    
    return back_substitution(U, y)

#Solving the equations for the Prandtl case

K = 70
alpha = np.longdouble(0.1) 
visc = np.longdouble(5)     
diff = np.longdouble(5)     
N = np.longdouble(0.01)    
L = np.longdouble(1000)


subdivisions = 100
tick = 10
points = np.arange(0,L/2+tick,tick)

def H(y):
    return np.longdouble(( 200 * (1 + np.cos(2 * np.longdouble(np.pi) * y/L)) ))
    #return 0
#    return 700 * 2 * abs(y) / L

def Bsfc(y):
    #return ( 0.1 * (1 + np.cos(2 * np.pi * y/L)) )
    return np.longdouble(0.1)
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
    
        R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

        Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
        
        S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
        S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
        
        phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
        Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
        
        m1 = np.longdouble(- np.sqrt(S1 + S2))
        m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
        m3 = np.clongdouble(m2.conjugate())
        
        
        def f1r(y):
            return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f1i(y):
            return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma1 = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
        def f2r(y):
            return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f2i(y):
            return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma2 = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
        
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
                equation2.append(np.longdouble(0))
                equation2.append(np.longdouble(0))
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
        return np.longdouble((Bsfc(y) * np.cos(2 * q * np.longdouble(np.pi) * y / L) ).real)
    def f4i(y):
        return np.clongdouble((Bsfc(y) * np.cos(2 * q * np.longdouble(np.pi) * y / L) ).imag)
    b.append(np.clongdouble(2/L * (quad(f4r,0,L/2,limit=subdivisions)[0] + quad(f4i,0,L/2,limit=subdivisions)[0]*1j)))
    
    
    if q != 0:
        final_system.append(equation2)
        b.append(np.longdouble(0))
    
    
    final_system.append(equation3)
    b.append(np.longdouble(0))


final_system = np.array(final_system)
b=np.array(b)


#Normal solver 
# X = np.linalg.solve(final_system,b)
    
# print (np.allclose(final_system @ X, b), '-----------------')

#LU solver 1
#LU, P = lu_factor(final_system)
#X = lu_solve((LU,P),b) 
#
#print (np.allclose(final_system @ X, b))

#LU solver 2
# P, Ls, U = scipy.linalg.lu(final_system)

# Bl = np.linalg.inv(P) @ b 

# Z = np.linalg.solve(Ls,Bl)

# X = np.linalg.solve(U,Z)

# print (np.allclose(final_system @ X, b))

#LU solver 3
X = lu_solve(final_system, b)

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
magic_value = 50
z = np.arange(0,2002,2) 
y = np.arange(-L,L+10,10) 
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
 
#Buoyancy for the prandtl case   
Bp = Bsfc(Y) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.cos(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )


    
##Plotting the buoyancy
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('Buoyancy',y=1.05)
plt.contourf(Y,Z,B,np.arange(-0.1,0.105,0.005),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
plt.colorbar()
plt.xlabel("Y [m]")
plt.ylabel("Z [m]")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'buoyancy.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
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
        if abs(Z[k][t] - H(Y[k][t])) < 0.1:
            if V[k][t] > 0.1:
                print (V[k][t],'fudeu geral -------------------------------------------------')
#            print (V[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '-----------------------------------------------------------------------------' )


##Plotting the V wind
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('V Wind')
plt.contourf(Y,Z,V,np.arange(-7,7.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
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
    E = 0
    for k in range(-K,K+1):
    
        R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

        Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
        
        S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
        S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
        
        phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
        Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
        
        m1 = np.longdouble(- np.sqrt(S1 + S2))
        m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
        m3 = np.clongdouble(m2.conjugate())
        
        
        def f1r(y):
            return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f1i(y):
            return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma1 = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
        def f2r(y):
            return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).real)
        def f2i(y):
            return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.pi * y / L) ).imag)
        gamma2 = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
        
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

#U for prandtl case:
Up = -Bsfc(Y)/N * np.sqrt(diff/visc) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.sin(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )


#Plotting the U wind
fig = plt.figure(figsize=(10,10)) # create a figure
plt.rcParams.update({'font.size':16})
plt.title('U Wind')
plt.contourf(Y,Z,U,np.arange(-25,26,1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Uwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


#Now, we need to get the values of the W wind 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
W = np.ones_like(Y)*[0]


for k in range(-K,K+1):
    
    R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

    Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
    
    S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
    S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
    
    phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
    Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
    
    m1 = np.longdouble(- np.sqrt(S1 + S2))
    m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
    m3 = np.clongdouble(m2.conjugate())
    
    if k != 0:
        W = W + np.cos(alpha)/visc * 4*(k)**2 *np.longdouble(np.pi)**2 / L**2 * (Ak[Aki.index(k)]*np.exp(m1*Z)/(m1**4) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**4)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**4)   ) * np.exp(2j * (k) * np.longdouble(np.pi) * Y / L) +  np.tan(alpha)*(Eq[Eqi.index(k)] * np.cos(2 * (k) * np.longdouble(np.pi) * Y / L) )
    else:
        W = W + np.tan(alpha)*(Eq[Eqi.index(k)] * np.cos(2 * (k) * np.longdouble(np.pi) * Y / L) )
        


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
plt.contourf(Y,Z,W,np.arange(-6,6.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Wwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



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
plt.contourf(Y,Z,P,np.arange(-15,15.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar(label='hPa')
plt.xlabel("Y axis")
plt.ylabel("Height")
#plt.xlim([-10000,10000])
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'pressure.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
plt.show()
plt.close()



#Calculating the streamlines 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
psi = np.ones_like(Y)*[0]


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
        psi = psi - np.tan(alpha) * 1j*L* Eq[Eqi.index(k)] / (2*k*np.pi) * np.exp(2j * (k) * np.pi * Y / L) + Eq[Eqi.index(0)] * Y * np.tan(alpha) - np.cos(alpha)/visc * 2j*k*np.pi/L *  ( Ak[Aki.index(k)]*np.exp(m1*Z)/(m1**4) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**4)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**4) ) * np.exp(2j * (k) * np.pi * Y / L)
    

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


V=V.real
W=W.real

##Plotting the streamlines
fig = plt.figure(figsize=(20,20)) 
plt.rcParams.update({'font.size':16})
plt.rcParams['contour.negative_linestyle'] = 'solid'
#plt.title('Streamfunction')
fig.add_subplot(1,1,1)
#plt.contourf(Y,Z,psi,np.arange(-300,305,5),cmap='seismic')
CS = plt.contour(Y,Z,psi,50,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Y,Z,U,np.arange(-100000,110000,10000),cmap='seismic')
#plt.contourf(Y,Z,psi,cmap='seismic')
jk = 20
#q = plt.quiver(Y[::jk,::jk],Z[::jk,::jk],V[::jk,::jk],W[::jk,::jk],scale=50,angles="xy")
#plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
#plt.streamplot(Y,Z,V,W,density = 3,arrowstyle='->',arrowsize = 1.5)
plt.xlabel("Y axis")
plt.ylabel("Height")
#plt.xlim([-5000,5000])
plt.xlim([-L,L])
plt.ylim([0,1500])
plt.subplots_adjust(bottom=0.07, top=0.99, hspace=0.1,right=0.99,left=0.05)
nameoffigure = 'streamlines.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()

#Calculating the U* and plotting it  
Ustar = np.ones_like(y)*[0]
for k in range(-K,K+1):
    Ustar = Ustar +  Eq[Eqi.index(k)] * np.cos(2 * (k) * np.pi * y / L)/np.cos(alpha)


##Plotting U* and streamlines
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


fig.add_subplot(2,1,1)
#plt.title('U star plot')
plt.plot(y,Ustar)
plt.xlabel("Y axis $(m)$")
plt.ylabel("U$^{\u2605}_\infty$ ($ms^{-1}$)")
#plt.xlim([-5000,5000])
plt.xlim([-L,L])
#plt.ylim([-6,2])
plt.ylim([-18,10])
plt.grid('True')

fig.add_subplot(2,1,2)
#CS = plt.contour(Y,Z,psi,50,colors='k')
CS = plt.contour(Y,Z,psi,30,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Y,Z,U,np.arange(-100000,110000,10000),cmap='seismic')
#plt.contourf(Y,Z,psi,cmap='seismic')
jk = 20
#q = plt.quiver(Y[::jk,::jk],Z[::jk,::jk],V[::jk,::jk],W[::jk,::jk],scale=50,angles="xy")
#plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
#plt.streamplot(Y,Z,V,W,density = 3,arrowstyle='->',arrowsize = 1.5)
plt.xlabel("Y axis $(m)$")
plt.ylabel("Height $(m)$")
#plt.xlim([-5000,5000])
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Ustar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


#Testing whether V integral is constant with height
Vintegral = []
Hmax = 400
for k in range(0,len(V)):
    if Z[k][0] > Hmax:
        Vsum = 0
        for t in range(0,len(V[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Vsum = Vsum + V[k][t]*abs(Y[0][0]-Y[0][1])
        Vintegral.append(Vsum)
        #print (Z[k][0])
print ("Max and min values of Vintegral", np.amax( Vintegral), np.amin( Vintegral))

jato=10
for k in range(0,len(Vintegral)):
    print(Vintegral[k].real,400+jato,"m" )
    jato = jato + 10
    
    
#Testing whether W integral is constant with height
Wintegral = []
for k in range(0,len(W)):
    if Z[k][0] > Hmax:
        Wsum = 0
        for t in range(0,len(W[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Wsum = Wsum + W[k][t]*abs(Y[0][0]-Y[0][1])
        Wintegral.append(Wsum)
        #print (Z[k][0])
print ("Max and min values of Wintegral", np.amax( Wintegral), np.amin( Wintegral))

jato=10
for k in range(0,len(Wintegral)):
    print(Wintegral[k].real,400+jato,"m" )
    jato = jato + 10
    

#Testing whether Ustar integral is constant with height
Ustsum = 0
for t in range(0,len(y)):
    if (y[t] >= -L/2 and y[t] <= L/2):
        Ustsum = Ustsum + Ustar[t]*abs(y[0]-y[1])
print ("Ustar integral value", Ustsum)


#Getting the U integral
Uintegral = []
integralZ = []
for k in range(0,len(U)):
    if Z[k][0] > Hmax:
        Usum = 0
        for t in range(0,len(U[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Usum = Usum + U[k][t]*abs(Y[0][0]-Y[0][1])
        Uintegral.append(Usum)
        integralZ.append(z[k])
        #print (Z[k][0])
print ("Max and min values of Uintegral", np.amax( Uintegral), np.amin( Uintegral))

#Getting the B integral
Bintegral = []
for k in range(0,len(B)):
    if Z[k][0] > Hmax:
        Bsum = 0
        for t in range(0,len(B[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Bsum = Bsum + B[k][t]*abs(Y[0][0]-Y[0][1])
        Bintegral.append(Bsum)
        #print (Z[k][0])
print ("Max and min values of Bintegral", np.amax( Bintegral), np.amin( Bintegral))

jato=10
for k in range(0,len(Vintegral)):
    print(Vintegral[k].real,400+jato,"m" )
    jato = jato + 10

#Testing wether C is a constant
Uintegral = np.array(Uintegral)
integralZ = np.array(integralZ)
Bintegral = np.array(Bintegral)
delta = (4*visc*diff)**(1/4) / np.sqrt(N * np.sin(alpha))
Cconst = np.sqrt( np.exp(2 * integralZ/delta) * (Bintegral**2 + (visc/diff) * N**2 * Uintegral**2 ) )

#Plotting the prandtl solution
Bpplot = np.array(Bp).T
Upplot = np.array(Up).T
fig,ax1=plt.subplots()
#plt.xticks(maxtime, time2plot, rotation='vertical')
plt.rcParams.update({'font.size':20})
fsize = 20

ax1.set_ylim([-0.02,0.1])
plt.xlabel('Z [m]',name='Arial',size=fsize)
plt.ylabel(r'B [$\rmms^{-2}$]',name='Arial',size=fsize)
plt.plot(z,Bpplot[0][:],linewidth=3,color='b')
ax1.tick_params('both', length=10, width=1, which='major')
ax2=ax1.twinx()
ax2.set_ylim([-4,4])
plt.plot(z,Upplot[0][:],linewidth=3,color='r')
plt.ylabel(r'U [$\rmms^{-1}$]',name='Arial',size=fsize)
ax1.set_xlim([0,1000])
ax1.set_xticks(np.arange(0,1100,100))
ax2.tick_params('both', length=10, width=1, which='major')

#Plotting the theta angle
Ustar_real = U * np.cos(alpha) + W * np.sin(alpha)
Wstar_real = - U * np.sin(alpha) + W * np.cos(alpha)
theta = np.arccos( abs(Ustar_real / np.sqrt(U**2 + W**2 )))
fraction = Ustar_real / np.sqrt(U**2 + W**2 )
theta = theta * 180 / np.pi


fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('Theta')
plt.contourf(Y,Z,theta,np.arange(0,101,1),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='Degrees')
plt.xlabel("Y axis")
plt.ylabel("Height")
#plt.xlim([-L,L])
#plt.ylim([0,1500])
nameoffigure = 'theta.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()

#Plotting the Ustar_real
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
plt.title('Ustar')
plt.contourf(Y,Z,Ustar_real,np.arange(-25,26,1),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
#plt.xlim([-L,L])
#plt.ylim([0,1500])


#Plotting the theta angle vs height at y=0
thetaplot = np.array(theta).T
fig,ax = plt.subplots(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})
ax.plot(thetaplot[50][1:],z[1:])
plt.xlabel("Theta")
plt.ylabel("Height (m)")
#plt.xlim([0,5])
ax.set_ylim([0,1500])
ax.set_xticks(np.arange(0,5.5,0.5))
ax.set_yticks(np.arange(0,1600,100))
ax.tick_params('both', length=10, width=1, which='major')


ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_tick_params(which='minor', bottom=False)











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

def lu(A):
    
    #Get the number of rows
    n = A.shape[0]
    
    U = A.copy()
    L = np.eye(n, dtype=np.longdouble)
    
    #Loop over rows
    for i in range(n):
            
        #Eliminate entries below i with row operations 
        #on U and reverse the row operations to 
        #manipulate L
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]
        
    return L, U


def forward_substitution(L, b):
    
    #Get number of rows
    n = L.shape[0]
    
    #Allocating space for the solution vector
    y = np.zeros_like(b, dtype=np.longdouble);
    
    #Here we perform the forward-substitution.  
    #Initializing  with the first row.
    y[0] = b[0] / L[0, 0]
    
    #Looping over rows in reverse (from the bottom  up),
    #starting with the second to last row, because  the 
    #last row solve was completed in the last step.
    for i in range(1, n):
        y[i] = (b[i] - np.dot(L[i,:i], y[:i])) / L[i,i]
        
    return y



def back_substitution(U, y):
    
    #Number of rows
    n = U.shape[0]
    
    #Allocating space for the solution vector
    x = np.zeros_like(y, dtype=np.longdouble);

    #Here we perform the back-substitution.  
    #Initializing with the last row.
    x[-1] = y[-1] / U[-1, -1]
    
    #Looping over rows in reverse (from the bottom up), 
    #starting with the second to last row, because the 
    #last row solve was completed in the last step.
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]
        
    return x


def lu_solve(A, b):
    
    L, U = lu(A)
    
    y = forward_substitution(L, b)
    
    return back_substitution(U, y)

#Solving the equations for the Prandtl case

K = 70
alpha = np.longdouble(0.1) 
visc = np.longdouble(5)     
diff = np.longdouble(5)     
N = np.longdouble(0.01)    
L = np.longdouble(5000)


subdivisions = 100 

tick = 10
points = np.arange(0,L/2+tick,tick)

def H(y):
    return np.longdouble(( 200 * (1 + np.cos(2 * np.longdouble(np.pi) * y/L)) ))   
    #return 0
#    return 700 * 2 * abs(y) / L

def Bsfc(y):
    return np.longdouble(0.1) * np.sin(2 * np.longdouble(np.pi) * y/L)
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
    
        R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

        Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
        
        S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
        S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
        
        phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
        Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
        
        m1 = np.longdouble(- np.sqrt(S1 + S2))
        m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
        m3 = np.clongdouble(m2.conjugate())
        
        
        def f1r(y):
            return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f1i(y):
            return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma1 = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
        def f2r(y):
            return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f2i(y):
            return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma2 = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
        
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
                equation2.append(np.longdouble(0))
                equation2.append(np.longdouble(0))
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
        return np.longdouble((Bsfc(y) * np.sin(2 * q * np.pi * y / L) ).real)
    def f4i(y):
        return np.clongdouble((Bsfc(y) * np.sin(2 * q * np.pi * y / L) ).imag)
    b.append(np.clongdouble(-2/L * (quad(f4r,0,L/2,limit=subdivisions)[0] + quad(f4i,0,L/2,limit=subdivisions)[0]*1j)))
    
    
    if q != 0:
        final_system.append(equation2)
        b.append(np.longdouble(0))
    
    
    final_system.append(equation3)
    b.append(np.longdouble(0))


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
# P, Ls, U = scipy.linalg.lu(final_system)

# Bl = np.linalg.inv(P) @ b 

# Z = np.linalg.solve(Ls,Bl)

# X = np.linalg.solve(U,Z)

# print (np.allclose(final_system @ X, b))

#LU solver 3
X = lu_solve(final_system, b)

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
z = np.arange(0,2010,10) 
y = np.arange(-L,L+10,10) 
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
plt.contourf(Y,Z,B,np.arange(-0.2,0.21,0.01),cmap='seismic')
#plt.contourf(Y,Z,B,cmap='seismic')
plt.colorbar(label='1/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])  
nameoffigure = 'buoyancy.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


#Getting the value of the V wind
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
V = np.ones_like(Y)*[0]


for k in range(-K,K+1):
    
    R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

    Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
        
    S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
    S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
        
    phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
    Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
        
    m1 = np.longdouble(- np.sqrt(S1 + S2))
    m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
    m3 = np.clongdouble(m2.conjugate())
    
    def f1r(y):
        return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f1i(y):
        return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma1k = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
    def f2r(y):
        return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f2i(y):
        return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma2k = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
    
    
    
    if k != 0:
        V = V + np.cos(alpha)/visc * 2j*k*np.longdouble(np.pi)/L *  ( 1j * Ek[Eki.index(k)]*np.exp(m1*Z)/(m1**3) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**3)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**3) ) * np.exp(2j * (k) * np.longdouble(np.pi) * Y / L)
        V = V + np.cos(alpha) / visc * 2 * k * np.longdouble(np.pi) / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + 2 * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) 

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
plt.contourf(Y,Z,V,np.arange(-7,7.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,V,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
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
    E = 0
    for k in range(-K,K+1):
    
        R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

        Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
            
        S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
        S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
            
        phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
        Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
            
        m1 = np.longdouble(- np.sqrt(S1 + S2))
        m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
        m3 = np.clongdouble(m2.conjugate())
        
        
        def f1r(y):
            return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f1i(y):
            return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma1 = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
        def f2r(y):
            return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).real)
        def f2i(y):
            return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (q - k) * np.longdouble(np.pi) * y / L) ).imag)
        gamma2 = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
        
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
plt.contourf(Y,Z,U,np.arange(-25,26,1),cmap='seismic')
#plt.contourf(Y,Z,U,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Uwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



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
plt.contourf(Y,Z,W,np.arange(-6,6.2,0.2),cmap='seismic')
#plt.contourf(Y,Z,W,cmap='seismic')
plt.colorbar(label='m/s')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Wwind.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()





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
plt.contourf(Y,Z,P,np.arange(-15,15.5,0.5),cmap='seismic')
#plt.contourf(Y,Z,P,cmap='seismic')
plt.colorbar(label='hPa')
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'pressure.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()




#Calculating the streamlines 
#z = np.arange(0,2010,10) 
#y = np.arange(-5000,5050,50) 
Y,Z = np.meshgrid(y,z)
psi = np.ones_like(Y)*[0]


for k in range(-K,K+1):
    
    R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

    Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
        
    S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
    S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
        
    phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
    Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
        
    m1 = np.longdouble(- np.sqrt(S1 + S2))
    m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
    m3 = np.clongdouble(m2.conjugate())
    
    
    def f1r(y):
        return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f1i(y):
        return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma1k = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
    def f2r(y):
        return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f2i(y):
        return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma2k = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
    
    if k != 0:
        psi = psi - np.tan(alpha) * 1j*L* Eq[Eqi.index(k)] / (2*k*np.pi) * np.exp(2j * (k) * np.pi * Y / L)  - np.cos(alpha)/visc * 2j*k*np.pi/L *  ( 1j*Ek[Eki.index(k)]*np.exp(m1*Z)/(m1**4) + Ck[Cki.index(k)]*np.exp(m2*Z)/(m2**4)  + Dk[Dki.index(k)]*np.exp(m3*Z)/(m3**4) ) * np.exp(2j * (k) * np.pi * Y / L)
        psi = psi - ( np.cos(alpha) / visc * 2 * k * np.pi / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + 2 * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) ) * Z
psi = psi + Eq[Eqi.index(0)] * Y * np.tan(alpha)

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



V=V.real
W=W.real

##Plotting the streamlines
fig = plt.figure(figsize=(20,20)) 
plt.rcParams.update({'font.size':16})
plt.rcParams['contour.negative_linestyle'] = 'solid'
#plt.title('Streamfunction')
fig.add_subplot(1,1,1)
#plt.contourf(Y,Z,psi,np.arange(-300,305,5),cmap='seismic')
CS = plt.contour(Y,Z,psi,50,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Y,Z,U,np.arange(-100000,110000,10000),cmap='seismic')
#plt.contourf(Y,Z,psi,cmap='seismic')
jk = 20
#q = plt.quiver(Y[::jk,::jk],Z[::jk,::jk],V[::jk,::jk],W[::jk,::jk],scale=50,angles="xy")
#plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
#plt.streamplot(Y,Z,V,W,density = 3,arrowstyle='->',arrowsize = 1.5)
plt.xlabel("Y axis")
plt.ylabel("Height")
plt.xlim([-L,L])
plt.ylim([0,1500])
plt.show()
plt.subplots_adjust(bottom=0.07, top=0.99, hspace=0.1,right=0.99,left=0.05)
nameoffigure = 'streamlines.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()


#Calculating the U* and plotting it  
Ustar = np.ones_like(y)*[0]
for k in range(-K,K+1):
    Ustar = Ustar +  1j*Eq[Eqi.index(k)] * np.sin(2 * (k) * np.pi * y / L)/np.cos(alpha)


##Plotting U* and streamlines
fig = plt.figure(figsize=(10,10)) 
plt.rcParams.update({'font.size':16})


fig.add_subplot(2,1,1)
#plt.title('U star plot')
plt.plot(y,Ustar)
plt.xlabel("Y axis $(m)$")
plt.ylabel("U$^{\u2605}_\infty$ ($ms^{-1}$)")
#plt.xlim([-5000,5000])
plt.xlim([-L,L])
#plt.ylim([-6,2])
plt.ylim([-18,10])
plt.grid('True')

fig.add_subplot(2,1,2)
#CS = plt.contour(Y,Z,psi,50,colors='k')
CS = plt.contour(Y,Z,psi,30,colors='k')
#plt.clabel(CS, fontsize=9, inline=True)
#plt.colorbar(label='m/s')
plt.contourf(Y,Z,U,np.arange(-100000,110000,10000),cmap='seismic')
#plt.contourf(Y,Z,psi,cmap='seismic')
#jk = 20
#q = plt.quiver(Y[::jk,::jk],Z[::jk,::jk],V[::jk,::jk],W[::jk,::jk],scale=50,angles="xy")
#plt.quiverkey(q, 1.03, 1.03, 2, label='2m/s')
#plt.streamplot(Y,Z,V,W,density = 3,arrowstyle='->',arrowsize = 1.5)
plt.xlabel("Y axis $(m)$")
plt.ylabel("Height $(m)$")
#plt.xlim([-5000,5000])
plt.xlim([-L,L])
plt.ylim([0,1500])
nameoffigure = 'Ustar.png'
string_in_string = "{}".format(nameoffigure)
plt.savefig('/home/owner/Documents/katabatic_flows/output/'+string_in_string)
#plt.show()
plt.close()



#Getting Vinfinity
Vinf=0
for k in range(-K,K+1):
    
    R = np.longdouble(2 * N**2 * np.cos(alpha)**2 / (visc * diff) * (k * np.longdouble(np.pi) / L)**2)

    Q = np.longdouble(N**2 * np.sin(alpha)**2 / (3 * visc * diff))
    
    S1 = np.longdouble(abs(R + np.sqrt(Q**3 + R**2) )**(1/3))
    S2 = np.longdouble(- abs( np.sqrt(Q**3 + R**2) -R )**(1/3))
    
    phi = np.longdouble(np.sqrt(S1**2 + S2**2 - S1*S2))
    Lk = np.longdouble(np.arccos(- (S1 + S2)/ (2 * phi) ))
    
    m1 = np.longdouble(- np.sqrt(S1 + S2))
    m2 = np.clongdouble(- np.sqrt(phi) * np.exp(1j * Lk/2))
    m3 = np.clongdouble(m2.conjugate())
    
    
    def f1r(y):
        return np.longdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f1i(y):
        return np.clongdouble((np.exp(m1 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma1k = np.clongdouble(2/L * (quad(f1r,0,L/2,limit=subdivisions)[0] + quad(f1i,0,L/2,limit=subdivisions)[0]*1j))
        
    def f2r(y):
        return np.longdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).real)
    def f2i(y):
        return np.clongdouble((np.exp(m2 * H(y)) * np.cos(2 * (- k) * np.longdouble(np.pi) * y / L) ).imag)
    gamma2k = np.clongdouble(2/L * (quad(f2r,0,L/2,limit=subdivisions)[0] + quad(f2i,0,L/2,limit=subdivisions)[0]*1j))
    
    if k != 0:
        Vinf = Vinf + ( np.cos(alpha) / visc * 2 * k * np.pi / L * ( Ek[Eki.index(k)] * gamma1k / m1**3 + 2 * (Ck[Cki.index(k)] * gamma2k / m2**3).imag  ) )


#Testing whether V integral is constant with height
Vintegral = []
Hmax = 400
for k in range(0,len(V)):
    if Z[k][0] > Hmax:
        Vsum = 0
        for t in range(0,len(V[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Vsum = Vsum + V[k][t]*abs(Y[0][0]-Y[0][1])
        Vintegral.append(Vsum)
        #print (Z[k][0])
print ("Max and min values of Vintegral", np.amax( Vintegral), np.amin( Vintegral))

jato=10
for k in range(0,len(Vintegral)):
    print(Vintegral[k].real,400+jato,"m" )
    jato = jato + 10
    
    
Wintegral = []
for k in range(0,len(W)):
    if Z[k][0] > Hmax:
        Wsum = 0
        for t in range(0,len(W[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Wsum = Wsum + W[k][t]*abs(Y[0][0]-Y[0][1])
        Wintegral.append(Wsum)
        #print (Z[k][0])
print ("Max and min values of Wintegral", np.amax( Wintegral), np.amin( Wintegral))

jato=10
for k in range(0,len(Wintegral)):
    print(Wintegral[k].real,400+jato,"m" )
    jato = jato + 10
    
#Testing whether Ustar integral is constant with height
Ustsum = 0
for t in range(0,len(y)):
    if (y[t] >= -L/2 and y[t] <= L/2):
        Ustsum = Ustsum + Ustar[t]*abs(y[0]-y[1])
print ("Ustar integral value", Ustsum)


#Getting the U integral
Uintegral = []
integralZ = []
for k in range(0,len(U)):
    if Z[k][0] > Hmax:
        Usum = 0
        for t in range(0,len(U[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Usum = Usum + U[k][t]*abs(Y[0][0]-Y[0][1])
        Uintegral.append(Usum)
        integralZ.append(z[k])
        #print (Z[k][0])
print ("Max and min values of Uintegral", np.amax( Uintegral), np.amin( Uintegral))

#Getting the B integral
Bintegral = []
for k in range(0,len(B)):
    if Z[k][0] > Hmax:
        Bsum = 0
        for t in range(0,len(B[0])):
            if (Y[k][t] >= -L/2 and Y[k][t] <= L/2):
                Bsum = Bsum + B[k][t]*abs(Y[0][0]-Y[0][1])
        Bintegral.append(Bsum)
        #print (Z[k][0])
print ("Max and min values of Bintegral", np.amax( Bintegral), np.amin( Bintegral))


#Testing wether C is a constant
Uintegral = np.array(Uintegral)
integralZ = np.array(integralZ)
Bintegral = np.array(Bintegral)
delta = (4*visc*diff)**(1/4) / np.sqrt(N * np.sin(alpha))
Cconst = np.sqrt( np.exp(2 * integralZ/delta) * (Bintegral**2 + (visc/diff) * N**2 * Uintegral**2 ) )

 


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

    
    





