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


K = 100
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(5)     
N = mpf(0.01)    
L = mpf(5000)

subdivisions = 100

pizao = mp.pi(dps=dps_value)

# def H(y):
#     return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )


def H(y):
    return ( mpf(300) * (mpf(1) + mp.cos(mpf(2) * pizao * (y-L/mpf(2))/L)) )


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


with open('/home/gome0004/kataprogramms/oddvariables4/Ek.pickle', 'wb') as handle:
    pickle.dump(Ek, handle, protocol=pickle.HIGHEST_PROTOCOL)


    


with open('/home/gome0004/kataprogramms/oddvariables4/Eki.pickle', 'wb') as handle:
    pickle.dump(Eki, handle, protocol=pickle.HIGHEST_PROTOCOL)


    
    
    
with open('/home/gome0004/kataprogramms/oddvariables4/Ck.pickle', 'wb') as handle:
    pickle.dump(Ck, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
    
    
with open('/home/gome0004/kataprogramms/oddvariables4/Cki.pickle', 'wb') as handle:
    pickle.dump(Cki, handle, protocol=pickle.HIGHEST_PROTOCOL)





with open('/home/gome0004/kataprogramms/oddvariables4/Dk.pickle', 'wb') as handle:
    pickle.dump(Dk, handle, protocol=pickle.HIGHEST_PROTOCOL)


    

with open('/home/gome0004/kataprogramms/oddvariables4/Dki.pickle', 'wb') as handle:
    pickle.dump(Dki, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
#Getting the Buoyancy value

z = np.arange(0,2010,1) 
y = np.arange(-float(L),float(L)+5,5) 
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
            
            

with open('/home/gome0004/kataprogramms/oddvariables4/B.pickle', 'wb') as handle:
    pickle.dump(B, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

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
#            print (V[k][t], Z[k][t], H(Y
    
    
with open('/home/gome0004/kataprogramms/oddvariables4/V.pickle', 'wb') as handle:
    pickle.dump(V, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    


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
##            print (U[k][t], Z[k][t], H(Y[k][t]), Y[k][t], '----------------------------            



with open('/home/gome0004/kataprogramms/oddvariables4/U.pickle', 'wb') as handle:
    pickle.dump(U, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/gome0004/kataprogramms/oddvariables4/Eq.pickle', 'wb') as handle:
    pickle.dump(Eq, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
with open('/home/gome0004/kataprogramms/oddvariables4/Eqi.pickle', 'wb') as handle:
    pickle.dump(Eqi, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    

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
##            print (W[k][t], Z[k][t], H(






with open('/home/gome0004/kataprogramms/oddvariables4/W.pickle', 'wb') as handle:
    pickle.dump(W, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    
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
            
            
            

with open('/home/gome0004/kataprogramms/oddvariables4/P.pickle', 'wb') as handle:
    pickle.dump(P, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    
    
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
##            print (psi[k][t], Z[k][t], H(Y[k][t]




with open('/home/gome0004/kataprogramms/oddvariables4/psi.pickle', 'wb') as handle:
    pickle.dump(psi, handle, protocol=pickle.HIGHEST_PROTOCOL)
