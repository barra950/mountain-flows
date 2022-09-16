import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np

#Solving the equations for the Prandtl case
dps_value = 100

mp.dps = dps_value


K = 50
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

#Now, we start converting out of mpf to float type
for k in range(0,len(Ak)):
    Ak[k] = float(Ak[k].real) + float(Ak[k].imag)*1j 
            
for k in range(0,len(Ck)):
    Ck[k] = float(Ck[k].real) + float(Ck[k].imag)*1j
    
for k in range(0,len(Dk)):
    Dk[k] = float(Dk[k].real) + float(Dk[k].imag)*1j
            
  
Xpost = []
for k in range(0,len(solution[0])):
    Xpost.append(solution[0][k])

Xpost = np.array(Xpost)

for k in range(0,len(Xpost)):
    Xpost[k] = float(Xpost[k].real) + float(Xpost[k].imag)*1j
    
    
#Getting the Buoyancy value



alpha = np.longdouble(0.1) 
visc = np.longdouble(5)     
diff = np.longdouble(5)     
N = np.longdouble(0.01)    
L = np.longdouble(1000)

def H(y):
    return np.longdouble(( 300 * (1 + np.cos(2 * np.longdouble(np.pi) * y/L)) ))
    #return 0
#    return 700 * 2 * abs(y) / L

def Bsfc(y):
    #return ( 0.1 * (1 + np.cos(2 * np.pi * y/L)) )
    return np.longdouble(0.1)
    #return 0.2 * 2 * abs(y) / L


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
plt.contourf(Y,Z,B,np.arange(-0.1,0.11,0.01),cmap='seismic')
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
#Up = -Bsfc(Y)/N * np.sqrt(diff/visc) * np.exp(-Z * np.sqrt(N * np.sin(alpha) ) / (4*visc*diff)**(1/4) ) * np.sin(np.sqrt(N*np.sin(alpha)) /((4*visc*diff)**(1/4))*Z )


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

    
             
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
