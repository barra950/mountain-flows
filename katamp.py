import mpmath
from mpmath import *
from mpmath import mp
import matplotlib.pyplot as plt

#Solving the equations for the Prandtl case
dps_value = 100

mp.dps = dps_value


K = 50
alpha = mpf(0.1) 
visc = mpf(5)     
diff = mpf(5)     
N = mpf(0.01)    
L = mpf(5000)

pizao = mp.pi(dps=dps_value)

def H(y):
    return ( mpf(400) * (mpf(1) + mp.cos(mpf(2) * pizao * y/L)) )

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

solution = []
for k in X:
    solution.append(k)



#Getting the values for Ak, Ck and Dk
import numpy as np
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
L = np.longdouble(5000)


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
plt.contourf(Y,Z,B,np.arange(-0.1,0.12,0.02),cmap='seismic')
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

    
             
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
