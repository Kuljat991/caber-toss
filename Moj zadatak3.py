import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import pi , sin, cos, sqrt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.optimize import bisect
from matplotlib.animation import FuncAnimation

"""
                                     VIZUALIZACIJA PROBLEMA

         Faza 1:      Faza 2:                Faza 3:                           Faza 4:
          zalet       izbacaj                  Let                          Prizemljenje
   |                 |       |                                          |                     |
   |                 |       |                                          |                     |
   |                 |       |                                          |                     |
   |                 |       |                                          |                     |
   |                 |       |                                          |                     |

                                    B                               
           B                 B    /O                                
          /O               /O    / /                               
         / /              / /   / /                                  A
        / /              / /   / /                                O\
       / /              / /   /O/ S                               \ \
      /O/ S          S /O/   / /                                   \ \
     / /              / /   / /                                     \ \
    / /              / /   / /                                       \O\ S
   / /              / /    O/                                         \ \
   O/               O/        A                                        \ \  B         S        A
     A                A                                                 \ \___________________      
_________________________________________________________________________\O__________O________O____________

"""
m = 50.0
l = 5.0
zalet = 7.0
ax0 = 2.5
g = 9.81
kut = 75.5

# Pocetni uvjeti
x0 = 0
vx0 = 0
y0 = 0.5
vy0 = 0
phi0 = kut*pi/180
vphi0 = 0
#___________________________________________________________________
#ZALET

def F1(f, t):
    x, dx, y, dy, theta, dtheta = f
    d2theta = 3/(2*l)*(ax0*sin(f[4])-g*cos(f[4]))
    d2y     = 0.
    d2x     = ax0
    return [dx,
            d2x,
            dy,
            d2y,
            dtheta,
            d2theta]

tstop1 = sqrt(2*zalet/ax0)
nt = 100.0

T1 = np.linspace(0, tstop1, nt)

rez1 = odeint(F1, [x0, vx0, y0, vy0, phi0, vphi0], T1)

X_A1    = rez1[:,0]
VX_A1   = rez1[:,1]
Y_A1    = rez1[:,2]
VY_A1   = rez1[:,3]
PHI_1  = rez1[:,4]
VPHI_1 = rez1[:,5]

X_B1 = X_A1+l*cos(PHI_1)
Y_B1 = Y_A1+l*sin(PHI_1)

X_S1 = X_A1+l/2*cos(PHI_1)
Y_S1 = Y_A1+l/2*sin(PHI_1)

"""
for i in range(0,np.size(T1),2):   
    plt.plot([X_A1[i],X_B1[i]],[Y_A1[i],Y_B1[i]], 'b-', lw=2)
    plt.plot([X_S1[i]],[Y_S1[i]], 'yo')
    plt.plot([X_A1[i]],[Y_A1[i]], 'go')
    plt.plot([X_B1[i]],[Y_B1[i]], 'ro')
plt.show()
"""

#___________________________________________________________________ 
#IZBACAJ
#koordinatni sustav se prebacuje u tocku S
#Pocetni uvijeti faze izbacaja

A = 6000.
B = 2000.

X_S2_poc = X_A1[-1]+l/2*cos(PHI_1[-1])
Y_S2_poc = Y_A1[-1]+l/2*sin(PHI_1[-1])
VX_S2_poc   = VX_A1[-1] + VPHI_1[-1]**2*l/2*cos(-pi/2+PHI_1[-1])
VY_S2_poc   = VY_A1[-1] + VPHI_1[-1]**2*l/2*sin(-pi/2+PHI_1[-1])
PHI_S2_poc  = PHI_1[-1]
VPHI_S2_poc = VPHI_1[-1]

tstop2 = tstop1 + 3.0
nt = 100.0

T2_tmp = np.linspace(tstop1, tstop2, nt)



def F2(h, t):
    x, dx, y, dy, theta, dtheta = h
    d2theta = ((A/m)*cos(theta)+(l/2)*((dtheta)**2)*sin(theta)*cos(theta))/((l/2*(cos(theta))**2)-2*l/3)
    d2y     = A/m-g-l/2*d2theta*cos(theta)+l/2*dtheta**2*sin(theta)
    d2x     = B/m+l/2*d2theta*sin(theta)+l/2*dtheta**2*sin(theta)
    return [dx,
            d2x,
            dy,
            d2y,
            dtheta,
            d2theta]


rez2 = odeint(F2, [X_S2_poc, VX_S2_poc, Y_S2_poc, VY_S2_poc, PHI_S2_poc, VPHI_S2_poc], T2_tmp)

X_S2_tmp     = rez2[:,0]
VX_S2_tmp    = rez2[:,1]
Y_S2_tmp     = rez2[:,2]
VY_S2_tmp    = rez2[:,3]
PHI_2_tmp   = rez2[:,4]
VPHI_2_tmp  = rez2[:,5]
#tt2 = ?

x2_interp = interp1d(T2_tmp, X_S2_tmp)
def x2_limit(t):
    return x2_interp(t) - (0.5 + X_S1[-1])
tt_2 = fsolve(x2_limit, T1[-1]) 

T2     = np.linspace(T1[-1], tt_2, 50)
X_S2   = np.interp(T2, T2_tmp, X_S2_tmp)
VX_S2  = np.interp(T2, T2_tmp, VX_S2_tmp)
Y_S2   = np.interp(T2, T2_tmp, Y_S2_tmp)
VY_S2  = np.interp(T2, T2_tmp, VY_S2_tmp)
PHI_2  = np.interp(T2, T2_tmp, PHI_2_tmp)
VPHI_2 = np.interp(T2, T2_tmp, VPHI_2_tmp)

X_A2 = X_S2-l/2*cos(PHI_2)
Y_A2 = Y_S2-l/2*sin(PHI_2)

X_B2 = X_S2+l/2*cos(PHI_2)
Y_B2 = Y_S2+l/2*sin(PHI_2)

"""
for i in range(0,np.size(T2),2):   
    plt.plot([X_A2[i],X_B2[i]],[Y_A2[i],Y_B2[i]], 'b-', lw=2)
    plt.plot([X_S2[i]],[Y_S2[i]], 'yo')
    plt.plot([X_A2[i]],[Y_A2[i]], 'go')
    plt.plot([X_B2[i]],[Y_B2[i]], 'ro')
plt.show()
"""

#___________________________________________________________________ 
#Let
#Pocetni uvijeti faze leta

X_S3_poc    = X_S2[-1]
Y_S3_poc    = Y_S2[-1]
VX_S3_poc   = VX_S2[-1]
VY_S3_poc   = VY_S2[-1]
PHI_S3_poc  = PHI_2[-1]
VPHI_S3_poc = VPHI_2[-1]

tstop3 = T2[-1] + 5.0
nt = 100.0

T3_tmp = np.linspace(T2[-1], tstop3, nt)

def F3(h, t):
    x, dx, y, dy, theta, dtheta = h
    d2theta = 0
    d2y     = -g
    d2x     = 0
    return [dx,
            d2x,
            dy,
            d2y,
            dtheta,
            d2theta]


rez3 = odeint(F3, [X_S3_poc, VX_S3_poc, Y_S3_poc, VY_S3_poc, PHI_S3_poc, VPHI_S3_poc], T3_tmp)

X_S3_tmp   = rez3[:,0]
VX_S3_tmp  = rez3[:,1]
Y_S3_tmp   = rez3[:,2]
VY_S3_tmp  = rez3[:,3]
PHI_3_tmp  = rez3[:,4]
VPHI_3_tmp = rez3[:,5]

Y_B3_tmp = Y_S3_tmp + l/2*sin(PHI_3_tmp)

y3b_interp = interp1d(T3_tmp, Y_B3_tmp)

def y3b_limit(t):
    if t < T3_tmp[0]:
        return y3b_interp(T3_tmp[0])
    if t > T3_tmp[-1]:
        return y3b_interp(T3_tmp[-1])
    return y3b_interp(t) 

tt_3 = bisect(y3b_limit, T3_tmp[0],20) 

T3     = np.linspace(T2[-1], tt_3, 50)
X_S3   = np.interp(T3, T3_tmp, X_S3_tmp)
VX_S3  = np.interp(T3, T3_tmp, VX_S3_tmp)
Y_S3   = np.interp(T3, T3_tmp, Y_S3_tmp)
VY_S3  = np.interp(T3, T3_tmp, VY_S3_tmp)
PHI_3  = np.interp(T3, T3_tmp, PHI_3_tmp)
VPHI_3 = np.interp(T3, T3_tmp, VPHI_3_tmp)



"""
T_all = np.append(T1, T2)
T_all = np.append(T_all, T3)
X_all = np.append(X1, X2)
X_all = np.append(X_all, X3)
Y_all = np.append(Y1, Y2)
Y_all = np.append(Y_all, Y3)
PHI_all = np.append(PHI1, PHI2)
PHI_all = np.append(PHI_all, PHI3)
"""

X_A3 = X_S3 - l/2*cos(PHI_3)
Y_A3 = Y_S3 - l/2*sin(PHI_3)

X_B3 = X_S3 + l/2*cos(PHI_3)
Y_B3 = Y_S3 + l/2*sin(PHI_3)

"""
for i in range(0,np.size(T3),2):   
    plt.plot([X_A3[i],X_B3[i]],[Y_A3[i],Y_B3[i]], 'b-', lw=2)
    plt.plot([X_S3[i]],[Y_S3[i]], 'yo')
    plt.plot([X_A3[i]],[Y_A3[i]], 'go')
    plt.plot([X_B3[i]],[Y_B3[i]], 'ro')
    plt.axis('equal')
plt.show()
"""
#___________________________________________________________________ 
#Prizemljenje
#Pocetni uvijeti faze prizemljenja

PHI_4_poc  = PHI_3[-1] + pi
VPHI_4_poc = VPHI_3[-1]

tstop4 = T3[-1] + 5
nt = 100.0

T4_tmp = np.linspace(T3[-1], tstop4, nt)

def F4(h, t):
    theta, dtheta = h
    d2theta = -3*g*cos(theta)/2*l
    return [dtheta,
            d2theta]


rez4 = odeint (F4,[PHI_4_poc, VPHI_4_poc], T4_tmp)

PHI_4_tmp   = rez4[:,0]
VPHI_4_tmp  = rez4[:,1]

noviPHI = np.rad2deg(PHI_4_tmp)

Y_A4_tmp = l*sin(PHI_4_tmp)

y4a_interp = interp1d(T4_tmp,Y_A4_tmp)

for i in T4_tmp:
    if (y4a_interp(i) < 0):
        t_bisect = i;
        break;


tt_4 = bisect(y4a_interp, T4_tmp[0], t_bisect)


T4    = np.linspace(T3[-1], tt_4, 100)
PHI4  = np.interp(T4, T3_tmp, PHI_4_tmp)

T4     = np.linspace(T3[-1], tt_4, 50)
PHI_4  = np.interp(T4, T4_tmp, PHI_4_tmp)
VPHI_4 = np.interp(T4, T4_tmp, VPHI_4_tmp)

X_B4 = np.zeros(50) + X_B3[-1]
Y_B4 = np.zeros(50)

X_A4 = X_B4 + l*cos(PHI_4)
Y_A4 = l*sin(PHI_4)

"""
for i in range(0,np.size(T4),2):   
    plt.plot([X_A4[i],X_B4[i]],[Y_A4[i],Y_B4[i]], 'b-', lw=2)
    plt.plot([X_A4[i]],[Y_A4[i]], 'go')
    plt.plot([X_B4[i]],[Y_B4[i]], 'ro')
    plt.axis('equal')
plt.show()
"""

T_all = T1
T_all = np.append(T_all,T2)
T_all = np.append(T_all,T3)
T_all = np.append(T_all,T4)

X_A_all = X_A1
X_A_all = np.append(X_A_all,X_A2)
X_A_all = np.append(X_A_all,X_A3)
X_A_all = np.append(X_A_all,X_A4)

Y_A_all = Y_A1
Y_A_all = np.append(Y_A_all,Y_A2)
Y_A_all = np.append(Y_A_all,Y_A3)
Y_A_all = np.append(Y_A_all,Y_A4)


X_B_all = X_B1
X_B_all = np.append(X_B_all,X_B2)
X_B_all = np.append(X_B_all,X_B3)
X_B_all = np.append(X_B_all,X_B4)

Y_B_all = Y_B1
Y_B_all = np.append(Y_B_all,Y_B2)
Y_B_all = np.append(Y_B_all,Y_B3)
Y_B_all = np.append(Y_B_all,Y_B4)

"""
for i in range(0,np.size(T_all),2):   
    plt.plot([X_A_all[i],X_B_all[i]],[Y_A_all[i],Y_B_all[i]], 'b-', lw=2)
    plt.plot([X_A_all[i]],[Y_A_all[i]], 'go')
    plt.plot([X_B_all[i]],[Y_B_all[i]], 'ro')
plt.show()
"""

Tfinal = T_all[-1]
T_final = np.linspace(0,Tfinal,100)
X_A_final = np.interp(T_final, T_all, X_A_all)
Y_A_final = np.interp(T_final, T_all, Y_A_all)

X_B_final = np.interp(T_final, T_all, X_B_all)
Y_B_final = np.interp(T_final, T_all, Y_B_all)


#plt.figure()
maxY = np.max(Y_B_final)+2
maxX = np.max([X_B_final[-1],X_A_final[-1]])+2
print (maxY)
for i in range(0,np.size(T_final)):
#    print "Frame", i
#    plt.clf()
    plt.axis('scaled')  
    plt.xlim(0,maxX)
    plt.ylim(0,maxY)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('a = '+ str(ax0) + ", " + "A =" + str(A))
    plt.plot([X_A_final[i],X_B_final[i]],[Y_A_final[i],Y_B_final[i]], 'b-', lw=2)
    plt.plot([X_A_final[i]],[Y_A_final[i]], 'go')
    plt.plot([X_B_final[i]],[Y_B_final[i]], 'ro')
#    plt.savefig('frame%03d.png'%i)
plt.show()








