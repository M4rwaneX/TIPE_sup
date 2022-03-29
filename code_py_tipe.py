import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


## oscillateurs couplés
# point fixe - ressort1 - masse1-ressort2-masse2
tf = 50        # s : temps final
N = 10000        # nombre de points
h = tf/(N-1)    # pas temporel
Fm=1 #

m1= 106  #kg
k1=90   # N/m

m2= 1   # kg
k2 = 0.1 # N/m
W1=np.sqrt(k1/m1)
W2=np.sqrt(k2/m2)
print('T1 =',2*3.14/W1)
print('T2 =',2*3.14/W2)
# on ajoute une force de frottement fluide f = -m.v/tau sur la masse m au bout du second ressort
tau1= 10
0 # frottement faible sur la masse 1
tau2 = 5 # constante de temps d'amortissement caractéristique du frottement


## accélérations : fonction dérivée

# équations :
# m1*a1=-k1*x1+k2*(x2-x1)-m1/tau1*v1 +Fm*np.cos(W1*t)
# et m2*a2=-k2*(x2-x1)-m2/tau2*v2

# Z =[x1,v1,x2,v2]
def fp(t,Z) :
    return Z[1],-k1/m1*Z[0]+k2/m1*(Z[2]-Z[0])-Z[1]/tau1+(Fm*np.cos(W1*t))/m1,Z[3],-k2/m2*(Z[2]-Z[0])-Z[3]/tau2+(Fm*np.cos(W1*t))/m2
    # dans l'ordre : les dérivées de x1(t), x1(t), x2(t), v2(t)

## résolution par solve_ivp :

Z0 =[0.1,0,0,0]
# on prendra implicitement x2(0) = 0 et v1(0)=v2(0) = 0 et x1(0) différent de 0
solution=solve_ivp(fp,[0,tf],Z0,t_eval=np.linspace(0,tf,N))


t_list=solution.t
x1_list=solution.y[0]
x2_list=solution.y[2]
v1_list=solution.y[1]
v2_list=solution.y[3]

'''plt.subplot(1,2,2)'''
plt.plot(t_list,x1_list,'b',label=("x1(t)"))
plt.plot(t_list,x2_list,'r',label=("x2(t)"))
plt.xlabel('t')
plt.grid()
plt.legend()

plt.show()
