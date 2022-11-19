
import numpy as np
 
s_0 = 999999
i_0 = 1
r_0 = 0
t_0 = 0
h = 1
contador = 1
beta = float(input("Beta: "))
gamma = float(input("Gamma: "))
N = 1000000

print("Sistema de ecs diferenciales del modelo SIR")

def funcions(t,s,i,r):
    return (-beta/N)*s*i

def funcioni(t,s,i,r):
    return (beta/N)*s*i-gamma*i

def funcionr(t, s,i,r):
    return gamma*i

interval = int(input("Interval: [days] "))

while(contador<=interval):
    k_1 = funcions(t_0,s_0,i_0,r_0)
    m_1 = funcioni(t_0,s_0,i_0,r_0)
    l_1 = funcionr(t_0,s_0,i_0,r_0)

    k_2 = funcions(t_0 + h*0.5, s_0 + h*.5*k_1, i_0 + h*0.5*m_1, r_0 + h*.5*l_1)
    m_2 = funcioni(t_0 + h*0.5, s_0 + h*.5*k_1, i_0 + h*0.5*m_1, r_0 + h*.5*l_1)
    l_2 = funcionr(t_0 + h*0.5, s_0 + h*.5*k_1, i_0 + h*0.5*m_1, r_0 + h*.5*l_1)

    k_3 = funcions(t_0 + h*0.5, s_0 + h*.5*k_2, i_0 + h*0.5*m_2, r_0 + h*.5*l_2)
    m_3 = funcioni(t_0 + h*0.5, s_0 + h*.5*k_2, i_0 + h*0.5*m_2, r_0 + h*.5*l_2)
    l_3 = funcionr(t_0 + h*0.5, s_0 + h*.5*k_2, i_0 + h*0.5*m_2, r_0 + h*.5*l_2)

    k_4 = funcions(t_0 + h*0.5, s_0 + h*k_3, i_0 + h*m_3, r_0 + h*l_3)
    m_4 = funcioni(t_0 + h*0.5, s_0 + h*k_3, i_0 + h*m_3, r_0 + h*l_3)
    l_4 = funcionr(t_0 + h*0.5, s_0 + h*k_3, i_0 + h*m_3, r_0 + h*l_3)

    s_0 = s_0 + (h/6)*(k_1+2*k_2+2*k_3+k_4)
    i_0 = i_0 + (h/6)*(m_1+2*m_2+2*m_3+m_4)
    r_0 = r_0 + (h/6)*(l_1+2*l_2+2*l_3+l_4)

    t_0 = t_0 + h
    contador = contador+1
    print("día :", t_0, "personas susceptibles:",s_0)
    print("día :",t_0, "personas infectados:",i_0)
    print("día: ", t_0, "personas recuperados:",r_0)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def f(s,t):
  n = s[0]
  c = s[1]
  j = s[2]
  dsdt = (-beta/N)*n*c
  didt = (beta/N)*n*c-gamma*c
  drdt = gamma*c
  return [dsdt, didt,drdt]

t = np.linspace(0,interval)
s0=[999999,1,0]

s = odeint(f,s0,t)

plt.plot(t,s[:,0],'r--', linewidth=2.0)
plt.plot(t,s[:,1],'b-', linewidth=2.0)
plt.plot(t,s[:,2],'g-', linewidth=2.0)
plt.xlabel("t")
plt.ylabel("Modelo SIR")
plt.legend(["S","I","R"])
plt.show()
