import numpy as np
import matplotlib.pyplot as plt

n = 1000000 # Población total

v_0 = 0 # Proporción de la población inicialmente vacunada

#### Valores iniciales ####
s_0 = (n - 1) * (1 - v_0) 
i_0 = 1
r_0 = (n - 1) * v_0
t_0 = 0

h = 1
contador = 1

#### Parámetros ####
beta = float(input("Beta: "))
gamma = float(input("Gamma: "))

b = 1/(70*365) # Natalidad
mu = 1/(70*365) # Mortalidad

r0 = beta/gamma

#### Campaña de vacunación ####
v = 0 # Proporción de vacunados (por día)


print("---- Sistema de ecuaciones diferenciales del modelo SIR ----")

def funcions(t,s,i,r):
    return (-beta/n)*s*i + b*n - mu*s - v*s

def funcioni(t,s,i,r):
    return (beta/n)*s*i - gamma*i - mu*i

def funcionr(t,s,i,r):
    return gamma*i - mu*r + v*s


def reff(s):
    return r0/n * s

iterations = int(input("Iteraciones: "))

s = [s_0]
i = [i_0]
r = [r_0]

while(contador <= iterations):
    n = s_0 + i_0 + r_0 # Actualizar población total

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

    s.append(s_0)
    i.append(i_0)
    r.append(r_0)

    t_0 = t_0 + h
    contador = contador + 1
    print("día:", t_0, "\n\tpersonas susceptibles:", s_0)
    print("\tpersonas infectadas:", i_0)
    print("\tpersonas recuperadas:",r_0)


#### Gráfico ####
t = np.linspace(0, t_0, iterations + 1)

plt.plot(t, s,'r--', linewidth=2.0)
plt.plot(t, i,'b-', linewidth=2.0)
plt.plot(t, r,'g-', linewidth=2.0)
plt.xlabel("t")
plt.ylabel("Modelo SIR")
plt.legend(["S","I","R"])
plt.show()

