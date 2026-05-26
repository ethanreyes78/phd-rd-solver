import numpy as np
#Define the Parameters
Lambda = .01 #population influx rate
alpha = .2 #environment-to-human transmission rate,
beta = .1 #human to human transmission
mu = 5*np.e**(-5) #natural death rate
omega = .01 #disease induced death rate
gamma = .2 #disease recovery rate
delta = .1 #removal rate of the pathogen from the aquatic environment
zi = 0.3 #shedding rate of infected hosts
g = 0.3 #is the intrinsic growth rate of the waterborne pathogen
K = 2.0 #carrying capacity
d = 1.0 #pathogen diffusion rate
v = 2.0 #pathogen advection rate that represents the speed of the river flow
#define diffusion rates of susceptible infected
c1 = 1.0
c2 = 2.0
c3 = 1.5
#Compute first reation term
def Reaction1(S, I, B):
return Lambda - alpha*S*B - beta*S*I - mu*S
#Compute second reaction term
def Reaction2(S, I, B):
return alpha *S*B + beta*S*I - (mu + omega + gamma)*I
#compute third reaction term
def Reaction3(I, R):
return gamma*I - mu*R
#compute fourth reaction term
def Reaction4(B, I):
return g*B * (1 - B/K) + zi*I - delta*B
def FiniteDiffSolver(L, T, nx, nt):
#initialize variables
dx = L / (nx - 1)
dt = T / nt
dS = np.zeros(nx)
dI = np.zeros(nx)
dR = np.zeros(nx)
dB = np.zeros(nx)
S = np.zeros(nx)
B = np.zeros(nx)
I = np.zeros(nx)
R = np.zeros(nx)
S[:] = 100.0 # Initial susceptible population
I[100] = 1.0 # Introduce 1 infected person in the middle of the domain

1

#iterate over the time nodes
for n in range(nt):
# enforce BCs
S[0] = S[1]
S[-1] = S[-2]
I[0] = I[1]
I[-1] = I[-2]
R[0] = R[1]
R[-1] = R[-2]
#special boundary conditions for B
B[0] = B[1] / (1 + dx*v/d)
B[-1] = B[-2] / (1 - dx*v/d)
dS[:] = 0.0
dI[:] = 0.0
dR[:] = 0.0
dB[:] = 0.0
#iterate over the spatial nodes
for i in range(1, nx-1):
#compute diffusion term + reaction term
dS[i] = c1 * (S[i+1] - 2*S[i] + S[i-1]) / (dx**2) + \
Reaction1(S[i], I[i], B[i])
dI[i] = c2 * (I[i+1] - 2*I[i] + I[i-1]) / (dx**2) + \
Reaction2(S[i], I[i], B[i])
dR[i] = c3 * (R[i+1] - 2*R[i] + R[i-1]) / (dx**2) + \
Reaction3(I[i], R[i])
dB[i] = d * (B[i+1] - 2*B[i] + B[i-1]) / (dx**2) - \
v * (B[i+1] - B[i-1]) / (2 * dx) + Reaction4(B[i], I[i])
#time update
S += dt * dS
I += dt * dI
R += dt * dR
B += dt * dB
return S, I, R, B

sol = FiniteDiffSolver(10.0, 1.0, 200, 5000)
print("S = " , sol[0])
print()
print("I = " , sol[1])
print()
print("R = " , sol[2])
print()
print("B = " , sol[3])

2

print()