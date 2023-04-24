import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scipint

# Lorenz System function definition
def Lorenz(y, t, rho, b):
    
    dydt = [10.0*(y[1]-y[0]), y[0]*(rho - y[2]) - y[1], y[0]*y[1]-8.0/3*y[2], 0]

    return dydt

# Load and view saved models
#saved_model = tf.keras.models.load_model('Henon_models/Henon_step=3_b=0_01')
saved_model = tf.keras.models.load_model('Lorenz_models/LorenzTRNDALL_step=3_rho=10LossN10000')
saved_model.summary()

M = 201

xpred = np.zeros((M,4))

rho = 35;
x0 = [1,2,3, rho]

xpred[0] = x0 #inital conditions
for m in range(1,M):
    xpred[m] = saved_model(xpred[m-1:m,:])

DataEndPred = 150

ti = 0
tfin = 500 #was 500
N = 10001
t = np.linspace(ti,tfin,N)
dt = t[2]-t[1]
print("dt = ", dt)

sol_rho_cust = scipint.odeint(Lorenz,x0, t, args=(rho,0))

fig = plt.figure()
plt.plot(sol_rho_cust[:,0],sol_rho_cust[:,2],'k') 

fig = plt.figure()
plt.plot(xpred[:DataEndPred,0], xpred[:DataEndPred,2],'b--o')
plt.plot(sol_rho_cust[:DataEndPred,0], sol_rho_cust[:DataEndPred,2], 'r.') 
plt.title('Forecasting Lorenz xz with Neural Networks', fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$z$', fontsize = 20)

fig = plt.figure()
plt.plot(xpred[:DataEndPred,0],'b--o')
plt.plot(sol_rho_cust[:DataEndPred,0], 'r.') 
plt.title('Forecasting Lorenz x-t with Neural Networks', fontsize = 20)
plt.xlabel('$n$', fontsize = 20)
plt.ylabel('$x$', fontsize = 20)

fig = plt.figure()
plt.plot(xpred[:DataEndPred,2],'b--o')
plt.plot(sol_rho_cust[:DataEndPred,2], 'r.') 
plt.title('Forecasting Lorenz z-t with Neural Networks', fontsize = 20)
plt.xlabel('$n$', fontsize = 20)
plt.ylabel('$z$', fontsize = 20)


plt.show(block=True)


# Save data as .mat file
import scipy.io
scipy.io.savemat('Lorenz_rho=35_Model10.mat', dict(xpred35 = xpred, sol_35 = sol_rho_cust[:M,:]))