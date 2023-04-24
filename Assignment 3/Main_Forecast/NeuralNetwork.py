#Forecasting the Lorenz system with Neural Network 
#Code adapted from https://github.com/jbramburger/DataDrivenDynSyst Learning Dynamics with Neural Networks\Forecast.ipynb

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

#------ Begin Generating Lorenz Data

# Initializations

ti = 0
tfin = 500
N = 5001
sol_rhoall = np.zeros((3*N,4))
t = np.linspace(ti,tfin,N)
dt = t[2]-t[1]
print("dt = ", dt)

# Lorenz parameter #1
rho = 10.0
x0 = [1,2,3, rho]
sol_rho10 = scipint.odeint(Lorenz,x0, t, args=(rho,0))
sol_rhoall[0:N,:] = sol_rho10

fig = plt.figure()
plt.plot(sol_rho10[:,0],sol_rho10[:,2],'k') 
plt.title('The Lorenz Attractor with parameter rho = 10.0', fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$z$', fontsize = 20)

# Lorenz parameter #2
rho = 28.0;
x0 = [1,2,3, rho]
sol_rho28 = scipint.odeint(Lorenz,x0, t, args=(rho,0))
sol_rhoall[N:2*N,:] = sol_rho28

fig = plt.figure()
plt.plot(sol_rho28[:,0],sol_rho28[:,2],'k') 
plt.title('The Lorenz Attractor rho = 28.0', fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$z$', fontsize = 20)

# Lorenz parameter #3
rho = 40.0;
x0 = [1,2,3, rho]
sol_rho40 = scipint.odeint(Lorenz,x0, t, args=(rho,0))
sol_rhoall[2*N:3*N,:] = sol_rho40

fig = plt.figure()
plt.plot(sol_rho40[:,0],sol_rho40[:,2],'k') 
plt.title('The Lorenz Attractor rho = 40.0', fontsize = 20)
plt.xlabel('$x$', fontsize = 20)
plt.ylabel('$z$', fontsize = 20)

#------ End Generating Lorenz Data


#------ Begin Forecastiong using neural networks with the Lorenz Data

forward_iters = 10 # Number of forward iterations
xnforward = [] #initialize matrix for training data


ENDrow = np.size(sol_rhoall,0)
for j in range(forward_iters):
    xnforward.append(sol_rhoall[j:ENDrow-forward_iters+j])

input("Press Enter to continue...")

# Initializes the neural network model
def init_model(num_hidden_layers = 10, num_neurons_per_layer = 100):

    model = tf.keras.Sequential()

    # Input is (x,y,z,rho)
    model.add(tf.keras.Input(4))

    for _ in range(num_hidden_layers):
        #adds the number of layer at each _ hidden layers
        model.add(tf.keras.layers.Dense(num_neurons_per_layer,
            activation=tf.keras.activations.get('relu'), 
            kernel_initializer='glorot_normal'))

        #model.add(tf.keras.layers.Dropout(0.3))

    # Output is (x,y,z,rho)
    model.add(tf.keras.layers.Dense(4))

    return model

def compute_loss(model, xnforward, steps):

    loss = 0

    for s in range(steps):

        if s == 0:
            xpred = model(xnforward[0])
        else:
            xpred = model(xpred) # x_(n+1) = model(x_n)

        xnp1 = xnforward[s+1] # Gets the next true point to compare with the prediction

        loss += tf.reduce_mean(tf.square(xpred-xnp1))/steps

    return loss

def get_grad(model, xnforward, steps):

    with tf.GradientTape(persistent=True) as tape:
        # This tape is for derivatives with respect to trainable vriables.
        tape.watch(model.trainable_variables)
        loss = compute_loss(model, xnforward, steps)

    g = tape.gradient(loss, model.trainable_variables)
    del tape

    return loss, g

# get neural network model
num_hidden_layers = 10
num_neurons_per_layer = 100;
model = init_model(num_hidden_layers,num_neurons_per_layer) 

# Learning rate chosen as increasing steps
lr = tf.keras.optimizers.schedules.PiecewiseConstantDecay([1000,4000,8000,15000,20000], [1e-3,5e-4,1e-4,5e-5,2.5e-5,1e-5])

optim = tf.keras.optimizers.Adam(learning_rate=lr)


# add time function from the time package
from time import time

steps = 3


@tf.function
def train_step():
    # Compute current loss and gradient w.r.t. parameters.
    loss, grad_theta = get_grad(model, xnforward, steps)

    # Perform gradient descent step
    optim.apply_gradients(zip(grad_theta, model.trainable_variables))

    return loss

# Number of training epochs
N_training = 30000
Loss_hist = [] # Matrix to collect losses

# Start timer
t0 = time()


# Train the data
for i in range(N_training+1):
    loss = train_step()

    Loss_hist.append(loss.numpy())

    if i%50 == 0:
        print('It {:05d}: loss = {:10.8e}'.format(i,loss))

# Print overal computation time
CompTime = time()-t0
print('\nComputation time:{} seconds'.format(CompTime))



# Use Trained Model to Forecast
M = 201

xpred = np.zeros((M,4))


rho = 10; # <--- Change parameter
xpred[0] = [1, 2, 3, rho] #inital conditions
for m in range(1,M):
    xpred[m] = model(xpred[m-1:m,:])


if rho == 28:
    DataStart = N
elif rho == 40:
    DataStart = 2*N
elif rho == 10:
    DataStart = 0
else:
    x0 = [1,2,3, rho]
    sol_rho_cust = scipint.odeint(Lorenz,x0, t, args=(rho,0))

    fig = plt.figure()
    plt.plot(sol_rho_cust[:,0],sol_rho_cust[:,2],'k') 
    plt.title('The Lorenz Attractor with parameter rho = other', fontsize = 20)
    plt.xlabel('$x$', fontsize = 20)
    plt.ylabel('$z$', fontsize = 20)

DataEndPred = 150
DataStart = 0
DataEnd = DataStart+DataEndPred



if rho == 10 or rho == 28 or  rho == 40:

    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,0],'b--o')
    plt.plot(sol_rhoall[DataStart:DataEnd,0], 'r.')
    plt.title('Forecasting Lorenz xt with Neural Networks', fontsize = 20)
    plt.xlabel('$t$', fontsize = 20)
    plt.ylabel('$x$', fontsize = 20)

    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,2],'b--o')
    plt.plot(sol_rhoall[DataStart:DataEnd,2], 'r.')
    plt.title('Forecasting Lorenz zt with Neural Networks', fontsize = 20)
    plt.xlabel('$t$', fontsize = 20)
    plt.ylabel('$z$', fontsize = 20)

    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,0], xpred[:DataEndPred,2],'b--o')
    plt.plot(sol_rhoall[DataStart:DataEnd,0],sol_rhoall[DataStart:DataEnd,2], 'r.')
    plt.title('Forecasting Lorenz xz with Neural Networks', fontsize = 20)
    plt.xlabel('$x$', fontsize = 20)
    plt.ylabel('$z$', fontsize = 20)


else:
    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,0],'b--o')
    plt.plot(sol_rho_cust[:DataEndPred,0], 'r.')
    plt.title('Forecasting Lorenz xt with Neural Networks', fontsize = 20)
    plt.xlabel('$t$', fontsize = 20)
    plt.ylabel('$x$', fontsize = 20)

    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,2],'b--o')
    plt.plot(sol_rho_cust[:DataEndPred,2], 'r.')
    plt.title('Forecasting Lorenz zt with Neural Networks', fontsize = 20)
    plt.xlabel('$t$', fontsize = 20)
    plt.ylabel('$z$', fontsize = 20)

    fig = plt.figure()
    plt.plot(xpred[:DataEndPred,0], xpred[:DataEndPred,2],'b--o')
    plt.plot(sol_rho_cust[:DataEndPred,0],sol_rho_cust[:DataEndPred,2], 'r.')
    plt.title('Forecasting Lorenz xz with Neural Networks', fontsize = 20)
    plt.xlabel('$x$', fontsize = 20)
    plt.ylabel('$z$', fontsize = 20)


#Show all the plots
plt.show(block=True)

#Save model
model.save('Lorenz_models/LorenzTRNDALL_step=3_rho=10largedt')

# Save data as .mat file
import scipy.io

Param = [dt, N, num_hidden_layers, num_neurons_per_layer, steps, CompTime]
if rho == 10 or rho == 28 or  rho == 40:
    scipy.io.savemat('LorenzTRNDALL_step=3_rho=10largedt.mat', dict(xpred = xpred, xtrue = sol_rhoall[DataStart:DataStart+M,:], Param = Param, rho = rho, loss = Loss_hist))
else:
    scipy.io.savemat('LorenzTRNDALL_step=3_rho=17LossN10000.mat', dict(xpred = xpred, xtrue = sol_rho_cust[:M,:], Param = Param, rho = rho, loss = Loss_hist))