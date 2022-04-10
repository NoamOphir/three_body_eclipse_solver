#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 12:35:16 2022

@author: noamophir
"""

#%% Packages 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#%% Algorithm Parameters
global G
global n
# global T_int

G     = (6.67408e-11 * (1/1.495978707)**3 
         *1e-33 * 1.98847e30 * (3600*24)**2)        ## Gravitational constant in units of au^3 sm^-1 days^-2

n     = 3                                           ## number of bodies
# year  = 115.5                                       ## [days]
# T_int = int(year*24*60)                             ## number of time intervals

# r     = np.zeros((n, T_int, 3))                     ## position
# v     = np.zeros((n, T_int, 3))                     ## velocity

def A(i, t, r):
    a_ = np.zeros(3)
    for j in range(n):
        if j != i:
            a_ += - G * m[j] * (r[i,t,:]-r[j,t,:]) / (np.linalg.norm(r[i,t,:]-r[j,t,:])**3)
    return a_

def rot_x(vector, theta):
    rot_mat = np.array([[1, 0                , 0             ],
                        [0, np.cos(theta)    , -np.sin(theta)],
                        [0, np.sin(theta)    , np.cos(theta) ]])
    return np.dot(rot_mat,vector)

# #%% Initial conditions 2-bodies
# global T_int

# year  = 120                                          ## [days]
# T_int = int(year*24*60)                             ## number of time intervals

# r     = np.zeros((n, T_int, 3))                     ## position
# v     = np.zeros((n, T_int, 3))                     ## velocity

# r[0,0,:] = np.array([0,0,0])
# r[1,0,:] = np.array([1,0,0])

# v[0,0,:] = np.array([0,0,0])
# v[1,0,:] = np.array([0,(2*np.pi)/365,0])

# a        = np.zeros((n,T_int,3))

# for i in range(n):
#     a[i,0,:] = A(i, 0, r)
# #%% Initial conditions 3-bodies: Binary planet and a star
# global m
# m1    = 1                                           ## masses in units of solar mass [sm]
# # m2    = 1/333000                                  #(Earth is 333000 times ligher then the sun)
# m2    = 0.001
# # m3    = 3
# m     = np.array([m1, m2, m2])                      ## masses array
# radii = np.array([0.0046524726,0,0])                ## object radii

# global T_int

# year  = 120                                          ## [days]
# T_int = int(year*24*60)                             ## number of time intervals

# r     = np.zeros((n, T_int, 3))                     ## position
# v     = np.zeros((n, T_int, 3))                     ## velocity

# r[0,0,:] = np.array([0,0,0])
# r[1,0,:] = np.array([1-0.01,0,0])
# r[2,0,:] = np.array([1+0.01,0,0])

# v[0,0,:] = np.array([0,0,0])
# v[1,0,:] = np.array([0,(2*np.pi)/365 - np.sqrt(G*m2/0.04),0])
# v[2,0,:] = np.array([0,(2*np.pi)/365 + np.sqrt(G*m2/0.04),0])

# a        = np.zeros((n,T_int,3))

# for i in range(n):
#     a[i,0,:] = A(i, 0, r)

#%% Initial conditions 3-bodies: Binary star and a planet
## Geometry: the Binary star moves on the x-y plan 
##           the plant       moves on the plan with angle theta with respect to the x-y plan  

global m

m1    = 5                                           ## masses in units of solar mass [solar mass]
# m2    = 1/333000                                  ## (Earth is 333000 times ligher then the sun)
m2    = 0.001                                       ## Jupiter-size planet [solar mass]
# m3    = 3
m     = np.array([m1, m1, m2])                      ## masses array [solar mass]
radii = np.array([0.005,0.005,0])                   ## object radii [au] (solar radii = 0.0046524726 au)

d1    = 0.01                                        ## distance m1 and center of rotation [au]
d2    = d1                                          ## distance m2 and center of rotation [au]
d3    = 1                                           ## distance m3 and center of rotation [au]

v1    = np.sqrt(G * m[1] * d1 / (d1 + d2)**2)       ## speed of m1 [au/days]
v2    = np.sqrt(G * m[0] * d2 / (d1 + d2)**2)       ## speed of m2 [au/days]
v3    = np.sqrt(G * (m[1]+m[0]) / d3)               ## speed of m3 [au/days]

P1    = 2 * np.pi * d1 / v1                         ## Time period between the binary objects  [days]
P2    = 2 * np.pi * d3 / v3                         ## Time period between the planet and the biany  [days]

global T_int

year  = P2                                          ## [days]
T_int = int(year*24*60)                             ## number of time intervals

r     = np.zeros((n, T_int, 3))                     ## position
v     = np.zeros((n, T_int, 3))                     ## velocity

theta = 90                                          ## angle (degrees) between x-y plan to the 
                                                    ## planet movement plan

r[0,0,:] = np.array([d1,0,0])                       ## inital position of m1 [au]
r[1,0,:] = np.array([-d2,0,0])                      ## inital position of m2 [au]
r[2,0,:] = np.array([d3,0,0])                       ## inital position of m3 [au]          

v[0,0,:] = np.array([0,v1,0])                       ## inital velocity of m1 [au/days]
v[1,0,:] = np.array([0,-v2,0])                      ## inital velocity of m2 [au/days]
v[2,0,:] = rot_x(np.array([0,v3,0]),theta*np.pi/180)## inital velocity of m3 [au/days]

a        = np.zeros((n,T_int,3))

for i in range(n):
    a[i,0,:] = A(i, 0, r)


print('Time period between the two stars is ' + str(round(P1*24,4)) +
      ' hours \nTime period of the planet is %.2f days' %P2)

#%% Functions
def Energy(r,v):
    E_kin = 0
    V     = 0
    for i in range(n):
        E_kin += 0.5 * m[i] * np.dot(v[i],v[i])
        for j in range(i-1):
            V += - G * m[i] * m[j] / np.linalg.norm(r[i,:]-r[j,:])
    
    return E_kin + V

def Angular_momentum(r,v):
    L = 0
    for i in range(n):
        L += np.cross(r[i,:],m[i]*v[i,:])
    return L

def transient(r, planet_indx, sun_indx):
    nu = 0
    
    for t in range(T_int):
        if (abs(r[planet_indx,t,0]-r[sun_indx,t,0])<radii[sun_indx]) and (r[planet_indx,t,1]<r[sun_indx,t,1]):
            nu += 1
    return nu

def transient_stat(r, v, radii, planet_indx, star_indx, theta, phi):
    theta       = np.pi * theta / 180
    phi         = np.pi * phi   / 180
    
    com_dist    = np.zeros(T_int)                      ### COM dist. at transiant
    area        = np.zeros(T_int)                      ### intersection area
    
    for t in range(T_int):
        r_rel = r[planet_indx,t,:]-r[star_indx,t,:]
        
        ## The viewer is in +infty on the z axis
        if (np.sqrt(r_rel[0]**2 + r_rel[1]**2) < radii[star_indx] and 
            (r[planet_indx,t,2]>r[star_indx,t,2])):
            
            com_dist[t] = np.sqrt(r_rel[0]**2 + r_rel[1] **2) / radii[star_indx]
        
        d     = np.sqrt(r_rel[0]**2 +r_rel[1]**2)
        rs    = radii[star_indx]
        rp    = radii[planet_indx]
        ds    = (rs**2-rp**2+d**2)/(2*d)
        dp    = d - ds
        
        if d < rs - rp and (r[planet_indx,t,2]>r[star_indx,t,2]):
            area[t] = 1
        elif d < rs + rp and (r[planet_indx,t,2]>r[star_indx,t,2]):
            area[t] = (rs**2 * np.arccos(ds/rs) - ds * np.sqrt(rs**2 - ds**2)
                       +rp**2 * np.arccos(dp/rp) - dp * np.sqrt(rp**2 - dp**2)) / (np.pi * rp**2)
        
    return com_dist , area

def transient_stat_x(r, v, radii, planet_indx, star_indx):
    
    com_dist    = np.zeros(T_int)                      ### COM dist. at transiant
    area        = np.zeros(T_int)                      ### intersection area
    
    for t in range(T_int):
        r_rel = r[planet_indx,t,:]-r[star_indx,t,:]
        
        ## The viewer is in +infty on the z axis
        if (np.sqrt(r_rel[1]**2 + r_rel[2]**2) < radii[star_indx] and 
            (r[planet_indx,t,0]>r[star_indx,t,0])):
            
            com_dist[t] = np.sqrt(r_rel[2]**2 + r_rel[1] **2) / radii[star_indx]
        
        d     = np.sqrt(r_rel[2]**2 +r_rel[1]**2)
        rs    = radii[star_indx]
        rp    = radii[planet_indx]
        ds    = (rs**2-rp**2+d**2)/(2*d)
        dp    = d - ds
        
        if d < rs - rp and (r[planet_indx,t,0]>r[star_indx,t,0]):
            area[t] = 1
        elif d < rs + rp and (r[planet_indx,t,0]>r[star_indx,t,0]):
            area[t] = (rs**2 * np.arccos(ds/rs) - ds * np.sqrt(rs**2 - ds**2)
                       +rp**2 * np.arccos(dp/rp) - dp * np.sqrt(rp**2 - dp**2)) / (np.pi * rp**2)
        
    return com_dist , area

def transient_counting_z(r):
    nu   = 0
    flag = 0
    
    for i in [0,1]:
        for t in range(T_int):   ### SHOULD BE A PERIOD TIME
            r_rel = r[2,t,:]-r[i,t,:]
            d     = np.sqrt(r_rel[0]**2 +r_rel[1]**2)
            rs    = radii[i]
            rp    = radii[2]
            
            if d < rs - rp and r[2,t,2] > r[i,t,2] and flag == 0:
                flag = 1
                nu   += 1
                
            elif d > rs - rp:
                flag = 0
    return nu

def transient_counting_x(r):
    nu   = 0
    flag = 0
    
    for i in [0,1]:
        for t in range(T_int):   ### SHOULD BE A PERIOD TIME
            r_rel = r[2,t,:]-r[i,t,:]
            d     = np.sqrt(r_rel[1]**2 +r_rel[2]**2)
            rs    = radii[i]
            rp    = radii[2]
            
            if d < rs - rp and r[2,t,0] > r[i,t,0] and flag == 0:
                flag = 1
                nu   += 1
                # print(i, t)                

            elif d > rs - rp:
                flag = 0
    return nu

def transient_counting(r, theta, phi):
    sint = np.sin(theta)
    cost = np.cos(theta)
    sinp = np.sin(phi)
    cosp = np.cos(phi)

    R    = np.array([[cost * cosp,  -cost * sinp,   sint],
                     [sint       , cost         , 0     ],
                     [-cosp * sint, sint * sinp  , cost]])       ## Rotation matrix
    
    nu   = 0
    flag = 0
    
    for i in [0,1]:
        for t in range(T_int):   ### SHOULD BE A PERIOD TIME
            rp_rot = np.dot(R, r[2,t,:])
            rs_rot = np.dot(R, r[i,t,:])
            
            r_rel  = rp_rot - rs_rot
            
            # r_rel  = np.dot(R, r[2,t,:]-r[i,t,:])
            d      = np.sqrt(r_rel[0]**2 +r_rel[1]**2)
            rs     = radii[i]
            rp     = radii[2]
            
            if d < rs - rp and rp_rot[2] > rs_rot[2] and flag == 0:
                flag = 1
                nu   += 1
                
            elif d > rs - rp:
                flag = 0
    return nu


def vec_par(r, v, radii, planet_indx, star_indx, theta, phi): ### ??
    theta = np.pi * theta / 180
    phi   = np.pi * phi   / 180
    
    b     = np.zeros((T_int,2))
    
    for t in range(T_int):
        # v_rel = v[planet_indx,t,:]-v[star_indx,t,:]
        r_rel = r[planet_indx,t,:]-r[star_indx,t,:]
        ## The viewer is in +infty on the x axis
        if np.linalg.norm(r_rel[1:2])<radii[star_indx] and (r[planet_indx,t,0]>r[star_indx,t,0]):
            b[t,:] = r_rel[1:2] / radii[star_indx]
            # print(r[planet_indx,t,:])
            # print(r[star_indx,t,:])
    return b
#%% Leapfrog
energy  = np.zeros(T_int-1)
ang_mom = np.zeros((T_int-1,3))
delta_t = 1/24/60                                    ## time interval = 1 min
for t in range(T_int-1):   
    for i in range(n):
        r[i,t+1,:]    = r[i,t,:] + v[i,t,:] * delta_t + 0.5 * a[i,t,:] * delta_t**2
        
    for i in range(n):        
        a[i,t+1,:]    = A(i, t+1, r)  
        v[i,t+1,:]    = v[i,t,:] + 0.5 * (a[i,t,:] + a[i,t+1,:]) * delta_t
    
    energy[t]    = Energy(r[:,t,:],v[:,t,:])
    ang_mom[t,:] = Angular_momentum(r[:,t,:], v[:,t,:])
    
    # print('\r', i, end='')
    print('\r Loading : %.2f %%' % (100*t/(T_int-1)) ,end=' ')
#%% Sanity check
fig, ax = plt.subplots(2,1,figsize=(7,8))
ax[0].set_xlabel("Time [days]")
ax[0].set_ylabel(r"energy of the system $[sm \cdot au^2 \cdot days^{-2}]$")
ax[0].set_ylim([min(energy)*0.99,max(energy)*1.01])

ax[0].plot(np.linspace(0,delta_t*T_int,T_int-1),energy)
ax[0].text(0.1*delta_t*T_int ,min(energy)*0.995,
            r"$\left< E \right>=$" + str(round(np.mean(energy),4)) +
            r"$\pm$" + str(round(np.std(energy),6)) + r"$\;[sm \cdot au^2 \cdot days{^-2}]$" +
            "\nError=" + str(round(100*np.std(energy)/np.mean(energy),4)) + "%")

ax[1].set_xlabel("Time [days]")
ax[1].set_ylabel(r"Angular momentum $L_z [sm \cdot au^2 \cdot days^{-1}]$")
ax[1].set_ylim([min(ang_mom[:,2])*0.99,max(ang_mom[:,2])*1.01])

ax[1].plot(np.linspace(0,delta_t*T_int,T_int-1),ang_mom[:,2])
ax[1].text(0.1*delta_t*T_int ,min(ang_mom[:,2])*0.995,
            r"$\left< L_z \right>=$" + str(round(np.mean(ang_mom[:,2]),6)) +
            r"$\pm$" + str(round(np.std(ang_mom[:,2]),6)) + r"$\;[sm \cdot au^2 \cdot days{^-1}]$" +
            "\nError=" + str(round(100*np.std(ang_mom[:,2])/np.mean(ang_mom[:,2]),6)) + "%")

plt.show()
#%% Animation
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
    
for t in range(T_int-1):   
    ax.cla()
    ax.set_xlim([-1.5,1.5])
    ax.set_ylim([-1.5,1.5])
    ax.set_zlim([-1.5,1.5])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')
    ax.set_zlabel('z [au]')
    
    ax.scatter(r[0,t,0],r[0,t,1],r[0,t,2],s=50, c='r')
    ax.scatter(r[1,t,0],r[1,t,1],r[1,t,2],s=50, c='b')
    ax.scatter(r[2,t,0],r[2,t,1],r[2,t,2],s=5,  c='g')


    ax.set_title("frame {} days".format(round(t/24,1)))
    # Note that using time.sleep does *not* work here!
    plt.pause(0.01)

#%% Graphing trajectories
fig, ax = plt.subplots(2,2,figsize=(10,10))

ax[0,0].set_xlim([-1.25*d3,1.25*d3])
ax[0,0].set_ylim([-1.25*d3,1.25*d3])
ax[0,0].set_xlabel("x [au]")
ax[0,0].set_ylabel("y [au]")

ax[0,0].scatter(r[0,:,0],r[0,:,1],s=10, alpha=0.5, label='star 1')
ax[0,0].scatter(r[1,:,0],r[1,:,1],s=10, alpha=0.5, label='star 2')
ax[0,0].plot(r[2,:,0],r[2,:,1], label='planet')

ax[0,0].grid()
ax[0,0].legend(loc=1)

ax[0,1].set_xlim([-1.25*d3,1.25*d3])
ax[0,1].set_ylim([-1.25*d3,1.25*d3])
ax[0,1].set_xlabel("y [au]")
ax[0,1].set_ylabel("z [au]")

ax[0,1].scatter(r[0,:,1],r[0,:,2],s=10, alpha=0.5, label='star 1')
ax[0,1].scatter(r[1,:,1],r[1,:,2],s=10, alpha=0.5, label='star 2')
ax[0,1].plot(r[2,:,1],r[2,:,2], label='planet')

ax[0,1].grid()
ax[0,1].legend(loc=1)

ax[1,0].set_xlim([-1.25*d3,1.25*d3])
ax[1,0].set_ylim([-1.25*d3,1.25*d3])
ax[1,0].set_xlabel("x [au]")
ax[1,0].set_ylabel("z [au]")

ax[1,0].scatter(r[0,:,0],r[0,:,2],s=10, alpha=0.5, label='star 1')
ax[1,0].scatter(r[1,:,0],r[1,:,2],s=10, alpha=0.5, label='star 2')
ax[1,0].plot(r[2,:,0],r[2,:,2], label='planet')

ax[1,0].grid()
ax[1,0].legend(loc=1)

plt.show()
#%% Ploting transient map
y1 = transient(r, 2, 1)
y2 = transient(r, 2, 0)

plt.plot(np.linspace(0,delta_t*T_int,T_int),y1)
plt.plot(np.linspace(0,delta_t*T_int,T_int),y2)

plt.xlabel("time [days]")

plt.show()

#%% Ploting impact distance from COM map
b1, a1 = transient_stat_x(r, v, radii, 2, 0)
b2, a2 = transient_stat_x(r, v, radii, 2, 1)

# plt.plot(np.linspace(0,delta_t*T_int,T_int),b1, label=r"$d_{planet-star1} \; [R_{\star}]$")
# plt.plot(np.linspace(0,delta_t*T_int,T_int),b2, label=r"$d_{planet-star2} \; [R_{\star}]$")

plt.plot(np.linspace(0,delta_t*T_int,T_int),a1, label=r"$A_{planet-star1} \; [\pi\cdot r_p^2]$")
plt.plot(np.linspace(0,delta_t*T_int,T_int),a2, label=r"$A_{planet-star2} \; [\pi\cdot r_p^2]$")

# ## Close-up
# time_left  = 28.25
# time_right = 30
# plt.plot(np.linspace(time_left,time_right,len(a1[int(time_left/delta_t):int(time_right/delta_t)])),
#          a1[int(time_left/delta_t):int(time_right/delta_t)],
#          label=r"$A_{planet-star1} \; [\pi\cdot r_p^2]$")
# plt.plot(np.linspace(time_left,time_right,len(a1[int(time_left/delta_t):int(time_right/delta_t)])),
#          a2[int(time_left/delta_t):int(time_right/delta_t)],
#          label=r"$A_{planet-star2} \; [\pi\cdot r_p^2]$")

plt.xlabel("time [days]")
# plt.ylabel(r"distance from centers of star and planet $[R_{\star}]$")
plt.legend()

plt.show()
#%% Animation of passing over the star
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot()
    
for t in range(int(T_int*0.6315),int(T_int*0.6338)):   
    ax.cla()
    ax.set_xlim([-0.25,0.25])
    ax.set_ylim([-0.25,0.25])
    # ax.set_zlim([-1.5,1.5])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')
    # ax.set_zlabel('z [au]')
    
    ax.scatter(r[0,t,1],r[0,t,2],s=50, c='r')
    ax.scatter(r[1,t,1],r[1,t,2],s=50, c='b')
    ax.scatter(r[2,t,1],r[2,t,2],s=5,  c='g')


    ax.set_title("frame {} days".format(round(t/24,1)))
    # Note that using time.sleep does *not* work here!
    plt.pause(0.01)
    # print("hello")
#%% Uniform dist. over a sphere

partition   = 1000
# Phi         = np.linspace(0,0.5*np.pi, partition)
# Theta       = np.linspace(0, np.pi,partition)
# Theta       = np.arccos(np.linspace(-1,1,partition))

# spacial_vec   = np.zeros((partition**2, 3))

# for i in range(partition):
#     for ii in range(partition):
#         spacial_vec[i*partition +ii,:]  = [np.sin(Theta[i])*np.cos(Phi[ii]),
#                                             np.sin(Theta[i])*np.sin(Phi[ii]),
#                                             np.cos(Theta[i])]

Phi         = np.random.rand(partition) * 2 * np.pi
Theta       = np.arccos(np.random.rand(partition) * 2 -1)
       
spacial_vec   = np.zeros((partition, 3))

for i in range(partition):
    spacial_vec[i,:] = [np.sin(Theta[i])*np.cos(Phi[i]),
                        np.sin(Theta[i])*np.sin(Phi[i]),
                        np.cos(Theta[i])]
#%% Ploting the distrabution
fig         = plt.figure()
ax          = fig.add_subplot(111, projection='3d')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])

ax.scatter(spacial_vec[:,0],spacial_vec[:,1],spacial_vec[:,2],s=2, alpha=0.5)

plt.show()

fig, ax = plt.subplots(1,figsize=(10,10))

# color = np.random.random(partition)
# ax.scatter(Phi, Theta, alpha=0.5, c=color)

ax.scatter(Phi, Theta, alpha=0.5)

ax.set_xlabel(r"$\phi$", fontsize=16)
ax.set_ylabel(r"$\theta$", fontsize=16)

ax.set_yticks([0,np.pi/2, np.pi])
ax.set_yticklabels([0,r"$\pi/2$",r"$\pi$"])

ax.set_xticks([0,np.pi/2, np.pi, np.pi*3/2, 2*np.pi])
ax.set_xticklabels([0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])

ax.grid()

plt.show()
#%% Stat. diffrent angles
trans = np.zeros(partition)

for i in range(partition):
    trans[i] = transient_counting(r, Theta[i], Phi[i])
    
    print('\r Loading: %.2f %%' %(i*100/partition), end='')

#%% Ploting the Stat. diffrent angles
fig, ax = plt.subplots(1,figsize=(10,10))

map1 = ax.scatter(Phi, Theta, c=trans, cmap='viridis')
# map1 = ax.imshow(np.stack([Phi, Theta]),cmap='viridis')
fig.colorbar(map1)
# plt.colorbar()

ax.set_xlabel(r"$\phi$", fontsize=16)
ax.set_ylabel(r"$\theta$", fontsize=16)

ax.set_yticks([0,np.pi/2, np.pi])
ax.set_yticklabels([0,r"$\pi/2$",r"$\pi$"])

ax.set_xticks([0,np.pi/2, np.pi, np.pi*3/2, 2*np.pi])
ax.set_xticklabels([0,r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])

ax.grid()

plt.show()

plt.hist(trans, bins=11)
plt.text(1,800,r"$N_{\geq 4} =$"+
         " %.f ( %.2f %% )\n" %(sum(i >= 4 for i in trans),(sum(i >= 4 for i in trans)/partition*100)) +
         r"$N_{\geq 3} =$"+
                  " %.f ( %.2f %% )" %(sum(i >= 3 for i in trans),(sum(i >= 3 for i in trans)/partition*100)),
         fontsize=16)

plt.xlabel("Number of eclipses in a single solid angle observation on a time period")

plt.show()
x = sum(i >= 4 for i in trans)

print(str(x) + " out of " + str(partition))