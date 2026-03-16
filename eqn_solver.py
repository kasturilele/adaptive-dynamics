#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 14:01:53 2025

@author: kasturilele
"""
#%%
#import all modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve 
import csv


#%%
# selection gradient from hitchhiker guide paper
def inequality_function(x, y):
    return (y - x) - (np.exp(y) - np.exp(x)) / 10

def tradeoff(r):
    return 0.001 * r


#tradeoff function from boots(1999) + various other tradeoff functions
def bt(b):
    return -0.0000206 * np.exp(3.2*b) #exponential 
    
    #ignore all the ones below this
    #return -0.00001 #constant
    #return -(0.0004*b) #linear 
    #return 1/(-2*(b + 0.12)) + 0.002 #convex boots
    #return 1/(-1*(b - 1.2)) + 0.01 #concave boots
    #return -0.0001 * b*b #square
    #return -(b*b*b) #square but more complex
    #return -0.00001 * b /(1 - b) #to get linear r-K tradeoff
    #return -0.00001 /(1-b) #inverse
    #return -0.00001 * b / ((b - 0.5)*(b-0.5) + 0.01) #inverse parabola tradeoff between r and K
    #return -0.00001 * b / np.log(10*b + 1) #to get log relationship between r and K
    #return -0.00001 * np.log(2-b) #to get log relationship between r and a
    #return -0.00001/(3*(b + 0.00001))  #modified convex boots
    #return -0.00001/(-1*(b - 1.2)) #modified concave boots


def bt1(b):
    return -0.0000206 * np.exp(3.2*b) #exponential 
    #return -0.00001 /(1-b) #inverse

def bt2(b):
    return -0.0000108* np.exp(3.5*b) #exponential 
    #return -0.00001 /(1.5-b) #inverse



# my selection gradient
def selection_gradient1(x, y):
    return y - (x * bt1(y)/bt1(x))
def selection_gradient2(x, y):
    return y - (x * bt2(y)/bt2(x))


# selection gradient without constraints
def sel_no_con(x, y):
    return 1 - (y/x)

#
def boots_function(x, y):
    return bt(y) - bt(x) - ((y-x)*(bt(x) - 0.003*1/x)/(0.003 + x))

def Kfun(r,a):
    return -r/a

# my selection gradient
def selection_gradient_paired(x, y):
    K1 = -1.5e-05
    K2 = 9e-10
    a22 = -0.00005
    return y + (bt(y)*(K1 - a22*x)/(a22*bt(x) - K2)) + ((K2*x - K1*bt(x))/(a22*bt(x) - K2))

# paired selection gradient where both species can have residents and mutants
def sel_grad_paired_new(x1, y1, x2, y2):
    a12 = -2.5e-05
    a21 = -2.5e-05
    den = bt2(x2)*bt1(x1) - a12*a21
    s1m = y1 + bt1(y1)*(a12*x2 - bt2(x2)*x1)/den + (a12*a21*x1 - bt1(x1)*a12*x2)/den
    s2m = y2 + bt2(y2)*(a21*x1 - bt1(x1)*x2)/den + (a12*a21*x2 - bt2(x2)*a21*x1)/den
    return s1m, s2m

def Kpaired(r1,r2,a11,a22,a12,a21):
    K1 = (r2*a12 - r1*a22)/(a11*a22 - a12*a21)
    K2 = (r1*a21 - r2*a11)/(a11*a22 - a12*a21)
    
    return K1,K2


e = 1 #endpoint for all 3 graphs

# Create a grid of x and y values
x_min, x_max = 0, e
y_min, y_max = 0, e
num_points = 200

x = np.linspace(x_min, x_max, num_points)
y = np.linspace(y_min, y_max, num_points)
X, Y = np.meshgrid(x, y)

# Evaluate the inequality at each point in the grid
Z1s = selection_gradient1(X, Y)
Z2s = selection_gradient2(X, Y)


# Create the plot
plt.figure(figsize=(6,6))

# Plot the region where the inequality is satisfied (Z > 0)
plt.contourf(X, Y, Z1s, levels=[0, Z1s.max()], colors=['lightblue'], extend='max')

# Optionally, plot the boundary where the expression equals 0
plt.contour(X, Y, Z1s, levels=[0], colors='red', linestyles='--')

plt.xlabel('resident')
plt.ylabel('mutant')
plt.title('Region where mutants can invade resident population (blue) for species 1')
plt.grid(True)

plt.show()

# Create the plot
plt.figure(figsize=(6,6))

# Plot the region where the inequality is satisfied (Z > 0)
plt.contourf(X, Y, Z2s, levels=[0, Z2s.max()], colors=['lightblue'], extend='max')

# Optionally, plot the boundary where the expression equals 0
plt.contour(X, Y, Z2s, levels=[0], colors='red', linestyles='--')

plt.xlabel('resident')
plt.ylabel('mutant')
plt.title('Region where mutants can invade resident population (blue) for species 2')
plt.grid(True)

plt.show()


#make plot of tradeoff function (between r and alpha)

xt = np.linspace(0, e, 100) # define the x space
yt1 = bt1(xt)
yt2 = bt2(xt)
 
#Create the plot
plt.plot(xt, yt1, color='#AF0040')
plt.plot(xt, yt2, color='#1E88E5')

#Add labels and title (optional)
plt.xlabel("r")
plt.ylabel("alpha")
plt.title("tradeoff between r and alpha")
plt.grid(True)
plt.show()

#make plot of tradeoff function (between r and K, knowing that K = r/alpha)

yk1 = Kfun(xt, yt1)
yk2 = Kfun(xt, yt2)

#Create the plot
plt.plot(xt, yk1, color='#AF0040')
plt.plot(xt, yk2, color='#1E88E5')

#Add labels and title (optional)
plt.xlabel("r")
plt.ylabel("K")
plt.title("tradeoff between r and K")
plt.grid(True)
plt.show()

# Evaluate the inequality at each point in the grid (for pairwise function)
#two z variables, one for each species
Z1,Z2 = sel_grad_paired_new(X, Y, X, Y)
#Z = boots_function(X, Y)


# Create the plot
plt.figure(figsize=(6,6))

# Plot the region where the inequality is satisfied (Z > 0)
plt.contourf(X, Y, Z1, levels=[0, Z1.max()], colors=['lightblue'], extend='max')

# Optionally, plot the boundary where the expression equals 0
plt.contour(X, Y, Z1, levels=[0], colors='orange', linestyles='--')

plt.xlabel('resident')
plt.ylabel('mutant')
plt.title('Region where mutants can invade resident population (blue) for species 1')
plt.grid(True)

plt.show()

# Create the plot
plt.figure(figsize=(6,6))

# Plot the region where the inequality is satisfied (Z > 0)
plt.contourf(X, Y, Z2, levels=[0, Z2.max()], colors=['lightgreen'], extend='max')

# Optionally, plot the boundary where the expression equals 0
plt.contour(X, Y, Z2, levels=[0], colors='red', linestyles='--')

plt.xlabel('resident')
plt.ylabel('mutant')
plt.title('Region where mutants can invade resident population (green) for species 2')
plt.grid(True)

plt.show()

#make plot of tradeoff function (between r and alpha)

x1t = np.linspace(0, e, 100) # define the x space
y1t = bt1(x1t)
x2t = np.linspace(0, e, 100) # define the x space
y2t = bt2(x2t)
 
#Create the plot
plt.plot(x1t, y1t, color='#AF0040')
plt.plot(x2t, y2t, color='#1E88E5')

#Add labels and title (optional)
plt.xlabel("r")
plt.ylabel("alpha")
plt.title("tradeoff between r and alpha (paired)")
plt.grid(True)
plt.show()

#make plot of tradeoff function (between r and K, knowing that K = r/alpha)

y1k,y2k = Kpaired(x1t, x2t, y1t, y2t, -2.5e-05, -2.5e-05)

#Create the plot
plt.plot(x1t, y1k, color='#AF0040')
plt.plot(x2t, y2k, color='#1E88E5')

#Add labels and title (optional)
plt.xlabel("r")
plt.ylabel("K")
plt.title("tradeoff between r and K (paired)")
plt.grid(True)
plt.show()

# =============================================================================
# #plot to see if protected dimorphisms are possible
# plt.figure(figsize=(6, 6))

# # Plot the region where the inequality is satisfied (Z > 0)
# plt.contourf(X, Y, Z, levels=[0, Z.max()], colors=['lightblue'], extend='max')
# plt.contour(X, Y, Z, levels=[0], colors='red', linestyles='--')

# # Plot the inverse of the region where the inequality is satisfied (Z > 0)
# plt.contourf(Y, X, Z, levels=[0, Z.max()], colors=['lightgreen'], extend='max')

# plt.xlabel('resident')
# plt.ylabel('mutant')
# plt.title('Region where mutants can invade resident population (blue), mirrored version (green)')
# plt.grid(True)

# plt.show()
# =============================================================================



#%%

#sanity check: numerically solving one-species fitness function
def ss_eqn1(var):
    ys = var
    xs = ys 
    eps = 1e-6
    
    #finite difference partial derivative
    dsdy1 = (selection_gradient1(xs, ys+eps) - selection_gradient1(xs, ys-eps)) / (2*eps)
    
    return dsdy1

def ss_eqn2(var):
    ys = var
    xs = ys 
    eps = 1e-6
    
    #finite difference partial derivative
    dsdy2 = (selection_gradient2(xs, ys+eps) - selection_gradient2(xs, ys-eps)) / (2*eps)
    
    return dsdy2

solsing1 = fsolve(ss_eqn1, [0.0])
solsing2 = fsolve(ss_eqn2, [0.0])

print(f"Equilibrium for single species: y1 = {solsing1[0]:.6f}, y2 = {solsing2[0]:.6f}")


#numerically solve the two-species fitness functions
def grad_equilibrium(vars):
    y1, y2 = vars
    x1, x2 = y1, y2  # or your fixed parameter values
    eps = 1e-6
    # finite-difference partial derivatives
    ds1_dy1 = (sel_grad_paired_new(x1, y1+eps, x2, y2)[0] - sel_grad_paired_new(x1, y1-eps, x2, y2)[0]) / (2*eps)
    ds2_dy2 = (sel_grad_paired_new(x1, y1, x2, y2+eps)[1] - sel_grad_paired_new(x1, y1, x2, y2-eps)[1]) / (2*eps)
    return [ds1_dy1, ds2_dy2]

sol = fsolve(grad_equilibrium, [0.0, 0.0])
sola = bt1(sol[0]), bt2(sol[1])

print(f"Equilibrium for two species: r1 = {sol[0]:.6f}, r2 = {sol[1]:.6f}")
print(f"Equilibrium for two species: a11 = {sola[0]:.7f}, a22 = {sola[1]:.7f}")



#%%

#solving the equilibrium equation over many possible values of constants in tradeoff function

def td1(b,c1,c2):
    return c1 * np.exp(c2*b) #exponential 

def td2(b,c1,c2):
    return c1 * np.exp(c2*b) #exponential 

def Kpaired2(r1,r2,a11,a22,a12,a21):
    K1 = (r2*a12 - r1*a22)/(a11*a22 - a12*a21)
    K2 = (r1*a21 - r2*a11)/(a11*a22 - a12*a21)
    
    return K1,K2

#single species selection gradient
def sel_grad(x, y, c1, c2):
    return y - (x * td1(y, c1, c2)/td1(x, c1, c2))

# paired selection gradient where both species can have residents and mutants
def sel_grad_paired_var(x1, y1, x2, y2, c1, c2, c3, c4):
    a12 = -2.5e-05
    a21 = -2.5e-05
    den = td2(x2,c3,c4)*td1(x1,c1,c2) - a12*a21
    s1m = y1 + td1(y1,c1,c2)*(a12*x2 - td2(x2,c3,c4)*x1)/den + (a12*a21*x1 - td1(x1,c1,c2)*a12*x2)/den
    s2m = y2 + td2(y2,c3,c4)*(a21*x1 - td1(x1,c1,c2)*x2)/den + (a12*a21*x2 - td2(x2,c3,c4)*a21*x1)/den
    return s1m, s2m

#function that gives predictions of coexistence
def pred_coexist(r1,r2,a11,a22,a12,a21):
    
    out = False
    # compute values of rho and f2/f1
    temp_rho = (a12 / a11) * (a21 / a22)
    temp1_f2f1 = (a11 / a22) * (a12 / a21)
    temp_f2f1 = np.sqrt(abs(temp1_f2f1)) * (r2 / r1)

    rho = np.sqrt(abs(temp_rho))
    f2f1 = 1 / temp_f2f1
    
    if (f2f1 > rho) and (f2f1 < 1 / rho):
        out = True
    
    return out

#solve the one-species fitness function (write same function twice as constants used are different for the two species)
def eq_1sp1(var):
    ys = var
    xs = ys 
    eps = 1e-6
    
    #finite difference partial derivative
    dsdy = (sel_grad(xs, ys+eps, con1, con2) - sel_grad(xs, ys-eps, con1, con2)) / 2*eps
    
    return dsdy

def eq_1sp2(var):
    ys = var
    xs = ys 
    eps = 1e-6
    
    #finite difference partial derivative
    dsdy = (sel_grad(xs, ys+eps, con3, con4) - sel_grad(xs, ys-eps, con3, con4)) / 2*eps
    
    return dsdy

#numerically solve the two-species fitness functions
def eq_2sp(vars):
    y1, y2 = vars
    x1, x2 = y1, y2  # or your fixed parameter values
    eps = 1e-6
    # finite-difference partial derivatives
    ds1_dy1 = (sel_grad_paired_var(x1, y1+eps, x2, y2, con1, con2, con3, con4)[0] - sel_grad_paired_var(x1, y1-eps, x2, y2, con1, con2, con3, con4)[0]) / (2*eps)
    ds2_dy2 = (sel_grad_paired_var(x1, y1, x2, y2+eps, con1, con2, con3, con4)[1] - sel_grad_paired_var(x1, y1, x2, y2-eps, con1, con2, con3, con4)[1]) / (2*eps)
    return [ds1_dy1, ds2_dy2]

c1s_list = np.linspace(-0.000001, -0.00005, 50)
c2s_list = np.linspace(0.5,5,50)
output_writeall = []
a12 = -2.5e-05
a21 = -2.5e-05

#fix two of the constants, vary the others


for i in range(0,len(c1s_list)):
    for j in range(0,len(c2s_list)):
        #define the constants outside the function
        output_temp = []
        con1, con2 = c1s_list[i], c2s_list[j]
        con3, con4 = -0.0000108, 3.5
        rpreds1 = fsolve(eq_1sp1, 0.0)
        rpreds2 = fsolve(eq_1sp2, 0.0)
        rpred1, rpred2 = fsolve(eq_2sp, [0.0, 0.0])
        print(f"Equilibrium for two species: y1 = {rpred1:.6f}, y2 = {rpred2:.6f}")
        apred11 = td1(rpred1, con1, con2)
        apred22 = td2(rpred2, con3, con4)
        apreds1 = td1(rpreds1, con1, con2)
        apreds2 = td1(rpreds2, con3, con4)
        outcomesing = pred_coexist(rpreds1, rpreds2, apreds1, apreds2, a12,a21)
        outcomepair = pred_coexist(rpred1, rpred2, apred11, apred22, a12,a21)
        Neq = Kpaired2(rpred1, rpred2, apred11, apred22, a12,a21)
        #print(outcome, Neq)
        #print()
        output_temp = [con1, con2, con3, con4, rpreds1, rpreds2, apreds1, apreds2, rpred1, rpred2, apred11, apred22, a12, a21, Neq[0], Neq[1], outcomesing, outcomepair]
        output_writeall.append(output_temp)

with open('/Users/kasturilele/Documents/SLiM/tradeoff_fun_ESS1.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f, dialect = 'excel')
    rowNames = ['C1','C2','C3','C4','r1_sing','r2_sing','a1_sing','a2_sing','r1','r2','a11','a22','a12','a21','N1_equilibrium','N2_equilibrium','outcome_single','outcome_pair']
      
    write.writerow(rowNames)
    write.writerows(output_writeall)

#c1s_list = np.linspace(-0.00001, -0.00003, 4)
#c2s_list = np.linspace(3,5,8)

# aijs1 = np.linspace(-1e-5,-1.05e-4,20)
# aijs2 = np.linspace(-1e-5,-1.05e-4,20)
# output_writeall = []

# for i in range(0,len(aijs1)):
#     for j in range(0,len(aijs2)):
#         #define the constants outside the function
#         output_temp = []
#         con1, con2, con3, con4 = -0.0000206, 3.2, -0.0000108, 3.5
#         rpreds1 = fsolve(eq_1sp1, 0.0)
#         rpreds2 = fsolve(eq_1sp2, 0.0)
#         rpred1, rpred2 = fsolve(eq_2sp, [0.0, 0.0])
#         a12 = aijs1[i]
#         a21 = aijs2[j]
#         #print(f"Equilibrium for two species: y1 = {rpred1:.6f}, y2 = {rpred2:.6f}")
#         apred11 = td1(rpred1, con1, con2)
#         apred22 = td2(rpred2, con3, con4)
#         apreds1 = td1(rpreds1, con1, con2)
#         apreds2 = td1(rpreds2, con3, con4)
#         outcomesing = pred_coexist(rpreds1, rpreds2, apreds1, apreds2, a12,a21)
#         outcomepair = pred_coexist(rpred1, rpred2, apred11, apred22, a12,a21)
#         Neq = Kpaired2(rpred1, rpred2, apred11, apred22, a12,a21)
#         #print(outcome, Neq)
#         #print()
#         output_temp = [con1, con2, con3, con4, rpreds1[0], rpreds2[0], apreds1[0], apreds2[0], rpred1, rpred2, apred11, apred22, a12, a21, Neq[0], Neq[1], outcomesing, outcomepair]
#         output_writeall.append(output_temp)

# with open('tradeoff_fun_ESS_aijs.csv', 'w') as f:
      
#     # using csv.writer method from CSV package
#     write = csv.writer(f, dialect = 'excel')
#     rowNames = ['C1','C2','C3','C4','r1_sing','r2_sing','a1_sing','a2_sing','r1','r2','a11','a22','a12','a21','N1_equilibrium','N2_equilibrium','outcome_single','outcome_pair']
      
#     write.writerow(rowNames)
#     write.writerows(output_writeall)

#%%

#testing gamma distributions for SLiM simulation



import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma


# Define the parameters of the Gamma distribution
# 'a' is the shape parameter (alpha or k)
# 'scale' is the scale parameter (theta or 1/beta)
# Note: scipy uses 'a' for shape and 'scale' for scale, where scale = 1/rate (beta)


a = 0.5 # Alpha 
mu = 1.0  # specify mean instead of scale because SLiM does that 

# Generate x values for the plot
# The range should be appropriate for the chosen parameters
x = np.linspace(0, 10, 500)

# Calculate the Probability Density Function (PDF)
pdf = gamma.pdf(x, a=a, scale=mu/a)

# Plot the Gamma distribution
plt.figure(figsize=(8, 5))
plt.plot(x, pdf, 'b-', lw=2, label=f'Gamma PDF (mean={mu}, scale={a})')
plt.title('Gamma Distribution')
plt.xlabel('x')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True)
plt.show()


#%%


