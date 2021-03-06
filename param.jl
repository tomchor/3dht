import Distributed

#------
# General
const Ndim = 3
Nprocs = Distributed.nprocs()
#------

#------
# Domain setup and and IC
U_dim=.005 # m/s
L_dim=512e-4 # m
N=2^5
ν_dim=1.5e-5 # m^2/s
Prod_dim = .5e-3 # m^2/s^3
#------

#------
L_scale=L_dim # m
U_scale=0.05 # m/s
#------

#------
L=L_dim/L_scale
ν=ν_dim/(L_scale*U_scale)
Prod=Prod_dim*L_scale/U_scale^3
k_peak = (2*π)/(0.8*L) # ~ From Pope
#------

#------
# Time setup
dt = .001
dt_out = 0.1
T_f = 10.
T_f = 0.5
out_T = 0:dt_out:T_f
Nt = trunc(Int,T_f/dt)
out_n = trunc.(Int,out_T./dt)
rms=0.001
#------

IC="iso"

#------
# To convert between energies
const Cf = 0.002943727133966569
#------
