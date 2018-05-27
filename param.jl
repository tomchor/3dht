#------
# General
const Ndim = 3
Nprocs = nprocs()
#------

#------
# Domain setup and and IC
L=1 # m
N=2^5
ν=1.5e-5 # m^2/s
#ν=0.0
#------
#------
# Time setup
dt = .0002
T_f = 5.
out_T = 0:.1:T_f
Nt = trunc(Int,T_f/dt)
out_n = trunc.(Int,out_T./dt)
rms=0.01
#------

IC="jet"

