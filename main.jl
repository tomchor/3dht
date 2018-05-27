import PyCall
@PyCall.pyimport numpy as np
@PyCall.pyimport xarray as xr

include("param.jl")
include("aux_schemes.jl")

println("N = $N")
Nx=N
Ny=N
Nz=N
Lx=L
Ly=L
Lz=L

#------
# Define grid
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;
x_center = linspace(dx/2, Lx-dx/2, Nx)
y_center = linspace(dy/2, Ly-dy/2, Ny)
z_center = linspace(dz/2, Lz-dz/2, Nz)
println("von Neumann: ", ν*dt*(dx^(-2) + dy^(-2) + dz^(-2)))
#------

#------
# Include FFT functions
include("variables.jl")
break
include("fft_mod.jl")
#------

#------
# IC
U = zeros(4, Ndim, Nx, Ny)
include("IC.jl")
#------

#------
# Start time integration
include("init_calc.jl")
#------

KE0, CFL0 = check_CFL(Uh[end,:,:,:], dx=dx, dy=dy, norm=[1, 1], dt=dt, reset=true);
function run_sim(Uh, NL)
    #-----
    # Print IC
    plot_1axis(vort(Uh0), join(["figs/zeta_",@sprintf("%06.02f", 0),".png"]), vm=100)
    plot_2axes(U0[1,:,:], U0[2,:,:], join(["figs/uv_",@sprintf("%06.02f", 0),".png"]), sym=true)
    #-----

    println("Starting loop ...")
    n=1
    # solver loop
    tic()
    for (jt, t) in enumerate(3*dt:dt:T_f)
        jt_tot=jt-1
        println("jt_tot ", jt_tot, " time= ", t)

        Uh, NL = advance_AB3_2D(Uh, NL, dt=dt)

        if jt_tot in out_n
            U = irplan*Uh[end,:,:,:]
            plot_1axis(vort(Uh[end,:,:,:]), join(["figs/zeta_",@sprintf("%06.02f", out_T[n]),".png"]), vm=100)
            plot_2axes(U[1,:,:], U[2,:,:], join(["figs/uv_",@sprintf("%06.02f", out_T[n]),".png"]), sym=true)
            n=n+1
        end

        check_CFL(Uh[end,:,:,:], dx=dx, dy=dy, norm=[KE0, 1], dt=dt, reset=false);
        if maximum(abs.(Uh))>2e10
            println("Blew up! Stopped execution at t = ",t)
            break
        end

    end
    toc();
    return Uh, NL
end
