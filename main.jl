import Printf
import PyCall
@PyCall.pyimport numpy as np
@PyCall.pyimport xarray as xr

include("param.jl")
include("aux_schemes.jl")

println("N = $N")
Nx=Ny=Nz=N
Lx=Ly=Lz=L

#------
# Define grid
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;
x_center = range(dx/2, stop=Lx-dx/2, step=dx)
y_center = range(dy/2, stop=Ly-dy/2, step=dy)
z_center = range(dz/2, stop=Lz-dz/2, step=dz)
println("von Neumann: ", ν*dt*(dx^(-2) + dy^(-2) + dz^(-2)))
#------

#------
# Include FFT functions
include("variables.jl")
include("fft_mod.jl")
#------

#------
# IC
include("IC.jl")
#------

#------
# Start time integration
include("init_calc.jl")
#------

KE0, CFL0 = check_CFL(Uh[end,:,:,:,:], dx=dx, dy=dy, norm=[1, 1], dt=dt, reset=true);
include("io.jl")
function run_sim(Uh, NL)
    #-----
    # Print IC
    plot_3axes(U0[1,:,Int(end//2),:], U0[2,:,Int(end//2),:], U0[3,:,Int(end//2),:], join(["output/uvw_",Printf.@sprintf("%06.02f", 0),".png"]), sym=true, vm=.4)
    #-----

    println("Starting loop ...")
    n=1
    # solver loop
    @time for (jt, t) in enumerate(3*dt:dt:T_f)
        jt_tot=jt-1
        println("jt_tot ", jt_tot, " time= ", t)

        Uh, NL = advance_AB3_3D(Uh, NL, dt=dt)

        if jt_tot in out_n
            U = irplan*Uh[end,:,:,:,:]
            DS["U"][:values] .= U[1,:,:,:]; DS["V"][:values]=U[2,:,:,:]; DS["W"][:values]=U[3,:,:,:];
            (DS*U_scale)[:transpose]()[:to_netcdf](join(["output/uvw_",Printf.@sprintf("%06.0f", 1e2*out_T[n]),".nc"]))
            plot_3axes(U[1,:,Int(end//2),:], U[2,:,Int(end//2),:], U[3,:,1,:], join(["output/uvw_",Printf.@sprintf("%06.02f", out_T[n]),".png"]), sym=true, vm=.4)
            n=n+1
        end

        check_CFL(Uh[end,:,:,:,:], dx=dx, dy=dy, norm=[KE0, 1], dt=dt, reset=false);
        if maximum(abs.(Uh))>2e10
            println("Blew up! Stopped execution at t = ",t)
            break
        end
        #GC.gc(); GC.gc(); GC.gc()
    end
    return Uh, NL, DS
end
