include("aux_mod.jl")

U0 = zeros(Ndim, Nx, Ny)
if IC=="jet"

    if Ndim==2
        u0 = [ zonal_jet(y; wid=.25, a=2*π) for x in x_center, y in y_center ]
        v0 = [ 0.0 for x in x_center, y in y_center ]

        U0[1,:,:] = u0
        U0[2,:,:] = v0
    else
        println("IC not coded for 3D")
    end

    #----
    # Add noise
    U0+=randn(size(U0))*rms
    #----

    #------
    # Get 
    Uh0 = rplan*U0
    #------

elseif IC=="blob"
    if Ndim==2
        u0 = [ gaussian(x, y) for x in x_center, y in y_center ]
        v0 = [ gaussian(x, y) for x in x_center, y in y_center ]

        U0[1,:,:] = u0
        U0[2,:,:] = v0
    else
        println("IC not coded for 3D")
    end

    #----
    # Add noise
    U0+=randn(size(U0))*0.0
    #----

    #------
    # Get 
    Uh0 = rplan*U0
    #------

elseif IC=="sine"
    if Ndim==2
        u0 = [ sin(2*π*x/Lx) for x in x_center, y in y_center ]
        v0 = [ 0.0 for x in x_center, y in y_center ]

        U0[1,:,:] = u0
        U0[2,:,:] = v0
    else
        println("IC not coded for 3D")
    end

    #----
    # Add noise
    U0+=randn(size(U0))*0.0
    #----

    #------
    # Get 
    Uh0 = rplan*U0
    #------

else
    println("IC not coded")
end


Uh0 = rm_div(Uh0)
println("MAX(∇⋅U) = ", maximum(∇(Uh0)))
