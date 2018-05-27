

function advance_AB3_3D(Uh::Array, NL::Array; dt::Float64=1)
    #------
    # Calculate derivatives [1,2,:,:] is y-deriv of u-component
    dUhdx = np.stack([Kx.*im.*Uh[end,:,:,:], Ky.*im.*Uh[end,:,:,:]], axis=1)
    #------

    #------
    # Get nonlinear term and calculate increment
    NL[end,:,:,:] = - get_nonlinear(Uh[end,:,:,:], dUhdx)
    NL_aux[1,:,:] = (1 - kxkx_k2).*NL[end,1,:,:] - kxky_k2.*NL[end,2,:,:]
    NL_aux[2,:,:] = - kxky_k2.*NL[end,1,:,:] + (1 - kyky_k2).*NL[end,2,:,:]
    NL[end,:,:,:] = NL_aux
    Q = sum(AB3_coeffs.*NL, 1)[1,:,:,:]
    #------

    #----
    # Shift the positions of previous time steps on last moment
    Uh = np.roll(Uh, -1, axis=0)
    NL = np.roll(NL, -1, axis=0)
    NL[end,:,:,:] = 0.0 # Not needed, but just to avoid mistakes
    #----

    #----
    # Calculate current time based on last
    Uh[end,:,:,:] = E.^(-1).*Uh[end-1,:,:,:] + dt*Q
    #----

    return Uh, NL
end


function get_nonlinear(vel::Array, dveldx::Array)
    NLx =   adv(vel[1,:,:,:], dveldx[1,1,:,:,:], apad, bpad, aphys, bphys, phys) + 
            adv(vel[2,:,:,:], dveldx[1,2,:,:,:], apad, bpad, aphys, bphys, phys) + 
            adv(vel[3,:,:,:], dveldx[1,3,:,:,:], apad, bpad, aphys, bphys, phys) 

    NLy =   adv(vel[1,:,:,:], dveldx[2,1,:,:,:], apad, bpad, aphys, bphys, phys) + 
            adv(vel[2,:,:,:], dveldx[2,2,:,:,:], apad, bpad, aphys, bphys, phys) +
            adv(vel[3,:,:,:], dveldx[2,3,:,:,:], apad, bpad, aphys, bphys, phys) 

    NLz =   adv(vel[1,:,:,:], dveldx[3,1,:,:,:], apad, bpad, aphys, bphys, phys) + 
            adv(vel[2,:,:,:], dveldx[3,2,:,:,:], apad, bpad, aphys, bphys, phys) +
            adv(vel[3,:,:,:], dveldx[3,3,:,:,:], apad, bpad, aphys, bphys, phys) 
    NL = np.stack([NLx, NLy, NLz], axis=0)
    return NL
end


