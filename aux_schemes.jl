
function advance_AB3_2D(Uh::Array, NL::Array; dt::Float64=1)
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
    #              U             dUdx                    V                dUdy
    #vel = rplan*ones(U0) # This was used for debugging
    NLx = adv(vel[1,:,:], dveldx[1,1,:,:], apad, bpad) + adv(vel[2,:,:], dveldx[1,2,:,:], apad, bpad) 
    #              U             dVdx                    V                dVdy
    NLy = adv(vel[1,:,:], dveldx[2,1,:,:], apad, bpad) + adv(vel[2,:,:], dveldx[2,2,:,:], apad, bpad) 
    NL = np.stack([NLx, NLy], axis=0)
    return NL
end


