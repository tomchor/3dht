

function advance_AB3_3D(Uh::Array, NL::Array; dt::Float64=1)
    #------
    # Calculate derivatives [1,2,:,:] is y-deriv of u-component
    dUhdx[:,1,:,:,:] .= Kx.*im.*Uh[end,:,:,:,:]; 
    dUhdx[:,2,:,:,:] .= Ky.*im.*Uh[end,:,:,:,:]; 
    dUhdx[:,3,:,:,:] .= Kz.*im.*Uh[end,:,:,:,:];
    #------

    #------
    # Get nonlinear term and calculate increment
    #NL[end,:,:,:,:] = -get_nonlinear(Uh[end,:,:,:,:], dUhdx) # Improve this
    NL[end,:,:,:,:] .= -get_nonlinear_aliased(Uh[end,:,:,:,:], dUhdx) # Improve this
    NL_aux[1,:,:,:] .= (1 .- kxkx_k2).*NL[end,1,:,:,:] .-       kxky_k2.*NL[end,2,:,:,:]       .- kxkz_k2.*NL[end,3,:,:,:]
    NL_aux[2,:,:,:] .=     - kxky_k2.*NL[end,1,:,:,:] .+ (1 .- kyky_k2).*NL[end,2,:,:,:]       .- kykz_k2.*NL[end,3,:,:,:]
    NL_aux[3,:,:,:] .=     - kxkz_k2.*NL[end,1,:,:,:]       .- kykz_k2.*NL[end,2,:,:,:] .+ (1 .- kzkz_k2).*NL[end,3,:,:,:]
    NL[end,:,:,:,:] .= NL_aux
    KE = Cf*sum(f_k.*Uh.*conj(Uh), dims=(2,3,4,5))*(dkx^Ndim)/(N^6)
    Q = sum(AB3_coeffs.*(NL .+ Prod*f_k.*Uh./(2*KE)), dims=1)[1,:,:,:,:]
    #------

    #----
    # Shift the positions of previous time steps on last moment
    Uh = np.roll(Uh, -1, axis=0)
    NL = np.roll(NL, -1, axis=0)
    NL[end,:,:,:,:] .= 0.0 # Not needed, but just to avoid mistakes
    #----

    #----
    # Calculate current time based on last
    Uh[end,:,:,:,:] .= E.^(-1).*Uh[end-1,:,:,:,:] .+ dt*Q
    #----

    return Uh, NL
end


function get_nonlinear(vel::Array, dveldx::Array)
    NL_aux2[1,:,:,:] .= adv(vel[1,:,:,:], dveldx[1,1,:,:,:], apad, bpad, aphys, bphys, phys) .+ 
                        adv(vel[2,:,:,:], dveldx[1,2,:,:,:], apad, bpad, aphys, bphys, phys) .+ 
                        adv(vel[3,:,:,:], dveldx[1,3,:,:,:], apad, bpad, aphys, bphys, phys) 

    NL_aux2[2,:,:,:] .= adv(vel[1,:,:,:], dveldx[2,1,:,:,:], apad, bpad, aphys, bphys, phys) .+ 
                        adv(vel[2,:,:,:], dveldx[2,2,:,:,:], apad, bpad, aphys, bphys, phys) .+
                        adv(vel[3,:,:,:], dveldx[2,3,:,:,:], apad, bpad, aphys, bphys, phys) 

    NL_aux2[3,:,:,:] .= adv(vel[1,:,:,:], dveldx[3,1,:,:,:], apad, bpad, aphys, bphys, phys) .+ 
                        adv(vel[2,:,:,:], dveldx[3,2,:,:,:], apad, bpad, aphys, bphys, phys) .+
                        adv(vel[3,:,:,:], dveldx[3,3,:,:,:], apad, bpad, aphys, bphys, phys) 
    return NL_aux2
end

function get_nonlinear_aliased(vel::Array, dveldx::Array)
    NL_aux2[1,:,:,:] .= adv_aliased(vel[1,:,:,:], dveldx[1,1,:,:,:]) .+ 
                        adv_aliased(vel[2,:,:,:], dveldx[1,2,:,:,:]) .+ 
                        adv_aliased(vel[3,:,:,:], dveldx[1,3,:,:,:]) 

    NL_aux2[2,:,:,:] .= adv_aliased(vel[1,:,:,:], dveldx[2,1,:,:,:]) .+ 
                        adv_aliased(vel[2,:,:,:], dveldx[2,2,:,:,:]) .+
                        adv_aliased(vel[3,:,:,:], dveldx[2,3,:,:,:]) 

    NL_aux2[3,:,:,:] .= adv_aliased(vel[1,:,:,:], dveldx[3,1,:,:,:]) .+ 
                        adv_aliased(vel[2,:,:,:], dveldx[3,2,:,:,:]) .+
                        adv_aliased(vel[3,:,:,:], dveldx[3,3,:,:,:]) 
    return NL_aux2
end


