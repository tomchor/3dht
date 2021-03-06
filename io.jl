using PyPlot; plt = PyPlot
plt.ioff()

function plot_3axes(A, B, C, name; sym=false, vm=1)
    vmin=-vm
    vmax=vm
    A=collect(A)
    B=collect(B)
    C=collect(C)
    fig, axes = plt.subplots(ncols=3, figsize=(15,5))
    #ds["U"][:isel](y=1)[:plot](ax=axes[1])
    ax=axes[1]
    ma=ax[:pcolormesh](Array(transpose(A)), vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    ax=axes[2]
    ma=ax[:pcolormesh](Array(transpose(B)), vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    ax=axes[3]
    ma=ax[:pcolormesh](Array(transpose(C)), vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    fig[:tight_layout]()
    fig[:savefig](name)
    plt.close("all")
end 


function plot_2axes(A, B, name; sym=false)
    vmin=-1
    vmax=1
    A=collect(A)
    B=collect(B)
    fig, axes = plt.subplots(ncols=2, figsize=(10,5))
    ax=axes[1]
    ma=ax[:pcolormesh](Array(transpose(A)), vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    ax=axes[2]
    ma=ax[:pcolormesh](Array(transpose(B)), vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    fig[:tight_layout]()
    fig[:savefig](name)
    plt.close("all")
end 


function plot_1axis(A, name; vm=1)
    println("plotting")
    A=collect(A)
    fig, axes = plt.subplots(ncols=1, figsize=(6,6))
    ax=axes
    ma=ax[:imshow](Array(transpose(A)), vmin=-vm, vmax=vm, interpolation="nearest", origin="lower", cmap="seismic")
    plt.colorbar(ma, ax=ax)
    fig[:tight_layout]()
    fig[:savefig](name)
    plt.close("all")
end 


function write_netcdfs(n)
    DS["U"][:values]=U[1,:,:,:]; DS["V"][:values]=U[2,:,:,:]; DS["W"][:values]=U[3,:,:,:];
    (DS*U_scale)[:transpose]()[:to_netcdf](join(["output/uvw_",Printf.@sprintf("%06.0f", 1e2*out_T[n]),".nc"]))
    DSf["Uh"][:values]=Uh[end,1,:,:,:]; DSf["Vh"][:values]=Uh[end,2,:,:,:]; DSf["Wh"][:values]=Uh[end,3,:,:,:];
    (DSf*U_scale)[:transpose]()[:to_netcdf](join(["output/uvwh_",Printf.@sprintf("%06.0f", 1e2*out_T[n]),".nc"]))
end

