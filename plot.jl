
function plot_2axes(A, B, name; sym=false)
    vmin=-1
    vmax=1
    A=collect(A)
    B=collect(B)
    fig, axes = plt.subplots(ncols=2, figsize=(10,5))
    ax=axes[1]
    ma=ax[:pcolormesh](A.', vmin=vmin, vmax=vmax)
    plt.colorbar(ma, ax=ax)
    ax=axes[2]
    ma=ax[:pcolormesh](B.', vmin=vmin, vmax=vmax)
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
    ma=ax[:imshow](A.', vmin=-vm, vmax=vm, interpolation="nearest", origin="lower", cmap="seismic")
    plt.colorbar(ma, ax=ax)
    fig[:tight_layout]()
    fig[:savefig](name)
    plt.close("all")
end 

