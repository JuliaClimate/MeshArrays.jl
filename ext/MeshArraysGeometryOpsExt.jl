
############################################################
#                                                          #
#   Geometric Operations on Polygons                       #
#                                                          #
############################################################

module MeshArraysGeometryOpsExt

    import GeometryOps as GO
    import MeshArrays: within_pol, GI

    function within_pol(pol; ID=0)
        name,geom=if ID==0
            "unknown",pol
        else
            pol[ID].name,pol[ID].geometry
        end
        rule = (x,y) -> GO.within(GI.Point(x,y), geom)
        (name,rule)
    end

	import MeshArrays: interpolation_setup, Interpolate
	import MeshArrays: GridSpec, GridLoad, GridLoadVar, land_mask
	
    function get_samples(nsamples=10000)	
        x=-180.0 .+360*rand(nsamples)
        y=-90.0 .+180*rand(nsamples)
        x,y
    end

    """
        example()

    ```
    using MeshArrays, GeometryOps, JLD2, DataDeps
    MAGO = Base.get_extension(MeshArrays, :MeshArraysGeometryOpsExt);
    output = MAGO.example()

    using CairoMakie

    function plot_poly_samples(Dint,pol0,pol1,x,y)
        heatmap(-179.75:0.5:179.75,-89.75:0.5:89.75,Dint,colormap=:reds)
        plot!(pol0,color=:yellow)
        plot!(pol1[1].geometry,color=:black)
        scatter!(x,y,color=:white,markersize=3)
        current_figure()
    end

    fig1=plot_poly_samples(Dint,pol0,pol1,x[ii],y[ii])	
    fig2=heatmap(μ*msk)
    ```
    """
    function example()
        λ=interpolation_setup()
        γ=GridSpec(ID=:LLC90)
        Γ=GridLoad(γ;option="light")
        D=Γ.Depth
        
        Dint=reshape(Interpolate(D,λ.f,λ.i,λ.j,λ.w),size(λ.lon))
        Dint[Dint.==0.0].=NaN

    	pol0=GO.polygonize(-179.75:0.5:179.75,-89.75:0.5:89.75, 5000 .< Dint)

    #	rule = (x,y) -> GO.within(GI.Point(x,y), geom)
        pol1=[(geometry=pol0.geom[1],name="test")]
        name,rule=within_pol(pol1,ID=1)

        msk=0.0 .+ 1.0*rule.(Γ.XC,Γ.YC)
        hFacC=GridLoadVar("hFacC",γ)
        μ=0.0 .+ 1.0*land_mask(hFacC[:,1])

        x,y=get_samples()
        ii=findall(rule.(x,y))
        return x,y,ii,pol0,pol1,Dint,μ*msk
    end

end #module MeshArraysGeometryOpsExt