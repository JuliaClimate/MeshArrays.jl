
############################################################
#                                                          #
#   Reading & Converting Polygons from Shapefile           #
#                                                          #
############################################################

module MeshArraysShapefileExt

    using Shapefile
    import MeshArrays: read_shp, to_polyarray

    """
        read_shp(fil; format=1)

    Call `Shapefile.Table` and `Shapefile.shapes`
    and return `polyarray` (default).
    
    - format 1 (default) : return `polyarray`
    - format 0.2 : return vector of `name,geometry`` named tuples.
    - format 0.1 : convert to `GI.coordinates`.

    ```
    using MeshArrays, DataDeps, Shapefile
    fil=MeshArrays.Dataset("countries_shp1",do_read=false)
    pol=MeshArrays.read_shp(fil,format=0.2)

    using CairoMakie
    MeshArraysMakieExt = Base.get_extension(MeshArrays, :MeshArraysMakieExt)
    lines(MeshArraysMakieExt.pol_to_Makie(pol))
    ```
    """
    function read_shp(fil; format=1)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
        if format==1
            np=length(table.NAME)
            tmp1=[(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
            to_polyarray(tmp1)
        elseif format==0.2
            np=length(table.NAME)
            [(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
        elseif format==0.1
            [Shapefile.GI.coordinates(a) for a in geoms]
        else
            error("unknown format")
        end
    end

end #module MeshArraysShapefileExt