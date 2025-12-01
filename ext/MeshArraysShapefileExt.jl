
############################################################
#                                                          #
#   Reading & Converting Polygons from Shapefile           #
#                                                          #
############################################################

module MeshArraysShapefileExt

    using Shapefile
    import MeshArrays: read_shp

    """
        read_shp(fil; format=1)

    - call `Shapefile.Table` then `Shapefile.shapes`.
    - format 1 : return vector of `name,geometry`` named tuples.
    - format 2 : convert to `GI.coordinates`.

    ```
    using MeshArrays, GeoJSON, DataDeps, Shapefile, CairoMakie
    fil=MeshArrays.Dataset("oceans_geojson1",do_read=false)
    pol=MeshArrays.read_shp(fil)
    ```
    """
    function read_shp(fil; format=1)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
        if format==1
            np=length(table.NAME)
            [(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
        elseif format==2
            [Shapefile.GI.coordinates(a) for a in geoms]
        else
            error("unknown format")
        end
    end

end #module MeshArraysShapefileExt