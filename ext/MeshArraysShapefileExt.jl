
############################################################
#                                                          #
#   Reading & Converting Polygons from Shapefile           #
#                                                          #
############################################################

module MeshArraysShapefileExt

    using Shapefile
    import MeshArrays: read_shp

    """
        read_shp(fil)

    - call `Shapefile.Table` then `Shapefile.shapes`.
    - convert to `GI.coordinates`.
    """
    function read_shp(fil; native=false)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
        if native
            np=length(table.NAME)
            [(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
        else
            [Shapefile.GI.coordinates(a) for a in geoms]
        end
    end

end #module MeshArraysShapefileExt