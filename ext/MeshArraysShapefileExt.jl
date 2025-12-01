
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
    function read_shp(fil)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
        [Shapefile.GI.coordinates(a) for a in geoms]
    end

end #module MeshArraysShapefileExt