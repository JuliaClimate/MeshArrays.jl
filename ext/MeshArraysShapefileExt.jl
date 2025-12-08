
############################################################
#                                                          #
#   Reading & Converting Polygons from Shapefile           #
#                                                          #
############################################################

module MeshArraysShapefileExt

    using Shapefile
    import MeshArrays: read_shp, to_polyarray

    function read_shp(fil; format=:polyarray)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
        if format==:polyarray
            np=length(table.NAME)
            tmp1=[(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
            to_polyarray(tmp1)
        elseif format==:Shapefile
            np=length(table.NAME)
            [(name=table.NAME[i],geometry=geoms[i]) for i in 1:np]
        elseif format==:coord
            [Shapefile.GI.coordinates(a) for a in geoms]
        else
            error("unknown format")
        end
    end

end #module MeshArraysShapefileExt