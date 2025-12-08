
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

end #module MeshArraysGeometryOpsExt