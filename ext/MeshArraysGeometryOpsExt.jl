
############################################################
#                                                          #
#   Geometric Operations on Polygons                       #
#                                                          #
############################################################

module MeshArraysGeometryOpsExt

    import GeometryOps as GO
    import MeshArrays: within_pol, GI

    """
        within_pol(pol; ID=1)

    Generate a `name,rule` pair to test if location 
        `lon,lat` is within polygon `pol[ID].geometry`.

    ```
    using MeshArrays, DataDeps, GeoJSON, GeometryOps
    pol=MeshArrays.Dataset("oceans_geojson1")
    name,rule=MeshArrays.within_pol(pol; ID=11)
    name,rule(-30,40)
    ```
    """
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