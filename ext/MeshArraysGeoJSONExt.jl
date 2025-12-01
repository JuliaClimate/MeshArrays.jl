
############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON             #
#                                                          #
############################################################

module MeshArraysGeoJSONExt

    using GeoJSON
    import MeshArrays: read_json

    """
        read_json(fil; native=false)

    - call `GeoJSON.read`
    - format 1 : return FeatureCollection (vector of name,geom pairs)
    - format 2 : convert to `GI.coordinates`.

    ```
    import MeshArrays, DataDeps, GeoJSON
    pol=MeshArrays.Dataset("oceans_geojson1")

    ii=11
    nam=pol[ii].name
    geom=pol[ii].geometry

    import GeometryOps as GO, GeoInterface as GI
    rule_pol = (x,y) -> GO.within(GI.Point(x,y), geom)
    rule_pol_vec = (x,y) -> rule_pol.(x,y)

    np=1000; lo=-180 .+360*rand(np); la=-90 .+180*rand(np);
    np_in=sum(rule_pol.(lo,la)); [np_in np-np_in nam]
    ```
    """
    function read_json(fil; format=1)
        if format==1
            GeoJSON.read(fil)
        elseif format==2
            jsonbytes = read(fil)
            geoms = GeoJSON.read(jsonbytes)
            #GeoJSON.GI.coordinates(geoms)
            [GeoJSON.GI.coordinates(a) for a in geoms]
        else
            error("unknown format")
        end
    end

end #module MeshArraysGeoJSONExt