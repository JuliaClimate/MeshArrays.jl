
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

    - (default) call `GeoJSON.read` then convert to `GI.coordinates`.
    - if `native` then skip conversion, just call `GeoJSON.read`.

    ```
    import MeshArrays, DataDeps, GeoJSON
    pol=MeshArrays.Dataset("oceans_geojson1",do_read=false)
    pol=MeshArrays.read_json(pol,native=true)

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
    function read_json(fil; native=false)
        if native
            GeoJSON.read(fil)
        else
            jsonbytes = read(fil)
            geoms = GeoJSON.read(jsonbytes)
            #GeoJSON.GI.coordinates(geoms)
            [GeoJSON.GI.coordinates(a) for a in geoms]
        end
    end

end #module MeshArraysGeoJSONExt