
############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON             #
#                                                          #
############################################################

module MeshArraysGeoJSONExt

    using GeoJSON
    import MeshArrays: read_json

    function read_json(fil)
        jsonbytes = read(fil)
        geoms = GeoJSON.read(jsonbytes)
        #GeoJSON.GI.coordinates(geoms)
        [GeoJSON.GI.coordinates(a) for a in geoms]
    end

end #module MeshArraysGeoJSONExt