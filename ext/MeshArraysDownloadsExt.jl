module MeshArraysDownloadsExt

    using Downloads

    import MeshArrays: download_polygons, interpolation_setup, unzip

    """
        download_polygons(ID::String)

    Download data and unzip (if needed) to `tempdir()`. Only works for predefined `ID`:

    - `ne_110m_admin_0_countries.shp`
    - `countries.geojson`
    """
    function download_polygons(ID::String)
        pth=tempdir()
        unzipfil="" #if provided then need to unzip + return this file name
        if ID=="ne_110m_admin_0_countries.shp"
            url="https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
            fil=joinpath(pth,"ne_110m_admin_0_countries.zip")
            unzipfil=joinpath(pth,"ne_110m_admin_0_countries.shp")
        elseif ID=="countries.geojson"
            fil=joinpath(pth,"countries.geojson")
            url = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/countries.geojson"
        else
            println("unknown file")
            fil="unknown"
            url="unknown"
        end
        !isfile(fil) ? Downloads.download(url,fil) : nothing
        if !isempty(unzipfil)
            unzip(fil)
            fil=unzipfil
        end
        fil
    end
    
    ##

    function interpolation_setup(;path=tempdir())
        fil=joinpath(tempdir(),"interp_coeffs_halfdeg.jld2")
        if !isfile(fil)
            url="https://zenodo.org/record/5784905/files/interp_coeffs_halfdeg.jld2"
            Downloads.download(url, fil)
        end
        fil
    end
    
end


