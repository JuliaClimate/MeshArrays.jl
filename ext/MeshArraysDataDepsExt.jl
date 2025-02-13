
module MeshArraysDataDepsExt

   using DataDeps, MeshArrays
   import MeshArrays: mydatadep
  

   __init__() = begin
      register(DataDep("countries_shp1",
         "Shapefile countries example, from https://www.naturalearthdata.com",
         ["https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"],
         "0f243aeac8ac6cf26f0417285b0bd33ac47f1b5bdb719fd3e0df37d03ea37110",post_fetch_method=unpack))
      register(DataDep("countries_geojson1",
         "GeoJSON countries example, from http://www.publicamundi.eu/",
         ["https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/countries.geojson"],
         "7ff2917a7a8b38e6392b0a27183d8f28ae9db1d81679526413a29b44fb7c3830"))
      register(DataDep("basemap_jpg1",
         "JPG basemap example, from https://upload.wikimedia.org/wikipedia/commons",
         ["https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"],
         "062582f953750c05c12514b9ff498b8db5b2c68e1e0563b0bea2f3385be9f5db"))
      register(DataDep("interp_halfdeg",
         "Interpolation coefficients example, from LLC90 grid to half degree grid, created with MeshArrays.jl",
         ["https://zenodo.org/record/5784905/files/interp_coeffs_halfdeg.jld2"],
         "a0e2806ad2a88c0acb5d3baef942fcbcf69fe27310487894cd38f6c782f55e99"))
      end

 
   """
       mydatadep(nam="countries_shp1")
       
   Download data dependency with predefined name; currently :

   ```
   "countries_shp1"
   "countries_geojson1"
   "basemap_jpg1"
   "interp_halfdeg"
   ```    
   """
   mydatadep(nam="countries_shp1") = begin
      withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
         if nam=="countries_shp1"
            datadep"countries_shp1"
         elseif nam=="countries_geojson1"
            datadep"countries_geojson1"
         elseif nam=="basemap_jpg1"
            datadep"basemap_jpg1"
         elseif nam=="interp_halfdeg"
            datadep"interp_halfdeg"
         else
            error("unknown data dependency")
         end
      end
   end

end
