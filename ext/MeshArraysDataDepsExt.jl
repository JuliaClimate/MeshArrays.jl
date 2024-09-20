
module MeshArraysDataDepsExt

   using DataDeps, MeshArrays
   import MeshArrays: mydatadep
  

   __init__() = begin
      register(DataDep("countries_shp1",
         "Shapefile countries example, from https://www.naturalearthdata.com",
         ["https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"],
         post_fetch_method=unpack))
      register(DataDep("countries_geojson1",
         "GeoJSON countries example, from http://www.publicamundi.eu/",
         ["https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/countries.geojson"]))
      register(DataDep("basemap_jpg1",
         "JPG basemap example, from https://upload.wikimedia.org/wikipedia/commons",
         ["https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"]))
      register(DataDep("interp_halfdeg",
         "Interpolation coefficients example, from LLC90 grid to half degree grid, created with MeshArrays.jl",
         ["https://zenodo.org/record/5784905/files/interp_coeffs_halfdeg.jld2"]))
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
