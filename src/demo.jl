
module demo

    import MeshArrays: read_JLD2, write_JLD2, Transect, rotate_points, rotate_XCYC
    import MeshArrays: edge_mask, MskToTab, shorter_paths!
    import MeshArrays: GRID_LLC90, GridSpec, MeshArray
    import MeshArrays: download_file, unzip

    """
        ocean_sections()

    Ocean sections are defined by longitude, latitude pairs.

    ```
    sec=demo.ocean_sections()
    ```
    """
    function ocean_sections()
        lonPairs=[]    
        latPairs=[]    
        namPairs=[]    
        
        push!(lonPairs,[-173 -164]); push!(latPairs,[65.5 65.5]); push!(namPairs,"Bering Strait");
        push!(lonPairs,[-5 -5]); push!(latPairs,[34 40]); push!(namPairs,"Gibraltar");
        push!(lonPairs,[-81 -77]); push!(latPairs,[28 26]); push!(namPairs,"Florida Strait");
        push!(lonPairs,[-81 -79]); push!(latPairs,[28 22]); push!(namPairs,"Florida Strait W1");
        push!(lonPairs,[-76 -76]); push!(latPairs,[21 8]); push!(namPairs,"Florida Strait S1");
        push!(lonPairs,[-77 -77]); push!(latPairs,[26 24]); push!(namPairs,"Florida Strait E1");
        push!(lonPairs,[-77 -77]); push!(latPairs,[24 22]); push!(namPairs,"Florida Strait E2");
        push!(lonPairs,[-65 -50]); push!(latPairs,[66 66]); push!(namPairs,"Davis Strait");
        push!(lonPairs,[-35 -20]); push!(latPairs,[67 65]); push!(namPairs,"Denmark Strait");
        push!(lonPairs,[-16 -7]); push!(latPairs,[65 62.5]); push!(namPairs,"Iceland Faroe");
        push!(lonPairs,[-6.5 -4]); push!(latPairs,[62.5 57]); push!(namPairs,"Faroe Scotland");
        push!(lonPairs,[-4 8]); push!(latPairs,[57 62]); push!(namPairs,"Scotland Norway");
        push!(lonPairs,[-68 -63]); push!(latPairs,[-54 -66]); push!(namPairs,"Drake Passage");
        push!(lonPairs,[103 103]); push!(latPairs,[4 -1]); push!(namPairs,"Indonesia W1");
        push!(lonPairs,[104 109]); push!(latPairs,[-3 -8]); push!(namPairs,"Indonesia W2");
        push!(lonPairs,[113 118]); push!(latPairs,[-8.5 -8.5]); push!(namPairs,"Indonesia W3");
        push!(lonPairs,[118 127 ]); push!(latPairs,[-8.5 -15]); push!(namPairs,"Indonesia W4");
        push!(lonPairs,[127 127]); push!(latPairs,[-25 -68]); push!(namPairs,"Australia Antarctica");
        push!(lonPairs,[38 46]); push!(latPairs,[-10 -22]); push!(namPairs,"Madagascar Channel");
        push!(lonPairs,[46 46]); push!(latPairs,[-22 -69]); push!(namPairs,"Madagascar Antarctica");
        push!(lonPairs,[20 20]); push!(latPairs,[-30 -69.5]); push!(namPairs,"South Africa Antarctica");
        push!(lonPairs,[-76 -72]); push!(latPairs,[21 18.5]); push!(namPairs,"Florida Strait E3");
        push!(lonPairs,[-72 -72]); push!(latPairs,[18.5 10]); push!(namPairs,"Florida Strait E4");

        (name=namPairs,lon=lonPairs,lat=latPairs)
    end

    """
        ocean_sections(Γ)

    Ocean sections are defined by longitude, latitude pairs.

    ```
    sec,pth_sec=demo.ocean_sections(Γ)
    ```
    """
    function ocean_sections(Γ)

        sec=ocean_sections()
        pth_sec=joinpath(tempdir(),"MeshArray_sections")
        !isdir(pth_sec) ? mkdir(pth_sec) : nothing
        
        for ii in 1:length(sec.name)
            lons=Float64.(sec.lon[ii])
            lats=Float64.(sec.lat[ii])
            name=sec.name[ii]
            Trsct=Transect(name,lons,lats,Γ)
            fil_Trsct=joinpath(pth_sec,"$(Trsct.name).jld2")
            !isfile(fil_Trsct) ? write_JLD2(fil_Trsct,tabC=Trsct.tabC,tabW=Trsct.tabW,tabS=Trsct.tabS) : nothing
        end

        sec,pth_sec
    end

    """
        one_section(Γ,lons,lats)

    One section gets defined by longitude pair, and latitude pair.

    ```
    my_section=demo.one_section(Γ,lons,lats)
    ```
    """
    function one_section(Γ,lons,lats)
        name="An Example Grid Path"
    
        x0,y0,z0,R=rotate_points(lons,lats)
        x,y,z=rotate_XCYC(Γ,R)
        mskCint=1.0*(z.>0)
        mskCedge,mskWedge,mskSedge=edge_mask(mskCint)
        
        #for plotting
        tabC=MskToTab(mskCedge)
        locCwhole=( lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
                    lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])
    
        mskCedge,mskWedge,mskSedge=shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))
    
        #for plotting
        tabC=MskToTab(mskCedge)
        
        ( 	lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
            lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])
    end

    ##

    """
        ocean_basins()

    Map of ocean basins on ECCO grid.

    ```
    basins=demo.ocean_basins()
    ```
    """
    function ocean_basins()        
        γ=GridSpec("LatLonCap",GRID_LLC90)
        fil_basin=joinpath(GRID_LLC90,"v4_basin.bin")
        basins=(map=read(fil_basin,MeshArray(γ,Float32)),
            name=["Pacific","Atlantic","indian","Arctic","Bering Sea","South China Sea","Gulf of Mexico",
            "Okhotsk Sea","Hudson Bay","Mediterranean Sea","Java Sea","North Sea","Japan Sea",
            "Timor Sea","East China Sea","Red Sea","Gulf","Baffin Bay","GIN Seas","Barents Sea"])
    end

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
        !isfile(fil) ? download_file(url,fil) : nothing
        if !isempty(unzipfil)
            unzip(fil)
            fil=unzipfil
        end
        fil
    end

    """
        get_basemap()

        Download (if needed) and read "Blue_Marble.jpg".
    """
    function get_basemap()
        dx=0.1
        lat=[j for i=-0.05:dx:359.95, j=-89.95:dx:89.95]; 
        lon=[i for i=-0.05:dx:359.95, j=-89.95:dx:89.95];
    
        earth_jpg=joinpath(tempdir(),"Blue_Marble.jpg")
        url="https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"
        !isfile(earth_jpg) ? MeshArrays.download_file(url,earth_jpg) : nothing
    
        earth_img=read_JLD2(earth_jpg)
        earth_img=reverse(permutedims(earth_img),dims=2)
        earth_img=circshift(earth_img,(1800,0))
    
        lon,lat,earth_img
    end
   

end
