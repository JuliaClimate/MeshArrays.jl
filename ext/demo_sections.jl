
module demo_sections

    ## the following does not use JLD2; should move to src?

    using MeshArrays
    import MeshArrays: ocean_sections, one_section, interpolation_setup

    """
        ocean_sections()

    Ocean sections are defined by longitude, latitude pairs.

    ```
    sec=ocean_sections()
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
    sec,pth_sec=ocean_sections(Γ)
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
            !isfile(fil_Trsct) ? MeshArrays.write_JLD2(fil_Trsct,tabC=Trsct.tabC,tabW=Trsct.tabW,tabS=Trsct.tabS) : nothing
        end

        sec,pth_sec
    end

    """
        one_section(Γ,lons,lats)

    One section gets defined by longitude pair, and latitude pair.

    ```
    my_section=one_section(Γ,lons,lats)
    ```
    """
    function one_section(Γ,lons,lats)
        name="An Example Grid Path"
    
        x0,y0,z0,R=MeshArrays.rotate_points(lons,lats)
        x,y,z=MeshArrays.rotate_XCYC(Γ,R)
        mskCint=1.0*(z.>0)
        mskCedge,mskWedge,mskSedge=MeshArrays.edge_mask(mskCint)
        
        #for plotting
        tabC=MeshArrays.MskToTab(mskCedge)
        locCwhole=( lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
                    lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])
    
        mskCedge,mskWedge,mskSedge=MeshArrays.shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))
    
        #for plotting
        tabC=MeshArrays.MskToTab(mskCedge)
        
        ( 	lon=[Γ.XC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)],
            lat=[Γ.YC[tabC[i,1]][tabC[i,2],tabC[i,3]] for i in 1:size(tabC,1)])
    end

    function interpolation_setup(fil::String)
        λ = MeshArrays.read_JLD2(fil)
        λ = MeshArrays.Dict_to_NamedTuple(λ)
    end

end
