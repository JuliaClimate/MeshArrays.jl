
############################################################
#                                                          #
#     Splitting Of LineString / coast lines                #
#                                                          #
############################################################

#`LineSplitting` module defines a `split` method for vectors of `LineString` objects.
#
#This is needed to fix e.g. coast line displays when lon_0 is not 0 but cutting polygons at lon_0+-180.

module LineSplitting
	import GeometryBasics.LineString
	import Observables.Observable
	import Base.split
	
	function regroup(tmp::Vector)
		coastlines_custom=LineString[]
		println(typeof(coastlines_custom))
		for ii in 1:length(tmp)
			push!(coastlines_custom,tmp[ii][:]...)
		end
		coastlines_custom
	end
	
	function split(tmp::Vector{<:LineString},lon0::Real)
		[split(a,lon0) for a in tmp]
	end
	
    #function split(tmp::Vector{Vector{Float64}},lon0::Float64)
    function split(tmp::LineString,lon0::Real)
        lon0<0.0 ? lon1=lon0+180 : lon1=lon0-180 
		np=length(tmp)
		tmp2=fill(0,np)
		for p in 1:np
			tmp1=tmp[p]
			tmp2[p]=maximum( [(tmp1[1][1]<=lon1)+2*(tmp1[2][1]>=lon1) , (tmp1[2][1]<=lon1)+2*(tmp1[1][1]>=lon1)] )
		end
		if sum(tmp2.==3)==0
			[tmp]
		else
            #println("splitting here")
			jj=[0;findall(tmp2.==3)...;np+1]
			[LineString(tmp[jj[ii]+1:jj[ii+1]-1]) for ii in 1:length(jj)-1]
		end
	end

	split(tmp::Vector{<:LineString},dest::Observable) = tmp

	function split(tmp::Vector{<:LineString},dest::String)
		if occursin("+lon_0",dest)
			tmp1=split(dest)
			tmp2=findall(occursin.(Ref("+lon_0"),tmp1))[1]
			lon_0=parse(Float64,split(tmp1[tmp2],"=")[2])
			regroup(split(tmp,lon_0))
		else
			tmp
		end
	end

end

############################################################
#                                                          #
#   Reading & Converting Polygons from GeoJSON/Shapefile   #
#                                                          #
############################################################

#`LineSplitting` module defines a `split` method for vectors of `LineString` objects.
#
#This is needed to fix e.g. coast line displays when lon_0 is not 0 but cutting polygons at lon_0+-180.

module PolygonReading

    using Downloads, ZipFile, GeoJSON, Shapefile
    using GeoInterface, GeometryBasics

    ## read data from file

    function read_json(fil)
        jsonbytes = read(fil)
        GeoJSON.read(jsonbytes)
    end

    function read_shp(fil)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
    end

    function process_json(ID="countries.geojson")
        fil=download_data_if_needed(ID)
        tmp2=read_json(fil)
        #tmp2=GeoInterface.coordinates(tmp2)
        tmp2=[GeoInterface.coordinates(a) for a in tmp2]

        tmp22=Vector{Point2{Float64}}[]
        for l1 in tmp2
            if isa(l1[1][1],Vector{Float64})
                push!(tmp22,geo2basic(l1[1]))
            else
                for l2 in l1
                    push!(tmp22,geo2basic(l2[1]))
                end
            end
        end
        tmp22
    end

    function process_shp(ID="ne_110m_admin_0_countries.shp")
        fil=download_data_if_needed(ID)
        tmp2=read_shp(fil)
        tmp2=[GeoInterface.coordinates(a) for a in tmp2]

        tmp22=Vector{Point2{Float64}}[]
        [[[push!(tmp22,geo2basic(l3)) for l3 in l2] for l2 in l1] for l1 in tmp2]
        tmp22
    end

    ## convert to GeometryBasics

    to_point2(a::Vector{<: T}) where T = Point2{T}(a[1], a[2])
    to_point2(a::AbstractVector{T}) where T <: Number = Point2{T}(a[1], a[2])

    """
        function geo2basic(vector::AbstractVector{<:AbstractVector})

    Source : @SimonDanisch , https://github.com/MakieOrg/GeoMakie.jl/pull/125
    """
    function geo2basic(vector::AbstractVector{<:AbstractVector})
        if isempty(vector)
            return Point{2, Float64}[]
        else
            # GeoJSON strips the eltype so we need to inspect the elements
            x = first(vector)
            if x isa AbstractVector && length(x) == 2 && x[1] isa Real
                return to_point2.(vector)
            elseif x isa AbstractVector && eltype(x) <: AbstractVector
                linestrings = map(x-> to_point2.(x), vector)
                return GeometryBasics.Polygon(linestrings[1], linestrings[2:end])
            else
                error("Unsupported eltype: $(x)")
            end
        end
    end

    # Download data if needed

    function download_data_if_needed(ID::String)
        pth=joinpath(tempdir(),"MeshArrays_Polygons")
        !ispath(pth) ? mkdir(pth) : nothing
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
    
    """
    function unzip(file,exdir="")
        
    Source : @sylvaticus, https://discourse.julialang.org/t/how-to-extract-a-file-in-a-zip-archive-without-using-os-specific-tools/34585/5
    """
    function unzip(file,exdir="")
        fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
        basePath = dirname(fileFullPath)
        outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(),exdir)))
        isdir(outPath) ? "" : mkdir(outPath)
        zarchive = ZipFile.Reader(fileFullPath)
        for f in zarchive.files
            fullFilePath = joinpath(outPath,f.name)
            if (endswith(f.name,"/") || endswith(f.name,"\\"))
                mkdir(fullFilePath)
            else
                write(fullFilePath, read(f))
            end
        end
        close(zarchive)
    end

end #module PolygonReading