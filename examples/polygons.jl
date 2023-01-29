
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
	
	function split(tmp::Vector,lon0=-160.0)
		[split(a,lon0) for a in tmp]
	end
	
	function split(tmp::LineString,lon0=-160.0)
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
			jj=[0;findall(tmp2.==3)...;np+1]
			[LineString(tmp[jj[ii]+1:jj[ii+1]-1]) for ii in 1:length(jj)-1]
		end
	end

	split(tmp::Vector,dest::Observable) = tmp

	function split(tmp::Vector,dest::String)
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

    using Downloads, GeoJSON, GeoInterface, Shapefile, GeometryBasics

    ## read data from file

    #fil=joinpath(tempdir(),"countries.geojson")
    #url = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/countries.geojson"
    #Downloads.download(url,fil)
    function get_land_geo_json()
        fil=joinpath(tempdir(),"countries.geojson")
        jsonbytes = read(fil)
        GeoJSON.read(jsonbytes)
    end

    #url="https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
    #fil=joinpath(tempdir(),"ne_110m_admin_0_countries.shp")
    #Downloads.download(url,fil)
    function get_land_geo_shp(fil)
        table = Shapefile.Table(fil)
        geoms = Shapefile.shapes(table)
    end

    ## convert to GeometryBasics

    to_point2(a::Vector{<: T}) where T = Point2{T}(a[1], a[2])
    to_point2(a::AbstractVector{T}) where T <: Number = Point2{T}(a[1], a[2])

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

    function download_tmp()
        url="https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
        fil=joinpath(tempdir(),"ne_110m_admin_0_countries.shp")
        !isfile(fil) ? Downloads.download(url,fil) : nothing
        fil
    end
    
    function read_tmp()
        #tmp1=get_land_geo_json()
        #tmp1a=GeoInterface.coordinates(tmp1[1])[1]
        #tmp1b=geo2basic(tmp1a)

        fil=download_tmp()
        tmp2=get_land_geo_shp(fil)
#        tmp2a=GeoInterface.coordinates(tmp2[1])[1]
#        [geo2basic(i) for i in tmp2a]
        tmp2=[GeoInterface.coordinates(a) for a in tmp2]
        tmp22=Vector{Point2{Float64}}[]
        [[[push!(tmp22,geo2basic(l3)) for l3 in l2] for l2 in l1] for l1 in tmp2]
        tmp22

end

end #module PolygonReading