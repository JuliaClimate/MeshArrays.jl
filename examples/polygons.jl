
using Downloads, GeoJSON, GeoInterface, Shapefile, GLMakie

## read data from file

#fil=joinpath(tempdir(),"countries.geojson")
#url = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/countries.geojson"
#Downloads.download(url,fil)
function get_land_geo_json()
    fil=joinpath(tempdir(),"countries.geojson")
    jsonbytes = read(fil)
    GeoJSON.read(jsonbytes)
end

tmp1=get_land_geo_json()
tmp1a=GeoInterface.coordinates(tmp1[1])[1]

#url="https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
#fil=joinpath(tempdir(),"ne_110m_admin_0_countries.shp")
#Downloads.download(url,fil)
function get_land_geo_shp()
    fil=joinpath(tempdir(),"ne_110m_admin_0_countries.shp")
    table = Shapefile.Table(fil)
    geoms = Shapefile.shapes(table)
end

tmp2=get_land_geo_shp()
tmp2a=GeoInterface.coordinates(tmp2[1])[1][1]

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

tmp1b=geo2basic(tmp1a)
tmp2b=geo2basic(tmp2a)

## Plot with Makie

lines(tmp1b,color=:black)
lines!(tmp2b,color=:red)
