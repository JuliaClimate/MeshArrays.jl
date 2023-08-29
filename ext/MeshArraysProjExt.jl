module MeshArraysProjExt

using Proj

"""
    Proj.Transformation(;MA_preset=1,lon0=0.0)

Presets for Proj in MeshArrays.    
"""
function Proj.Transformation(;MA_preset=1,lon0=0.0)
    source="+proj=longlat +datum=WGS84"
    if MA_preset==2 
        dest="+proj=eqearth +lon_0=$(lon0) +lat_1=0.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80"
    elseif MA_preset==1
        dest="+proj=wintri +lon_0=$(lon0) "
    elseif MA_preset==3
        dest="+proj=longlat +datum=WGS84 +lon_0=$(lon0)"
    end
    Proj.Transformation(source,dest, always_xy=true) 
end

end

