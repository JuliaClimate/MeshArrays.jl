

using NCDatasets, GLMakie

file="ocean_daily.static.nc"
ds=Dataset(file)

#list below should be double checked (done quickly, likely to have errors)
variable_pairs=(
(:Coriolis, :Coriolis),
(:areacello, :RAC),
(:areacello_bu, :RAZ),
(:areacello_cu, :RAW),
(:areacello_cv, :RAS),
(:deptho, :Depth),
(:dxCu, :DXW),
(:dxCv, :DXS),
(:dxt, :DXC),
(:dyCu, :DYW),
(:dyCv, :DYS),
(:dyt, :DYC),
(:geolat, :YC),
(:geolat_c, :YC),
(:geolat_u, :YW),
(:geolat_v, :YS),
(:geolon, :XC),
(:geolon_c, :XC),
(:geolon_u, :XW),
(:geolon_v, :XS),
(:sftof, :SeaAreaFraction),
(:wet, :hFacC),
(:wet_c, :hFacC),
(:wet_u, :hFacW),
(:wet_v, :hFacS),
(:xh, :XC1d),
(:xq, :XZ1d),
(:yh, :YC1d),
(:yq, :XZ1d),
)

G=NamedTuple[]

preproc(x) = [(ismissing(i) ? NaN : i) for i in x]
  
for i in 1:length(variable_pairs)
    (n_in,n_out)=variable_pairs[i]
    tmp=try
        MeshArray(preproc(ds[n_in][:,:]))
    catch
#        println("skip $(n_in)")
        ds[n_in][:]
    end
    push!(G,NamedTuple{(n_out,)}((tmp,)))
end

G=merge(G...)
heatmap(G.Depth)

