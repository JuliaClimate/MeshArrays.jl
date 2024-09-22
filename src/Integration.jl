
module Integration

using Distributed, SharedArrays, Glob
import MeshArrays: read_JLD2, write_JLD2
import MeshArrays: GridLoad, GridLoadVar, GridSpec
import MeshArrays: demo, MeshArray, gridmask

##

example(;option=:loops,regions=:dlat_10,depths=[(0,7000)]) = begin
  g=GridSpec(ID=:LLC90); Γ=GridLoad(g)
  G=(hFacC=GridLoadVar("hFacC",g),RF=GridLoadVar("RF",g),
     RC=GridLoadVar("RC",g),RAC=GridLoadVar("RAC",g),
     DRF=GridLoadVar("DRF",g))
  G=merge(Γ,G)

  M=define_sums(option=option, regions=regions, grid=G, depths=depths)

  diags=joinpath(pwd(),"diags")
  files=glob("state_3d_set1*.data",diags)

  G,M,files
end

##

DEPTHS=[(0,100),(100,200),(200,300),(300,400),(400,500),(500,600),(600,700),
        (700,800),(800,1000),(1000,1200),(1200,1400),(1400,1700),(1700,2000),
        (2000,2500),(2500,3000),(3000,4000),(4000,5000),(5000,7000)]

lon180(x)=Float64(x>180.0 ? x-360.0 : x)

"""
    define_regions(;option=:global,grid::NamedTuple)

Define regional integration mask (one value period region).
"""
function define_regions(;option=:global,grid::NamedTuple)
 if option==:basins
  demo.ocean_basins()
 elseif option==:global
  mask=1.0*(grid.hFacC[:,1].>0)
  (mask=mask,name=["Global"])
 elseif option==:dlat_10
  mask=1.0*(grid.hFacC[:,1].>0)
  la=grid.YC
  lats=[-90 ; -75:10:75 ; 90]
  nl=length(lats)-1
  name=[Symbol("lat_$(lats[l])_to_$(lats[l+1])") for l in 1:nl]
  [mask[findall((mask.>0)*(la.>=lats[l])*(la.<lats[l+1]))].=l for l in 1:nl]
  (mask=mask,name=name)
 elseif isa(option,Tuple)
  dlo=option[1]; dla=option[2]
  mask=1.0*(grid.hFacC[:,1].>0)

  lo=lon180.(grid.XC)
  lons=collect(-180:dlo:180)
  nlo=length(lons)-1
  la=grid.YC
  lats=[-90 ; -75:dla:75 ; 90]
  nla=length(lats)-1

  name=Symbol[]
  for i_a in 1:nla
   for i_o in 1:nlo
    if sum((mask.>0)*(la.>=lats[i_a])*(la.<lats[i_a+1])*
           (lo.>=lons[i_o])*(lo.<lons[i_o+1]))>0
    t_a="$(lats[i_a])Nto$(lats[i_a+1])N"
    t_o="$(lons[i_o])Eto$(lons[i_o+1])E"
    push!(name,Symbol(t_a*"_"*t_o))
    mask[findall((mask.>0)*(la.>=lats[i_a])*(la.<lats[i_a+1])
	 *(lo.>=lons[i_o])*(lo.<lons[i_o+1]))].=length(name)
    end
   end
  end
  (mask=mask,name=name)
 else
  error("unknown option")
 end
end

layer_mask(dF,d0,d1)=begin
  md=0*dF[1:end-1]
  for k in 1:length(dF)-1
    dF0=abs(dF[k])
    dF1=abs(dF[k+1])
    md[k]=max(0.0,min(d1,dF1)-max(d0,dF0))/(dF1-dF0)
  end
  md
end

"""
    define_sums(;option=:loops, grid::NamedTuple, regions=:global, depths=[(0,7000)])

Define regional integration function for each basin and depth range.
"""
function define_sums(;option=:loops, grid::NamedTuple, regions=:global, depths=[(0,7000)])
  dep=(isa(depths,Tuple) ? [depths] : depths)
  nd=length(dep)
  rgns=define_regions(option=regions,grid=grid) 
  nb=length(rgns.name)
  allones=1.0 .+0*grid.hFacC
  nr=length(grid.RC)

  xymsk(b) = 1.0*(rgns.mask.==b)
  func_h(X,b)=sum(xymsk(b)*X)
  tmp2d=MeshArray(grid.XC.grid,Float32)

  zmsk(d0,d1,k) = layer_mask(grid.RF,d0,d1)[k] 
  func(X,b,d0,d1)=sum([sum(xymsk(b)*zmsk(d0,d1,k)*X[:,k]*
         grid.DRF[k]*grid.hFacC[:,k]*grid.RAC) for k in 1:nr])
  tmp3d=MeshArray(grid.XC.grid,Float32,nr)

  function func_v(X,d0,d1)
    tmp2d.=0.0
    for k in 1:nr
      tmp2d.+=zmsk(d0,d1,k)*X[:,k]*grid.DRF[k]*grid.hFacC[:,k]*grid.RAC
    end
    tmp2d
  end

  if option==:streamlined_loop
  #ocn_surf=[sum(xymsk(b)*(grid.hFacC[:,1].>0)*grid.RAC) for b in 1:nb]
  BX=(name=String[],volsum=Function[],volume=Float64[],
       ocn_surf=Float64[],tmp2d=tmp2d,tmp3d=tmp3d)
  for b in 1:nb
   for d in 1:nd
    (d0,d1)=dep[d]
    n=string(rgns.name[b])*"_dep_$(d0)_to_$(d1)"
    @inline f=X->func(X,b,d0,d1)
    v=f(allones)
    push!(BX.name,n)
    push!(BX.volsum,f)
    push!(BX.volume,v)
    #push!(BX.ocn_surf,ocn_surf[b])
   end
  end
  end

  BXh=(name=String[],hsum=Function[],tmp2d=tmp2d,tmp3d=tmp3d)
  for b in 1:nb
    n=string(rgns.name[b])
    @inline f=X->func_h(X,b)
    push!(BXh.name,n)
    push!(BXh.hsum,f)
  end

  BXv=(name=String[],vint=Function[],tmp2d=tmp2d,tmp3d=tmp3d)
  for d in 1:nd
    (d0,d1)=dep[d]
    n="$(d0)-$(d1)m"
    @inline f=X->func_v(X,d0,d1)
    push!(BXv.name,n)
    push!(BXv.vint,f)
  end

  if option==:loops
    gridmask(rgns.mask,BXh.name,depths,BXh.hsum,BXv.vint,tmp2d,tmp3d)
  elseif option==:streamlined_loop
    gridmask(rgns.mask,BX.name,depths,BX.volsum,[],tmp2d,tmp3d)
  else
    error("unknown option")
  end
end

#nonan(x)=[(isnan(y) ? 0.0 : y) for y in x]

"""
    loops(mask::gridmask; files=String[], var=:THETA, rd=read)

```
begin
  @everywhere using MeshArrays, MITgcm
  @everywhere rd(F,var,tim,tmp)=read(read_mdsio(F,var),tmp)
  @everywhere G,M,files=Integration.example()
    #,regions=(30,10),depths=Integration.DEPTHS)
end;

H=Integration.loops(M,files=files,rd=rd)
# Hbis=Integration.streamlined_loop(M,files=files,rd=rd)
```

and to save results:

```
using JLD2; output_path="test.jld2"
jldsave(output_path; depths=M.depths, integral=H, volume=vol, name=M.names)
```

where vol is calculated as follows:

```
#option=:loops
M.tmp2d.=M.v_int[1](allones)
vol=[b(M.tmp2d) for b in M.h_sum]
```

or 

```
#option=:streamlined_loop
allones=1.0 .+0*G.hFacC
vol=[b(allones) for b in M.h_sum]
```
"""
function loops(mask::gridmask; files=String[], var=:THETA, rd=read)
  nt=length(files)
  nh=length(mask.names)
  nv=length(mask.depths)
  BA=SharedArray{Float64}(nh,nv,nt)
  @sync @distributed for t in 1:nt
    mod(t,10)==0 ? println(t) : nothing
    F=files[t]
    ext=split(F,".")[end]
    mask.tmp3d.=rd(F,var,t,mask.tmp3d)
    for layer in 1:nv
       mask.tmp2d.=mask.v_int[layer](mask.tmp3d)
       BA[:,layer,t]=[b(mask.tmp2d) for b in mask.h_sum]
    end
    GC.gc()
  end
  Array(BA)
end

"""
    streamlined_loop(mask::gridmask; files=String[], var=:THETA, rd=read)

Alternate approach to loops, where loops are streamlined in a single dimension.
"""
function streamlined_loop(mask::gridmask; files=String[], var=:THETA, rd=read)
  nt=length(files)
  nb=length(mask.names)
  BA=SharedArray{Float64}(nb,nt)
  @sync @distributed for t in 1:nt
    mod(t,10)==0 ? println(t) : nothing
    F=files[t]
    ext=split(F,".")[end]
    mask.tmp3d.=rd(F,var,t,mask.tmp3d)
    BA[:,t]=[b(mask.tmp3d) for b in mask.h_sum]
    GC.gc()
  end
  BA
end

##

volumes(M::gridmask,G::NamedTuple)=begin
  allones=1.0 .+0*G.hFacC
  vol=zeros(length(M.names),length(M.depths))
  for j in 1:length(M.depths)
    M.tmp2d.=M.v_int[j](allones)
    vol[:,j]=[b(M.tmp2d) for b in M.h_sum]
  end
  vol
end

end

