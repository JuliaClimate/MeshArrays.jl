
module ERA5_to_ECCO4

using MeshArrays, Interpolations, NCDatasets, DataFrames, Statistics, FortranFiles, Dates

## list of variables etc

list_in=["dlw","dsw","pres","rain","d2m","tmp2m_degC","u10m","ustr","v10m","vstr","wspeed"];
list_ds=["msdwlwrf","msdwswrf","sp","tp","d2m","t2m","u10","metss","v10","mntss","..."]

offset=zeros(12)
offset[6]=-273.15
factor=ones(12)
factor[1]=-1.0
factor[2]=-1.0
factor[4]=1/3600
factor[7]=1.0
factor[8]=-1.0
factor[9]=1.0
factor[10]=-1.0

## compute specific humidity

Rdry=287.0597 ; Rvap=461.5250 ; a1=611.21 ; a3=17.502 ; a4=32.19 ; T0=273.16
#Calculation of E saturation water vapour from Teten's formula
E(dtas)=a1*exp(a3*(dtas-T0)/(dtas-a4))
#Calculation of saturation specific humidity at 2m qsat  (equal to huss)
qsat(ps,E)=(Rdry/Rvap)*E/(ps-((1-Rdry/Rvap)*E))

## interpolation and bin average code

function setup_arrays_ERA5(Γ)
    lon=collect(0.0:0.25:359.75)
    lat=collect(90.0:-0.25:-90.0)
    ds_knn=calc_knn(lon,lat,Γ)
    A,B=qt_arrays(lon,lat,Γ)
    x=zeros(length(lon),length(lat))
    (lon=lon,lat=lat,knn=ds_knn,A=A,B=B,x=x)
end

function calc_knn(q_lon,q_lat,Γ)
    lo=[q_lon[i] for i in 1:length(q_lon), j in 1:length(q_lat)]
    la=[q_lat[j] for i in 1:length(q_lon), j in 1:length(q_lat)]
    (f,i,j,c)=knn(Γ.XC,Γ.YC,lo[:],la[:])
    c
end

q_fx(lon)=mod(lon,360.0)*4+1
q_fy(lat)=(-lat+90.0)*4+1

function qt_arrays(q_lon,q_lat,Γ)
  ni=length(q_lon)
  nj=length(q_lat)
  A=zeros(ni+1,nj)
  B=MeshArray(Γ.XC.grid)
  (A,B)
end

function calc_interp(a,A,B,Γ)
    A[1:end-1,:].=a
    A[end,:].=a[1,:]
    itp = interpolate(A, BSpline(Linear()));
    for i in eachindex(B)
        B[i][:].=[itp(q_fx(Γ.XC[i][ij]),q_fy(Γ.YC[i][ij])) for ij in eachindex(B[i])]
    end
end

storage_arrays(nx,ny,nt)=( y=NaN*zeros(nx,ny),n=NaN*zeros(nx,ny),z=NaN*zeros(nx,ny,nt),zi=NaN*zeros(nx,ny,nt) )

storage_arrays(γ::gcmgrid,L::NamedTuple) = begin
  (nx,ny)=γ.ioSize
  (nnx,nny)=(length(L.lon),length(L.lat))

  buffer2d=Array{Union{Missing, Float64}}(undef, nnx, nny)
  buff2d  =Array{Union{Missing, Float64}}(undef, nnx, nny)
  A=Array{Float32}(undef, nnx, nny)

  ( y=NaN*zeros(nx,ny),yi=NaN*zeros(nx,ny),n=NaN*zeros(nx,ny) , 
    buffer2d=buffer2d,buff2d=buff2d,A=A )
end

## average to 3h

function ave3h(a,buffer2d,buff2d,y,d,r,rec0; path="ERA5/", variable=1) 
    a.=0.0
    for h in 1:3
      rec=(d-1)*24+(r-1)*3+h
      mo=findall(rec0.<rec)[end]
      mo<10 ? mon="0$(mo)" : mon="$(mo)"

      fil=joinpath(path,"$(y)","ERA5_$(y)_$(mon).nc")
      ds=Dataset(fil)
      jj=rec-rec0[mo]
      ndims(ds["d2m"])==4 ? ii=(:,:,1,jj) : ii=(:,:,jj)
      iii=(:,:,2,jj)
      if variable==11
        buffer2d.=ds["u10"][ii...]
        ismissing(buffer2d[1]) ? buffer2d.=ds["u10"][iii...] : nothing
        a.+= (buffer2d.^2)
        buffer2d.=ds["v10"][ii...]
        ismissing(buffer2d[1]) ? buffer2d.=ds["v10"][iii...] : nothing
        a.+= (buffer2d.^2)
      else
        w=list_ds[variable]
        buffer2d.=ds[w][ii...]
        ismissing(buffer2d[1]) ? buffer2d.=ds[w][iii...] : nothing
        a.+= factor[variable] .*(buffer2d .+offset[variable])
      end
      close(ds)
    end
    a.=a/3
end

function loop_over_years(variable=1,path="ERA5/")
  y0=1941
  ny=83
  ndmax=366

  γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
  Γ=GridLoad(γ)

  path_out=joinpath(path,"ERA5_llc90")

  L=setup_arrays_ERA5(Γ)
  df=DataFrame(:index=>L.knn[:],:x=>L.x[:])
  S=storage_arrays(γ,L)

  for y in y0.+(0:ny-1)
    nd=min(Day(DateTime(y+1,1)-DateTime(y,1)).value,ndmax)
    fil_out=joinpath(path_out,"llc90_ERA5_$(list_in[variable])_$y")
    println(fil_out)
    if !isfile(fil_out)
      g=FortranFiles.FortranFile(fil_out,"w",access="direct",recl=90*1170*4,convert="big-endian")
      dt=Dates.Day.([[DateTime(y,m+1,1)-DateTime(y,m,1) for m in 1:11];DateTime(y+1,1,1)-DateTime(y,12,1)])
      rec0=[0;cumsum([a.value for a in dt])*24]

      for d in 1:nd
        mod(d,10)==0 ? println([y d]) : nothing
        for r in 1:8
          ave3h(S.A,S.buffer2d,S.buff2d,y,d,r,rec0; path=path, variable=variable)
          variable==11 ? S.A.=sqrt.(S.A) : nothing        

          calc_interp(S.A,L.A,L.B,Γ)
          S.yi.=write(L.B)

          df.x.=S.A[:]
          gdf=groupby(df,:index)
          df2=combine(gdf, :x => mean, nrow)

          S.y.=NaN
          [S.y[df2.index[i]]=df2.x_mean[i] for i in 1:size(df2,1)]
          [S.y[i]=(isnan(S.y[i]) ? S.yi[i] : S.y[i]) for i in eachindex(S.y)]

          write(g,rec=r+(d-1)*8,Float32.(write(L.B)))
        end
      end

      close(g)
    else
      println("skipping : already done")
    end
  end
  S
end

end #module ERA5interp

