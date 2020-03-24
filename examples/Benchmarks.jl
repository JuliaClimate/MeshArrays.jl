# -*- coding: utf-8 -*-
# # Benchmarking `MeshArrays` performance
#
# - Started in 2018/09 (`gaelforget`)
# - Reformated in 2019/08 (`gaelforget`)
# - Added to repo in 2020/01 (`gaelforget`)

# **2020/01/25 using julia 1.3**
#
# _note the speedup in `GridLoad` & `LatitudeCircles` & `demo3` benchmarks_
#
# ```
# gradient          263.030 μs (243 allocations: 3.69 MiB)
# smooth             24.818 ms (11314 allocations: 238.19 MiB)
#
# GridLoad           75.318 ms (7913 allocations: 663.22 MiB)
# read_bin          177.079 μs (65 allocations: 2.01 MiB)
# LatitudeCircles   603.469 ms (260018 allocations: 3.17 GiB)
# ThroughFlow loop   52.728 ms (752206 allocations: 23.48 MiB)
# demo3             643.774 ms (1012062 allocations: 3.19 GiB)
# ```
#
# **2019/08/08 using julia 1.1**
#
# ```
# gradient          299.661 μs (318 allocations: 3.69 MiB)
# smooth             34.006 ms (15572 allocations: 238.32 MiB)
#
# GridLoad          337.105 ms (1960 allocations: 530.35 MiB)
# read_bin          191.141 μs (61 allocations: 1.61 MiB)
# LatitudeCircles     2.904 s (1511608 allocations: 4.42 GiB)
# ThroughFlow loop   65.804 ms (752743 allocations: 24.46 MiB)
# demo3               2.581 s (2264189 allocations: 4.44 GiB)
# ```

using BenchmarkTools
using MeshArrays
p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

cd("/Users/gforget/mywork/data/")

# ## Benchmark 1
#
# This uses `demo2 / smooth` as done in `2018/09` to test exchanges and array operations.

γ,Γ=GridOfOnes("CubeSphere",6,100);
(Rini,Rend,DXCsm,DYCsm)=demo2(Γ);
@btime (dFLDdx, dFLDdy)=gradient(Rini,Γ);
@btime Rend=smooth(Rini,DXCsm,DYCsm,Γ);

# **2019/08/08**
# ```
# 299.661 μs (318 allocations: 3.69 MiB)
#  34.006 ms (15572 allocations: 238.32 MiB)
# ```

# ## Benchmark 2
#
# This uses `demo3` to test indexing type operations.

mygrid=GridSpec("LatLonCap","GRID_LLC90/")
GridVariables=GridLoad(mygrid)
TrspX=mygrid.read(mygrid.path*"TrspX.bin",MeshArray(mygrid,Float32))
TrspY=mygrid.read(mygrid.path*"TrspY.bin",MeshArray(mygrid,Float32))
TauX=mygrid.read(mygrid.path*"TauX.bin",MeshArray(mygrid,Float32))
TauY=mygrid.read(mygrid.path*"TauY.bin",MeshArray(mygrid,Float32))
SSH=mygrid.read(mygrid.path*"SSH.bin",MeshArray(mygrid,Float32))
(UV, LC, Tr)=demo3(TrspX,TrspY,GridVariables);

@btime GridVariables=GridLoad(mygrid)
@btime mygrid.read(mygrid.path*"TrspX.bin",MeshArray(mygrid,Float32))
@btime LC=LatitudeCircles(-89.0:89.0,GridVariables)
@btime for i=1:length(LC); ThroughFlow(UV,LC[i],GridVariables); end
@btime (UV, LC, Tr)=demo3(TrspX,TrspY,GridVariables);

# **2019/08/08**
# ```
# 257.936 ms (1960 allocations: 530.35 MiB)
# 177.847 μs (61 allocations: 1.61 MiB)
#   2.297 s  (1511608 allocations: 4.42 GiB)
#  60.967 ms (752743 allocations: 24.46 MiB)
#   2.305 s  (2264189 allocations: 4.44 GiB)
# ```
