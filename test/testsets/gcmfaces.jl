@testset "gcmfaces type:" begin
    for ID in (:PeriodicDomain, :CS32, :LLC270)
        γ=GridSpec(ID=ID)
        MeshArrays.gcmfaces(γ)
        MeshArrays.gcmfaces(γ,Float32)
        tmp=MeshArrays.gcmfaces(γ,Float32,3)
        tmp[3,1,2]=1.0
        view(tmp,:,1,1:2)
        MeshArrays.gcmfaces(γ,tmp.f)
        MeshArrays.fsize(tmp)
        MeshArrays.fsize(tmp,2)
        MeshArrays.fsize(tmp.f)
        MeshArrays.fsize(tmp.f,2)
        size(tmp)
        size(tmp,3)
        tmp1=similar(tmp)
        2 .*tmp1
        findall(tmp.>0)
        MeshArrays.nFacesEtc(tmp)
        @suppress show(tmp)

        x=tmp[1:10,1,1:2]; y=x[2]; x[3]=1.0
        view(x,1:3,:,1)
        MeshArrays.gcmsubset(γ,x.f,x.fSize,x.aSize,x.i,x.iSize)
        MeshArrays.fsize(x.f)
        MeshArrays.fsize(x.f,2)
        size(x)
        size(x,3)

        MeshArrays.fijind(tmp,10)
        MeshArrays.fsize(tmp)
        @test isa(tmp,MeshArrays.gcmfaces)

        tmp=MeshArray(γ)
        tmp1=findall(tmp.>0)
        tmp[tmp1]
        tmp[tmp1].=1.0
        size(tmp1)
        tmp1[2]
        view(tmp1,:)
        @suppress show(tmp1)
        similar(tmp1)
        #tmp[tmp1]

        @suppress show(tmp)
        MeshArrays.getindexetc(tmp,2)
        MeshArray(γ,tmp.f,meta=tmp.meta)
        MeshArray(γ,meta=tmp.meta)
        MeshArray(γ,Float32,3,4)
        MeshArray(γ,Float32,tmp.fSize,tmp.fIndex,2,3)

        tmp1=MeshArray(γ,Float32,3)
        MeshArray(γ,tmp1.f,meta=tmp1.meta)
        MeshArrays.getindexetc(tmp1,2,1)
    end
end
