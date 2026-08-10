@testset "MeshArrays tests:" begin
    for nTopo=1:4
        if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
        elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
        elseif nTopo==3; grTopo="PeriodicChannel"; nFaces=1; N=400;
        elseif nTopo==4; grTopo="PeriodicDomain"; nFaces=1; N=400;
        end;
        Npt=nFaces*N*N
        γ=MeshArrays.GridSpec_ones(grTopo,nP=N)
        Γ=MeshArrays.GridLoad_ones(γ;option="full")
        @test γ.class == grTopo
        Rini= 0.; Rend= 0.;
        (Rini,Rend,DXCsm,DYCsm)=demo2(Γ);
        @test isa(Rend,MeshArray)
        @test sum(isfinite.(Rend)) == Npt
        Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
        Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
        #println([Sini Send])
        @test isapprox(Sini,1.000; atol=1e-2)
        @test isapprox(Send,0.093; atol=1e-2)
        (dRdx,dRdy)=gradient(Rend,Γ)
        MeshArrays.exchange_main(Rend,2)
        MeshArrays.exchange_main(dRdx,dRdy,1)
        (dRdx_e,dRdy_e)=MeshArrays.exchange_main(dRdx,dRdy,1)
        @test isa(dRdx_e,MeshArray_wh)
    end
end

@testset "exchange! in-place tests:" begin
    for nTopo=1:4
        if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
        elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
        elseif nTopo==3; grTopo="PeriodicChannel"; nFaces=1; N=400;
        elseif nTopo==4; grTopo="PeriodicDomain"; nFaces=1; N=400;
        end;
        γ=MeshArrays.GridSpec_ones(grTopo,nP=N)
        Γ=MeshArrays.GridLoad_ones(γ;option="full")
        (_,Rend,_,_)=demo2(Γ)
        (dRdx,dRdy)=gradient(Rend,Γ)

        # scalar: exchange! result must match exchange
        ref_s  = exchange(Rend)
        buf_s  = exchange_alloc(Rend)
        exchange!(buf_s, Rend)
        @test isa(buf_s, MeshArray_wh)
        for i in eachindex(buf_s.MA.f)
            @test buf_s.MA.f[i] ≈ ref_s.MA.f[i]
        end

        # vector: exchange! result must match exchange
        (ref_u, ref_v) = exchange(dRdx, dRdy)
        (buf_u, buf_v) = exchange_alloc(dRdx, dRdy)
        exchange!(buf_u, buf_v, dRdx, dRdy)
        @test isa(buf_u, MeshArray_wh)
        @test isa(buf_v, MeshArray_wh)
        for i in eachindex(buf_u.MA.f)
            @test buf_u.MA.f[i] ≈ ref_u.MA.f[i]
            @test buf_v.MA.f[i] ≈ ref_v.MA.f[i]
        end

        # idempotency: calling exchange! twice gives the same result
        exchange!(buf_s, Rend)
        for i in eachindex(buf_s.MA.f)
            @test buf_s.MA.f[i] ≈ ref_s.MA.f[i]
        end
    end
end

