
begin
    using Drifters

    s=size(G.DXC[1])
    u=MeshArray(fill(1.0,s))*G.DXC
    v=MeshArray(fill(0.5,s))*G.DYC

    u=u/G.DXC; u[findall(isnan.(u))].=0;
    v=v/G.DYC; v[findall(isnan.(v))].=0;
#    (u,v)=exchange(u,v,1)
    if !isdefined(Main,:func)
        func=(u -> MeshArrays.update_location_dpdo!(u,Γ))
    end

    F=FlowFields(u,u,v,v,(0.,10.),func)

	x=s[1]*(0.4 .+ 0.2*rand(100))
	y=s[2]*(0.4 .+ 0.2*rand(100))
    I=Individuals(F,x,y,ones(size(x)))

    solve!(I)

    p=1:100;
    scatter(I.🔴.x[p],I.🔴.y[p],color=:blue)
    p=101:200;
    scatter!(I.🔴.x[p],I.🔴.y[p],color=:red)

    current_figure()
end

