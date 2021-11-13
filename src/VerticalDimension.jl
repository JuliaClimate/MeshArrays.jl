"""
    isosurface(θ,T,Γ)

Depth of isosurface θ=T
"""
function isosurface(θ,T,Γ)
    d=NaN*similar(θ[:,1])
    nr=size(θ,2)
    for j=1:size(d,1)
        for k=1:nr-1
            i=findall(isnan.(d[j]).&(θ[j,k].>T).&(θ[j,k+1].<=T))
            a=(θ[j,k][i] .- T)./(θ[j,k][i] .- θ[j,k+1][i])
            d[j][i]=(1 .- a).*Γ.RC[k] + a.*Γ.RC[k+1]
            i=findall(isnan.(d[j]).&(θ[j,k].<=T).&(θ[j,k+1].>T))
            a=(θ[j,k+1][i] .- T)./(θ[j,k+1][i] .- θ[j,k][i])
            d[j][i]=(1 .- a).*Γ.RC[k+1] + a.*Γ.RC[k]
        end
    end
    return d
end
