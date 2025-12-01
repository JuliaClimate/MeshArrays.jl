
module MeshArraysJLD2Ext

    using MeshArrays, JLD2
    import MeshArrays: read_jld2, write_jld2

    """
       write_jld2(file;kwargs...)

    Call `JLD2.jldsave(file;kwargs...)`
    """
    write_jld2(fil;kwargs...) = jldsave(fil;kwargs...)

    """
       read_jld2(file)

    Call `FileIO.load(file)`. Works for not just `.jld2`, but also `.jpg`, etc files.
    """
    read_jld2(a) = load(a)

end
