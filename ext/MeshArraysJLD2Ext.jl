
module MeshArraysJLD2Ext

    using MeshArrays, JLD2
    import MeshArrays: read_JLD2, write_JLD2

    """
       write_JLD2(file;kwargs...)

    Call `JLD2.jldsave(file;kwargs...)`
    """
    write_JLD2(fil;kwargs...) = jldsave(fil;kwargs...)

    """
       read_JLD2(file)

    Call `FileIO.load(file)`. Works `.jld2`, `.jpg`, etc files.
    """
    read_JLD2(a) = load(a)

end
