
module MeshArraysJLD2Ext

    using MeshArrays, JLD2
    import MeshArrays: read_JLD2, write_JLD2

    write_JLD2(a...) = jldsave(a...)
    read_JLD2(a) = load(a)

end
