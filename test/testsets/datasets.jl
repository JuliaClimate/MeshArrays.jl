@testset "Dataset error handling:" begin
    # Unknown names must throw, not silently return the name string
    @test_throws Exception MeshArrays.Dataset("nonexistent_dataset")
    @test_throws Exception MeshArrays.Dataset("nonexistent_dataset", do_read=false)
end
