module MeshArraysDownloadsExt

    using Downloads

    import MeshArrays: download_file
    
    download_file(url,file) = Downloads.download(url,file)
    
end


