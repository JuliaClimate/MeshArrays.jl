
# Tutorials

## Basics

[The basics tutorial](basics.html) ([code link](https://raw.githubusercontent.com/JuliaClimate/MeshArrays.jl/master/examples/basics.jl)) illustrates how the `MeshArrays.jl` data structures let us write generic code readily applicable to whole families of grids. 

It focuses on a global workflow (smoothing) that requires communication across the entire gridded domain -- a key feature provided by `MeshArrays.jl`. 
The same workflow is repeated three times, for different grid configurations commonly used in numerical models.


Grid scale noise           |  Smoothed noise
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/118325229-2d883d80-b4d1-11eb-953b-ddbb11bcfe1b.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/118325093-f31ea080-b4d0-11eb-8c6e-8cd0cc2cc255.png)

## Geography

[The geography tutorial](geography.html) ([code link](https://raw.githubusercontent.com/JuliaClimate/MeshArrays.jl/master/examples/geography.jl)) deals with interpolation, projection, and vizualization of gridded fields in geographic coordinates. 


Ocean Depth Map          |  Great Circle Path
:------------------------------:|:---------------------------------:
![raw](https://user-images.githubusercontent.com/20276764/144878637-1412679c-f1e6-4491-a8f1-43d729aa224d.png)  |  ![smooth](https://user-images.githubusercontent.com/20276764/144878668-5b681d5e-79b1-45e0-99d0-f80d2afeba8c.png)


## Vector Fields

[The vectors tutorial](vectors.html) ([code link](https://raw.githubusercontent.com/JuliaClimate/MeshArrays.jl/master/examples/vectors.jl)) illustrates how `MeshArrays.jl` represents gridded, vector fields. This enables
 analyses of various quantities, like heat, flow with oceanic currents and atmospheric winds within the climate system. 

Horizontal Streamfunction          |  Vertical Streamfunction
:------------------------------:|:---------------------------------:
![raw](https://github.com/JuliaClimate/GlobalOceanNotebooks/raw/master/OceanTransports/Streamfunction.png)  |  ![smooth](https://github.com/JuliaClimate/GlobalOceanNotebooks/raw/master/OceanTransports/MOC.png)
