# ShadyRay
Path tracing engine for evaluating a complex light sources algorithm

Done as part of undergraduate thesis work on Rendering Complex Light Sources at MAHE, Manipal, India.

## Features
* Forwards path tracing algorithm with random light sampling.
* Various material models including microfacet reflection and transmission (Beckmann distribution), Lambert, and Ashikhmin-Shirley for coated surfaces.
* Lua scripting for scene description.
* Particle tracing for complex light sources, with object point sampling via barycentric coordinates.
* Volume clustering and lookup for point sources.
* Parallelization via Intel TBB.
* Fast ray-triangle intersection with masking via Intel Embree 2.x.
* Point and general emissive surface lights.
* Textures with normal mapping.
