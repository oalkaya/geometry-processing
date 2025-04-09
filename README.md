# Geometry Processing Assignments

This repository contains a collection of assignments exploring key concepts in geometry processing, from mesh manipulation to surface reconstruction and character deformation.

## Exercise Overview

| Task | Description |
|------|-------------|
| **Assignment 1: Mesh "Hello World"** | Work with the provided code framework to perform basic mesh operations: compute mesh connectivity, surface normals, connected components, and implement a basic subdivision scheme (√3-subdivision). Familiarizes users with `libigl` and the mesh viewer. |
| **Assignment 2: Implicit Surface Reconstruction** | Construct an implicit function from point clouds using MLS (Moving Least Squares), sample it over a volumetric grid, and extract a mesh using the marching cubes algorithm. Tasks include constraint creation, spatial indexing, MLS interpolation, and normal estimation. |
| **Assignment 3: Normals and Curvature (Optional)** | Explore and implement several discrete differential geometry concepts. Compute different types of normals, estimate curvatures, and implement both explicit/implicit Laplacian smoothing and bilateral denoising. |
| **Assignment 4: Mesh Parameterization** | Implement four mesh parameterization techniques—Spring energy, Harmonic (Dirichlet) energy, LSCM, and ARAP. Solve sparse linear systems under both fixed and free boundary conditions and visualize the resulting distortion. |
| **Assignment 5: Detail-Preserving Mesh Editing** | Perform interactive mesh deformation using a two-level multiresolution approach. Includes Laplacian smoothing, deformation transfer, and transferring high-frequency surface details while maintaining interactivity. |
| **Assignment 6: Skeletal Animation and Harmonic Skinning** | Implement a full skeletal-based deformation pipeline: handle selection, harmonic skinning weight computation, Linear Blend Skinning (LBS), Dual Quaternion Skinning (DQS), and Per-face Skinning with Poisson Stitching. Extend to context-aware deformation using example poses. |