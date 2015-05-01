# CaTBioMechanics
Computational and Theoretical Bio-mechanics
Created on 05/01/2015 by Walter Kou, Northwestern University

# What is this
This is a 3D finite-element based model for bio-mechanics,
including the model for material (elasticity) and neural-controlled activation (active fiber).
# Implementation?
c++, MPI (serial mesh, but parallel solver, currently)
build on the top of library: Libmesh (http://libmesh.github.io/)

# Version 0: implicit 3D tube model
  Geometry: 3D tube
  material: passive + active
  sovler: implicit FE + petsc
