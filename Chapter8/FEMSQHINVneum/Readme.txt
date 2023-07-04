------------
file: README
------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose:
%
%    piecewise linear finite elements with uniform refinement and 
%    linear multilevel iterative solution of the discretized equations
%
% details:
%
%    this code solves the following linear elliptic equation:
%
%       - \nabla \cdot (a \nabla u) + b u = f in Omega in R^2
%
%    boundary conditions:
%                               u = g   on \Gamma_D 
%              a \nabla u \cdot n = h   on \Gamma_N
% 
%    where \partial \Omega = \Gamma = \Gamma_D \cup \Gamma_N.
%
%    for implementation details of the FEM discretization, refer to:
%
%        "Finite Element Solution of Elliptic Boundary Value Problems"
%            by Axelsson and Barker
%
%    and to the lecture notes of M. Holst.
%
%    for details of the multilevel numerical method employed, refer to
%    the lecture notes of M. Holst.
%
% author:   michael holst 101095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file description:
-----------------

    README      ==> this file
    main.m      ==> main driver; defines mesh, refines, discretizes, solves.
    mesh.m      ==> mesh and domain definition routine
    mesh_dir.m  ==> example "mesh.m" file for purely dirichlet bc's
    mesh_neu.m  ==> example "mesh.m" file for mix of dirichlet and neumann bc's
    refin.m     ==> uniform mesh refinement routine
    edgtp.m     ==> a mesh search algorithm used by the refinement procedure
    isnod.m     ==> a mesh search algorithm used by the refinement procedure
    mastr.m     ==> master element information definitions
    assem.m     ==> assembly routine (currently a copy of "assem_f.m")
    assem_a.m   ==> matrix assembly following very closely Axelsson and Barker
    assem_f.m   ==> much faster matrix assembly routine with loops rearranged
    aa.m        ==> function "a" in above equation
    bb.m        ==> function "b" in above equation
    ff.m        ==> function "f" in above equation
    gg.m        ==> function "g" in above equation
    hh.m        ==> function "h" in above equation
    uu.m        ==> function "u" in above equation (analytical solution)
    draw.m      ==> draws a finite element mesh
    drawf.m     ==> draws a function defined over a finite element mesh
    mypaus.m    ==> pauses and displays a message
  

