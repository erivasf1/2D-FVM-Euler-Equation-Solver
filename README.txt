//1D FVM Data Allocation

1) Mesh:
  1.1) in user-provided mesh specs the following are given:
    1.1a) x-range (xmin and xmax)
    1.1b) cellnum (number of cells)
  1.2) the following is outputted from GenerateMesh fcn.:
    1.2a) xcoords: list of xcoord locations for every CELL FACE including xmin and xmax (pretty much a linspace of xmin,xmax,face# which = cellnum+1)

2) Data indexing
  2.1) Before adding ghost cells (SetBoundaryConditions fcn.)
    2.1a) cell_list[0,....,cellnum-1] -- only consisting of interior nodes
    2.1b) identifying corresponding cell faces: cell_i = cell_list[i]; right face = xcoords[i+1] & left face = xcoords[i]
  2.2) After adding 4 ghost cells (2 per side of nozzle)
    2.2a) cell_list[0,....,cellnum+3] -- consisting of interior + ghost cells 
    2.2b) identifying corresponding cell faces: cell_i = cell_list[i]; right face = xcoords[i-1] & left face = xcoords[i-2]

CURRENTLY:
Copies labeled "TEST" are copies of its corresponding file with the main nuance being the absence of pointers. Try getting this to run first (hopefully converge) and then use pointers

THINGS LEARNED:
1) Residuals really are just checks on how well the cell(or domain) satisfies the governing equations (i.e. conservation laws). If the residual is high then solver will want to "greatly" change the values in the cell in hopes of it generating a lesser residual
