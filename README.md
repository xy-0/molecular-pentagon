# Python Script for Calculating Aspect Ratios (*R*)

This code is designed to process files containing coordinates of the vertices of Wigner-Seitz (WS) cells and calculate the *R* values for all cells in a crystal lattice.

**Parameters:**
- `grid_number`: The number of points used to compute the moment of inertia tensor.
- `input_dir`: The directory containing all files ([Filename].txt) representing the vertices of the convex WS cell.
- `output_dir`: The directory where output files will be saved.

**Outputs:**
- Coordinates of all internal points ([Filename].txt).
- A .xyz file for molecular dynamics simulations ([Filename].xyz).
- A file containing all the calculated results (output.txt).

**Example Files:**
An example of input files (vertices of WS cells corresponding to the two Wyckoff positions in the A15 phase) and output results are provided in the 'examples' folder.

Code contributed by X.Y. and X.K.
