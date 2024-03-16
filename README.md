# Vector-Field-Synthesis

Vector field synthesis algorithms

Research with Professor Guoning Chen and Graduate Student Nguyen Phan under the Cahsi LREU

Algorithm essentially formulates a linear system from the input and then solves it.
Current algorithm is able to reconstruct a sparse vector field.
The next iteration will also take into account the Jacobian of the field.
Subsequently we will tackle the problem of uncstructured meshes and non-linear inputs such as magnitudes.

Tests use matplotlib to display vector fields.
Vector field data is stored in ply files.
Numpy takes care of the heavy linear algebra calculations
