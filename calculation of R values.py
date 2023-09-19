# coding: utf-8
'''
    Computes the aspect ratio of a convex polyhedron.
    Calculates the eigenvalues and eigenvectors of the moment of inertia tensor for a polyhedron.

    Parameters:
    grid_number: the number of points used to compute the moment of inertia tensor.
    input_dir: The local directory containing the files ([Filename].txt) representing the vertices of the convex polyhedron.
    output_dir: The local directory where output files will be saved.

    Returns:
    The coordinates of all internal points ([Filename].txt).
    An .xyz file suitable for molecular dynamics simulations ([Filename].xyz).
    An output.txt file containing all the calculated results.
'''

import ast
import os
import numpy as np
import MDAnalysis as mda
from scipy.spatial import Delaunay

def write_xyz(coordinates, output_file, comment="Coordinates of carbon atoms"):
    """
    Write coordinates to an .xyz file format.

    :param coordinates: List of tuples containing x, y, z coordinates
    :param output_file: Name of the output .xyz file
    :param comment: Comment to be written on the second line of the file
    """
    with open(output_file, 'w') as f:
        # Number of atoms
        f.write(str(len(coordinates)) + "\n")
        # Comment line
        f.write(comment + "\n")
        # Coordinates
        for coord in coordinates:
            f.write(f"C {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

def get_inertia_eig(output_dir, filename):
    my_file = output_dir+'\\{:s}.txt'.format(filename)
    my_xyz = output_dir+'\\{:s}.xyz'.format(filename)
    my_out = output_dir+'\\output.txt'
    if not os.path.exists(my_xyz):
        with open(my_file, 'r') as file:
            data_string = file.read().strip()

        # Convert the string to a Python 2D list (matrix)
        data_matrix = ast.literal_eval(data_string)
        write_xyz(data_matrix, my_xyz)
    u = mda.Universe(my_xyz, in_memory=True)
    mol = u.select_atoms("all")

    mass = mol.total_mass()

    inertia = mol.moment_of_inertia(pbc=False)
    eigval, eigvec = np.linalg.eig(inertia)
    sorted_indices = np.argsort(eigval)
    # Sort eigenvalues
    sorted_eigenvalues = eigval[sorted_indices]/mass
    # Sort eigenvectors according to the sorted_indices
    sorted_eigenvectors = eigvec[:, sorted_indices]

    str1 = '{:s}: '.format(filename)
    strEval = " ".join(f"{item:.8g}" for item in sorted_eigenvalues) + ' '
    strEvec1 = " ".join(f"{item:.8g}" for item in sorted_eigenvectors[:, 0]) + ' '
    strEvec2 = " ".join(f"{item:.8g}" for item in sorted_eigenvectors[:, 1]) + ' '
    strEvec3 = " ".join(f"{item:.8g}" for item in sorted_eigenvectors[:, 2]) + ' '
    file_out = open(my_out, 'a')
    if (sorted_eigenvalues[1]-sorted_eigenvalues[0])>(sorted_eigenvalues[2]-sorted_eigenvalues[1]):
        aspect_ratio=sorted_eigenvalues[2]/sorted_eigenvalues[0]
    else:
        aspect_ratio = sorted_eigenvalues[0] / sorted_eigenvalues[2]
    file_out.write(filename+' '+"{:.8g}".format(aspect_ratio)+' '+strEval + strEvec1 + strEvec2 + strEvec3+'\n')
    file_out.close()

def is_point_inside_convex_polyhedron(point, polyhedron_vertices):
    """
    Check if a 3D point is inside a convex polyhedron.

    Parameters:
    point: A tuple representing the 3D coordinates of the test point (x, y, z).
    polyhedron_vertices: A list of tuples representing the 3D coordinates of the polyhedron's vertices.

    Returns:
    True if the point is inside the polyhedron, False otherwise.
    """
    polyhedron_points = np.array(polyhedron_vertices)
    triangulation = Delaunay(polyhedron_points)
    test_point = np.array(point)
    simplex_index = triangulation.find_simplex(test_point)
    return simplex_index != -1

def get_inside_points(vertex, grid_number):
    inside_points = []
    x_max = max(x for x, _, _ in vertex)
    x_min = min(x for x, _, _ in vertex)
    x_dim = x_max - x_min
    y_max = max(x for _, x, _ in vertex)
    y_min = min(x for _, x, _ in vertex)
    y_dim = y_max - y_min
    z_max = max(x for _, _, x in vertex)
    z_min = min(x for _, _, x in vertex)
    z_dim = z_max - z_min
    grid_size = (x_dim * y_dim * z_dim / (2 * grid_number)) ** (1 / 3)
    for x in np.arange(x_min, x_max + 2 * grid_size, grid_size).tolist():
        for y in np.arange(y_min, y_max + 2 * grid_size, grid_size).tolist():
            for z in np.arange(z_min, z_max + 2 * grid_size, grid_size).tolist():
                if is_point_inside_convex_polyhedron((x, y, z), vertex) == True:
                    inside_points.append((x, y, z))
    return inside_points

def get_aspect_ratio(grid_number, input_vertex_cor_dir, output_dir):
    file_out = open(output_dir + '\\output.txt', 'a')
    file_out.write('''In the 'output.txt' file, each column corresponds to the following:
1: Name of the polyhedron.
2: Aspect ratio (R) of the polyhedron
3-5: Eigenvalues of the moment of inertia tensor in ascending order.
6-8: Eigenvector corresponding to the minimal eigenvalue of the moment of inertia tensor.
9-11: Eigenvector corresponding to the intermediate eigenvalue of the moment of inertia tensor.
12-14: Eigenvector corresponding to the maximal eigenvalue of the moment of inertia tensor.
##########################################\n''')
    file_out.close()
    file_list = os.listdir(input_vertex_cor_dir)
    for file in file_list:
        file_out_path = output_dir + '\\' + file[:-4] + '.txt'
        file_out = open(file_out_path, 'a')
        vertex2 = []
        with open(input_vertex_cor_dir + '\\' + file, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline()
                if not lines:
                    break
                vertex1 = lines.split()
                vertex3 = []
                for i in vertex1:
                    if '\n' not in i:
                        vertex3.append(float(i))
                    else:
                        vertex3.append(float(i.replace('\n', '')))
                vertex2.append(vertex3)
            inside_points = get_inside_points(vertex2, grid_number)
            file_out.write(str(inside_points) + '\n')
        file_out.close()
        get_inertia_eig(output_dir, file[:-4])

if __name__ == "__main__":
    #Roughly define how many internal points are considered to calculate the moment of inertia tensor
    grid_number=100000
    #Define the folder of input files
    input_dir=r'/home/input'
    #Define the folder of ouput files
    output_dir=r'/home/output'
    get_aspect_ratio(grid_number,input_dir,output_dir)





