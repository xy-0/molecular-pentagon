import numpy as np

def d2r(a): # degree to radian
    return a * np.pi / 180

long_short_ratio = []
file = r'C:\0_Work_MIT\Xianyou Project\Voronoi Statistics\Rls for Xianyou\A15-p1-2.vol' # Reading a output file of Voro++
with open(file, 'r') as file_to_read:
    while True:
        line = file_to_read.readline()
        if not line:
            break
        coordinates_string = line.split(' ')
        coordinates_vector = []
        for i in coordinates_string:
            if i != '':
                coordinates_vertice = i[i.index('(') + 1:i.index(')')].split(',')
                for j in coordinates_vertice:
                    if j.find('e') != -1:
                        coordinates_vertice[coordinates_vertice.index(j)] = float(0)
                    else:
                        coordinates_vertice[coordinates_vertice.index(j)] = float(j)
                coordinates_vector.append(np.mat(coordinates_vertice).T)

        projection_longest=[0,0,0]
        for theta in range(0, 359, 5):  # Low-precision search (5 degree interval) of long axis of the current Wigner-Seitz cell in a spherical coordinate system
            for phi in range(0, 179, 5):
                vector_current = np.mat([np.sin(d2r(theta)) * np.cos(d2r(phi)),
                                      np.sin(d2r(phi)) * np.sin(d2r(theta)),
                                      np.cos(d2r(theta))])
                projection_sum = 0
                projection_array = []
                for i in coordinates_vector:

                    projection_current = float(vector_current * i)
                    projection_array.append(projection_current)
                projection_array.sort()
                projection_sum = abs(projection_array[-1]) + abs(projection_array[0]) # The longest length of projection that all vertices could have along the direction of vector_current
                if projection_sum > projection_longest[2]:
                    projection_longest=[theta, phi, projection_sum]

        for theta in np.arange(projection_longest[0] - 5, projection_longest[0] + 6, 1):  # High-precision search (1 degree interval) of long axis of the current Wigner-Seitz cell in a spherical coordinate system
            for phi in np.arange(projection_longest[1] - 5, projection_longest[1] + 6, 1):
                vector_current = np.mat([np.sin(d2r(theta)) * np.cos(d2r(phi)),
                                      np.sin(d2r(phi)) * np.sin(d2r(theta)),
                                      np.cos(d2r(theta))])
                projection_sum = 0
                projection_array = []
                for i in coordinates_vector:
                    projection_current = float(vector_current * i)
                    projection_array.append(projection_current)
                projection_array.sort()
                projection_sum = abs(projection_array[-1]) + abs(projection_array[0])
                if projection_sum > projection_longest[2]:
                    projection_longest=[theta, phi, projection_sum]

        vector_longest = np.mat([np.sin(d2r(projection_longest[0])) * np.cos(d2r(projection_longest[1])),
                                 np.sin(d2r(projection_longest[1])) * np.sin(d2r(projection_longest[0])),
                                 np.cos(d2r(projection_longest[0]))])

        vector_a = np.cross(vector_longest, [1, 0,
                                              0])  # Looking for an arbitray combination of linearly independent vector_a and vector_b that is in the normal plane of vector-longest
        if np.linalg.norm(vector_a) != 0:
            vector_a /=  np.linalg.norm(vector_a)
            vector_b = np.cross(vector_longest, vector_a)
            vector_b /= np.linalg.norm(vector_b)
        else:
            vector_a = np.mat([0, 1, 0])
            vector_b = np.mat([0, 0, 1])
        projection_sum = 0
        for alpha in range(0, 180): # Calculate the short axis on the plane defined by vector_a and vector_b
            vector_short = [np.sin(d2r(alpha)) * vector_a[0][0] + np.cos(d2r(alpha)) * vector_b[0][0],
                            np.sin(d2r(alpha)) * vector_a[0][1] + np.cos(d2r(alpha)) * vector_b[0][1],
                            np.sin(d2r(alpha)) * vector_a[0][2] + np.cos(d2r(alpha)) * vector_b[0][2]]
            projection_array = []
            for i in coordinates_vector:
                projection_short = float(vector_short * i)
                projection_array.append(projection_short)
            projection_array.sort()
            projection_sum += abs(projection_array[-1]) + abs(projection_array[0]) #
        short_axis=projection_sum/180 # The length of short axis is the average value of vertices' longest projection length on the plane defined by vector_a and vector_b
        long_short_ratio.append(projection_longest[2] / short_axis)
print(long_short_ratio)