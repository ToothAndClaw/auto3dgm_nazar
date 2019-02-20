from auto3dgm.mesh.meshfactory import MeshFactory
from numpy import array, zeros

def print_polyfaces(polydata):
    cell_n = polydata.GetNumberOfCells()
    point_n = polydata.GetCell(1).GetNumberOfPoints()
    faces_new = zeros((cell_n, point_n), dtype=int)
    for i in range(cell_n):
        print('cell/polygon ' + str(i))
        for j in range(point_n):
            print(polydata.GetCell(i).GetPointId(j))

vertices = array([[0, 0, 0],[1, 0, 0],[1, 1, 0],[0, 1, 0],[2, 2, 1]], dtype=int)
faces = array([
	[0, 1, 2],
	[1, 2, 3],
	[1, 3, 4]], dtype=int)

mesh = MeshFactory.mesh_from_data(vertices, faces, False)

print(mesh)
print(mesh.vertices)
print(mesh.faces)

"""
Print output should look like this:

<auto3dgm.mesh.mesh.Mesh object at 0x1048962b0>
[[0 0 0]
 [1 0 0]
 [1 1 0]
 [0 1 0]
 [2 2 1]]
[[0 1 2]
 [1 3 4]
 [1 2 4]]
"""

#Another test for the off files:
#Setup: There is a valid off file tooth.off
a=MeshFactory.mesh_from_file("tooth.off")
"""
<class 'vtkCommonDataModelPython.vtkPolyData'>
"""

a.vertices
"""
array([[ 23.1074,  12.3061,  44.2893],
       [ 23.1281,  12.3142,  44.2809],
       [ 23.1233,  12.296 ,  44.2963],
       ..., 
       [ 22.2795,  14.2197,  47.1709],
       [ 22.2679,  14.236 ,  47.1686],
       [ 22.232 ,  14.2798,  47.163 ]])
"""

# This fails currently, but is not a show stopper
a.faces
"""
array([], shape=(0, 3), dtype=float64)
"""
