import numpy as np
from vtk.util.numpy_support import vtk_to_numpy # calling vtk_to_numpy doesn't work

class Mesh:
    #params: self, and a VTK Object called vtk_mesh
    def __init__(self, vtk_mesh):
        center = vtk_mesh.vtkCenterOfMass()
        center.SetInputData(vtk_mesh.GetPolyData())
        center.SetUseScalarAsWeights(False)
        center.Update()
        self.centerpoint = center.GetCenter()

        transform = vtk_mesh.vtkTransform()
        transform.translate(-center[0], -center[1], -center[2])
        transformt = vtk_mesh.vtkTransformPolyDataFilter()
        transformt.SetInputData(vtk_mesh.GetPolyData())
        transformt.SetTransform(transform)
        transformt.Update()

        self.polydata = transformt.GetOutput()

    ''' we don't have to use getters like this, but using it in this case bc 1) 
    attributes are now pointers and don't have to be updated, 2) can't be 
    overwritten arbitrarily but we could create a setter so that any 
    update goes into the polydata object '''
    @property 
    def vertices(self):
        return vtk_to_numpy(self.polydata.GetPoints().GetData())

    ''' I wrote ugly code here, but vertex IDs are better than raw per-face point 
    coordinates and VTK does not make getting this easy '''
    @property
    def faces(self):
        cell_n = self.polydata.GetNumberOfCells()
        point_n = self.polydata.GetCell(0).GetNumberOfPoints()
        faces = np.zeros((cell_n, point_n), dtype=int)
        for i in range(cell_n):
            for j in range(point_n):
                faces[i, j] = self.polydata.GetCell(i).GetPointId(j)
        return faces