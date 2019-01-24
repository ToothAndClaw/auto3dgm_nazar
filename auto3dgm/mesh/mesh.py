import numpy as np
from vtk.util.numpy_support import vtk_to_numpy # calling vtk_to_numpy doesn't work
import vtk
#import vtkCenterOfMass

class Mesh:
    #params: self, and a VTK Object called vtk_mesh
    def __init__(self, vtk_mesh):
        center = vtk.vtkCenterOfMass()
        print(type(vtk_mesh))
        center.SetInputData(vtk_mesh)
        center.SetUseScalarsAsWeights(False)
        center.Update()
        self.centerpoint = center.GetCenter()

        transform = vtk.vtkTransform()
        transform.Translate(-self.centerpoint[0], -self.centerpoint[1], -self.centerpoint[2])
        transformt = vtk.vtkTransformPolyDataFilter()
        transformt.SetInputData(vtk_mesh)
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

    @property
    def centroid(self):
        #average of all vertices
        '''
        arrset = vtk_to_numpy(self.polydata.GetPoints().GetData())
        ret = [0, 0, 0]
        count = 0
        for point in arrset:
            ret += point
            count += 1
        return [val/count for val in ret]
        '''
        return self.centerpoint

    @property
    def centroid2(self):
        centerFilter = vtk.vtkCenterOfMass()
        centerFilter.SetInputData(self.polydata)
        centerFilter.SetUseScalarAsWeights(False)
        centerFilter.update()
        return centerFilter.GetCenter()

    @property
    def scale(self):
        return np.linalg.norm(vtk_to_numpy(self.polydata.GetPoints().GetData()), 'fro')
