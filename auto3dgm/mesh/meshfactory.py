from os.path import isfile, splitext

from auto3dgm.mesh.mesh import Mesh
from numpy import array, ndarray, newaxis, concatenate, empty
from vtk import vtkPLYReader,vtkOBJReader,vtkSTLReader,vtkPolyData, vtkPoints, vtkCellArray
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk, numpy_to_vtkIdTypeArray
from warnings import warn


class MeshFactory(object):
    """Stub of Mesh Factory class.
    Longer class information...
    Attributes:
        attr1: Information.
        attr2: Information.

    MeshFactory (Class with Static Methods)
        createFromFile: receives a string file location, creates a VTK mesh object, returns Mesh object
        createFromData: receives vertex and/or face arrays, creates a VTK mesh object, returns Mesh object
        Rationale: Using a factory decouples mesh creation logic from VTK,
        enabling possibly complex loading behavior and even using non-VTK Mesh classes in the future
        Notes: Should be capable of creating meshes from files with formats .ply, .obj, .stl, and .off
    """

    @staticmethod
    def mesh_from_file(file_path):
        """Returns a VTK PolyData object from a mesh.
        TODO ing file types: off. (The PyMesh Package might help.)"""
        allowed_filetypes = ['.ply', '.obj','.stl']

        if isfile(file_path) and splitext(file_path)[1] in allowed_filetypes:
            if splitext(file_path)[1] == '.ply':
                reader = vtkPLYReader()

            elif splitext(file_path)[1] == '.obj':
                reader = vtkOBJReader()

            elif splitext(file_path)[1] == '.stl':
                reader = vtkSTLReader()

            reader.SetFileName(file_path)
            reader.Update()
            
            polydata = reader.GetOutput()
            if isinstance(polydata, vtkPolyData):
                return Mesh(polydata)
            else:
                msg = 'VTK reader output type expected {}, but got {}'.format(
                    'vtkCommonDataModelPython.vtkPolyData', type(polydata))
                raise TypeError(msg)
        else:
            msg = 'File {} not present or not allowed filetype: {}'.format(
                file_path, ', '.join(allowed_filetypes))
            raise OSError(msg)

    @staticmethod
    def mesh_from_data(vertices, faces=empty([0,0]), deep=True):
        """Returns a VTK PolyData object from vertex and face ndarrays"""
        polydata = vtkPolyData()

        # vertices
        points = vtkPoints()
        points.SetData(numpy_to_vtk(vertices, deep=deep))
        polydata.SetPoints(points)

        # faces
        if isinstance(faces, ndarray) and faces.ndim == 2 and faces.shape[1] == 3:
            faces = concatenate((array([3, 3, 3])[:, newaxis], faces), axis=1)
            cells = vtkCellArray()
            nf = faces.shape[0]
            vtk_id_array = numpy_to_vtkIdTypeArray(faces.ravel(), deep=deep)
            cells.SetCells(nf, vtk_id_array)
            polydata.SetPolys(cells)

        return Mesh(polydata)
