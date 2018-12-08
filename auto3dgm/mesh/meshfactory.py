from os.path import isfile, splitext

from mesh import Mesh
from vtk import vtkPLYReader, vtkPolyData
from warnings import warn

class MeshFactory(object):
    """Stub of Mesh Factory class.

    Longer class information...

    Attributes:
        attr1: Information.
        attr2: Information.
    """

    @staticmethod
    def mesh_from_file(file_path):
        """Returns a VTK PolyData object from a mesh.

        Currently only works with PLY, TODO add other file types."""
        allowed_filetypes = ['.ply']

        if isfile(file_path) and splitext(file_path)[1] in allowed_filetypes:
            reader = vtkPLYReader()
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