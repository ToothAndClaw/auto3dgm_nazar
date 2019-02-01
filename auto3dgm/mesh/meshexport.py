from auto3dgm.mesh.meshfactory import MeshFactory
import vtk


"Static methods:writeToFile: Receives string file location, write options (i.e., file format), "
"and Mesh object; writes a mesh file"

" Important Comment" 

"1) Please change the Output file name if need be "

class MeshExport():

      @staticmethod
      def file_output(file_path):
           polydata_mesh=meshfactory.MeshFactory.mesh_from_file(file_path)
           #######    PLY    ########################           
           plywriter = vtk.vtkPLYWriter()
           plywriter.SetFileName(file_path+"_mesh_exported.ply")
           plywriter.SetInputConnection(polydata_mesh.GetOutputPort())
           plywriter.SetFileTypeToBinary()
           plywriter.Write()
            #######    STL    ########################           
           stlwriter = vtk.vtkSTLWriter()
           stlwriter.SetFileName(file_path+"_mesh_exported.stl")
           stlwriter.SetInputConnection(polydata_mesh.GetOutputPort())
           stlwriter.SetFileTypeToBinary()
           stlwriter.Write()       
            #######    OBJ    ########################           
           
           objwriter = vtk.vtkOBJWriter()           
           objwriter.SetFileName(file_path+"_mesh_exported.obj")
           objwriter.SetInputConnection(polydata_mesh.GetOutputPort())
           objwriter.SetFileTypeToBinary()
           objwriter.Write()  
           







