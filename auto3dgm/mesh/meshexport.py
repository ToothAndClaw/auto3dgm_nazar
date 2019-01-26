from os.path import isfile, splitext

from mesh import Mesh
from auto3dgm.mesh.meshfactory import MeshFactory

from vtk import vtkSTLWriter, vtkPolyData, vtkPLYWriter, vtkOBJExporter
from warnings import warn

import sys
import vtk
import os


class MeshExport(object):
    #params: self, and a VTK Object called vtk_mesh
    @staticmethod
    def mesh_from_file(file_path,outputname):
           vtk_mesh=MeshFactory.mesh_from_file(file_path)
           poly_data=vtk_mesh.polydata
           vtkfilename="output.vtk"
           writer = vtk.vtkPolyDataWriter();
           writer.SetFileName(vtkfilename);
           writer.SetInputData(poly_data)
           writer.Write()
 ###########################################################################          
           
           reader = vtk.vtkPolyDataReader()
           reader.SetFileName(vtkfilename)

           stlWriter = vtkSTLWriter(); 
           name_stl=os.path.join(outputname,".stl")
           stlWriter.SetFileName(name_stl); 
           stlWriter.SetInputConnection(reader.GetOutputPort()); 
           stlWriter.Write(); 

           plyWriter = vtkPLYWriter(); 
           name_ply=os.path.join(outputname,".ply")
           plyWriter.SetFileName(name_ply); 
           plyWriter.SetInputConnection(reader.GetOutputPort()); 
           plyWriter.Write(); 

 ##############################################################################         
           name_off=os.path.join(outputname,".off")
           f=open(name_off, "a+")
           f.write("OFF\n")


# Read the VTK file (scalars)
           reader = vtk.vtkPolyDataReader()
           reader.SetFileName(vtkfilename)
           reader.Update()
           polydata = reader.GetOutput()
           CellArray = polydata.GetPolys()
           Polygons = CellArray.GetData()

# of vertices
           for j in range(polydata.GetNumberOfPoints()): 
               j += 1 # i=i+1
   

# of vertices in the polygon... assuming they are same for all   
           bb=Polygons.GetValue(0)   

# number of faces 
           b=[]
           for i_v in xrange(0,  Polygons.GetNumberOfTuples()):
                long_list=Polygons.GetValue(i_v)
                b.append(long_list)
                i_v += 1 
           ver_num=i_v/(1+bb) 
   
# writes in the off file the number of vertices , faces 
           f.write(str(j))
           f.write(" ")
           f.write(str(ver_num))
           f.write("\n")



# writes the coordinates of the vertices
           for i in range(polydata.GetNumberOfPoints()): 
              p = polydata.GetPoints().GetPoint(i) #p is a tuple with the x,y & z
              char_rem = ['(', ',', ')']
              q=str(p).translate(None, ''.join(char_rem))
              f.write(q)
              f.write("\n")
   
   
#writes the vertices indices for the faces
           for i in range(i_v):
    #print b[i]
             f.write(str(b[i]))
             f.write(" ")
             i+=1
             if (i) % (1+bb) == 0:
               f.write("\n")
               continue
           f.close()






