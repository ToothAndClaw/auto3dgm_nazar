import vtk
#import numpy as np

output_file_name = 'plate.off'
input_mesh = 'plate.vtk'


# wants to append the file
f=open(output_file_name, "a+")
f.write("OFF\n")


# Read the VTK file (scalars)
reader = vtk.vtkPolyDataReader()
reader.SetFileName(input_mesh)
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
print i_v   
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





# edges
#edges = vtk.vtkExtractEdges()
#edges.SetInput(polydata)
#edge_mapper = vtk.vtkPolyDataMapper()
#edge_mapper.SetInput(edges.GetOutput())

f.close()


#writer = vtk.vtkPolyDataWriter()
#writer.Write()