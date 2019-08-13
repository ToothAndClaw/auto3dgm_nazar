import vtk

class MeshExport(object):

    def file_output(fp,mesh,format='ply'):
        tmp=mesh.name.split("/")
        t=len(tmp)
        name=tmp[t-1].split(".")[0]
        if format=='ply':
            writer=vtk.vtkPLYWriter()
        if format=='stl':
            writer=vtk.vtkPLYWriter()
        if format=='obj':
            writer=vtk.vtkOBJWriter()
        writer.SetInputData(mesh.polydata)
        writer.SetFileName(fp+name+"."+format)
        writer.Write()    
        return(True)
