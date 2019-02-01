import glob
from auto3dgm.mesh.meshfactory import MeshFactory 
from auto3dgm.dataset.datasetcollection import DatasetCollection
class DatasetFactory(object):
    #params: directory string
    #params: filetype
    @staticmethod
    def ds_from_dir(directorystring,ftype='ply'):
        path=directorystring
        searchstring=path+'*'+'.'+ftype
        files=glob.glob(searchstring)
        if not files:
            msg = 'No .'+ftype+' files were found in '+path
            raise OSError(msg)
        meshes=[]
        mesheslist=[]
        for file in files: 
            meshes.append(MeshFactory.mesh_from_file(file))
        mesheslist.append(meshes)
        return(DatasetCollection(datasets=mesheslist))
    #params: a list of files
    @staticmethod
    def ds_from_filelist(filelist):
        meshes=[]
        mesheslist=[]
        for file in files: 
            meshes.append(MeshFactory.mesh_from_file(file))
        mesheslist.append(meshes)
        return(DatasetCollection(datasets=mesheslist))




