import glob
from auto3dgm.mesh.meshfactory import MeshFactory 
from auto3dgm.dataset.dataset import Dataset
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
        for file in files: 
            meshes.append(MeshFactory.mesh_from_file(file))
        return(Dataset(meshes=meshes))
    #params: a list of files
    @staticmethod
    def ds_from_filelist(filelist):
        meshes=[]
        for file in files: 
            meshes.append(MeshFactory.mesh_from_file(file))
        return(Dataset(meshes=meshes))
