from auto3dgm.dataset.datasetfactory import DatasetFactory
from auto3dgm.mesh.subsample import Subsample

mesh_dir = "/Users/jmw110/sandbox/auto3dgm_test/"

dc = DatasetFactory.ds_from_dir(mesh_dir)
print(dc)
print(dc.datasets)
print(dc.analysis_sets) # Why is the initial dataset copied as an analysis_set?

ss_points = [100, 200]
ss = Subsample(dc.datasets[0], ss_points, 'FPS') # this absolutely does not work
