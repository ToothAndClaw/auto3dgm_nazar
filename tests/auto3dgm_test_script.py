import auto3dgm

dataset_coll = auto3dgm.dataset.datasetfactory.DatasetFactory.ds_from_dir("fixtures/sample")
orig_meshes = dataset_coll.datasets[0]

ss = auto3dgm.mesh.subsample.Subsample(pointNumber=[10, 20], meshes=orig_meshes)
ss_res = ss.ret

low_res_meshes = []
for name, mesh in ss_res[10]['output']['output'].items():
	mesh.name = name
	low_res_meshes.append(mesh)

high_res_meshes = []
for name, mesh in ss_res[20]['output']['output'].items():
	mesh.name = name
	high_res_meshes.append(mesh)

m1 = low_res_meshes[0]

imp.reload(auto3dgm.analysis.correspondence)
corr = auto3dgm.analysis.correspondence.Correspondence(meshes=low_res_meshes)

imp.reload(auto3dgm.analysis.correspondence)
corr_high = auto3dgm.analysis.correspondence.Correspondence(meshes=high_res_meshes, initial_alignment=ga)

corr2 = auto3dgm.analysis.correspondence.Correspondence(meshes=low_res_meshes)
