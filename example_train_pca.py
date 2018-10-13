import os
import numpy as np

from sfeal import core as sf

reload(sf)

sfmesh = sf.MESH()
sfmodel = sf.SSM()

""" Some configurations. Modify as required. """
config = dict()
config["root dir"] = "/hpc/mosa004/Lung/Data"  # specify where the root directory for lung meshes are
config["study"] = "Human_Lung_Atlas"  # specify the study
config["path"] = os.path.join(config["root dir"], config["study"])
config["volume"] = "TLC"  # specify the imaged volume
config["fitted_mesh_dir"] = "Lung/SurfaceFEMesh"
config["subjects"] = ["P2BRP042-H18", "P2BRP031-H684", "P2BRP030-H682"]  # subjects going into the PCA
config["lung"] = "LR"  # can be L or R (left or right lung only) or LR (both lungs)
config["morphic original mesh path"] = "morphic_original"
config["morphic aligned mesh path"] = "morphic_aligned"
config["morphic mesh name"] = "Lung_fitted.mesh"
config["reference lung"] = "P2BRP030-H682"
config["number of modes"] = "full"


def _get_mesh():
    """
    Gets the fitted meshes and convert them to morphic meshes.

    :return:
    meshes list containing morphic mesh paths.
    morphic reference_mesh used in _align() for Procrustes registration.

    """
    path = config["path"]
    volume = config["volume"]
    fitted_mesh_dir = config["fitted_mesh_dir"]
    subjects = config["subjects"]
    meshes = []

    for sub in subjects:
        full_path = os.path.join(path, sub, volume, fitted_mesh_dir)
        output_path = sfmesh.convert_cm_mesh(full_path, lung=config["lung"])

        if not os.path.exists(full_path + '/' + config["morphic original mesh path"] + '/' + config["morphic mesh name"]):
            sfmesh.generate_mesh(output_path, lung=config["lung"], save=True)

        for mesh_file in os.listdir(output_path):
            if mesh_file.endswith(config["morphic mesh name"]):
                if sub == config["reference lung"]:
                    reference_mesh = os.path.join(output_path, mesh_file)
                meshes.append(os.path.join(output_path, mesh_file))

    return meshes, reference_mesh


def _align(m, r):
    """
    Method calling the align_mesh() method from SFEAL.core to Procrustes align all the meshes to the reference mesh
    specified above.

    :param m: meshes list from _get_mesh()
    :param r: reference_mesh from _get_mesh()
    :return: aligned_mesh_objects from SFEAL.core.align_mesh()
    """
    meshes = m
    reference_mesh = r

    aligned_mesh_objects = []

    for mesh in meshes:
        d, Z, tform, m = sfmesh.align_mesh(reference_mesh, mesh, scaling=True, reflection='best')
        aligned_mesh_objects.append(m)

    return aligned_mesh_objects


def _prepare_sfeal():
    """
    Preapring the meshes for PCA.

    :return: Prepared SFEAL model object.
    """
    path = config["path"]
    volume = config["volume"]
    fitted_mesh_dir = config["fitted_mesh_dir"]
    subjects = config["subjects"]
    aligned_meshes = []

    for sub in subjects:
        full_path = os.path.join(path, sub, volume, fitted_mesh_dir, config["morphic aligned mesh path"])
        for mesh_file in os.listdir(full_path):
            if mesh_file.endswith(config["morphic mesh name"]):
                aligned_meshes.append(os.path.join(full_path, mesh_file))

    for mesh in range(len(aligned_meshes)):
        sfmodel.add_mesh(aligned_meshes[mesh])

    return sfmodel, aligned_meshes


def _get_score(sfeal_model, aligned_mesh_names):
    """
    Calculates mesh PCA scores and stores them in a csv file 'scores.csv' in the SFEAL module directory

    :param sfeal_model:
    :param aligned_mesh_names:
    :return:
    """
    sf = sfeal_model
    subject_names = sorted(aligned_mesh_names)
    m_distance = list()
    score_array = np.chararray((len(subject_names), sf.num_modes + 1), itemsize=25)
    for i in range(len(subject_names)):
        score, ratio = sf.get_score(subject_names[i])
        mah_distance = sf.get_mahalanobis()
        m_distance.append(mah_distance)
        score_array[i][0] = subject_names[i].split('/')[6]
        for j in range(sf.num_modes):
            score_array[i][j+1] = score['MODE   '+str(j+1)+' SCORE']

    np.savetxt('scores.csv', score_array, delimiter=',', fmt='%s')
    return score_array, m_distance


def main():
    meshes, ref_mesh = _get_mesh()
    aligned_mesh_objs = _align(meshes, ref_mesh)
    number_of_subjects_in_pca = len(config["subjects"])

    """ Create SFEAL and perform PCA """
    sf, aligned_mesh_names = _prepare_sfeal()
    number_of_modes = number_of_subjects_in_pca - 1 if config["number of modes"] == "full" else config["number of modes"]
    pmesh, _ = sf.pca_train(num_modes=number_of_modes)
    sf.save_mesh_id()

    scores, mahalanobis_distance = _get_score(sf, aligned_mesh_names)

    return sf, pmesh


if __name__ == "__main__":
    sf, pmesh = main()
