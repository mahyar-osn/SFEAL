import os
import numpy as np

from sfeal import core as sf

from subjects import subjects as list_of_subjects

from mlr import EiMLR, EeMLR

# reload(sf)

sfmesh = sf.MESH()
sfmodel = sf.SSM()

""" Some configurations. Modify as required. """
config = dict()
# config["root dir"] = "/hpc/mosa004/Lung/Data"  # specify where the root directory for lung meshes are
# config["study"] = "Human_IPF"  # specify the study
# config["path"] = os.path.join(config["root dir"], config["study"])
# config["volume"] = "TLC"  # specify the imaged volume
# config["scale"] = True  # whether lung size is normalised
# config["fitted_mesh_dir"] = "Lung/SurfaceFEMesh/KatherineMesh"
# # config["subjects"] = list_of_subjects  # subjects going into the PCA
# config["subjects"] = os.listdir(config["path"])  # subjects going into the PCA
# config["lung"] = "LR"  # can be L or R (left or right lung only) or LR (both lungs)
# config["morphic original mesh path"] = "morphic_original"
# config["morphic aligned mesh path"] = "morphic_aligned"
# config["morphic mesh name"] = "Lung_fitted.mesh"
# # config["reference lung"] = "IPF613"
# config["subjects for pca"] = []
# config["number of modes"] = 5
config["root dir"] = "/hpc/mosa004/SFEAL/data"  # specify where the root directory for lung meshes are
config["study"] = "HLA_HA"  # specify the study
config["path"] = os.path.join(config["root dir"], config["study"])
config["volume"] = "Insp"  # specify the imaged volume
config["scale"] = True  # whether lung size is normalised
config["fitted_mesh_dir"] = "original/cmiss"
config["subjects"] = list_of_subjects  # subjects going into the PCA
# config["subjects"] = os.listdir(config["path"])  # subjects going into the PCA
config["lung"] = "LR"  # can be L or R (left or right lung only) or LR (both lungs)
config["morphic original mesh path"] = "morphic_original"
config["morphic aligned mesh path"] = "morphic_aligned"
config["morphic mesh name"] = "Lung_fitted.mesh"
config["reference lung"] = "AGING025"
config["subjects for pca"] = []
config["number of modes"] = 5


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
    not_exist_meshes = []

    for sub in subjects:

        full_path = os.path.join(path, sub, volume, fitted_mesh_dir)
        output_path = sfmesh.convert_cm_mesh(full_path, lung=config["lung"])

        if output_path is None:
            not_exist_meshes.append(sub)
            continue

        config["subjects for pca"].append(sub)
        if not os.path.exists(full_path + '/' + config["morphic original mesh path"] + '/' + config["morphic mesh name"]):
            sfmesh.generate_mesh(output_path, lung=config["lung"], save=True)

        for mesh_file in os.listdir(output_path):
            if mesh_file.endswith(config["morphic mesh name"]):
                # if sub == config["reference lung"]:
                #     reference_mesh = '/hpc/mosa004/SFEAL/data/HLA_HA/AGING025/Insp/original/cmiss/morphic_original/Lung_fitted.mesh'
                meshes.append(os.path.join(output_path, mesh_file))

    reference_mesh = '/hpc/mosa004/SFEAL/data/HLA_HA/AGING025/Insp/original/cmiss/morphic_original/Lung_fitted.mesh'

    return meshes, reference_mesh


def _align(m, r, scaling=False):
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
        d, Z, tform, m = sfmesh.align_mesh(reference_mesh, mesh, scaling=scaling, reflection='best')
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
    subjects = config["subjects for pca"]
    aligned_meshes = []

    for sub in subjects:
        full_path = os.path.join(path, sub, volume, fitted_mesh_dir, config["morphic aligned mesh path"])
        if not os.path.exists(full_path):
            print('morphic_aligned directory does not exist for subject {}. Skipping...'.format(sub))
            continue
        for mesh_file in os.listdir(full_path):
            if mesh_file.endswith(config["morphic mesh name"]):
                aligned_meshes.append(os.path.join(full_path, mesh_file))

    for mesh in range(len(aligned_meshes)):
        sfmodel.add_mesh(aligned_meshes[mesh])

    return sfmodel, aligned_meshes


def _get_score(sfeal_model, mesh, aligned_mesh_names):
    """
    Calculates mesh PCA scores and stores them in a csv file 'scores.csv' in the SFEAL module directory

    :param sfeal_model:
    :param aligned_mesh_names:
    :return:
    """
    pmesh = mesh
    sf = sfeal_model
    subject_names = sorted(aligned_mesh_names)
    m_distance = list()
    score_array = np.chararray((len(subject_names), sf.get_number_of_modes() + 1), itemsize=25)
    for i in range(len(subject_names)):
        score, ratio = sf.get_score(subject_names[i])
        mah_distance = sf.get_mahalanobis()
        m_distance.append(mah_distance)
        score_array[i][0] = subject_names[i].split('/')[6]
        for j in range(sf.num_modes):
            score_array[i][j+1] = score['MODE   '+str(j+1)+' SCORE']

    np.savetxt('aging_scores.csv', score_array, delimiter=',', fmt='%s')
    return m_distance


def _read_file():
    import pandas as pd
    f = '/people/mosa004/Desktop/predict_subjects_for_Hari.csv'
    df = pd.read_csv(f, header=0, delimiter=',')
    return df


def main():
    meshes, ref_mesh = _get_mesh()
    aligned_mesh_objs = _align(meshes, ref_mesh, config["scale"])
    number_of_subjects_in_pca = len(config["subjects for pca"])

    """ Create SFEAL and perform PCA """
    sf, aligned_mesh_names = _prepare_sfeal()
    number_of_modes = number_of_subjects_in_pca - 1 if config["number of modes"] == "full" else config["number of modes"]
    pmesh, _ = sf.pca_train(num_modes=number_of_modes)
    sf.save_mesh_id()

    pr_path = "/hpc/mosa004/Lung/Data/Human_IPF"
    pr_path_1 = "/TLC/Lung/SurfaceFEMesh/KatherineMesh/morphic_aligned/Lung_fitted.mesh"

    projected_weights = dict()
    for project_subject in os.listdir(pr_path):
        pr_sub_path = pr_path + "/" + project_subject + pr_path_1
        if os.path.exists(pr_sub_path):
            print("Subject = {}".format(project_subject))
            w, r = sf.project_new_mesh(pr_sub_path)
            projected_weights[project_subject] = w

    weight_list = list()
    for subjects in projected_weights.keys():
        weight_list.append(subjects)
        for modes in sorted(projected_weights[subjects]):
            weight_list.append(projected_weights[subjects][modes])

    a = np.asarray(weight_list, dtype=object)
    b = a.reshape(30, 6)
    np.savetxt('/people/mosa004/Desktop/ipf_weights.txt', b, fmt='%s')

    # _get_score(sf, pmesh, aligned_mesh_names)

    # modes = [1, 2, 3, 4]
    # weights = [0, 0, 0, 0]
    #
    # for m in modes:
    #     for w in np.linspace(-2.5, 2.5, 2):
    #         if m == 1:
    #             weights[0] = w
    #         elif m == 2:
    #             weights[1] = w
    #         elif m == 3:
    #             weights[2] = w
    #         elif m == 4:
    #             weights[3] = w
    #
    #         print weights
    #
    #         if w == -2.5:
    #             config["export_name"] = 'mode_{}_N25_no_scale'.format(m)
    #         else:
    #             config["export_name"] = 'mode_{}_P25_no_scale'.format(m)
    #
    #         sf.export_to_cm(pmesh, weights, name=config["export_name"], lung='L', show_mesh=False)
    #         sf.export_to_cm(pmesh, weights, name=config["export_name"], lung='R', show_mesh=False)
    #
    #         weights = [0, 0, 0, 0]

    # config["export_name"] = 'mode_1_N25_no_scale'
    # sf.export_to_cm(pmesh, weights, name=config["export_name"], lung='L', show_mesh=False)

    #
    # pft_df = _read_file()
    # pfts = dict()
    # for subject, row in pft_df.iterrows():
    #     pfts['subject'] = row['subjects']
    #     pfts['age'] = row['age']
    #     pfts['bmi'] = row['bmi']
    #     pfts['fvc'] = row['fvc']
    #     pfts['dlco'] = row['dlco']
    #     pfts['tlc'] = row['tlc']
    #     pfts['rv'] = row['rv']
    #     pfts['frc'] = row['frc']
    #     pfts['fev1'] = row['fev1']
    #     pfts['vc'] = row['vc']
    #     pfts['pefr'] = row['pefr']
    #     pfts['rvtlc'] = row['rvtlc']
    #
    #     ei = EiMLR(pfts['age'], pfts['fvc'], pfts['dlco'], pfts['bmi'], pfts['rvtlc'], pfts['tlc'])
    #     ee = EeMLR(pfts['age'], pfts['frc'], pfts['rv'], pfts['rvtlc'], pfts['dlco'], pfts['fev1'], pfts['pefr'],
    #                pfts['vc'])
    #
    #     eim1 = ei.predict_m1()
    #     eim2 = ei.predict_m2()
    #     eim3 = ei.predict_m3()
    #
    #     eem1 = ee.predict_m1()
    #     eem2 = ee.predict_m2()
    #     eem3 = ee.predict_m3()
    #
    #     if config["volume"] == 'Insp':
    #         weights = [eim1, eim2, eim3]
    #         sf.export_to_cm(pmesh, weights, name=pfts['subject'], lung='L', show_mesh=False)
    #         sf.export_to_cm(pmesh, weights, name=pfts['subject'], lung='R', show_mesh=False)
    #     elif config["volume"] == 'Expn':
    #         weights = [eem1, eem2, eem3]
    #         sf.export_to_cm(pmesh, weights, name=pfts['subject'], lung='L', show_mesh=False)
    #         sf.export_to_cm(pmesh, weights, name=pfts['subject'], lung='R', show_mesh=False)

    # scores, mahalanobis_distance = _get_score(sf, aligned_mesh_names)

    return sf, pmesh


if __name__ == "__main__":
    sf, pmesh = main()
