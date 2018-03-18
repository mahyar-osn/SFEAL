import numpy
import morphic


class SSM(object):

    def __init__(self):

        self.X = []
        self.input_mesh = None
        self.groups = None
        self.mesh = None
        self.pcamesh = None
        self.pmesh = None
        self.score_0 = None
        self.z_score = []
        self.ratio = {}
        self.fname = None
        self.mean = None
        self.score_z = None
        self.SD = None
        self.nodes = None
        self.new_data = []
        self.weights = []
        self.lung = None

    def add_mesh(self, mesh, index=0):
        mesh = morphic.Mesh(str(mesh))
        if self.input_mesh is None:
            self.input_mesh = mesh
        if isinstance(mesh, str):
            mesh = morphic.Mesh(mesh)
        x = []
        if self.groups is None:
            for node in mesh.nodes:
                if not isinstance(node, morphic.mesher.DepNode):
                    x.extend(node.values.flatten().tolist())

        else:
            for node in mesh.nodes:
                if node.in_group(self.groups):
                    x.extend(node.values.flatten().tolist())

        self.X.append(x)

    def pca_setup(self, num_modes=2):
        from sklearn import decomposition
        self.X = numpy.array(self.X)
        self.num_modes = num_modes
        self.pca = decomposition.PCA(n_components=num_modes)
        self.pca.fit(self.X)
        self.mean = self.pca.mean_
        self.components = self.pca.components_.T
        self.variance = self.pca.explained_variance_
        self.generate_mesh()
        return self.mesh, self.X

    def generate_mesh(self):
        self.mesh = morphic.Mesh()
        weights = numpy.zeros(self.num_modes + 1)
        weights[0] = 1.0
        self.mesh.add_stdnode('weights', weights)
        variance = numpy.zeros(self.num_modes + 1)
        variance[0] = 1.0
        variance[1:] = numpy.sqrt(self.variance)
        self.mesh.add_stdnode('variance', variance)
        idx = 0
        if self.groups is None:
            for node in self.input_mesh.nodes:
                nsize = node.values.size
                x = self.get_pca_node_values(node, idx)
                self.mesh.add_pcanode(node.id, x, 'weights', 'variance', group='pca')
                idx += nsize

        else:
            for node in self.input_mesh.nodes:
                nsize = node.values.size
                if node.in_group(self.groups):
                    x = self.get_pca_node_values(node, idx)
                    self.mesh.add_pcanode(node.id, x, 'weights', 'variance', group='pca')
                    idx += nsize
                else:
                    if isinstance(node, morphic.mesher.StdNode):
                        self.mesh.add_stdnode(node.id, node.values)
                    elif isinstance(node, morphic.mesher.DepNode):
                        self.mesh.add_depnode(node.id, node.element, node.node,
                                              shape=node.shape, scale=node.scale)
                    if isinstance(node, morphic.mesher.PCANode):
                        raise Exception('Not implemented')

        for element in self.input_mesh.elements:
            self.mesh.add_element(element.id, element.basis, element.node_ids)

        self.mesh.generate()

    def get_pca_node_values(self, node, idx):
        nsize = node.values.size
        if len(node.shape) == 1:
            pca_node_shape = (node.shape[0], 1, self.num_modes)
            x = numpy.zeros((node.shape[0], 1, self.num_modes + 1))
            x[:, 0, 0] = self.mean[idx:idx + nsize].reshape(node.shape)
            x[:, :, 1:] = self.components[idx:idx + nsize, :].reshape(pca_node_shape)
            return x
        if len(node.shape) == 2:
            pca_node_shape = (node.shape[0], node.shape[1], self.num_modes)
            x = numpy.zeros((node.shape[0], node.shape[1], self.num_modes + 1))
            x[:, :, 0] = self.mean[idx:idx + nsize].reshape(node.shape)
            x[:, :, 1:] = self.components[idx:idx + nsize, :].reshape(pca_node_shape)
            return x
        print 'Cannot reshape this node when generating pca mesh'

    def calculate_score(self, mesh_file_names, mesh_file, save=False):
        if not self.new_data:
            pass
        else:
            self.new_data = []
        subject_name = mesh_file
        print '\n\t=========================================\n'
        print '\t   Please wait... \n'
        totalSubjects = 0
        subjectStore = []
        x = []
        for i in range(len(mesh_file_names)):
            single_mesh = mesh_file_names[i]
            mesh = morphic.Mesh(str(single_mesh))
            totalSubjects += 1
            subjectStore.append(mesh_file_names[i])
            nodes = mesh.get_nodes()
            size = len(nodes)
            for node in mesh.nodes:
                x.extend(node.values)
                X1 = numpy.asarray(x)

        X = X1.reshape((totalSubjects, size * 12))
        print '\t   Total number of subjects in pca = %d' % totalSubjects
        num_modes = totalSubjects - 1
        from sklearn import decomposition
        pca = decomposition.PCA(n_components=num_modes)
        pca.fit(X)
        pca_mean = pca.mean_
        pca_mean = pca_mean.reshape((1, size * 12))
        pca_components = pca.components_.T
        pca_variance = pca.explained_variance_
        pca_explained_variance = pca.explained_variance_ratio_
        self.ratio = {'MODE_{} RATIO'.format(m + 1): '{:.2f}'.format(float(pca_explained_variance[m]))
                      for m in range(len(pca_explained_variance))}
        dataset = dict()
        for i in range(len(subjectStore)):
            dataset.update({subjectStore[i]: i})

        count = len(pca_variance)
        mode_count = []
        for i in range(len(pca_variance)):
            mode_count.append(i + 1)

        print '\t   Total modes of variation = %d' % count
        print '\t   Projecting Subject: %s' % subject_name
        mode_scores = []
        for j in range(len(dataset)):
            subject = X[j] - pca_mean
            score = numpy.dot(subject, pca_components)
            mode_scores.append(score[0][0:count])

        self.SD = numpy.std(mode_scores, axis=0)
        self.mean = numpy.average(mode_scores, axis=0)
        number = dataset[subject_name]
        subject_0 = X[number] - pca_mean
        self.score_0 = numpy.dot(subject_0, pca_components)
        self.score_0 = self.score_0[0][0:count]
        self.score_z = self.convert_scores(self.score_0, self.SD, self.mean)

        self.score_z = {'MODE_{} SCORE'.format(m + 1): '{:.2f}'.format(float(self.score_z[m]))
                      for m in range(len(self.score_z))}

        if save:
            # save scores to a csv file

            print '\n\t   ______________\n'
            print '\t   SORRY!'
            print '\t   SAVE OPTION NOT IMPLEMENTED YET'
            print '\n\t   ______________\n'

            """
            TODO

            import os
            import csv
            if self.fname == None:
                print '\n\t   ______________\n'
                self.fname = raw_input('\t   PLEASE NAME YOUR FILE: ')
            else:
                output_dir = 'output/'
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                weights_file = output_dir+'%s.csv' % self.fname
                if os.path.isfile(weights_file):
                    print '\n\t   ______________\n'
                    print '\t   File already exists'
                    print '\t   Appending...\n'
                else:
                    weights_file = output_dir+'%s.csv' % self.fname
                    text_file = open(weights_file, "a")
                    text_file.write("Mesh")
                    text_file.close()
                    for i in range(len(mode_count)):
                        f = open(weights_file)
                        data = [item for item in csv.reader(f)]
                        f.close()
                        new_column = mode_count
                        self.new_data = []
                        for j, item in enumerate(data):
                            try:
                                item.append(new_column[i])
                            except IndexError, e:
                                item.append("placeholder")
                            self.new_data.append(item)
                        f = open(weights_file, 'w')
                        csv.writer(f).writerows(self.new_data)
                        f.close()
                if not self.new_data:
                    pass
                else:
                    self.new_data = []
                # saving = True
                # if saving:
                text_file = open(weights_file, "a")
                text_file.write("%s" % subject_name)
                text_file.close()

                # if os.path.isfile(weights_file):
                #     self.new_data = []
                #     new = []
                #     for i in range(len(mode_count)):
                #         f = open(weights_file)
                #         data = [item for item in csv.reader(f)]
                #         f.close()
                #
                #         data = [data[-1]]
                #         for j, item in enumerate(data):
                #             try:
                #                 item.append("%.2f" % self.score_z[i])
                #             except IndexError, e:
                #                 item.append("placeholder")
                #             self.new_data.append(item)
                #         f = open(weights_file, 'a')
                #         csv.writer(f).writerows(self.new_data)
                #         f.close()
                print "\t   Scores saved to %s \n" % weights_file
                print '\t   ______________'
                
            """
        print '\n\t=========================================\n'
        return self.score_z, self.ratio

    def convert_scores(self, scores, SD, mean):
        for i in range(len(scores)):
            self.z_score.append((scores[i] - mean[i]) / SD[i])

        return self.z_score

    def export_to_cm(self, pmesh, weights, name='default', lung='l', show_mesh=False):
        if not self.weights:
            pass
        else:
            self.weights = []

        if lung == 'R':
            self.lung = 'Right'
        elif lung == 'L':
            self.lung = 'Left'
        else:
            raise Exception("'lung' argument can ONLY be L OR R!")

        self.weights = weights
        self.pmesh = pmesh
        self.pmesh.nodes['weights'].values[1:] = 0  # reset weights to zero
        self.pmesh.nodes['weights'].values[0] = 1  # adding average
        for numMode in range(len(self.weights)):
            self.pmesh.nodes['weights'].values[numMode + 1] = self.weights[numMode]
            self.pmesh.update_pca_nodes()

        # saving
        import os
        import pandas as pd
        import subprocess
        from useful_files import nodes

        if self.nodes is not None:
            print "SHOULDN'T BE HERE!"
            self.nodes = None

        self.nodes = nodes.Nodes()

        output_dir = 'output/export_to_cm/%s' % name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        input_folder = '../useful_files'

        path_to_export_mesh = output_dir
        temp_file = '%s_reconstructed_temp.csv' % self.lung
        output_file = '%s_reconstructed' % self.lung
        save_temp_file = output_dir + '/%s' % temp_file
        save_output_file = output_dir + '/%s' % output_file
        ipnode_file = '%s_reconstructed' % self.lung
        path_to_ipnode_file = '%s/%s' % (output_dir, ipnode_file)

        path_to_com_file = os.path.join(os.path.dirname(__file__),
                                        input_folder, 'perl_com', 'ipnode2exnode.com')
        path_to_cmgui_file = os.path.join(os.path.dirname(__file__),
                                          input_folder, 'perl_com', 'cmgui.com')
        ip2ex_perl = os.path.join(os.path.dirname(__file__),
                                  input_folder, 'perl_com', 'ipnode2exnode.pl')
        ip2ex_cm = os.path.join(os.path.dirname(__file__),
                                input_folder, 'perl_com', 'ipnode2exnode')
        cmgui_file = os.path.join(os.path.dirname(__file__),
                                  input_folder, 'perl_com', 'cmgui')
        param_file = os.path.join(os.path.dirname(__file__),
                                  input_folder, 'perl_com', '3d_fitting')
        versions_file = os.path.join(os.path.dirname(__file__),
                                     input_folder, 'perl_com', 'versions')
        base_file = os.path.join(os.path.dirname(__file__),
                                 input_folder, 'perl_com', 'BiCubic_Surface_Unit')

        if self.lung == 'Right':
            node_file = 'nodes_%s.csv' % self.lung
            input_file = os.path.join(os.path.dirname(__file__),
                                      input_folder, node_file)
            elem_file = os.path.join(os.path.dirname(__file__),
                                     input_folder, 'perl_com', 'templateRight')

            nodes = self.nodes.set_nodes(lung='right')
            for node_number in nodes:
                node = self.pmesh.nodes[node_number]
                nodeValues = node.values
                with open(save_temp_file, 'a') as f:
                    numpy.savetxt(f, nodeValues)

            a = pd.read_csv(input_file)
            b = pd.read_csv(save_temp_file, delimiter=' ')
            result = pd.concat([a, b], axis=1)
            result.to_csv('%s.csv' % save_output_file, sep=' ', index=False)
            os.remove(save_temp_file)
            py2ip_right_file = os.path.join(os.path.dirname(__file__),
                                            input_folder, 'perl_com', 'py2ip_right.pl')

            subprocess.call(["perl", py2ip_right_file, "%s.csv" % save_output_file,
                             "%s.ipnode" % path_to_ipnode_file])

            with open(path_to_com_file, 'wb') as comfile:
                comfile.write(" set echo on;\n")
                comfile.write(" fem def param;r;{0}".format("%s;\n" % param_file))
                comfile.write(" fem def coor;r;{0}".format("%s;\n" % versions_file))
                comfile.write(" fem def base;r;{0}".format("%s;\n" % base_file))
                comfile.write(" fem def node;r;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem def elem;r;{0}".format("%s;\n" % elem_file))
                comfile.write(" fem export node;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem export elem;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem def node;w;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem quit;\n")

            if show_mesh:
                with open(path_to_cmgui_file, 'wb') as comfile:
                    comfile.write(" gfx read node {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx read elem {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx cre egroup fissure;\n")
                    comfile.write(" gfx mod egroup fissure add 51..62;\n")
                    comfile.write(
                        " gfx mod g_e {0} general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e {0} lines coordinate coordinates select_on material green selected_material default_selected;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e fissure general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n")
                    comfile.write(" gfx mod g_e fissure surfaces material tissue;\n")
                    comfile.write(" gfx edit scene;\n")
                    comfile.write(" gfx cre win;\n")
                show_cmgui = 'show'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui],
                                shell=False)

            else:
                show_cmgui = 'no'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui], shell=True)

            print "\n\t=========================================\n"
            print "\t   ALL MESH FILES EXPORTED TO:"
            print "\n\t\t   %s " % path_to_export_mesh
            print "\n\t=========================================\n"

            os.remove(save_output_file + '.csv')

        elif self.lung == 'Left':
            node_file = 'nodes_%s.csv' % self.lung
            input_file = os.path.join(os.path.dirname(__file__),
                                      input_folder, node_file)
            elem_file = os.path.join(os.path.dirname(__file__),
                                     input_folder, 'perl_com', 'templateLeft')

            nodes = self.nodes.set_nodes(lung='left')
            for node_number in nodes:
                node = self.pmesh.nodes[node_number]
                nodeValues = node.values
                with open(save_temp_file, 'a') as f:
                    numpy.savetxt(f, nodeValues)

            a = pd.read_csv(input_file)
            b = pd.read_csv(save_temp_file, delimiter=' ')
            result = pd.concat([a, b], axis=1)
            result.to_csv('%s.csv' % save_output_file, sep=' ', index=False)
            os.remove(save_temp_file)

            py2ip_left_file = os.path.join(os.path.dirname(__file__),
                                           input_folder, 'perl_com', 'py2ip_left.pl')
            subprocess.call(["perl", py2ip_left_file, "%s.csv" % save_output_file,
                             "%s.ipnode" % path_to_ipnode_file])

            with open(path_to_com_file, 'wb') as comfile:
                comfile.write(" set echo on;\n")
                comfile.write(" fem def param;r;{0}".format("%s;\n" % param_file))
                comfile.write(" fem def coor;r;{0}".format("%s;\n" % versions_file))
                comfile.write(" fem def base;r;{0}".format("%s;\n" % base_file))
                comfile.write(" fem def node;r;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem def elem;r;{0}".format("%s;\n" % elem_file))
                comfile.write(" fem export node;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem export elem;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem def node;w;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem quit;\n")

            if show_mesh:
                with open(path_to_cmgui_file, 'wb') as comfile:
                    comfile.write(" gfx read node {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx read elem {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx cre egroup fissure;\n")
                    comfile.write(" gfx mod egroup fissure add 111..118;\n")
                    comfile.write(
                        " gfx mod g_e {0} general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e {0} lines coordinate coordinates select_on material green selected_material default_selected;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e fissure general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n")
                    comfile.write(" gfx mod g_e fissure surfaces material tissue;\n")
                    comfile.write(" gfx edit scene;\n")
                    comfile.write(" gfx cre win;\n")
                show_cmgui = 'show'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui],
                                shell=False)

            else:
                show_cmgui = 'no'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui], shell=True)

            print "\n\t=========================================\n"
            print "\t   ALL MESH FILES EXPORTED TO:"
            print "\n\t\t   %s " % path_to_export_mesh
            print "\n\t=========================================\n"

            os.remove(save_output_file + '.csv')

        return None

    def project_new_mesh(self, mesh_file_names, mesh_file):
        if not self.new_data:
            pass
        else:
            self.new_data = []
        subject_name = mesh_file
        print '\n\t=========================================\n'
        print '\t   Please wait... \n'
        totalSubjects = 0
        subjectStore = []
        x = []
        y = []
        for i in range(len(mesh_file_names)):
            single_mesh = mesh_file_names[i]
            mesh = morphic.Mesh(str(single_mesh))
            totalSubjects += 1
            subjectStore.append(mesh_file_names[i])
            nodes = mesh.get_nodes()
            size = len(nodes)
            for node in mesh.nodes:
                x.extend(node.values)
                X1 = numpy.asarray(x)

        X = X1.reshape((totalSubjects, size * 12))
        print '\t   Total number of subjects in pca = %d' % totalSubjects
        num_modes = totalSubjects - 1
        from sklearn import decomposition
        pca = decomposition.PCA(n_components=num_modes)
        pca.fit(X)
        pca_mean = pca.mean_
        pca_mean = pca_mean.reshape((1, size * 12))
        pca_components = pca.components_.T
        pca_variance = pca.explained_variance_
        pca_explained_variance = pca.explained_variance_ratio_
        self.ratio = {'MODE_{} RATIO'.format(m + 1): '{:.2f}'.format(float(pca_explained_variance[m]))
                      for m in range(len(pca_explained_variance))}
        dataset = dict()
        for i in range(len(subjectStore)):
            dataset.update({subjectStore[i]: i})

        count = len(pca_variance)
        mode_count = []
        for i in range(len(pca_variance)):
            mode_count.append(i + 1)

        print '\t   Total modes of variation = %d' % count
        print '\t   Projecting Subject: %s' % subject_name

        mode_scores = []
        for j in range(len(dataset)):
            subject = X[j] - pca_mean
            score = numpy.dot(subject, pca_components)
            mode_scores.append(score[0][0:count])

        if self.SD is not None:
            self.SD = None
        if self.mean is not None:
            self.mean = None

        self.SD = numpy.std(mode_scores, axis=0)
        self.mean = numpy.average(mode_scores, axis=0)

        project_mesh_path = mesh_file
        project_mesh = morphic.Mesh(project_mesh_path)
        for node in project_mesh.nodes:
            y.extend(node.values)
            Y = numpy.asarray(y)

        Y = Y.reshape((size * 12))
        subject_0 = Y - pca_mean

        if self.score_0 is not None:
            self.score_0 = None
        if self.score_z is not None:
            self.score_z = None

        self.score_0 = numpy.dot(subject_0, pca_components)
        self.score_0 = self.score_0[0][0:count]
        self.score_z = self.convert_scores(self.score_0, self.SD, self.mean)

        self.score_z = {'MODE_{} SCORE'.format(m + 1): '{:.2f}'.format(float(self.score_z[m]))
                      for m in range(len(self.score_z))}

        print '\n\t=========================================\n'
        return self.score_z, self.ratio

class MESH(object):

