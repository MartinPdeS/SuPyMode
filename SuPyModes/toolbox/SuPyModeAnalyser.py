#-------------------------Importations------------------------------------------
""" standard imports """
import numpy as np
import copy

""" package imports """
global config, args
import config.config as config
import config.args as args
#-------------------------Importations------------------------------------------


class ModeTracker(object):
    """ This class analyse two set of modes in order to compare them whith
    their overlap and then arrange their label in order to make them
    coincide.

    """

    def __init__(self, data, iter, Fields):



        self._data = data

        self.iter = iter

        self.Fields = Fields


    def compare_label(self):
        """
        The method construct the dictionnay with the previous and new labels
        added to the overlap values.

        returns:
            : param  self.New_Fields : Fields with corrected labels.
            : type self.New_Fields : dict

        """

        for new_label in range(config.N_solution):
            rec = []
            for old_label in range(config.N_solution):

                overlap = np.abs( np.sum(self.Fields.loc[(self.iter-1,new_label)].values * self.Fields.loc[(self.iter-2,old_label)].values ) )

                self._data.loc[(self.iter - 1, old_label), ('overlap', new_label)] = overlap

                rec.append(overlap)

            if np.max(rec) < 0.5:
                print('WARNING: Bad field concordance [{0:.3f}], ITR evolution is probably too fast'.format(np.max(rec)))

            best_fit_label = np.argmax(rec)

            if best_fit_label != new_label:
                self.replace_best_fit(new_label, best_fit_label)



        return self._data, self.Fields


    def replace_best_fit(self, new_label, best_fit_label):
        """
        The method find the best fit for the new mode labels (nomenclature
        start with "_")

        """

        if not args.quiet:

            print('Mode swapping happening')

        old_row = copy.copy( self._data.loc[(self.iter-1,best_fit_label)] )

        new_row = copy.copy( self._data.loc[(self.iter-1,new_label)] )

        self._data.loc[(self.iter-1,best_fit_label)] = new_row

        self._data.loc[(self.iter-1,new_label)] = old_row

        old_field = copy.copy(self.Fields.loc[(self.iter-2,best_fit_label)])

        new_field = copy.copy( self.Fields.loc[ (self.iter-1, new_label) ] )

        self.Fields.loc[ (self.iter-1, new_label) ] = old_field

        self.Fields.loc[ (self.iter-2,best_fit_label) ] = new_field




class SuPyAnalyser(object):
    "NOT USED"

    def __init__(self):
        pass



    def make_circle(self, m, rayon):
        """
        arguments:

            : param m:
            : type m: .

            : param rayon:
            : type rayon: .

        returns:
            : param :
            : type :



        calls:
            :call1:

        """

        d = 2*m + 1

        rayon_x, rayon_y = d/2, d/2

        x, y = np.indices((d, d))

        return (np.abs(np.hypot(rayon_x - x, rayon_y - y)-rayon) < 0.5).astype(int)


    def rotate_bound(self, image, angle):
        """
        arguments:

            : param image:
            : type image:

            :param angle:
            :type angle:


        returns:
            : param :
            : type :



        calls:
            :call1:

        """


        image = self.vector[image]

        (h, w) = image.shape[:2]

        (cX, cY) = (w // 2, h // 2)

        M = cv2.getRotationMatrix2D((cX, cY), -angle, 1.0)

        cos = np.abs(M[0, 0])

        sin = np.abs(M[0, 1])

        # compute the new bounding dimensions of the image
        nW = int((h * sin) + (w * cos))

        nH = int((h * cos) + (w * sin))

        # adjust the rotation matrix to take into account translation
        M[0, 2] += (nW / 2) - cX

        M[1, 2] += (nH / 2) - cY

        # perform the actual rotation and return the image
        return cv2.warpAffine(image, M, (nW, nH))


    def analyse_mode(self, mode = 'mode1'):
        """
        arguments:

            : param mode:
            : type mode:

        returns:
            : param :
            : type :


        calls:
            :call1:
        """

        vec = self.vector[mode]

        middle_index = self.m // 2

        line = vec[:,middle_index]

        edges2 = copy.copy(self.edges)

        for i in range(0, self.m):

            if edges2[0][i] > 0.3*self.m and edges2[0][i] < 0.7*self.m:

                edges2[0][i] = 0

                edges2[1][i] = 0

        a = np.nonzero(self.edges)

        parameter_0 = 0
        parameter_1 = 0

        print('LP%i%i'% (parameter_0, parameter_1))
        return 'LP%i%i'% (parameter_0, parameter_1)


    def analyse_symmetry(self, mode):
        """
        arguments:

            : param mode:
            : type mode: .

        calls:
            :call1:

        """

        a = self.make_circle(sel.m, self.m/2)

        plt.plot(a * self.vector['mode0'])


    def analyse_zeros(self, mode):
        """
        arguments:

            : param mode:
            : type mode:

        returns:
            : param zeros_parameter:
            : type zeros_parameter:

        calls:
            :call1:

        """

        zeros_parameter = 0

        vec = self.vector[mode]

        middle_index = self.m // 2

        line = vec[:,middle_index]

        low_border = 0.1 * self.m

        up_border = 0.9 * self.m

        for j in range(len(line)-1):

            if line[j]*line[j+1] < 0 and j > low_border and j < up_border:

                zeros_parameter += 1

        return zeros_parameter
