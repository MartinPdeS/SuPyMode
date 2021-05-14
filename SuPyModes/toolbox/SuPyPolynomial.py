#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Polynomial(object):


    def __init__(self, N, M):
        self.N = N
        self.M = M


    def polynomial_eq(self, degree):
        if degree == 2:
            return self.degree_2

        if degree == 3:
            return self.degree_3

        if degree == 4:
            return self.degree_4

        if degree == 5:
            return self.degree_5

        if degree == 6:
            return self.degree_6


    def degree_2(self, variables):

        VAR = variables
        EQ = []

        for j in range(len(variables)):
            EQ.append(
                      VAR[0]*self.M[j]**1 + \
                      VAR[1]*self.M[j]**0 - \
                      self.N[j]
                      )

        return EQ


    def degree_3(self, variables):

        VAR = variables
        EQ = []

        for j in range(len(variables)):
            EQ.append(
                      VAR[0]*self.M[j]**2 + \
                      VAR[1]*self.M[j]**1 + \
                      VAR[2]*self.M[j]**0  - \
                      self.N[j]
                      )

        return EQ


    def degree_4(self, variables):

        VAR = variables
        EQ = []

        for j in range(len(variables)):
            EQ.append(
                      VAR[0]*self.M[j]**3 + \
                      VAR[1]*self.M[j]**2 + \
                      VAR[2]*self.M[j]**1 + \
                      VAR[3]*self.M[j]**0 - \
                      self.N[j]
                      )

        return EQ


    def degree_5(self, variables):

        VAR = variables
        EQ = []

        for j in range(len(variables)):
            EQ.append(
                      VAR[0]*self.M[j]**4 + \
                      VAR[1]*self.M[j]**3 + \
                      VAR[2]*self.M[j]**2 + \
                      VAR[3]*self.M[j]**1 + \
                      VAR[4]*self.M[j]**0 - \
                      self.N[j]
                      )

        return EQ


    def degree_6(self, variables):

        VAR = variables
        EQ = []

        for j in range(len(variables)):
            EQ.append(
                      VAR[0]*self.M[j]**5 + \
                      VAR[1]*self.M[j]**4 + \
                      VAR[2]*self.M[j]**3 + \
                      VAR[3]*self.M[j]**2 + \
                      VAR[4]*self.M[j]**1 + \
                      VAR[5]*self.M[j]**0 - \
                      self.N[j]
                      )

        return EQ
