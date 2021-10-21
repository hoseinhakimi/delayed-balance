import numpy as np
import random


class TriangularBalance:
    def __init__(self, size, initialRatio, delay):
        # Size of the networ
        self.Size = size
        # how much it takes for agents to sense changes
        self.delay = delay
        # Caclulate the number of triangles connected to one link
        self.TrianglesOnLink = self.LinkTriangles()
        # Calculate the number of all triangles in the fully connected network
        self.Triangles = self.TriangleCount()
        # The portion of friendly link
        self.InitRatio = initialRatio
        # This will initialize the network based on given parameters
        self.NetworkInitiator()

    # region Initial Functions

    # Count all the triangles in network
    def TriangleCount(self):
        numinator = (self.Size) * (self.Size - 1) * (self.Size - 2)
        denominator = 6
        return(numinator/denominator)

    # Count the triangles connected to a link
    def LinkTriangles(self):
        trianglesOnLink = self.Size - 2
        return(trianglesOnLink)

    # This function initializes the network at the begining of the run
    def NetworkInitiator(self):
        # Generate a matrix with random elements between 0 and 1
        tempMatrix = np.random.rand(self.Size, self.Size)
        # Change elements smaller than friendship ratio to -1 (Will be changed to  1 soon)
        tempMatrix[tempMatrix < self.InitRatio] = -1
        # Change elements greater than friendship ratio to  1 (Will be changed to -1 soon)
        tempMatrix[tempMatrix >= self.InitRatio] = 1
        # The Created marix is not symmetric (and does not reperesnt a Graph) so we symmetrize it by this trick
        # Another point is that we change the sign of links here as we promised earlier
        adjMatrix = np.tril(-tempMatrix, -1) + np.tril(-tempMatrix, -1).T
        # Put this matrix to the class InitialNetwork and Network
        self.InitialNetwork = adjMatrix
        # To save the initial state we store the network in another variable
        self.Network = self.InitialNetwork
        # Birth matrix
        self.BirthTime = np.zeros((self.Size, self.Size))
        # Time of system
        self.SystemTime = 1
        # Calculate the Energy of our network
        self.Energy = self.NetworkEnergy()

    # This function calculates the total energy of the netwrok
    def NetworkEnergy(self):
        # Calculate the sum of product of all two pairs on each link
        netLen2Path = np.matmul(self.Network, self.Network)
        # Calculate the energy of all triangles on each link
        # (not exactly energy it needs a multiply by -1 to be energy)
        energyMat = np.multiply(self.Network, netLen2Path)
        # Every link is counted 2 times
        unnormalTotalEnergy = np.sum(energyMat) / 2
        # Every triangle is counted 3 times
        unnormalTotalEnergy = unnormalTotalEnergy / 3
        # We want energy to be between -1 to +1
        totalEnergy = float(-unnormalTotalEnergy) / self.Triangles
        return(totalEnergy)
    # endregion

    # region Dynamics
    # Calculation of the attribution of one link in the total energy of the network
    def LinkEnergy(self, adjTuple):
        # Get the link's sign
        linkSign = self.Network[adjTuple]
        # Adjacent links
        linkRow = self.effectiveAdj[adjTuple[0]]
        linkCol = self.effectiveAdj[adjTuple[1]]
        # The energy of triangles on the link
        linkEng = float(-1.0 * np.inner(linkRow, linkCol)
                        * linkSign) / self.Triangles
        return(linkEng)


    def LinkAuxEnergy(self, adjTuple):
         # Get the link's sign
        linkSign = self.auxNetwork[adjTuple]
        # Adjacent links
        linkRow = self.auxNetwork[adjTuple[0]]
        linkCol = self.auxNetwork[adjTuple[1]]
        # The energy of triangles on the link
        linkEng = float(-1.0 * np.inner(linkRow, linkCol) * linkSign) / self.Triangles
        return(linkEng)

    def LinkEnergyCheck(self, linkEnergy, linkSign):
        # This part is for manipulating E=0 case
        randomSign = random.sample((1, -1), 1)
        addedEnergy = randomSign[0] / (2 * self.Triangles)
        tempEnergy = linkEnergy + addedEnergy
        # **************************************
        # link's new sign and energy change
        tempSign = - np.sign(tempEnergy) * linkSign
        # *********************************
        return(tempSign)

    def forwardOneStep(self):
        lateChange = self.changeTimeLine[-(self.delay+1)]
        self.effectiveAdj[lateChange[0], lateChange[1]] -= lateChange[2]
        self.effectiveAdj[lateChange[1], lateChange[0]] -= lateChange[2]

    def makeFuture(self):
        self.auxNetwork = self.Network.copy()
        for Time in range(self.delay):
            self.LinkBaseDynamics()

    def backwardAdj(self):
        temp = np.copy(self.Network)
        for i in range(self.delay):
            item = self.changeTimeLine[::1][i]
            temp[item[0], item[1]] += item[2]
            temp[item[1], item[0]] += item[2]
        self.effectiveAdj = temp

    def realEnergyChange(self, link, newSign):
        # Get the link's sign
        linkSign = self.Network[link]
        # Adjacent links
        linkRow = self.Network[link[0]]
        linkCol = self.Network[link[1]]
        # The energy of triangles on the link
        linkEng = float(-1.0 * np.inner(linkRow, linkCol)) / self.Triangles
        change = (newSign - linkSign) * linkEng
        return change

    def LinkAuxDynamics(self):
        # choose a random link
        link        = tuple(random.sample(range(0, self.Size-1), 2))
        # get the sign
        linkSign    = self.auxNetwork[link]
        # link energy
        linkEnergy  = self.LinkAuxEnergy(link)
        # check if it will change
        engStat     = self.LinkEnergyCheck(linkEnergy, linkSign)
        # how much the sign should change
        signChange  = linkSign - engStat[1]
        # change system's energy and link's sign
        self.auxNetwork[link]             -= signChange
        self.auxNetwork[link[1]][link[0]] -= signChange

    def LinkDelayedDynamics(self):
        # choose a random link
        link = tuple(random.sample(range(0, self.Size-1), 2))
        # get the sign
        linkSign = self.Network[link]
        # link energy
        linkEnergy = self.LinkEnergy(link)
        # check if it will change
        engStat = self.LinkEnergyCheck(linkEnergy, linkSign)
        # how much the sign should change
        realEChange = self.realEnergyChange(link, engStat)
        signChange = linkSign - engStat
        # change system's energy and link's sign
        self.Energy += realEChange
        self.Network[link] = engStat
        self.Network[link[1]][link[0]] = engStat
        self.changeTimeLine.append([link[0], link[1], signChange])
        self.forwardOneStep()

    def initialRandomChange(self):
        for _ in range(self.delay):
            link = tuple(random.sample(range(0, self.Size-1), 2))
            value = np.random.choice([-2, 0, 2])
            self.changeTimeLine.append([link[0], link[1], value])

       

    def DelayedDynamics(self):
        # dynamics time length
        self.changeTimeLine = []
        self.initialRandomChange()
        self.backwardAdj()
        itterateLength = self.Size ** self.ItterateExp *2
        TimeLine = np.zeros((itterateLength, 2))
        for Time in range(itterateLength):
            TimeLine[Time, 0] = self.Network.mean()
            TimeLine[Time, 1] = self.Energy
            self.LinkDelayedDynamics()
        return(TimeLine)

    def TriadDynamics(self, itterateExp):
        self.ItterateExp = itterateExp
        self.TimeLine = self.DelayedDynamics()
    # endregion
