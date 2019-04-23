##########
# A Model for Teaching Diffusion of a Drug in a Hydrogel for C. elegans
# Author: Hamilton A. White
# Ph.D. Student: WPI/UMass Medical School Joint Program
# Date: 4/21/19
#
# Adapted From
#   Originally a temperature diffusion program developed by Garbel Nervadof [3/21/16, URI: https://www.codeproject.com/Articles/1087025/Using-Python-to-Solve-Computational-Physics-Proble, CPOL v. 1.02]
#
# Copyright:
#   (c) 2019 - present: Hamilton A. White, All Rights Reserved (except were applicable by law).
#   The original code was from Garbel Nervadof (see above).  Source code for his work can be found at the URL above, and is under a CPOL v. 1.02 license.
#   Modifications to the original work in creating this derivation are the responsibility of H.A.W.
#
# About the program:
#   In essence, data about the hydrogel is entered into the program, along with information as to the number of iterations to run and size of the output.  When
#   the program begins, it calculates the expected diffusion of materials into the hydrogel, and then outputs that data as both a plot and a data file.
#   We take and output information on the concentration at the following points around the worm (which we assume to be a cylinder, displayed in cross-section):
#
#           B
#       A   O   C
#           D
#
#   In this diagram, "O" represents the worm that is in the exact center of the cube-shaped hydrogel channel.
#
# WARNING:
#   This program outputs a lot of files with every run.  Please run each simulation in a unique folder with no other files to avoid losing data.
#   This program is also a work in progress, and is in the process of constant development.
#
#
# Important Papers for Experimental Data:
#   Group 2: https://doi.org/10.1371/journal.pone.0063283
#   Group 3: https://doi.org/10.1039/C6RA06130C
#   Group 4: https://doi.org/10.1038/s41598-017-13755-9
#   Group 5: https://doi.org/10.3390/ijms19051491
#
# Student Goal: Get to an 85% concentration of chemical at the worm surface in 200 iterations!
#
##########



# Requirements
import numpy as np
import matplotlib.pyplot as plt
import csv

# Hydrogel Parameters
Et = 0.2 # Porosity, a measure of the amount of open space within the hydrogel on a range of 0 --> 1.
drug_size = 4.5e-7 # Drug size in centimeters
pore_size = 9.5e-6 # Size of the pore, in centimeters. Citation: S. Lin, N. Sangaj, T. Razafiarison, C. Zhang, and S. Varghese, “Influence of Physical Properties of Biomaterials on Cellular Behavior,” Pharm Res, vol. 28, no. 6, pp. 1422–1430, Jun. 2011.
tor = 2.2 # Tortuosity: a measure of the amount of twisting in the path through a porous substrate.  This value was calculated from an average of the ulk tortuosity found in the literature.  Citation: P. Aggarwal et al., “Correlation of chromatographic performance with morphological features of organic polymer monoliths,” Journal of Chromatography A, vol. 1334, pp. 20–29, Mar. 2014.
diffusivity_NO = 2.6e-5 # cm^2*s^-1 Nitric Oxide
diffusivity_EtOH = 0.84e-5 # cm^2*s^-1 Ethanol
diffusivity_CO2 = 1.92e-5 # cm^2*s^-1 Carbon Dioxide
diffusivity_N2 = 1.88e-5 # cm^2*s^-1 Nitrogen
diffusivity_acetone = 1.16e-5 # cm^2*s^-1 Acetone
de = drug_size/pore_size # Constrictivity: a measure of the ratio of the diffusing particle size to the size of the pore.

#D = (diffusivity * Et * de) / (tor)


# Set maximum iterations and define all times that the system will be iterated through
end = 500
start = 0
d = 20
data = []
dict = {}

for n in range(start,end+1,d):
    iter= range(start,end+1,d)
    counts = len(iter)
    print("Please wait for a moment")
    
    # Inputs
    diff = 10
    maxIter = n
    lenX = lenY = 60 # This is a hard-coded value that shouldn't be changed, since the size of the worm in this program and the concentration data that the program outputs is dependent on this size.
    delta = 1

    # Define output sample positions, based on size of the output plot defined in the input section
    wormA_x = int(lenX/2 - 0.416*lenY)
    wormA_y = int(lenY/2)
    wormB_x = int(lenX/2)
    wormB_y = int(lenY/2 + 0.416*lenY)
    wormC_x = int(lenX/2 + 0.416*lenY)
    wormC_y = int(lenY/2)
    wormD_x = int(lenX/2)
    wormD_y = int(lenY/2 - 0.416*lenY)
    pos = np.matrix([[wormA_x,wormA_y],[wormB_x,wormB_y],[wormC_x,wormC_y],[wormD_x,wormD_y]])
    print("Positions A, B, C, and D measured are:" + str(pos))
    print("If the positions above are not: [[ 5 30], [30 54], [54 30], [30  5]], Then reset lenX = lenY statement in the code above to be a value of 60.")

    # Boundary condition
    Ttop = 100
    Tbottom = 0
    Tleft = 100
    Tright = 100

    # Initial guess of interior grid
    Tguess = 0

    # Set colour interpolation and colour map
    colorinterpolation = 100
    colourMap = plt.cm.jet #you can try: colourMap = plt.cm.coolwarm

    # Set meshgrid
    X, Y = np.meshgrid(np.arange(0, lenX), np.arange(0, lenY))

    # Set array size and set the interior value with Tguess
    T = np.empty((lenX, lenY))
    T.fill(Tguess)

    # Set Boundary condition
    T[(lenY-1):, :] = Ttop
    T[:1, :] = Tbottom
    T[:, (lenX-1):] = Tright
    T[:, :1] = Tleft

    # Iteration (We assume that the iteration is convergence in maxIter = 500)
    for iteration in range(0, maxIter):
        for i in range(1, lenX-1, delta):
            for j in range(1, lenY-1, delta):
                T[i, j] = 0.25 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1])
                

    print("Iteration #" + str(n) + " Finished")

    # Configure the contour
    plt.title("Contour of Temperature %s" % (n))
    plt.contourf(X, Y, T, colorinterpolation, cmap=colourMap)

    # Set Colorbar
    plt.colorbar()

    # Create circle representing C. elegans
    circle1=plt.Circle((lenX/2,lenY/2),0.416*lenY,color='black')
    concentrations = np.array([T[wormA_y,wormA_x],T[wormB_y,wormB_x],T[wormC_y,wormC_x],T[wormD_y,wormD_x]]).tolist() # This is backwards from your general perception of X & Y due to the order of concentration calculation above in T[i,j].
    print("The concentrations at A, B, C, and D are: " + str(concentrations))
    plt.gcf().gca().add_artist(circle1)
    plt.savefig('diffusion.t%s.jpg' % (n), transparent=True, dpi=600)

    # Show the result in the plot window
    plt.show(block=False)
    plt.pause(0.25)
    plt.close()
    print("")

    #Export data files
    data.append([n,concentrations])
    print("")

np.savetxt("diffdata.txt", data, fmt='%s', delimiter=",")
    

            

    
    



