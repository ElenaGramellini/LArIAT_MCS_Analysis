from ROOT import *
import os
import math
import argparse
import numpy as np
from array import array


# python loop over floats
def frange(start, stop, step):
    x = start
    while x < stop:
        yield x
        x += step


###################################################################
########     Take file name from the command line   ###############
######## The default is HighlandFormula_histo.root  ###############
###################################################################
parser = argparse.ArgumentParser()
parser.add_argument("fileName"     , nargs='?', default = "HighlandFormula_histo.root", type = str, help="insert fileName")
args     = parser.parse_args()
fileName = args.fileName


###################################################################
####################       Get TTree       ########################
###################################################################
inFile   = TFile.Open(fileName)
tTree    = inFile.Get("HighlandFormula/trackResTree")


### Define the momentum range      ###
### Momentum range: 400 - 1200 MeV ###
MomentumHistos_List = []
maxMomentum         = 1200. # MeV
minMomentum         = 400.  # MeV
binSize             = 50.   # MeV
nBinMom             = int((maxMomentum-minMomentum)/binSize)

# Define the histogram for the theta3D^2 distributions
for i in xrange(nBinMom):
    lowerBound = int(minMomentum + i*binSize)
    upperBound = int(minMomentum + (i+1)*binSize)
    theta3D_Temp = TH1D("theta_"+str(lowerBound)+"_"+str(upperBound)+"MeV", "theta_"+str(lowerBound)+"_"+str(upperBound)+"MeV", 300, 0., 0.03 )
    MomentumHistos_List.append(theta3D_Temp)



print "Entries: ", tTree.GetEntry()
sillyCount =0 
for entry in tTree:
     sillyCount += 1
     if not sillyCount % 10000:
          print sillyCount

     # Get the important variables from the ttree
     momentum   = entry.wcP        # This is the momentum at the WC4 
     eLoss      = entry.energyLoss # This is the energy loss between WC and the TPC Front Face
     angle3D    = entry.theta_3d
     angle3D_2  = angle3D*angle3D


     # TO DO HERE:
     # Calculate the momentum at the point where we are measuring the MCS



     # We need to use the calculation of the momentum to go on 
     # for now, I'm just using the WC mometum. 
     # Find the right histo to fill
     momentumBin = int(momentum/binSize - minMomentum/binSize)
     # Let's try to avoid stupid shit
     if momentumBin > 0 and momentumBin < nBinMom:
         MomentumHistos_List[momentumBin].Fill(angle3D_2)


#Define the money plot
highlandPlot = TH1D("HighlandFormula", "highlandFormula; Momentum [MeV/c]; #sigma_{MCS} [rad]", int(maxMomentum/binSize), 0., maxMomentum )

#Define the fitting function and the list of our parameters
myExpo = TF1("expo","expo",0,0.01)
sigma_List    = []
sigmaErr_List = []

# Loop on the theta3d^2 histograms. 
# For each momentum bin: 
#   1. fit an exponential 
#   2. retreive the slope parameter (alpha in the paper)
#   3. calculate sigma
#   4. calculate the error on sigma from the error on alpha
#   5. Fill the Highland formula plot
for h in MomentumHistos_List:
    sigma    = -1000.
    sigmaErr = -1000.

    h.Fit(myExpo)
    if myExpo.GetParameter(1):
        par1         = myExpo.GetParameter(1)
        par1Err      = myExpo.GetParError(1)
        sigma        = TMath.Sqrt(1./(-2.*par1))
        sigmaErr     = (sigma /2.) * (par1Err/par1)  # simple error propagation sigma = sqrt(  par1 / 2 )



    sigma_List   .append(sigma)
    sigmaErr_List.append(sigmaErr)

    # Estimate the rigth bin in highland formula
    titleStrings  =  (h.GetTitle()).split("_")
    fake_momentum = float(titleStrings[1]) + binSize/2.
    momentumBin   = int(fake_momentum/binSize) +1
    print titleStrings, fake_momentum, momentumBin, sigma
    highlandPlot.SetBinContent(momentumBin, sigma)
    highlandPlot.SetBinError(momentumBin, sigmaErr)



# Calculating the and storing the theoretical highland formula
S2 = 13.6;
c = 299792458;
epsilon = 0.038;
mass = 105.7; #(For charged Muon)                                                                                                                                                                

SigmaExp = []
PExp     = []
zero     = []

for momExp in frange(100,1200,0.1): 
    velocity_numerator   = (momExp/mass)*(momExp/mass);
    velocity_denomenator = 1 +  velocity_numerator /(c*c) ;
    velocity = velocity_numerator/velocity_denomenator;
    beta = velocity / c;
    Term1 = (S2) / (momExp * beta * c);
    Term3 = 1 + epsilon ;
    TestSigma = Term1*Term3
    SigmaExp.append(TestSigma);
    PExp.append(momExp);
    zero.append(0);

g4x      = array('f', PExp)
g4y      = array('f', SigmaExp )
g4exl    = array('f', zero)
g4exr    = array('f', zero)

nPoints=len(g4x)
gr      = TGraphErrors ( nPoints , g4x , g4y     , g4exl, g4exr )
gr.SetTitle("TheoreticalHF; Momentum [MeV/c]; Sigma [rad]")
gr . SetLineWidth(2) ;
gr . SetLineColor(kRed) ;
gr . SetFillColor(0)




# Save everything
outFile = TFile("MomVsTheta3d2_Highland.root","recreate")
outFile.cd()

for h in MomentumHistos_List:
    h.Write()

highlandPlot.Write()
gr.Write()
outFile.Write()
outFile.Close()


raw_input()  



