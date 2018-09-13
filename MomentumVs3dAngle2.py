from ROOT import *
import os
import math
import argparse
import numpy as np
from array import array



def frange(start, stop, step):
    x = start
    while x < stop:
        yield x
        x += step


###################################################################
####################       Get TTree       ########################
###################################################################
fileName = "HighlandFormula_histo.root"
inFile   = TFile.Open(fileName)
tTree    = inFile.Get("HighlandFormula/trackResTree")

# Quick plot of the 2D Mom Vs Angle distribution
cMomVsAngle = TCanvas("cMomVsAngle" ,"cMomVsAngle" ,0 ,0 ,600 ,600)
cMomVsAngle.cd()
tTree.Draw("wcP:theta_3d","theta_3d > 0 && theta_3d < 1","colz")


## Momentum range: 400 - 1200 MeV
MomentumHistos_List = []
maxMomentum         = 1200.
minMomentum         = 400.
binSize             = 50.
nBinMom             = int((maxMomentum-minMomentum)/binSize)

for i in xrange(nBinMom):
    lowerBound = int(minMomentum + i*binSize)
    upperBound = int(minMomentum + (i+1)*binSize)
    theta3D_Temp = TH1D("theta_"+str(lowerBound)+"_"+str(upperBound)+"MeV", "theta_"+str(lowerBound)+"_"+str(upperBound)+"MeV", 300, 0., 0.03 )
    MomentumHistos_List.append(theta3D_Temp)

#print MomentumHistos_List

print "Entries: ", tTree.GetEntry()
sillyCount =0 
for entry in tTree:
     sillyCount += 1
     if not sillyCount % 10000:
          print sillyCount

     momentum   = entry.wcP
     angle3D    = entry.theta_3d
     angle3D_2  = angle3D*angle3D


     # Find the right histo to fill
     momentumBin = int(momentum/binSize - minMomentum/binSize)
     # Let's try to avoid stupid shit
     if momentumBin > 0 and momentumBin < nBinMom:
         MomentumHistos_List[momentumBin].Fill(angle3D_2)



highlandPlot = TH1D("HighlandFormula", "highlandFormula", int(maxMomentum/binSize), 0., maxMomentum )


myExpo = TF1("expo","expo",0,0.01)
sigma_List    = []
sigmaErr_List = []

print
print
for h in MomentumHistos_List:
    sigma    = -1000.
    sigmaErr = -1000.

    h.Fit(myExpo)
    if myExpo.GetParameter(1):
        par1         = myExpo.GetParameter(1)
        par1Err      = myExpo.GetParError(1)
        sigma        = TMath.Sqrt(1./(-2.*par1))
        sigmaErr     = (sigma /2.) * (par1Err/par1)  # simple error propagation sigma = sqrt(  par1 / 2 )
#        print "expo par: ", myExpo.GetParameter(1)," +- ", myExpo.GetParError(1)
#        print "par1    : ", par1," +- ", par1Err


    sigma_List   .append(sigma)
    sigmaErr_List.append(sigmaErr)

    # Estimate the rigth bin in highland formula
    titleStrings  =  (h.GetTitle()).split("_")
    fake_momentum = float(titleStrings[1]) + binSize/2.
    momentumBin   = int(fake_momentum/binSize) +1
    print titleStrings, fake_momentum, momentumBin, sigma
    highlandPlot.SetBinContent(momentumBin, sigma)
    highlandPlot.SetBinError(momentumBin, sigmaErr)

print

for i in xrange(len(MomentumHistos_List)):
    print (MomentumHistos_List[i]).GetTitle(), (MomentumHistos_List[i]).GetEntries(), sigma_List[i], sigmaErr_List[i]

print
print

# Calculating the true highland formula
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
#gr . GetXaxis().SetRangeUser(0,1000)
#gr . GetYaxis().SetRangeUser(0,2.)
gr . SetLineWidth(2) ;
gr . SetLineColor(kRed) ;
gr . SetFillColor(0)





outFile = TFile("MomVsTheta3d2_Highland.root","recreate")
outFile.cd()

for h in MomentumHistos_List:
    h.Write()

highlandPlot.Write()
gr.Write()
outFile.Write()
outFile.Close()

#

'''
# Comparison between plots

cRaw.Divide(3,1)
pRaw1 = cRaw.cd(1)
pRaw1.SetGrid()
recoData_Int.Draw("pe")

pRaw2 = cRaw.cd(2)
pRaw2.SetGrid()
recoData_Inc.Draw("pe")
cRaw.Update()

pRaw3 = cRaw.cd(3)
pRaw3.SetGrid()
rawXS = recoData_Int.Clone("rawXS")
rawXS.Sumw2()
rawXS.Divide(recoData_Inc)
rawXS.Scale(101.10968)
rawXS.Draw("pe")

cRaw.Update()
'''

raw_input()  



