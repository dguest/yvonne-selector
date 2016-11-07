#!/usr/bin/env python

import ROOT as r
import numpy as np
import glob

from ROOT import gROOT, TCanvas, TF1, TChain

class FatJet:
    def __init__(self, tree, idx):
        # make an object representing the jet at index `idx` by pulling
        # interesting attributes out of the tree.
        self.m = tree.fjet_m[idx]
        self.pt = tree.fjet_pt[idx]
        self.pt_untrim = tree.ut_asso_fjet_pt[idx]
        self.tau21 = tree.fjet_tau21[idx]
        self.tau21_untrim = tree.ut_asso_fjet_tau21[idx]


def rootName (var):
    if var =='Z+js':
        return "*Zqq*root"
    elif var == "W+js":
        return "*Wqq*root"
    elif var == "Bkgnd":
        return "*JZ*root"

    # rewriting using a function:
    # 1. input: a variable called rootname : Z+js,W+js or Bkgnd.
    # 2. Read in the list of files with the designated keyword
    # 3. Create the TChain Tree and output that.

def treeCreation(var):
    rootFilesName = rootName(var)
    #Initialize the fileList
    fileList=[]
    #Creating a list with all the files with the keywords Wqq, Zqq or JZ.
    #fileList=glob.glob(rootFilesName)
    #Defining a TChain
    treeChain=TChain("zpj")
    # Adding rootfiles with the tree into the TChains
    treeChain.Add(rootFilesName)
    # printing an intermediate check statement
    print "Loading Trees of %r. i.e. rootfiles with the keyword %r" %(var,rootName(var))
    # printing the number of events in the tree.
    nEntries =treeChain.GetEntries()
    print "There are %d events in the tree" %nEntries
    # return the treeChain
    return (treeChain)



    # Initializating the histograms
    
    '''H_Tau21P_S1Z = r.TH1F("H_Tau21P_S1Z","Selection1: Tau21P dist;Tau21P;Event",50,-1,3)
    H_Tau21P_S1W = r.TH1F("H_Tau21P_S1W","Selection1: Tau21P dist;Tau21P;Event",50,-1,3)
    H_Tau21P_S1Bg = r.TH1F("H_Tau21P_S1Bg","Selection1: Tau21P dist;Tau21P;Event",50,-1,3)

    H_Tau21P_S2Z = r.TH1F("H_Tau21P_S2Z","Selection2: Tau21P dist;Tau21P;Event",50,-1,3)
    H_Tau21P_S2W = r.TH1F("H_Tau21P_S2W","Selection2: Tau21P dist;Tau21P;Event",50,-1,3)
    H_Tau21P_S2Bg = r.TH1F("H_Tau21P_S2Bg","Selection2: Tau21P dist;Tau21P;Event",50,-1,3)
'''
    # loop thru the Z+js events and do stuff
def selection(treeChain, selectionType):
    #Initialization of the lists:
    S1_jets=[]
    S2_jets=[]
    S1_bg=[]
    S2_bg=[]
    for i,evt in enumerate(treeChain):
        #if i%10==0:
        print "processing event %d" % i
        print "size of the treeChain: %d" %(treeChain.GetEntries())
        # here's one way to do it... //Ignored; using the class method
        # simultaneously loop thru the list of fatjets pT/mass/tau21:
        for (pt, mass, tau21) in zip(evt.fjet_pt, evt.fjet_m, evt.fjet_tau21):
            if pt<400e3: continue
            if mass < 20e3: continue
            print "Event %d: found some jet w/ like high pt and stuff: (%g, %g, %g)" % (i, pt, mass, tau21)

        # but IMO a clearer way is to build up a list of jet "objects".
        # see the very simple class definition at the top of this file.
        n_fatjets = len(evt.fjet_m)
        fatjets_all = [FatJet(treeChain, idx) for idx in xrange(n_fatjets)]

        # calculate the tau21' values for each fatjet
        for jet in fatjets_all:
            # note that in python we can arbitrarilly assign new attributes to an object
        # neglecting jets that has a jet.pt_untrimed = 0
            checkPT=False
            if jet.pt_untrim == 0.:
                checkPT=True
                continue
            else:
                jet.tau21_prime = jet.tau21 + 0.04 * np.log(jet.m**2 / jet.pt_untrim / 1e6)

        # now slim down the list of fatjets by making some selection requirements on pT and mass
        selected = fatjets_all
        selected = filter(lambda j: j.m>20e3, selected)
        selected = filter(lambda j: j.pt>200e3, selected)
        selected = filter(lambda j: j.pt_untrim!=0., selected)

        # printing out information
        print "(Event %d) selected %d of %d fatjets!" % (i, len(selected), n_fatjets)
        print "Length of fatjet_all: %r" %(len(fatjets_all))
        quitLoopingM2=False
        for iselected, jet in enumerate(selected):
            # print the selected jets and compare their tau21 and tau21' values
            print "\t%d: (pT, tau21, tau21') = (%d, %.2f, %.2f) \n" % (iselected, jet.pt/1e3, jet.tau21, jet.tau21_prime)
        #Method 1
        if selectionType ==1:
            if len(selected)>1:
                if (selected[0].tau21_prime< selected[1].tau21_prime):
                    S1_jets.append(selected[0])
                    print "method1:0th jet appended"
                else:
                #    S1_bg.append(selected[0])
                    S1_jets.append(selected[1])
                    print "length of selected:%r" %(len(selected))
                    for counter in range (2,len(selected)):
                        S1_bg.append(selected[counter])
                        print "counter counter :%r" %(counter)
                        print "method1:1st jet selected as signal" 
            if len(selected)==1:
                S1_jets.append(selected[0])

            

        #Method 2
        if (selectionType ==2) & (quitLoopingM2==False):
            if selected[iselected].tau21_prime<0.5:
                S2_jets.append(jet)
                quitLoopingM2=True;

            for k in range(0,len(selected)):
                if k!=iselected:
                    S2_bg.append(jet)
            print "method2: (%d)th  jet seleced as signal" %(iselected)
            break
         
    if selectionType==1:
        return S1_jets
    if selectionType==2:    
        return S2_jets

        #else: S2_bg.append(jet)

def makeHist(jets,title,var,setLogy):
    c = TCanvas( 'c' , title, 200,10,700,500)
    hist=[]

    i=0
    while i<len(jets):
        histName = "H"+var+ "_signalType"+sigType[i]
        titleName = title+ sigType[i]
        hist[i] = r.TH1F(histName,titleName,50,-1,3)
        hist[i].SetFillColor(histColor[i])

        for jet in jets[i]:
           # print "The selected jet of %r:: (PT, tau21, tau21P) = (%d, %2f %.2f) \n" % (title, jet[i].pt/1e3, jet[i].tau21, jet[i].tau21_prime)
            hist[i].Fill(jet[i].tau21_prime)
        i+=1

    if i==0:
        hist[i].Draw()
    else:
        hist[i].Draw("same")
    c.Update()
    pdfName= title+"pdf"
    if setLogy:
        c.SetLogy
        c.SaveAs(pdfName)

if __name__ == "__main__":
    from argparse import ArgumentParser
    # read the command line arguments
    parser = ArgumentParser(description="Loop through some ntuple and read stuff")
    parser.add_argument("--tree-name", default="zpj", help="The name of the tree to read.")
    parser.add_argument("input_file", help="the input ntuple")
    args = parser.parse_args()

    sigType=["Z+js","W+js","Bkgnd"]
    histColor=[7,8,6]

    selectedFJetsS1={}
    selectedFJetsS2={}
    tau21PHistS1={}
    tau21PHistS2={}

    print "the type of selectedFJetsS1 %r" %type(selectedFJetsS1)

    for i in sigType:
        rawTree = treeCreation(i)
        selectedFJetsS1[i]=selection(rawTree,1)
        selectedFJetsS2[i]=selection(rawTree,2)

    tau21PHistS1= makeHist(selectedFJetsS1,"Tau21P Distribution:Selection1","Tau21P", True)
    tau21PHistS2=makeHist(selectedFJetsS2,"Tau21P Distribution:Selection 2","Tau21P",True)
