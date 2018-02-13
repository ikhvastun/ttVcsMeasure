from ROOT import *
from array import array
#import ctypes

gStyle.SetPaintTextFormat("4.2f");
gStyle.SetOptStat(0)

c = TChain('FakeElectrons/fakeTree')
c.Add("/Users/illiakhvastunov/Desktop/CERN/MCsamples/80X/ttW_ReReco_NewMC_NewJEC/singleTop.root")
c.Add("/Users/illiakhvastunov/Desktop/CERN/MCsamples/80X/ttW_ReReco_NewMC_NewJEC/singleTopBar.root")

ptBorders = [30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 670.]
etaBins = [0, 0.8,1.6, 2.4]

histoB = TProfile2D("histoB","histoB", len(ptBorders)-1, array('d', ptBorders), len(etaBins)-1, array('d', etaBins))
histoC = TProfile2D("histoC","histoC", len(ptBorders)-1, array('d', ptBorders), len(etaBins)-1, array('d', etaBins))
histoL = TProfile2D("histoL","histoL", len(ptBorders)-1, array('d', ptBorders), len(etaBins)-1, array('d', etaBins))

c.Draw("_csv>0.5426:abs(_jetEta):_jetPt >> histoB", "Iteration$<_n_Jets&&_jetFlavour == 5")
c.Draw("_csv>0.5426:abs(_jetEta):_jetPt >> histoC", "Iteration$<_n_Jets&&_jetFlavour == 4")
c.Draw("_csv>0.5426:abs(_jetEta):_jetPt >> histoL", "Iteration$<_n_Jets&&_jetFlavour == 0")

#canvas = TCanvas("canvas","canvas", 1200, 400)
#canvas.Divide(3,1)
#canvas.cd(1)
#histoB.Draw('COLZtexte')
#canvas.cd(2)
#histoC.Draw('COLZtexte')
#canvas.cd(3)
#histoL.Draw('COLZtexte')

fileOutput = TFile("btagEff_Loose.root", "Recreate")
histoB.Write();
histoC.Write();
histoL.Write();
fileOutput.Close()

#print histo.GetBinContent(histo.FindBin(200., 2.))
#print histo.GetBinError(histo.FindBin(200., 2.))
