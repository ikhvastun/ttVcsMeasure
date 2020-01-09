from ROOT import TFile
from ROOT import gDirectory
from ROOT import TH1D
from ROOT import TH2F

fm16 = TFile("looseToTight_2016_m_3l.root", "READ")
fm17 = TFile("looseToTight_2017_m_3l.root", "READ")
fe16 = TFile("looseToTight_2016_e_3l.root", "READ")
fe17 = TFile("looseToTight_2017_e_3l.root", "READ")

histm16 = fm16.Get("EGamma_SF2D")
histm17 = fm17.Get("EGamma_SF2D")
histe16 = fe16.Get("EGamma_SF2D")
histe17 = fe17.Get("EGamma_SF2D")

histm16_dataEff = fm16.Get("EGamma_EffData2D")
histm17_dataEff = fm17.Get("EGamma_EffData2D")
histe16_dataEff = fe16.Get("EGamma_EffData2D")
histe17_dataEff = fe17.Get("EGamma_EffData2D")

histm16_mcEff = fm16.Get("EGamma_EffMC2D")
histm17_mcEff = fm17.Get("EGamma_EffMC2D")
histe16_mcEff = fe16.Get("EGamma_EffMC2D")
histe17_mcEff = fe17.Get("EGamma_EffMC2D")




histNewm16 = TH2F("EGamma_SF2D", "EGamma_SF2D", histm16.GetNbinsY(), 0.0, 120.0, histm16.GetNbinsX(), 0.0, 2.5)
histNewm17 = TH2F("EGamma_SF2D", "EGamma_SF2D", histm17.GetNbinsY(), 0.0, 120.0, histm17.GetNbinsX(), 0.0, 2.5)
histNewe16 = TH2F("EGamma_SF2D", "EGamma_SF2D", histe16.GetNbinsY(), 0.0, 120.0, histe16.GetNbinsX(), 0.0, 2.5)
histNewe17 = TH2F("EGamma_SF2D", "EGamma_SF2D", histe17.GetNbinsY(), 0.0, 120.0, histe17.GetNbinsX(), 0.0, 2.5)

histNewm16_unc = TH2F("EGamma_SF2D", "EGamma_SF2D", histm16.GetNbinsY(), 0.0, 120.0, histm16.GetNbinsX(), 0.0, 2.5)
histNewm17_unc = TH2F("EGamma_SF2D", "EGamma_SF2D", histm17.GetNbinsY(), 0.0, 120.0, histm17.GetNbinsX(), 0.0, 2.5)
histNewe16_unc = TH2F("EGamma_SF2D", "EGamma_SF2D", histe16.GetNbinsY(), 0.0, 120.0, histe16.GetNbinsX(), 0.0, 2.5)
histNewe17_unc = TH2F("EGamma_SF2D", "EGamma_SF2D", histe17.GetNbinsY(), 0.0, 120.0, histe17.GetNbinsX(), 0.0, 2.5)

histNewm16.Sumw2()
histNewm17.Sumw2()
histNewe16.Sumw2()
histNewe17.Sumw2()

foutput = TFile("test.root", "recreate")

for b in range(1,histNewm16.GetNbinsX()+1):
  for c in range(1,histNewm16.GetNbinsY()+1):
    histNewm16.SetBinContent(c,b, histm16.GetBinError( b, c) )
    histNewm17.SetBinContent(c,b, histm17.GetBinError( b, c) )
    histNewe16.SetBinContent(c,b, histe16.GetBinError( b, c) )
    histNewe17.SetBinContent(c,b, histe17.GetBinError( b, c) )

#canvas = TCanvas("canvas","canvas", 1200, 1200)
#canvas.Divide(2,2)
#canvas.cd(1)
#histoB.Draw('COLZtexte')
#canvas.cd(2)
#histoC.Draw('COLZtexte')
#canvas.cd(3)
#histoL.Draw('COLZtexte')
#
#histNewm16.Draw("colz")

histNewm16.Write()
histNewm17.Write()
histNewe16.Write()
histNewe17.Write()
f.Close()
foutput.Close()
