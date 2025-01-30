import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

rfile = rt.TFile("out_photon_studies_temp.root","read")

hEphotons = rfile.Get("hEphotons")
hEvCos = rfile.Get("hEvCos")
hEtot  = rfile.Get("hEtot")
hPtot  = rfile.Get("hPtot")
hPara  = rfile.Get("hPara")
hPara2 = rfile.Get("hPara2")

cEphotons = rt.TCanvas("cEphotons","",800,600)
hEphotons.Draw("colz")
cEphotons.Update()


cEvCos = rt.TCanvas("cEvCos","",800,600)
hEvCos.Draw("colz")
cEvCos.Update()

cEtot = rt.TCanvas("cEtot","",800,600)
hEtot.Draw()
cEtot.Update()

cPtot = rt.TCanvas("cPtot","",800,600)
hPtot.Draw()
print("Ptot integral: ",hPtot.Integral())
cPtot.Update()

cPara = rt.TCanvas("cPara","",800,600)
hPara.Draw("colz")
cPara.Update()

cPara2 = rt.TCanvas("cPara2","",800,600)
hPara2.Draw("colz")
cPara2.Update()


input()
