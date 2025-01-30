import os,sys
import ROOT as rt
rt.gStyle.SetPadBottomMargin(0.15)
rt.gStyle.SetOptStat(0)

r = rt.TFile("out_gamma_truth_temp.root","read")

hX = r.Get("hX_truephoton__bnbnu")
hY = r.Get("hY_truephoton__bnbnu")
hZ = r.Get("hZ_truephoton__bnbnu")

# X-position plot
cX = rt.TCanvas("cX","",800,600)
hX.GetXaxis().SetTitleSize(0.05)
hX.SetLineWidth(2)
hX.Draw("histE1")
tX1 = rt.TLine(10.0,0.0,10.0,2500.0)
tX2 = rt.TLine(246.0,0.0,246.0,2500.0)

tX1.SetLineColor(rt.kRed)
tX1.SetLineWidth(2)
tX1.SetLineStyle(2)
tX1.Draw()

tX2.SetLineColor(rt.kRed)
tX2.SetLineWidth(2)
tX2.SetLineStyle(2)
tX2.Draw()

tXt = rt.TText(245,2930,"Entries: %.1f"%(hX.Integral()))
tXt.Draw()

cX.Update()
cX.SaveAs("cVertX_ingv_anyprimaryphoton.png")

# Y-position plot
cY = rt.TCanvas("cY","",800,600)
hY.GetXaxis().SetTitleSize(0.05)
hY.SetLineWidth(2)
hY.Draw("histE1")
tY1 = rt.TLine(-106.5,0.0,-106.5,2700.0)
tY2 = rt.TLine(106.5,0.0,106.5,2700.0)

tY1.SetLineColor(rt.kRed)
tY1.SetLineWidth(2)
tY1.SetLineStyle(2)
tY1.Draw()

tY2.SetLineColor(rt.kRed)
tY2.SetLineWidth(2)
tY2.SetLineStyle(2)
tY2.Draw()

tYt = rt.TText(80,2980,"Entries: %.1f"%(hY.Integral()))
tYt.Draw()

cY.Update()
cY.SaveAs("cVertY_ingv_anyprimaryphoton.png")

# Z-position plot
cZ = rt.TCanvas("cZ","",800,600)
#hZ.Rebin(3)
hZ.GetXaxis().SetTitleSize(0.05)
hZ.SetLineWidth(2)
hZ.Draw("histE1")
tZ1 = rt.TLine(10.0,0.0,10.0,820.0)
tZ2 = rt.TLine(1026.0,0.0,1026.0,820.0)

tZ1.SetLineColor(rt.kRed)
tZ1.SetLineWidth(2)
tZ1.SetLineStyle(2)
tZ1.Draw()

tZ2.SetLineColor(rt.kRed)
tZ2.SetLineWidth(2)
tZ2.SetLineStyle(2)
tZ2.Draw()

tZt = rt.TText(740,870,"Entries: %.1f"%(hZ.Integral()))
tZt.Draw()

cZ.Update()
cZ.SaveAs("cVertZ_ingv_anyprimaryphoton.png")

input()
