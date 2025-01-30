import os,sys
import ROOT as rt

fbase = "out_temp_run3_CV_500k.root"
fcompare = "out_temp_edepfix.root"

rbase = rt.TFile(fbase)
rcompare = rt.TFile(fcompare)

varlist = ["hEnu__bnbnu"]

keep_v = []
for var in varlist:
    c = rt.TCanvas("c"+var,"var",800,600)
    hbase = rbase.Get(var)
    hcompare = rcompare.Get(var)

    hcompare.SetLineColor(rt.kRed)

    hbase.Draw("histE1")
    hcompare.Draw("histE1same")
    print("[",var,"] ---------------")
    print("  base: ",hbase.Integral())
    print("  compare: ",hcompare.Integral())
    c.Update()

    keep_v += [c,hbase,hcompare]

print("[enter] to end")
input()


