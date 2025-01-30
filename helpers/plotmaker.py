import os,sys
import ROOT as rt

class PlotMaker:

    def __init__(self,targetpot=6.6e20):
        # needs path to certain ntuples
        ntuple_folder = "~/working/data/ntuples/"
        #self.mcntuple_paths = {"bnbnu":ntuple_folder+"/dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV.root"}
        self.mcntuple_paths = {"bnbnu":ntuple_folder+"/ntuple_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV_5000files_20241202.root"}
        # need special extbnb ntuple
        self.extbnb = None
        self.hist_v = {}
        self.targetpot = targetpot

    def load_ntuples(self):
        self.tree = {}
        self.is_mc = {}
        self.pot = {}
        self.rfile = {}
        self.nentries = {}
        print("Load ntuples: ---------------------------")
        for name in self.mcntuple_paths:
            self.rfile[name] = rt.TFile( self.mcntuple_paths[name], "open" )
            self.tree[name] = self.rfile[name].Get("EventTree")
            self.is_mc[name] = True
            self.nentries[name] = self.tree[name].GetEntries()
            pottree = self.rfile[name].Get("potTree")
            pot = 0.0
            for i in range(pottree.GetEntries()):
                pottree.GetEntry(i)
                pot += pottree.totGoodPOT
            print("[",name,"] nentries=",self.nentries[name]," pot=",pot)
            self.pot[name] = pot

    def addHist( self, name, hist_template, fill_var_fn, selection_fn, get_weight_fn ):
        self.hist_v[name] = {"hist":hist_template,
                             "var_fn":fill_var_fn,
                             "select_fn":selection_fn,
                             "weight_fn":get_weight_fn}
        

    def runloop(self):
        """
        
        to get variables, 

        """
        
        for treename in self.tree:

            hist_v = {}
            for varname in self.hist_v:
                htemplate = self.hist_v[varname]["hist"]
                hname = htemplate.GetName() +"__"+treename
                hist = htemplate.Clone(hname)
                hist.Reset()
                hist_v[varname] = hist

            tree = self.tree[treename]
            for ientry in range(self.nentries[treename]):
                if ientry>0 and ientry%1000==0:
                    print("running entry[",ientry,"] of ",self.nentries[treename])
                tree.GetEntry(ientry)

                for varname in self.hist_v:
                    var_fn = self.hist_v[varname]["var_fn"]
                    select_fn = self.hist_v[varname]["select_fn"]
                    cutresult = select_fn( tree, self.is_mc[treename] )
                    if not cutresult:
                        continue
                    weight_fn = self.hist_v[varname]["weight_fn"]
                    out = var_fn( tree, self.is_mc[treename] )
                    w = weight_fn( tree, self.is_mc[treename] )
                    if out is not tuple:
                        hist_v[varname].Fill( out, w )
                    elif len(out)==1:
                        hist_v[varname].Fill( out[0], w )
                    else:
                        raise ValueError("cannot support output with ",len(out)," outputs")
            # scale
            for hname in hist_v:
                if self.is_mc[treename]:
                    scale = self.targetpot/self.pot[treename]
                    hist_v[hname].Scale(scale)
                hist_v[hname].Write()

                
        
