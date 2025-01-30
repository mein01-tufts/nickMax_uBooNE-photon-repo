import os,sys
import numpy as np
import ROOT as rt

ntuple_folder = "~/working/data/ntuples/"
rfile = rt.TFile(ntuple_folder+"/ntuple_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV_5000files_20241202.root","read")

tree = rfile.Get("EventTree")
nentries = tree.GetEntries()
pottree = rfile.Get("potTree")
pot = 0.0
for i in range(pottree.GetEntries()):
    pottree.GetEntry(i)
    pot += pottree.totGoodPOT
print("nentries=",nentries," pot=",pot)


# what are the questions we want to answer:
#  1. is there a correlation between the true energy and the energy deposited?
#  2. is there a tagged-gamma/vertex game we can play to verify our description of photon interactions and acceptances?
#      We can use 1 photon with nuclear activity (>=1 p) plus the pion invariant mass to define the expected path and energy
#        of the second photon.  We will have some rate finding the 2nd photon and some rate of missing it.
#        These are the 2gamma+Np and 1gamma+Np events. Assuming true NCpi0 production, this ratio depends on what physics?
#         - photon interaction physics/propagation, do we have the interaction length correct? Do we get the distribution trunk energy correct?
#         - detector physics: for a given energy deposition, do we get the response we expect?
#  3. we can parameterize the data using three variables:
#     Given the photon and vertex, we can define the expected cos theta, or zentih angle relative to the vertex-to-first-photon axis.
#     However, we do not know the azimuthal angle. We can calculate the expected energy and expected angle (on the same plane)
#     and caculate the deviations of the second photon to that. This is in some ways, equivalent to calculating the invariant mass.
#     But we also want to put the 'not found' events in the same plot in order to make a kind of ratio.
#     Maybe we define a threshold of 'expected photon found' and place the event in that bin. We have migratory issues here.
#     Maybe we fill the histogram with a prior distribution?
# 4. For demonstration of technique, I can change the 'interaction length' and weight it down.
#    I can also change the energy scale, either for total vis, or first edep, and quantify the size of effect I can tease apart. 

# we want to make the following plots:
# 1. for every true photon with an above threshold EDep, we calculate the "fromwall" distance
# 2. for photons originating inside the TPC, we can

def get_primary_photons(tree):

    nuvtx = np.array( (tree.trueVtxX, tree.trueVtxY, tree.trueVtxZ) )
    n_vis_photons = 0
    prim_photons_v = []
    

    for i in range(tree.nTrueSimParts):
        if tree.trueSimPartPDG[i]!=22:
            continue

        partvtx = np.array( (tree.trueSimPartX[i], tree.trueSimPartY[i], tree.trueSimPartZ[i]) )
        dvtx = nuvtx-partvtx
        dist = np.sqrt( (dvtx*dvtx).sum() )

        # check its a primary vertex
        if dist>0.5:
            continue

        # passes vertex check: so primary photon

        # get the energy deposit
        pixsum = np.zeros(3) 
        pixsum[0] = tree.trueSimPartPixelSumUplane[i]*0.016
        pixsum[1] = tree.trueSimPartPixelSumVplane[i]*0.016
        pixsum[2] = tree.trueSimPartPixelSumYplane[i]*0.016
        #print(pixsum)
        pixsum_max = np.max(pixsum)
        if pixsum_max>10.0:
            n_vis_photons += 1

        # get 4-momentum
        pmom4 = np.array( (tree.trueSimPartE[i], tree.trueSimPartPx[i], tree.trueSimPartPy[i], tree.trueSimPartPz[i]) )

        photon_info = {"pmom4":pmom4,
                       "tid":tree.trueSimPartTID[i],
                       "mid":tree.trueSimPartMID[i],
                       "start":nuvtx,
                       "pixsumMeV":pixsum,
                       "pixsumMax":pixsum_max}

        prim_photons_v.append( photon_info )

    return n_vis_photons,prim_photons_v
    
out = rt.TFile("out_photon_studies_temp.root","recreate")

hW = rt.TH1D("hW",";invariant mass (MeV^{2})",100,0,1000)
hEphotons = rt.TH2D("hEphotons",";highest photon energy (MeV);lowest photon energy (MeV)",50,0,500,50,0,500)
hEvCos = rt.TH2D("hEvCos",";photon energy (MeV);cos",100,0,1000,50,-1.0,1.0)
hPara  = rt.TH2D("hPara",";photon energy (MeV);cos",50,0,10,50,0.0,2.0)
hPara2 = rt.TH2D("hPara2",";photon energy (MeV);cos",100,0,1000,50,0.0,2.0)
hEtot = rt.TH1D("hEtot","",100,0,1000)
hPtot = rt.TH1D("hPtot","",100,0,1000)

for ientry in range(nentries):
    tree.GetEntry(ientry)

    if ientry>0 and ientry%10000==0:
        print("entry[",ientry,"] of ",nentries)

    # select true primary photons
    n_vis_photons,prim_photon_v = get_primary_photons(tree)

    if n_vis_photons==0:
        continue

    # get true vertex
    nuvtx = prim_photon_v[0]["start"]

    # pair up photons
    #W = (p+q)^2 = pp+qq+2p.q = 0 + 0 + 2*(Ep*Eq - EpEq*cos(theta_pq)) = EpEq*(1-cos(theta_pq))

    n = len(prim_photon_v)
    for i in range(n):
        iphoton = prim_photon_v[i]
        for j in range(i+1,n):
            jphoton = prim_photon_v[j]
            if iphoton["mid"]==jphoton["mid"]:
                # same mother, calc invariant mass
                Ei = iphoton["pmom4"][0]
                Ej = jphoton["pmom4"][0]
                EpEq = iphoton["pmom4"][0]*jphoton["pmom4"][0]
                cos_pq3 = np.sum(iphoton["pmom4"][1:]*jphoton["pmom4"][1:])/EpEq

                ptot = iphoton["pmom4"][1:]+jphoton["pmom4"][1:]
                ptot = np.sqrt( np.sum(ptot*ptot) )

                W2 = 2*EpEq*(1-cos_pq3)
                if W2>0:
                    W = np.sqrt(W2)
                else:
                    W = 0.0
                #print("---------------------------------")
                #print("pmom_i: ",iphoton["pmom4"])
                #print("pmom_j: ",jphoton["pmom4"])
                #print(W,EpEq,cos_pq3,tree.xsecWeight)
                hW.Fill(W,tree.xsecWeight)

                if iphoton["pmom4"][0]>jphoton["pmom4"][0]:
                    hEphotons.Fill( iphoton["pmom4"][0], jphoton["pmom4"][0], tree.xsecWeight )
                else:
                    hEphotons.Fill( jphoton["pmom4"][0], iphoton["pmom4"][0], tree.xsecWeight )

                hEvCos.Fill( iphoton["pmom4"][0], cos_pq3, tree.xsecWeight )
                hEvCos.Fill( jphoton["pmom4"][0], cos_pq3, tree.xsecWeight )

                hEtot.Fill( iphoton["pmom4"][0]+jphoton["pmom4"][0], tree.xsecWeight )
                hPtot.Fill( ptot, tree.xsecWeight )

                #print(2*Ei*Ej/135.0,"  ",1-cos_pq3)
                hPara.Fill( 2*Ei/135.0, 1-cos_pq3, 0.5*tree.xsecWeight )
                hPara.Fill( 2*Ej/135.0, 1-cos_pq3, 0.5*tree.xsecWeight )
                hPara2.Fill( 2*Ei*Ej/135.0, 1-cos_pq3, tree.xsecWeight )       


out.Write()
out.Close()
