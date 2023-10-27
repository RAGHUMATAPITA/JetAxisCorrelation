#include "call_libraries.h"

TH1D* hEvents = new TH1D("hEvents", "", 10, 0, 10);
TH1D* hpthat = new TH1D("hpthat", "", 200, 0, 1200.);
TH1D* hpthatW = new TH1D("hpthatW", "", 100, 0, 1.);
TH1D* hCent = new TH1D("hCent", "", 200, 0, 200);
TH1D* hZvtx = new TH1D("hZvtx", "", 200, -20, 20);

// for test
TH1D* hptreco = new TH1D("hptreco", "", 200., 0., 1000.);
TH1D* hptreco_Matched = new TH1D("hptreco_Matched", "", 200., 0., 1000.);
TH1D* hptreco_Matched_PM = new TH1D("hptreco_Matched_PM", "", 200., 0., 1000.);

TH1D* hptgen = new TH1D("hptgen", "", 200., 0., 1000.);
TH1D* hptgen_Matched = new TH1D("hptgen_Matched", "", 200., 0., 1000.);
TH1D* hptgen_Matched_PM = new TH1D("hptgen_Matched_PM", "", 200., 0., 1000.);

const int NCentbin = 5;
double CentbinEdge[NCentbin+1] = {0., 20., 60., 100., 160., 200};
TH1D* hCentBin = new TH1D("hCentBin", "", NCentbin, CentbinEdge);

int    bins4D_jet[4]   =   { 200  , 50   , 64           , NCentbin        };
double xmin4D_jet[4]   =   { 0.   , -2.5 , -TMath::Pi() , 0.              };
double xmax4D_jet[4]   =   { 1000.,  2.5 , TMath::Pi()  , (double)NCentbin};

// Jet histo

/*
TH2D* hqqbar_Scan_Gen_noW = new TH2D("hqqbar_Scan_Gen_noW", "", 15, 0., 15., 15, 0., 15.);
TH2D* hqqbar_Scan_Gen_W = new TH2D("hqqbar_Scan_Gen_W", "", 15, 0., 15., 15, 0., 15.);

TH2D* hqqbar_Scan_Reco_noW = new TH2D("hqqbar_Scan_Reco_noW", "", 15, 0., 15., 15, 0., 15.);
TH2D* hqqbar_Scan_Reco_W = new TH2D("hqqbar_Scan_Reco_W", "", 15, 0., 15., 15, 0., 15.);
*/

int bins_refparton[3]    = {16,  16,   11 };
double xmin_refparton[3] = {0.,  0.,   0. };
double xmax_refparton[3] = {16., 16.,  11.};

THnSparseD* hqqbar_Scan_Gen_noW = new THnSparseD("hqqbar_Scan_Gen_noW", "", 3, bins_refparton, xmin_refparton, xmax_refparton);
THnSparseD* hqqbar_Scan_Gen_W = new THnSparseD("hqqbar_Scan_Gen_W", "", 3, bins_refparton, xmin_refparton, xmax_refparton);

THnSparseD* hqqbar_Scan_Reco_noW = new THnSparseD("hqqbar_Scan_Reco_noW", "", 3, bins_refparton, xmin_refparton, xmax_refparton);
THnSparseD* hqqbar_Scan_Reco_W = new THnSparseD("hqqbar_Scan_Reco_W", "", 3, bins_refparton, xmin_refparton, xmax_refparton);

// MC reco or data
 //w/o pthat weight
THnSparseD * hJet_RawpT_Eta_Phi_ctbin_nopTCut_noW = new THnSparseD("hJet_RawpT_Eta_Phi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_noW = new THnSparseD("hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_Eta_Phi_ctbin_nopTCut_noW = new THnSparseD("hJet_CorrpT_Eta_Phi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_noW = new THnSparseD("hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

THnSparseD * hJet_RawpT_Eta_Phi_ctbin_pTCut_noW = new THnSparseD("hJet_RawpT_Eta_Phi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW = new THnSparseD("hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW = new THnSparseD("hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW = new THnSparseD("hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

 //w/o pthat weight
THnSparseD * hJet_RawpT_Eta_Phi_ctbin_nopTCut_W = new THnSparseD("hJet_RawpT_Eta_Phi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_W = new THnSparseD("hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_Eta_Phi_ctbin_nopTCut_W = new THnSparseD("hJet_CorrpT_Eta_Phi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_W = new THnSparseD("hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

THnSparseD * hJet_RawpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hJet_RawpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W = new THnSparseD("hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hJet_CorrpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W = new THnSparseD("hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

int    bins_pT_refparton[3]   =   { 200  , 8 , NCentbin        };
double xmin_pT_refparton[3]   =   { 0.   , 0., 0               };
double xmax_pT_refparton[3]   =   { 1000., 8., (double)NCentbin};

// Matched, unmatched
THnSparseD * hJet_CorrpT_refparton_ctbin_W = new THnSparseD("hJet_CorrpT_refparton_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_CorrpT_refpartonB_ctbin_W = new THnSparseD("hJet_CorrpT_refpartonB_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_CorrpT_refpartonB_ctbin_noW = new THnSparseD("hJet_CorrpT_refpartonB_ctbin_noW", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_MCorrpT_refparton_ctbin_W = new THnSparseD("hJet_MCorrpT_refparton_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_MCorrpT_refpartonB_ctbin_W = new THnSparseD("hJet_MCorrpT_refpartonB_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_UCorrpT_refparton_ctbin_W = new THnSparseD("hJet_UCorrpT_refparton_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

THnSparseD * hJet_UCorrpT_refpartonB_ctbin_W = new THnSparseD("hJet_UCorrpT_refpartonB_ctbin_W", "", 3, bins_pT_refparton, xmin_pT_refparton, xmax_pT_refparton);

// MC gen
// w/o pthat weight
THnSparseD * hJet_GenpT_Eta_Phi_ctbin_nopTCut_noW = new THnSparseD("hJet_GenpT_Eta_Phi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_noW = new THnSparseD("hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

THnSparseD * hJet_GenpT_Eta_Phi_ctbin_pTCut_noW = new THnSparseD("hJet_GenpT_Eta_Phi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW = new THnSparseD("hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// w/ pthat weight
THnSparseD * hJet_GenpT_Eta_Phi_ctbin_nopTCut_W = new THnSparseD("hJet_GenpT_Eta_Phi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_W = new THnSparseD("hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

THnSparseD * hJet_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hJet_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W = new THnSparseD("hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// leading and subleading jet                                               
// w/o weight                                                               

THnSparseD * hld_RawpT_Eta_Phi_ctbin_noW = new THnSparseD("hld_RawpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_RawpT_Eta_Phi_ctbin_noW = new THnSparseD("hsld_RawpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hld_RawpT_WTAEta_WTAPhi_ctbin_noW = new THnSparseD("hld_RawpT_WTAEta_WTAPhi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_RawpT_WTAEta_WTAPhi_ctbin_noW = new THnSparseD("hsld_RawpT_WTAEta_WTAPhi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
         
THnSparseD * hld_CorrpT_Eta_Phi_ctbin_noW = new THnSparseD("hld_CorrpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_CorrpT_Eta_Phi_ctbin_noW = new THnSparseD("hsld_CorrpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hld_CorrpT_WTAEta_WTAPhi_ctbin_noW = new THnSparseD("hld_CorrpT_WTAEta_WTAPhi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_CorrpT_WTAEta_WTAPhi_ctbin_noW = new THnSparseD("hsld_CorrpT_WTAEta_WTAPhi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// w/ weight
THnSparseD * hld_RawpT_Eta_Phi_ctbin_W = new THnSparseD("hld_RawpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_RawpT_Eta_Phi_ctbin_W = new THnSparseD("hsld_RawpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hld_RawpT_WTAEta_WTAPhi_ctbin_W = new THnSparseD("hld_RawpT_WTAEta_WTAPhi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_RawpT_WTAEta_WTAPhi_ctbin_W = new THnSparseD("hsld_RawpT_WTAEta_WTAPhi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

THnSparseD * hld_CorrpT_Eta_Phi_ctbin_W = new THnSparseD("hld_CorrpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_CorrpT_Eta_Phi_ctbin_W = new THnSparseD("hsld_CorrpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hld_CorrpT_WTAEta_WTAPhi_ctbin_W = new THnSparseD("hld_CorrpT_WTAEta_WTAPhi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_CorrpT_WTAEta_WTAPhi_ctbin_W = new THnSparseD("hsld_CorrpT_WTAEta_WTAPhi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);


//MC gen
// w/o weight
THnSparseD * hld_genpT_Eta_Phi_ctbin_noW = new THnSparseD("hld_genpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_genpT_Eta_Phi_ctbin_noW = new THnSparseD("hsld_genpT_Eta_Phi_ctbin_noW", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// w/ weight
THnSparseD * hld_genpT_Eta_Phi_ctbin_W = new THnSparseD("hld_genpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hsld_genpT_Eta_Phi_ctbin_W = new THnSparseD("hsld_genpT_Eta_Phi_ctbin_W", "", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int    refpTnopTCutbins2D_jet[2]   =   { 400   ,  NCentbin        };
double refpTnopTCutxmin2D_jet[2]   =   { -1000.,  0.              };
double refpTnopTCutxmax2D_jet[2]   =   { 1000. ,  (double)NCentbin};

int    refpTpTCut0bins2D_jet[2]   =   { 400   ,  NCentbin        };
double refpTpTCut0xmin2D_jet[2]   =   { -1000.,  0.              };
double refpTpTCut0xmax2D_jet[2]   =   { 1000. ,  (double)NCentbin};

int    bins2D_jet[2]   =   { 200   ,  NCentbin        };
double xmin2D_jet[2]   =   { 0.    ,  0.              };
double xmax2D_jet[2]   =   { 1000. ,  (double)NCentbin};

//ref_pt histo
//w/o pthat weight
THnSparseD * hJet_refpT_ctbin_pTCut_noW = new THnSparseD("hJet_refpT_ctbin_pTCut_noW", "", 2, bins2D_jet, xmin2D_jet, xmax2D_jet);
THnSparseD * hJet_refpT_ctbin_nopTCut_noW = new THnSparseD("hJet_refpT_ctbin_nopTCut_noW", "", 2, refpTnopTCutbins2D_jet, refpTnopTCutxmin2D_jet, refpTnopTCutxmax2D_jet);
THnSparseD * hJet_refpT_ctbin_pTCut0_noW = new THnSparseD("hJet_refpT_ctbin_pTCut0_noW", "", 2, refpTpTCut0bins2D_jet, refpTpTCut0xmin2D_jet, refpTpTCut0xmax2D_jet);

//w/ pthat weight
THnSparseD * hJet_refpT_ctbin_pTCut_W = new THnSparseD("hJet_refpT_ctbin_pTCut_W", "", 2, bins2D_jet, xmin2D_jet, xmax2D_jet);
THnSparseD * hJet_refpT_ctbin_nopTCut_W = new THnSparseD("hJet_refpT_ctbin_nopTCut_W", "", 2, refpTnopTCutbins2D_jet, refpTnopTCutxmin2D_jet, refpTnopTCutxmax2D_jet);
THnSparseD * hJet_refpT_ctbin_pTCut0_W = new THnSparseD("hJet_refpT_ctbin_pTCut0_W", "", 2, refpTpTCut0bins2D_jet, refpTpTCut0xmin2D_jet, refpTpTCut0xmax2D_jet);

// JES
int    bins4D_jes[4]   =   { 500 ,  50   , 8 , NCentbin        };
double xmin4D_jes[4]   =   { 0.  ,  0.   , 0., 0               };
double xmax4D_jes[4]   =   { 5.  ,  1000., 8., (double)NCentbin};

//w/o pthat weight
THnSparseD * hJes_rawpT_refpT_pTCut0_refparton_ctbin_noW = new THnSparseD("hJes_rawpT_refpT_pTCut0_refparton_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_rawpT_refpT_pTCut_refparton_ctbin_noW = new THnSparseD("hJes_rawpT_refpT_pTCut_refparton_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut0_refparton_ctbin_noW = new THnSparseD("hJes_CorrpT_refpT_pTCut0_refparton_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut_refparton_ctbin_noW = new THnSparseD("hJes_CorrpT_refpT_pTCut_refparton_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);

THnSparseD * hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_noW = new THnSparseD("hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW = new THnSparseD("hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW = new THnSparseD("hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW = new THnSparseD("hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);

//w/ pthat weight
THnSparseD * hJes_rawpT_refpT_pTCut0_refparton_ctbin_W = new THnSparseD("hJes_rawpT_refpT_pTCut0_refparton_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_rawpT_refpT_pTCut_refparton_ctbin_W = new THnSparseD("hJes_rawpT_refpT_pTCut_refparton_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut0_refparton_ctbin_W = new THnSparseD("hJes_CorrpT_refpT_pTCut0_refparton_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut_refparton_ctbin_W = new THnSparseD("hJes_CorrpT_refpT_pTCut_refparton_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);

THnSparseD * hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_W = new THnSparseD("hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W = new THnSparseD("hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_W = new THnSparseD("hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD * hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W = new THnSparseD("hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W", "", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);

//JER
int    bins4D_jer[4]   =   { 400 ,  50   , 8 , NCentbin        };
double xmin4D_jer[4]   =   { -2. ,  0.   , 0., 0               };
double xmax4D_jer[4]   =   { 2.  ,  1000., 8., (double)NCentbin};
//w/o pthat weight
THnSparseD * hJer_rawpT_refpT_pTCut0_refparton_ctbin_noW = new THnSparseD("hJer_rawpT_refpT_pTCut0_refparton_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_rawpT_refpT_pTCut_refparton_ctbin_noW = new THnSparseD("hJer_rawpT_refpT_pTCut_refparton_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut0_refparton_ctbin_noW = new THnSparseD("hJer_CorrpT_refpT_pTCut0_refparton_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut_refparton_ctbin_noW = new THnSparseD("hJer_CorrpT_refpT_pTCut_refparton_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);


THnSparseD * hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_noW = new THnSparseD("hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW = new THnSparseD("hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW = new THnSparseD("hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW = new THnSparseD("hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);

//w/ pthat weight
THnSparseD * hJer_rawpT_refpT_pTCut0_refparton_ctbin_W = new THnSparseD("hJer_rawpT_refpT_pTCut0_refparton_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_rawpT_refpT_pTCut_refparton_ctbin_W = new THnSparseD("hJer_rawpT_refpT_pTCut_refparton_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut0_refparton_ctbin_W = new THnSparseD("hJer_CorrpT_refpT_pTCut0_refparton_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut_refparton_ctbin_W = new THnSparseD("hJer_CorrpT_refpT_pTCut_refparton_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);


THnSparseD * hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_W = new THnSparseD("hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W = new THnSparseD("hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_W = new THnSparseD("hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD * hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W = new THnSparseD("hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W", "", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);


//correlation histograms
//MC gen
//const int DR_nbins = 120

const int DR_nbins = 18;

//double DR_bins_edge[DR_nbins+1] = {0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.16, 0.20, 0.25, 0.30, 0.4};
double DR_bins_edge[DR_nbins+1] = {0.0, 0.0035, 0.007, 0.011, 0.016, 0.022, 0.029, 0.037, 0.046, 0.056, 0.067, 0.080, 0.095, 0.115, 0.14, 0.17, 0.22, 0.30, 0.4};

const int jtpT_nbins = 11;
double jtpT_bins_edge[jtpT_nbins+1] = {40., 60., 80., 100., 120., 140., 160., 180., 200., 300., 500., 5020.};
TH1D* hJetpTBin = new TH1D("hJetpTBin", "", jtpT_nbins, jtpT_bins_edge);

const int ndim_DR = 4;

int    bins_njets[ndim_DR]   =   { 1        ,jtpT_nbins         ,8 , NCentbin        };
double xmin_njets[ndim_DR]   =   { 1.       ,0.                 ,0., 0.              };
double xmax_njets[ndim_DR]   =   { 2.       ,(double)jtpT_nbins ,8., (double)NCentbin};


int    bins_DR[ndim_DR]   =   { DR_nbins               ,jtpT_nbins         ,8 , NCentbin};
//double xmin_DR[ndim_DR]   =   { 0.        ,0.                 ,0., 0.              };
//double xmax_DR[ndim_DR]   =   { 0.4       ,(double)jtpT_nbins ,8., (double)NCentbin};
double xmin_DR[ndim_DR]   =   { DR_bins_edge[0]        ,0.                 ,0., 0.              };
double xmax_DR[ndim_DR]   =   { DR_bins_edge[DR_nbins] ,(double)jtpT_nbins ,8., (double)NCentbin};

int nEtaBins = 201;
int nPhiBins = 199;

int    bins_DEta[ndim_DR]   =   {nEtaBins    ,jtpT_nbins         ,8 , NCentbin};
double xmin_DEta[ndim_DR]   =   {-0.5        ,0.                 ,0., 0.              };
double xmax_DEta[ndim_DR]   =   { 0.5        ,(double)jtpT_nbins ,8., (double)NCentbin};

int    bins_DPhi[ndim_DR]   =   {nPhiBins       ,jtpT_nbins         ,8 , NCentbin};
double xmin_DPhi[ndim_DR]   =   {-0.5           ,0.                 ,0., 0.              };
double xmax_DPhi[ndim_DR]   =   { 0.5           ,(double)jtpT_nbins ,8., (double)NCentbin};

//binning gor quark and antiquark
int    bins_njets_PM[ndim_DR]   =   { 1        ,jtpT_nbins         ,16 , NCentbin        };
double xmin_njets_PM[ndim_DR]   =   { 1.       ,0.                 ,0.,  0.              };
double xmax_njets_PM[ndim_DR]   =   { 2.       ,(double)jtpT_nbins ,16., (double)NCentbin};

int    bins_DR_PM[ndim_DR]   =   { DR_nbins               ,jtpT_nbins         ,16 , NCentbin};
//double xmin_DR_PM[ndim_DR]   =   { 0.        ,0.                 ,0.,  0.              };
//double xmax_DR_PM[ndim_DR]   =   { 0.4       ,(double)jtpT_nbins ,16., (double)NCentbin};
double xmin_DR_PM[ndim_DR]   =   { DR_bins_edge[0]        ,0.                 ,0.,  0.              };
double xmax_DR_PM[ndim_DR]   =   { DR_bins_edge[DR_nbins] ,(double)jtpT_nbins ,16., (double)NCentbin};

//MC gen
// w/o pthat weight
THnSparseD* hNGenJets_noW = new THnSparseD("hNGenJets_noW", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_Gen_EWTA_noW = new THnSparseD("hDR_Gen_EWTA_noW", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_Gen_EWTA_noW = new THnSparseD("hDEta_Gen_EWTA_noW", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_Gen_EWTA_noW = new THnSparseD("hDPhi_Gen_EWTA_noW", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);
//TH1D* hGenCorrDR_EWTA_noW_[NCentbin][jtpT_nbins];
// w/ pthat weight
THnSparseD* hNGenJets_W = new THnSparseD("hNGenJets_W", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_Gen_EWTA_W = new THnSparseD("hDR_Gen_EWTA_W", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_Gen_EWTA_W = new THnSparseD("hDEta_Gen_EWTA_W", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_Gen_EWTA_W = new THnSparseD("hDPhi_Gen_EWTA_W", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

// MC ref gen matched
// w/o pthat weight
THnSparseD* hNGenJets_Matched_noW = new THnSparseD("hNGenJets_Matched_noW", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_Gen_EWTA_Matched_noW = new THnSparseD("hDR_Gen_EWTA_Matched_noW", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_Gen_EWTA_Matched_noW = new THnSparseD("hDEta_Gen_EWTA_Matched_noW", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_Gen_EWTA_Matched_noW = new THnSparseD("hDPhi_Gen_EWTA_Matched_noW", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

// w/ pthat weight
THnSparseD* hNGenJets_Matched_W = new THnSparseD("hNGenJets_Matched_W", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_Gen_EWTA_Matched_W = new THnSparseD("hDR_Gen_EWTA_Matched_W", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_Gen_EWTA_Matched_W = new THnSparseD("hDEta_Gen_EWTA_Matched_W", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_Gen_EWTA_Matched_W = new THnSparseD("hDPhi_Gen_EWTA_Matched_W", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

THnSparseD* hNGenJets_Matched_PM_W = new THnSparseD("hNGenJets_Matched_PM_W", "", ndim_DR, bins_njets_PM, xmin_njets_PM, xmax_njets_PM);
THnSparseD* hDR_Gen_EWTA_Matched_PM_W = new THnSparseD("hDR_Gen_EWTA_Matched_PM_W", "", ndim_DR, bins_DR_PM, xmin_DR_PM, xmax_DR_PM);

// MC reco or data
// w/o pthat weight 
THnSparseD* hNJets_noW = new THnSparseD("hNJets_noW", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_EWTA_noW = new THnSparseD("hDR_EWTA_noW", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_EWTA_noW = new THnSparseD("hDEta_EWTA_noW", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_EWTA_noW = new THnSparseD("hDPhi_EWTA_noW", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

// w/ pthat weight 
THnSparseD* hNJets_W = new THnSparseD("hNJets_W", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_EWTA_W = new THnSparseD("hDR_EWTA_W", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_EWTA_W = new THnSparseD("hDEta_EWTA_W", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_EWTA_W = new THnSparseD("hDPhi_EWTA_W", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

//reco matched 
//w/o pthat weight
THnSparseD* hNJets_Matched_noW = new THnSparseD("hNJets_Matched_noW", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_EWTA_Matched_noW = new THnSparseD("hDR_EWTA_Matched_noW", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_EWTA_Matched_noW = new THnSparseD("hDEta_EWTA_Matched_noW", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_EWTA_Matched_noW = new THnSparseD("hDPhi_EWTA_Matched_noW", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

// w/ pthat weight
THnSparseD* hNJets_Matched_W = new THnSparseD("hNJets_Matched_W", "", ndim_DR, bins_njets, xmin_njets, xmax_njets);
THnSparseD* hDR_EWTA_Matched_W = new THnSparseD("hDR_EWTA_Matched_W", "", ndim_DR, bins_DR, xmin_DR, xmax_DR);
THnSparseD* hDEta_EWTA_Matched_W = new THnSparseD("hDEta_EWTA_Matched_W", "", ndim_DR, bins_DEta, xmin_DEta, xmax_DEta);
THnSparseD* hDPhi_EWTA_Matched_W = new THnSparseD("hDPhi_EWTA_Matched_W", "", ndim_DR, bins_DPhi, xmin_DPhi, xmax_DPhi);

THnSparseD* hNJets_Matched_PM_W = new THnSparseD("hNJets_Matched_PM_W", "", ndim_DR, bins_njets_PM, xmin_njets_PM, xmax_njets_PM);
THnSparseD* hDR_EWTA_Matched_PM_W = new THnSparseD("hDR_EWTA_Matched_PM_W", "", ndim_DR, bins_DR_PM, xmin_DR_PM, xmax_DR_PM);

void sumw2()
{
  //Event histo
  hpthat->Sumw2();
  hpthatW->Sumw2();
  hEvents->Sumw2();
  hCent->Sumw2();
  hZvtx->Sumw2();

  // for test
  hptreco->Sumw2();
  hptreco_Matched->Sumw2();
  hptreco_Matched_PM->Sumw2();

  hptgen->Sumw2();
  hptgen_Matched->Sumw2();
  hptgen_Matched_PM->Sumw2();

  // Jet histo 
  hqqbar_Scan_Gen_noW->Sumw2();
  hqqbar_Scan_Gen_W->Sumw2();

  hqqbar_Scan_Reco_noW->Sumw2();
  hqqbar_Scan_Reco_W->Sumw2();

  //MC reco
  //w/o pthat weight
  hJet_RawpT_Eta_Phi_ctbin_nopTCut_noW->Sumw2();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Sumw2();
  hJet_CorrpT_Eta_Phi_ctbin_nopTCut_noW->Sumw2();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Sumw2();
  hJet_RawpT_Eta_Phi_ctbin_pTCut_noW->Sumw2();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Sumw2();
  hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW->Sumw2();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Sumw2();

  //w/ pthat weight
  hJet_RawpT_Eta_Phi_ctbin_nopTCut_W->Sumw2();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Sumw2();
  hJet_CorrpT_Eta_Phi_ctbin_nopTCut_W->Sumw2();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Sumw2();
  hJet_RawpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W->Sumw2();
  hJet_CorrpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W->Sumw2();

  // matched, unmatched, total 
  hJet_CorrpT_refparton_ctbin_W->Sumw2();
  hJet_CorrpT_refpartonB_ctbin_W->Sumw2();
  hJet_CorrpT_refpartonB_ctbin_noW->Sumw2();

  hJet_MCorrpT_refparton_ctbin_W->Sumw2();
  hJet_MCorrpT_refpartonB_ctbin_W->Sumw2();

  hJet_UCorrpT_refparton_ctbin_W->Sumw2();
  hJet_UCorrpT_refpartonB_ctbin_W->Sumw2();

  //MC gen
  //w/o pthat weight
  hJet_GenpT_Eta_Phi_ctbin_nopTCut_noW->Sumw2();
  hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Sumw2();
  hJet_GenpT_Eta_Phi_ctbin_pTCut_noW->Sumw2();
  hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Sumw2();

  //w/ pthat weight
  hJet_GenpT_Eta_Phi_ctbin_nopTCut_W->Sumw2();
  hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Sumw2();
  hJet_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W->Sumw2();

  // leading and subleading histo
  //w/o pthat weight     
  hld_RawpT_Eta_Phi_ctbin_noW->Sumw2();
  hsld_RawpT_Eta_Phi_ctbin_noW->Sumw2();
  hld_CorrpT_Eta_Phi_ctbin_noW->Sumw2();
  hsld_CorrpT_Eta_Phi_ctbin_noW->Sumw2();

  // w/ weight
  hld_RawpT_Eta_Phi_ctbin_W->Sumw2();
  hsld_RawpT_Eta_Phi_ctbin_W->Sumw2();
  hld_CorrpT_Eta_Phi_ctbin_W->Sumw2();
  hsld_CorrpT_Eta_Phi_ctbin_W->Sumw2();

  //MC gen
  // w/o weight
  hld_genpT_Eta_Phi_ctbin_noW->Sumw2();
  hsld_genpT_Eta_Phi_ctbin_noW->Sumw2();

  // w/ weight
  hld_genpT_Eta_Phi_ctbin_W->Sumw2();
  hsld_genpT_Eta_Phi_ctbin_W->Sumw2();

  //ref_pt histo
  //w/o pthat weight
  hJet_refpT_ctbin_nopTCut_noW->Sumw2();
  hJet_refpT_ctbin_pTCut_noW->Sumw2();
  hJet_refpT_ctbin_pTCut0_noW->Sumw2();

  //w/ pthat weight
  hJet_refpT_ctbin_nopTCut_W->Sumw2();
  hJet_refpT_ctbin_pTCut_W->Sumw2();
  hJet_refpT_ctbin_pTCut0_W->Sumw2();

  //JES
  //w/o pthat weight   
  hJes_rawpT_refpT_pTCut0_refparton_ctbin_noW->Sumw2();
  hJes_rawpT_refpT_pTCut_refparton_ctbin_noW->Sumw2();
  hJes_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Sumw2();
  hJes_CorrpT_refpT_pTCut_refparton_ctbin_noW->Sumw2();

  hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Sumw2();
  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Sumw2();
  hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Sumw2();
  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Sumw2();
  //w/ pthat weight
  hJes_rawpT_refpT_pTCut0_refparton_ctbin_W->Sumw2();
  hJes_rawpT_refpT_pTCut_refparton_ctbin_W->Sumw2();
  hJes_CorrpT_refpT_pTCut0_refparton_ctbin_W->Sumw2();
  hJes_CorrpT_refpT_pTCut_refparton_ctbin_W->Sumw2();

  hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Sumw2();
  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W->Sumw2();
  hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Sumw2();
  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Sumw2();
  
  //JER
  //w/o pthat weight
  hJer_rawpT_refpT_pTCut0_refparton_ctbin_noW->Sumw2();
  hJer_rawpT_refpT_pTCut_refparton_ctbin_noW->Sumw2();
  hJer_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Sumw2();
  hJer_CorrpT_refpT_pTCut_refparton_ctbin_noW->Sumw2();

  hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Sumw2();
  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Sumw2();
  hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Sumw2();
  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Sumw2();

  //w/ pthat weight
  hJer_rawpT_refpT_pTCut0_refparton_ctbin_W->Sumw2();
  hJer_rawpT_refpT_pTCut_refparton_ctbin_W->Sumw2();
  hJer_CorrpT_refpT_pTCut0_refparton_ctbin_W->Sumw2();
  hJer_CorrpT_refpT_pTCut_refparton_ctbin_W->Sumw2();

  hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Sumw2();
  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W->Sumw2();
  hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Sumw2();
  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Sumw2();

  // correlation histogram
  //MC gen
  // w/o pthat weight

  hDR_Gen_EWTA_noW->GetAxis(0)->Set(DR_nbins, DR_bins_edge); // set variable bin

  hNGenJets_noW->Sumw2();
  hDR_Gen_EWTA_noW->Sumw2();
  hDEta_Gen_EWTA_noW->Sumw2();
  hDPhi_Gen_EWTA_noW->Sumw2();

  /*
  for(int ict = 0; ict < NCentbin; ict++)
    {
      for(int ipt = 0; ipt < jtpT_nbins; ipt++)
	{
	  hGenCorrDR_EWTA_noW_[ict][ipt] = new TH1D(Form("hGenCorrDR_EWTA_noW_%d_%d", ict, ipt), "", 80, 0., 0.2);
	  hGenCorrDR_EWTA_noW_[ict][ipt]->Sumw2();
	}
    }
  */


  // w/ pthat weight
  hDR_Gen_EWTA_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge); // set variable bin

  hNGenJets_W->Sumw2();
  hDR_Gen_EWTA_W->Sumw2();
  hDEta_Gen_EWTA_W->Sumw2();
  hDPhi_Gen_EWTA_W->Sumw2();

  //MC ref gen matched
  //w/o pthat weight 
  hDR_Gen_EWTA_Matched_noW->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin

  hNGenJets_Matched_noW->Sumw2();
  hDR_Gen_EWTA_Matched_noW->Sumw2();
  hDEta_Gen_EWTA_Matched_noW->Sumw2();
  hDPhi_Gen_EWTA_Matched_noW->Sumw2();

  //w/ pthat weight 

  hDR_Gen_EWTA_Matched_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin

  hNGenJets_Matched_W->Sumw2();
  hDR_Gen_EWTA_Matched_W->Sumw2();
  hDEta_Gen_EWTA_Matched_W->Sumw2();
  hDPhi_Gen_EWTA_Matched_W->Sumw2();

  hDR_Gen_EWTA_Matched_PM_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin 

  hNGenJets_Matched_PM_W->Sumw2();
  hDR_Gen_EWTA_Matched_PM_W->Sumw2();

  //MC reco or data
  // w/o pthat weight   

  hDR_EWTA_noW->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin 

  hNJets_noW->Sumw2();
  hDR_EWTA_noW->Sumw2();
  hDEta_EWTA_noW->Sumw2();
  hDPhi_EWTA_noW->Sumw2();

  // w/ pthat weight  

  hDR_EWTA_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin

  hNJets_W->Sumw2();
  hDR_EWTA_W->Sumw2();
  hDEta_EWTA_W->Sumw2();
  hDPhi_EWTA_W->Sumw2();

  //reco matched 
  //w/o pthat weight

  hDR_EWTA_Matched_noW->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin 

  hNJets_Matched_noW->Sumw2();
  hDR_EWTA_Matched_noW->Sumw2();
  hDEta_EWTA_Matched_noW->Sumw2();
  hDPhi_EWTA_Matched_noW->Sumw2();

  //w/ pthat weight

  hDR_EWTA_Matched_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin 

  hNJets_Matched_W->Sumw2();
  hDR_EWTA_Matched_W->Sumw2();
  hDEta_EWTA_Matched_W->Sumw2();
  hDPhi_EWTA_Matched_W->Sumw2();

  hDR_EWTA_Matched_PM_W->GetAxis(0)->Set(DR_nbins, DR_bins_edge);  // set variable bin 

  hNJets_Matched_PM_W->Sumw2();
  hDR_EWTA_Matched_PM_W->Sumw2();
}

void Write_Event_hist(bool is_MC)
{
  hEvents->Write();
  if(is_MC)
    {
      hpthat->Write();
      hpthatW->Write();
    }
  hCent->Write();
  hZvtx->Write();
}

void Write_Jet_QA_hist(bool is_MC) // Write QA histograms
{
  // for test
  hptreco->Write();
  if(is_MC)
    {
      hptreco_Matched->Write();
      hptreco_Matched_PM->Write();
      hptgen->Write();
      hptgen_Matched->Write();
      hptgen_Matched_PM->Write();
    }

  // Jet histo
  //w/o pthat weight
  hJet_RawpT_Eta_Phi_ctbin_nopTCut_noW->Write();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Write();
  hJet_CorrpT_Eta_Phi_ctbin_nopTCut_noW->Write();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Write();
  hJet_RawpT_Eta_Phi_ctbin_pTCut_noW->Write();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Write();
  hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW->Write();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Write();

  //w/ pthat weight
  hJet_RawpT_Eta_Phi_ctbin_nopTCut_W->Write();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Write();
  hJet_CorrpT_Eta_Phi_ctbin_nopTCut_W->Write();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Write();
  hJet_RawpT_Eta_Phi_ctbin_pTCut_W->Write();
  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W->Write();
  hJet_CorrpT_Eta_Phi_ctbin_pTCut_W->Write();
  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W->Write();
  
  // leading and subleading histo
  //w/o pthat weight     
  hld_RawpT_Eta_Phi_ctbin_noW->Write();
  hsld_RawpT_Eta_Phi_ctbin_noW->Write();
  hld_CorrpT_Eta_Phi_ctbin_noW->Write();
  hsld_CorrpT_Eta_Phi_ctbin_noW->Write();

  // w/ weight
  hld_RawpT_Eta_Phi_ctbin_W->Write();
  hsld_RawpT_Eta_Phi_ctbin_W->Write();
  hld_CorrpT_Eta_Phi_ctbin_W->Write();
  hsld_CorrpT_Eta_Phi_ctbin_W->Write();

  if(is_MC)
    {
      hqqbar_Scan_Gen_noW->Write();
      hqqbar_Scan_Gen_W->Write();

      hqqbar_Scan_Reco_noW->Write();
      hqqbar_Scan_Reco_W->Write();
      
      //ref_pt histo
      //w/o pthat weight
      hJet_refpT_ctbin_nopTCut_noW->Write();
      hJet_refpT_ctbin_pTCut_noW->Write();
      hJet_refpT_ctbin_pTCut0_noW->Write();
      //w/ pthat weight
      hJet_refpT_ctbin_nopTCut_W->Write();
      hJet_refpT_ctbin_pTCut_W->Write();
      hJet_refpT_ctbin_pTCut0_W->Write();

      // unmatched, matched and total
      hJet_CorrpT_refparton_ctbin_W->Write();
      hJet_CorrpT_refpartonB_ctbin_W->Write();
      hJet_CorrpT_refpartonB_ctbin_noW->Write();

      hJet_MCorrpT_refparton_ctbin_W->Write();
      hJet_MCorrpT_refpartonB_ctbin_W->Write();

      hJet_UCorrpT_refparton_ctbin_W->Write();
      hJet_UCorrpT_refpartonB_ctbin_W->Write();
      
      //MC gen Jet histo
      //w/o pthat weight
      hJet_GenpT_Eta_Phi_ctbin_nopTCut_noW->Write();
      hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Write();
      hJet_GenpT_Eta_Phi_ctbin_pTCut_noW->Write();
      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Write();
      
      //w/ pthat weight
      hJet_GenpT_Eta_Phi_ctbin_nopTCut_W->Write();
      hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Write();
      hJet_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W->Write();

      // leading and subleading histo
      //w/o pthat weight
      hld_genpT_Eta_Phi_ctbin_noW->Write();
      hsld_genpT_Eta_Phi_ctbin_noW->Write();

      //w/o pthat weight
      hld_genpT_Eta_Phi_ctbin_noW->Write();
      hsld_genpT_Eta_Phi_ctbin_noW->Write();
    }
}

void Write_JES_JER_hist(bool is_JES_JER) // Write JES nd JER histograms
{
  if(is_JES_JER)
    {
      //JES
      //w/o pthat weight
      hJes_rawpT_refpT_pTCut0_refparton_ctbin_noW->Write();
      hJes_rawpT_refpT_pTCut_refparton_ctbin_noW->Write();
      hJes_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Write();
      hJes_CorrpT_refpT_pTCut_refparton_ctbin_noW->Write();
      
      hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Write();
      hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Write();
      hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Write();
      hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Write();
      
      //w/ pthat weight
      hJes_rawpT_refpT_pTCut0_refparton_ctbin_W->Write();
      hJes_rawpT_refpT_pTCut_refparton_ctbin_W->Write();
      hJes_CorrpT_refpT_pTCut0_refparton_ctbin_W->Write();
      hJes_CorrpT_refpT_pTCut_refparton_ctbin_W->Write();
      
      hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Write();
      hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W->Write();
      hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Write();
      hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Write();
      
      //JER
      //w/o pthat weight
      hJer_rawpT_refpT_pTCut0_refparton_ctbin_noW->Write();
      hJer_rawpT_refpT_pTCut_refparton_ctbin_noW->Write();
      hJer_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Write();
      hJer_CorrpT_refpT_pTCut_refparton_ctbin_noW->Write();
      
      hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Write();
      hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Write();
      hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Write();
      hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Write();
      
      //w/ pthat weight
      hJer_rawpT_refpT_pTCut0_refparton_ctbin_W->Write();
      hJer_rawpT_refpT_pTCut_refparton_ctbin_W->Write();
      hJer_CorrpT_refpT_pTCut0_refparton_ctbin_W->Write();
      hJer_CorrpT_refpT_pTCut_refparton_ctbin_W->Write();
      
      hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Write();
      hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W->Write();
      hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Write();
      hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Write();
    }
}

void Write_Corr_hist(bool is_MC)
{

  // w/o pthat weight 
  hNJets_noW->Write();
  hDR_EWTA_noW->Write();
  hDEta_EWTA_noW->Write();
  hDPhi_EWTA_noW->Write();

  // w/ pthat weight  
  hNJets_W->Write();                                                                                                      
  hDR_EWTA_W->Write();
  hDEta_EWTA_W->Write();
  hDPhi_EWTA_W->Write();

  if(is_MC)
    {
      // reco matched 
      //w/o pthat weight
      hNJets_Matched_noW->Write();                                                                                                      
      hDR_EWTA_Matched_noW->Write();
      hDEta_EWTA_Matched_noW->Write();
      hDPhi_EWTA_Matched_noW->Write();

      //w/ pthat weight
      hNJets_Matched_W->Write();                                                                                                      
      hDR_EWTA_Matched_W->Write();
      hDEta_EWTA_Matched_W->Write();
      hDPhi_EWTA_Matched_W->Write();

      hNJets_Matched_PM_W->Write();
      hDR_EWTA_Matched_PM_W->Write();

      // MC gen
      // w/o pthat weight
      hNGenJets_noW->Write();
      hDR_Gen_EWTA_noW->Write();
      hDEta_Gen_EWTA_noW->Write();
      hDPhi_Gen_EWTA_noW->Write();
      /*
      for(int ict = 0; ict < NCentbin; ict++)
	{
	  for(int ipt = 0; ipt < jtpT_nbins; ipt++)
	    {
	      hGenCorrDR_EWTA_noW_[ict][ipt]->Write();
	    }
	}
      */
      // w/ pthat weight 
      hNGenJets_W->Write();
      hDR_Gen_EWTA_W->Write();
      hDEta_Gen_EWTA_W->Write();
      hDPhi_Gen_EWTA_W->Write();

      // MC ref gen matched
      // w/o pthat weight                                                                                                                    
      hNGenJets_Matched_noW->Write();
      hDR_Gen_EWTA_Matched_noW->Write();
      hDEta_Gen_EWTA_Matched_noW->Write();
      hDPhi_Gen_EWTA_Matched_noW->Write();

      // w/ pthat weight   
      hNGenJets_Matched_W->Write();
      hDR_Gen_EWTA_Matched_W->Write();
      hDEta_Gen_EWTA_Matched_W->Write();
      hDPhi_Gen_EWTA_Matched_W->Write();

      hNGenJets_Matched_PM_W->Write();
      hDR_Gen_EWTA_Matched_PM_W->Write();
    }
}
