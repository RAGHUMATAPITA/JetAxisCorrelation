#include "call_libraries.h"  // call libraries from ROOT and C++

// declare variables

//~~~~~~~~~~~~~~~~~~~~~ event quantities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
float vertexz; // event z vertex
int hiBin; // event centrality (used if use_centrality = true in input_variables.h)
UInt_t run;
UInt_t lumi;
ULong64_t event;
//Int_t event;

// events quantities from gen
float weight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)

// trigger quantities
int jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

//event filter
std::vector<int> event_filter_bool(10);

//~~~~~~~~~~~~~~~~~~~~~ jet and track quantities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// reco jets
const int nmax = 99999;
int nref;         // number of jets
float jtpt[nmax]; // jet pT
float jteta[nmax]; // jet eta
float jtphi[nmax]; // jet phi
float rawpt[nmax]; // jet pT without JEC
float trackMax[nmax]; // track maximum pT in a jet
float WTAeta[nmax]; // WTA jet eta
float WTAphi[nmax]; // WTA jet phi

// reco tracks
int ntrk;                 // number of track
float trkpt[nmax];       // track pT
float trketa[nmax];      // track eta
float trkphi[nmax];      // track phi
float trkpterr[nmax];    // track pT error (uncertainty)
float trkdcaxy[nmax];    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaz[nmax];     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaxyerr[nmax]; // track dxy error (uncertainty)
float trkdcazerr[nmax];  // track dxy error (uncertainty)
float trkchi2[nmax];     // track reconstruction chi2 of the fitting
float pfEcal[nmax];      // particle flow energy deposit in ECAL
float pfHcal[nmax];      // particle flow energy deposit in HCAL
float trkmva[nmax];      // track mva for each step
UChar_t trkalgo[nmax];       // track algorithm/step
UChar_t trkndof[nmax];       // track number of degrees of freedom in the fitting 
int trkcharge[nmax];     // track charge
UChar_t trknhits[nmax];      // number of hits in the tracker
UChar_t trknlayer[nmax];     // number of layers with measurement in the tracker
bool highpur[nmax];      // tracker steps MVA selection

// gen jets
int ngen;             // number of gen jets
float gen_jtpt[nmax];  // gen jet pT
float gen_jteta[nmax]; // gen jet eta
float gen_jtphi[nmax]; // gen jet phi
float gen_WTAeta[nmax]; // gen WTA jet eta
float gen_WTAphi[nmax]; // gen ETA jet phi

// matched jets
float refpt[nmax]; // jet pT matched with Gen pT
float refeta[nmax]; // jet eta matched with Gen eta
float refphi[nmax]; // jet phi matched with Gen phi
int refparton_flavor[nmax]; // jet flavor matched with Gen flavor
int refparton_flavorForB[nmax]; // jet flavor matched with Gen flavor
//float refparton_flavorForB[nmax]; //for pp (jtpartonflavor stored in float) 

// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;   // gen particle pid

void read_tree(TChain *tree, bool is_MC, TString jet_trigger, TString colliding_system, std::vector<TString> event_filterstr)
{
  tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files
  
  // enable branches of interest -> see definition of each variables above
  
  // event quantities
  
  tree->SetBranchStatus(Form("%s",jet_trigger.Data()), 1);
  tree->SetBranchAddress(Form("%s",jet_trigger.Data()), &jet_trigger_bit);
  
  tree->SetBranchStatus("vz", 1);
  tree->SetBranchAddress("vz", &vertexz);
  
  tree->SetBranchStatus("evt", 1);
  tree->SetBranchAddress("evt", &event);
  
  tree->SetBranchStatus("run", 1);
  tree->SetBranchAddress("run", &run);

  tree->SetBranchStatus("lumi", 1);
  tree->SetBranchAddress("lumi", &lumi);

  if(colliding_system=="PbPb")
    {
      tree->SetBranchStatus("hiBin", 1);
      tree->SetBranchAddress("hiBin", &hiBin);
    }
  
  for(int i = 0; i < (int) event_filterstr.size(); i++)
    {
      tree->SetBranchStatus(Form("%s",event_filterstr[i].Data()), 1);
      tree->SetBranchAddress(Form("%s",event_filterstr[i].Data()),&event_filter_bool[i]);
    }
  
  if(is_MC)
    {
      tree->SetBranchStatus("weight", 1);
      tree->SetBranchAddress("weight", &weight);
      
      tree->SetBranchStatus("pthat", 1); 
      tree->SetBranchAddress("pthat", &pthat);
    }

    // reco jet quantities
    tree->SetBranchStatus("nref", 1);
    tree->SetBranchAddress("nref", &nref);
    
    tree->SetBranchStatus("rawpt", 1);
    tree->SetBranchAddress("rawpt", &rawpt);
    
    tree->SetBranchStatus("trackMax", 1);
    tree->SetBranchAddress("trackMax", &trackMax);
    
    tree->SetBranchStatus("jtpt", 1);
    tree->SetBranchAddress("jtpt", &jtpt);
    
    tree->SetBranchStatus("jteta", 1);
    tree->SetBranchAddress("jteta", &jteta);
    
    tree->SetBranchStatus("jtphi", 1);
    tree->SetBranchAddress("jtphi", &jtphi);
    
    tree->SetBranchStatus("WTAeta", 1);
    tree->SetBranchAddress("WTAeta", &WTAeta);
    
    tree->SetBranchStatus("WTAphi", 1);
    tree->SetBranchAddress("WTAphi", &WTAphi);    
  
    if (is_MC)
      {
	// gen jet quantities
        tree->SetBranchStatus("ngen", 1);
	tree->SetBranchAddress("ngen", &ngen);
	
        tree->SetBranchStatus("genpt", 1);
	tree->SetBranchAddress("genpt", &gen_jtpt);
	
	tree->SetBranchStatus("geneta", 1);
	tree->SetBranchAddress("geneta", &gen_jteta);
	 
	tree->SetBranchStatus("genphi", 1);
	tree->SetBranchAddress("genphi", &gen_jtphi);
	
	tree->SetBranchStatus("WTAgeneta", 1);
	tree->SetBranchAddress("WTAgeneta", &gen_WTAeta);
	
	tree->SetBranchStatus("WTAgenphi", 1);
	tree->SetBranchAddress("WTAgenphi", &gen_WTAphi);

	//matched gen jet quantities
	tree->SetBranchStatus("refpt", 1);
        tree->SetBranchAddress("refpt", &refpt);
	
        tree->SetBranchStatus("refeta", 1);
        tree->SetBranchAddress("refeta", &refeta);
	
        tree->SetBranchStatus("refphi", 1);
	tree->SetBranchAddress("refphi", &refphi);

	if(colliding_system=="PbPb")
	  {
	    tree->SetBranchStatus("refparton_flavor", 1);
	    tree->SetBranchAddress("refparton_flavor", &refparton_flavor);

	    tree->SetBranchStatus("refparton_flavorForB", 1);
	    tree->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
	  }
	else if(colliding_system=="pp")
	  {
	    tree->SetBranchStatus("refparton_flavor", 1);
            tree->SetBranchAddress("refparton_flavor", &refparton_flavor);
	    
	    tree->SetBranchStatus("jtPartonFlavor", 1);
	    tree->SetBranchAddress("jtPartonFlavor", &refparton_flavorForB);
	  }

	// gen particle quantities
	tree->SetBranchStatus("pt", 1);
	tree->SetBranchAddress("pt", &gen_trkpt);
	
        tree->SetBranchStatus("eta", 1);
	tree->SetBranchAddress("eta", &gen_trketa);
	
        tree->SetBranchStatus("phi", 1);
	tree->SetBranchAddress("phi", &gen_trkphi);
	
        tree->SetBranchStatus("chg", 1);
	tree->SetBranchAddress("chg", &gen_trkchg);
	/*	
        tree->SetBranchStatus("pdg", 1);
	tree->SetBranchAddress("pdg", &gen_trkpid);
	
        tree->SetBranchStatus("sube", 1);
        tree->SetBranchAddress("sube", &gen_trksube);
	*/
      }


    // track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchAddress("nTrk", &ntrk);
    
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchAddress("trkPt", &trkpt);
    
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchAddress("trkEta", &trketa);
    
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchAddress("trkPhi", &trkphi);

    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchAddress("trkCharge", &trkcharge);    

    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchAddress("trkDxy1", &trkdcaxy);
    
    tree->SetBranchStatus("trkDxyError1", 1);
    tree->SetBranchAddress("trkDxyError1", &trkdcaxyerr);
    
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchAddress("trkDz1", &trkdcaz);
    
    tree->SetBranchStatus("trkDzError1", 1);
    tree->SetBranchAddress("trkDzError1", &trkdcazerr);
    
    tree->SetBranchStatus("trkChi2", 1);
    tree->SetBranchAddress("trkChi2", &trkchi2);
    
    tree->SetBranchStatus("trkNdof", 1);
    tree->SetBranchAddress("trkNdof", &trkndof);
    
    tree->SetBranchStatus("trkNHit", 1);
    tree->SetBranchAddress("trkNHit", &trknhits);    
    
    tree->SetBranchStatus("trkNlayer", 1);
    tree->SetBranchAddress("trkNlayer", &trknlayer);
    
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchAddress("highPurity", &highpur);
    
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    
    tree->SetBranchStatus("pfHcal", 1);
    tree->SetBranchAddress("pfHcal", &pfHcal);

    if(colliding_system=="PbPb")
      {
	//special for 2018 PbPb MC
        tree->SetBranchStatus("trkMVA", 1);
	tree->SetBranchAddress("trkMVA", &trkmva);
	
        tree->SetBranchStatus("trkAlgo", 1);
	tree->SetBranchAddress("trkAlgo", &trkalgo);
      }
}
