#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // input variables 
#include "histogram_definition.h" // define histograms
#include "read_tree.h" // read the TChains
#include "particleid.h"  // call for particle id
#include "vector_definition.h"  // define all the vectors
#include "JetCorrector.h" // reader for JEC
//#include "JetUncertainty.h" // reader for JEU
#include "function_defination.h" // function defined here
 
void Tree_Analyzer_MAOD(TString input_file, int itxtoutFile, TString out_file, TString colliding_system)
{
  std::cout<<"colliding_system is: "<<colliding_system.Data()<<std::endl;
  
  if(colliding_system == "pp")
    {
      jet_collection = "ak4PFJetAnalyzer";
      jet_trigger = "HLT_HIAK4CaloJet80_v1";
      JEC_file = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
      use_cent = false;
      doEventFilter = false;

      event_filter_str.resize(0);
      event_filter_str.push_back("pBeamScrapingFilter");
      event_filter_str.push_back("pPAprimaryVertexFilter");
      event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
    }
  
  std::cout<<"event_filter_str size is: "<<event_filter_str.size()<<std::endl;
  std::cout<<endl;

  sumw2();
  TH1::SetDefaultSumw2();

  TFile *xjetFrac = new TFile();
  TH1D* hXjetFrac_vs_genpT = new TH1D();

  if(is_MC && doEventFilter)
    {
      xjetFrac = TFile::Open(Form("ReweightFile/%s", xJetsw.Data()));
      hXjetFrac_vs_genpT = (TH1D*)xjetFrac->Get("hUM_ld_CorrpT_gt_noW");
    }


  //TFile *xjetFrac = TFile::Open("XjetFraction_vs_GenpT_NewNew_Bin10000_NoGenCut_Function.root");
  //TF1* fXjetFrac_vs_genpT = (TF1*)xjetFrac->Get("fgaus");


  vector<string> JECFiles;
  //JECFiles.push_back(Form("JEC_files/%s", JEC_file.c_str()));
  JECFiles.push_back(Form("JEC_files/%s", JEC_file.Data()));

  JetCorrector JEC(JECFiles);
  
  fstream openInputFile;
  openInputFile.open(Form("%s",input_file.Data()), ios::in);
  if(!openInputFile.is_open())
    {
      cout << "List of input files not founded!" << endl;
      return;
    }
  // Make a chain and a vector of file names
  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(openInputFile, file_chain))
    {
      file_name_vector.push_back(file_chain);
      //file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
    }
  openInputFile.close();

  // Read the trees to be added in the Chain
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *ski_tree = new TChain("skimanalysis/HltTree");
  TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data()));
  
  /*
  TChain *trk_tree = new TChain("ppTrack/trackTree");

  TChain *gen_tree = new TChain();
  if(is_MC)
    {
      gen_tree = new TChain("HiGenParticleAna/hi");
    }
  */
  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *testfile = TFile::Open(*listIterator);
      
      if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered))
	{
	  cout << "File: " << *listIterator << " failed to open" << endl;
	  continue;
	}
      
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;

      hlt_tree->Add(*listIterator);
      hea_tree->Add(*listIterator);
      ski_tree->Add(*listIterator);
      jet_tree->Add(*listIterator);


      /*    
      trk_tree->Add(*listIterator);

      if(is_MC)
	{
	  gen_tree->Add(*listIterator);
	}
      */
      
      //testfile->Close();
    }
  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(ski_tree);	
  hlt_tree->AddFriend(jet_tree);

  /*
  hlt_tree->AddFriend(trk_tree);
  if(is_MC)
    {
      hlt_tree->AddFriend(gen_tree);
    }
  */

  read_tree(hlt_tree, is_MC, jet_trigger.Data(), colliding_system.Data(), event_filter_str); // access the tree informations
  
  int nevents = hlt_tree->GetEntries(); // number of events
 
  cout << "Total number of events in those files: "<< nevents << endl;

  int count = 0;
  
  for(int i = 0; i < nevents; i++) //event loop start
  //for(int i = 0; i < 2000; i++) //event loop start
    {
      hlt_tree->GetEntry(i);

      if(i%10000 == 0)
	{
	  std::cout<<i<<"  events running"<<std::endl;
	}
      count = count+nref;
      
      if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue; // vertex cut

      hEvents->AddBinContent(1,1);

      int centbin;
      int ctbin = -1;
      
      if(use_cent)
	{
	  centbin = hiBin;
	  ctbin = hCentBin->FindBin(centbin) - 1;
	  if(centbin >= centCut) continue;
	}
      else // if you use multiplicity 
	{
	  centbin = -1;
	  ctbin = 1;
	  if(centbin >= mult_Cut) continue;
	}

      hEvents->AddBinContent(2,1);

      //std::cout<<"Centrality is: "<<centbin<<"  "<<ctbin<<std::endl;
      double ptHatw = 1.;
      double pTHat = 0.;

      //double testw = weight;
      if(is_MC)
	{
	  ptHatw = weight;
	  pTHat = pthat;
	  
	  if(pTHat == 0 || pTHat <= pthat_cut) continue;
	}
      
      //std::cout<<"before: "<<vertexz<<"   "<<jet_trigger_bit<<"  "<<event<<"  "<<run<<"   "<<centbin<<"  "<<ptHatw<<"  "<<pTHat<<std::endl;
      
      
      //std::cout<<"after: "<<vertexz<<"   "<<jet_trigger_bit<<"  "<<event<<"  "<<run<<"   "<<centbin<<"  "<<ptHatw<<"  "<<pTHat<<std::endl;
      
      hEvents->AddBinContent(3,1);
      
      
      /*
      for(int ii = 0; ii < (int) event_filter_str.size(); ii++)
	{
	  if (event_filter_bool[ii] != 1) continue;
	}
      */

      hEvents->AddBinContent(4,1);

      //if(jet_trigger_bit != 1) continue;

      hEvents->AddBinContent(5,1);

      if(nref <= 0) continue; // if there is no jets in an event

      hEvents->AddBinContent(6,1);

      // prepare for event filter
      double mum_ld_CorrpT = -999., mum_ld_CorrEta = -999., mum_ld_CorrPhi = -999.;// leading matched unmatched jet quantities   
      int mum_ld_index = -999, mum_ld_falvorB = -9999; double mum_ld_refpt = -9999.;

      
      std::vector<double> highest_GenpT;
      double high_GenpT = -999.;

      bool excludeEvents = false;
      double genweight = ptHatw;
      
      if(is_MC && doEventFilter)
        {
          for(int iref = 0; iref < (int)nref; iref++) // reco jet loop to calculate leading jet for the filter
            {
              if(trackMax[iref]/rawpt[iref] < 0.01)continue; // Cut for jets for very low maxium pT track                  
              if(trackMax[iref]/rawpt[iref] > 0.98)continue; // Cut for jets where all the pT is taken by one track        
	      
              JEC.SetJetPT(rawpt[iref]);
              JEC.SetJetEta(jteta[iref]);
              JEC.SetJetPhi(jtphi[iref]);
              float jet_ptcorr = JEC.GetCorrectedPT();
	      
	      int refpartonnn = -99;
              if(fabs(refparton_flavor[iref]) >= 1 && fabs(refparton_flavor[iref]) <= 6)
                {
                  refpartonnn = fabs(refparton_flavor[iref]);
                }
              else if(fabs(refparton_flavor[iref]) == 21)
                {
                  refpartonnn = 7;
                }
              else {refpartonnn = 0;}
	      
	      int refpartonnnB = -99;
              if(fabs(refparton_flavorForB[iref]) >= 1 && fabs(refparton_flavorForB[iref]) <= 6)
                {
                  refpartonnnB = fabs(refparton_flavorForB[iref]);
                }
              else if(fabs(refparton_flavorForB[iref]) == 21)
                {
                  refpartonnnB = 7;
                }
              else{refpartonnnB = 0;}
	    
	      if(refpartonnn == -99 || refpartonnnB == -99) continue;

	      find_leading_Jets(jet_ptcorr, jteta[iref], jtphi[iref], iref, refpt[iref], refpartonnnB, mum_ld_CorrpT, mum_ld_CorrEta,mum_ld_CorrPhi, mum_ld_index, mum_ld_refpt, mum_ld_falvorB);

	    } // reco jet loop end
			
	  for (int igen = 0; igen < (int)ngen; igen++) // gen jet loop to calculate leading jet for the filter
	    {
	      highest_GenpT.push_back(gen_jtpt[igen]);
	    }

	  std::sort(highest_GenpT.begin(), highest_GenpT.end()); // sort the vector
	  
	  if(highest_GenpT.size() >= 1)
	    {
	      high_GenpT = highest_GenpT[highest_GenpT.size()-1];
	    }

	  highest_GenpT.clear(); // clear the gen pT vector
	  	  
          if(mum_ld_index >= 0)
            {
              if((TMath::Abs(mum_ld_CorrEta) < jet_eta_max_cut ) && (mum_ld_CorrpT > leading_pT_min && mum_ld_CorrpT < leading_pT_max))
		{
                  if((int)refpt[mum_ld_index] == -999 && mum_ld_CorrpT > pTHat) // && (int)high_GenpT != -999)
                    {
                      excludeEvents = true;
                    }
                  else excludeEvents = false;
                }
            }
	  
          if(!excludeEvents && (int)high_GenpT != -999)
            {
              double genw = hXjetFrac_vs_genpT->GetBinContent(hXjetFrac_vs_genpT->FindBin(high_GenpT));    
              //double genw = fXjetFrac_vs_genpT->Eval(high_GenpT);
              genweight = (ptHatw)*(1./(1. - genw));
            }
        } // is_MC and doEventFilter condition
      
      //std::cout<<"excludeEvents: "<<(int)excludeEvents<<"  "<<genweight<<"  "<<mum_ld_CorrpT<<"  "<<pTHat<<std::endl;
      
      if(doEventFilter)
	{
	  if(excludeEvents) continue; // filter the events
	}

      hEvents->AddBinContent(7,1);

      hpthat->Fill(pTHat);
      hpthatW->Fill(ptHatw);

      // resize all the vectors used to store jet quantity
      //vectors for checking  
      FilteredJet_jtEta.resize(0);
      FilteredJet_jtPhi.resize(0);
      FilteredJet_WTAEta.resize(0);
      FilteredJet_WTAPhi.resize(0);

      // reco jets vector 
      FilteredJet_RawpT_Vec.resize(0);
      FilteredJet_CorrpT_Vec.resize(0);
      FilteredWTAJet_RawpT_Vec.resize(0);
      FilteredWTAJet_CorrpT_Vec.resize(0);

      refparton_Vec.resize(0);
      refpartonB_Vec.resize(0);
      
      //reco matched jets vector
      Matched_FilteredJet_CorrpT_Vec.resize(0);
      Matched_FilteredWTAJet_CorrpT_Vec.resize(0);

      Matched_refparton_Vec.resize(0);
      Matched_refpartonB_Vec.resize(0);

      Matched_refparton_PM_Vec.resize(0);
      Matched_refpartonB_PM_Vec.resize(0);

      //gen jets vector
      FilteredJet_GenpT_Vec.resize(0);      
      FilteredWTAJet_GenpT_Vec.resize(0);      

      Genrefparton_Vec.resize(0);
      GenrefpartonB_Vec.resize(0);

      //gen matched jets vector
      Matched_FilteredJet_GenpT_Vec.resize(0);
      Matched_FilteredWTAJet_GenpT_Vec.resize(0);

      Matched_Genrefparton_Vec.resize(0);
      Matched_GenrefpartonB_Vec.resize(0);

      Matched_Genrefparton_PM_Vec.resize(0);
      Matched_GenrefpartonB_PM_Vec.resize(0);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Jet loop start
      
      double ld_RawpT=-999., ld_RawEta=-999., ld_RawPhi=-999.; // leading Raw jet quantities                          
      double sld_RawpT=-999., sld_RawEta=-999, sld_RawPhi=-999; // subleading Raw jet quantities                       
      double ld_CorrpT=-999., ld_CorrEta=-999., ld_CorrPhi=-999.; // leading jet quantities                       
      double sld_CorrpT=-999., sld_CorrEta=-999, sld_CorrPhi=-999; // subleading jet quantities    
      double ld_RawWTApT=-999.,  ld_RawWTAEta=-999., ld_RawWTAPhi=-999.; // leading Raw WTA jet quantities             
      double sld_RawWTApT=-999., sld_RawWTAEta=-999, sld_RawWTAPhi=-999; // subleading Raw WTA jet quantities     
      double ld_CorrWTApT=-999.,  ld_CorrWTAEta=-999., ld_CorrWTAPhi=-999.; // leading Jt WTA jet quantities     
      double sld_CorrWTApT=-999., sld_CorrWTAEta=-999, sld_CorrWTAPhi=-999; // subleading Jt WTA jet quantities                            
    
      for (int j = 0; j < nref; j++) //Jet loop start
	{
	  if(trackMax[j]/rawpt[j] < 0.01)continue; // Cut for jets for very low maxium pT track
	  if(trackMax[j]/rawpt[j] > 0.98)continue; // Cut for jets where all the pT is taken by one track

	  JEC.SetJetPT(rawpt[j]); 
	  JEC.SetJetEta(jteta[j]); 
	  JEC.SetJetPhi(jtphi[j]);

	  float jet_pt_corr = JEC.GetCorrectedPT();

	  find_leading_subleading_Jets(rawpt[j], jteta[j], jtphi[j], ld_RawpT, ld_RawEta, ld_RawPhi, sld_RawpT, sld_RawEta, sld_RawPhi);
	  find_leading_subleading_Jets(jet_pt_corr, jteta[j], jtphi[j], ld_CorrpT, ld_CorrEta, ld_CorrPhi, sld_CorrpT, sld_CorrEta, sld_CorrPhi);

	  find_leading_subleading_Jets(rawpt[j], WTAeta[j], WTAphi[j], ld_RawWTApT, ld_RawWTAEta, ld_RawWTAPhi, sld_RawWTApT, sld_RawWTAEta, sld_RawWTAPhi);
	  find_leading_subleading_Jets(jet_pt_corr, WTAeta[j], WTAphi[j], ld_CorrWTApT, ld_CorrWTAEta, ld_CorrWTAPhi, sld_CorrWTApT, sld_CorrWTAEta, sld_CorrWTAPhi);

	  //std::cout<<"raw pT and corrected pT: "<<rawpt[j]<<"  "<<jtpt[j]<<"  "<<jet_pt_corr<<std::endl;

	  //fill the vectors for correlations

	  TVector3 FilteredJet_RawpT;
	  FilteredJet_RawpT.SetPtEtaPhi(rawpt[j], jteta[j], jtphi[j]);
	  FilteredJet_RawpT_Vec.push_back(FilteredJet_RawpT);
	  
	  TVector3 FilteredJet_CorrpT;
	  FilteredJet_CorrpT.SetPtEtaPhi(jet_pt_corr, jteta[j], jtphi[j]);
	  FilteredJet_CorrpT_Vec.push_back(FilteredJet_CorrpT);
	  
	  TVector3 FilteredWTAJet_RawpT;
	  FilteredWTAJet_RawpT.SetPtEtaPhi(rawpt[j], WTAeta[j], WTAphi[j]);
	  FilteredWTAJet_RawpT_Vec.push_back(FilteredWTAJet_RawpT);
	  
	  TVector3 FilteredWTAJet_CorrpT;
	  FilteredWTAJet_CorrpT.SetPtEtaPhi(jet_pt_corr, WTAeta[j], WTAphi[j]);
	  FilteredWTAJet_CorrpT_Vec.push_back(FilteredWTAJet_CorrpT);
	  
	  refparton_Vec.push_back(1);
	  refpartonB_Vec.push_back(1);
	  
	  //FilteredJet_jtEta.push_back(jteta[j]);
	  //FilteredJet_jtPhi.push_back(jtphi[j]);
	  
	  //FilteredJet_WTAEta.push_back(WTAeta[j]);
	  //FilteredJet_WTAPhi.push_back(WTAphi[j]);
	  
	  double Jet_RawpT_Eta_Phi_ctbin[4] = {rawpt[j], jteta[j], jtphi[j], (double)ctbin};
	  double Jet_RawpT_WTAEta_WTAPhi_ctbin[4] = {rawpt[j], WTAeta[j], WTAphi[j], (double)ctbin};
	  
	  double Jet_CorrpT_Eta_Phi_ctbin[4] = {jet_pt_corr, jteta[j], jtphi[j], (double)ctbin};
	  double Jet_CorrpT_WTAEta_WTAPhi_ctbin[4] = {jet_pt_corr, WTAeta[j], WTAphi[j], (double)ctbin};

	  
	  if(jteta[j] > jet_eta_min_cut && jteta[j] < jet_eta_max_cut) 
	    {
	      // w/o pthat weight 
	      hJet_RawpT_Eta_Phi_ctbin_nopTCut_noW->Fill(Jet_RawpT_Eta_Phi_ctbin);
	      hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin);
	      
	      hJet_CorrpT_Eta_Phi_ctbin_nopTCut_noW->Fill(Jet_CorrpT_Eta_Phi_ctbin);
	      hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin);

	      // w/ pthat weight 
	      hJet_RawpT_Eta_Phi_ctbin_nopTCut_W->Fill(Jet_RawpT_Eta_Phi_ctbin, ptHatw);
	      hJet_RawpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin, ptHatw);
	      
	      hJet_CorrpT_Eta_Phi_ctbin_nopTCut_W->Fill(Jet_CorrpT_Eta_Phi_ctbin, ptHatw);
	      hJet_CorrpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin, ptHatw);
	      
	      if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
		{
		  // w/o pthat weight 
		  
		  hJet_RawpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_RawpT_Eta_Phi_ctbin);
		  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin);
		  
		  hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_CorrpT_Eta_Phi_ctbin);
		  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin);

		  // w/ pthat weight 
		  hJet_RawpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_RawpT_Eta_Phi_ctbin, ptHatw);
		  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin, ptHatw);
		  
		  hJet_CorrpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_CorrpT_Eta_Phi_ctbin, ptHatw);
		  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin, ptHatw);
		  
		} // jet pT cut
	    } // jet eta cut
	
	  //if(rawpt[j] <= jet_pt_min_cut || rawpt[j] >= jet_pt_max_cut) continue; // Raw jet pt cut
	  //if(jtpt[j] <= jet_pt_min_cut || jtpt[j] >= jet_pt_max_cut) continue; // jet pt cut
	  //if(jteta[j] <= jet_eta_min_cut || jteta[j] >= jet_eta_max_cut) continue; // jet eta cut
	  //if(WTAeta[j] <= jet_eta_min_cut || WTAeta[j] >= jet_eta_max_cut) continue; // WTA jet eta cut

	  if(is_MC)
	    {
	      double Jet_refpT_ctbin[2] = {refpt[j], (double)ctbin};

	      hJet_refpT_ctbin_nopTCut_noW->Fill(Jet_refpT_ctbin);
	      hJet_refpT_ctbin_nopTCut_W->Fill(Jet_refpT_ctbin, ptHatw);
	      
	      if(refpt[j] > jet_pt_min_cut && refpt[j] < jet_pt_max_cut)
		{
		  hJet_refpT_ctbin_pTCut_noW->Fill(Jet_refpT_ctbin);
		  hJet_refpT_ctbin_pTCut_W->Fill(Jet_refpT_ctbin, ptHatw);
		}
	      
	      int refparton = -99, refparton_PM = -99; 
	      
	      if(fabs(refparton_flavor[j]) >= 1 && fabs(refparton_flavor[j]) <= 6)
		{
		  refparton = fabs(refparton_flavor[j]);

		  if(refparton_flavor[j] > 0)
		    {
		      refparton_PM = fabs(refparton_flavor[j]);
		    }
		  else if(refparton_flavor[j] < 0)
		    {
		      refparton_PM = fabs(refparton_flavor[j])+8;
		    }
		}
	      else if(fabs(refparton_flavor[j]) == 21)
		{
		  refparton = 7;
		  refparton_PM = 7;
		}
	      else
		{
		  refparton = 0;
		  refparton_PM = 0;
		}
	      
	      int refpartonB = -99, refpartonB_PM = -99;
	      if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6)
		{
		  refpartonB = fabs(refparton_flavorForB[j]);

		  if(refparton_flavorForB[j] > 0)
                    {
                      refpartonB_PM = fabs(refparton_flavorForB[j]);
		    }
                  else if(refparton_flavorForB[j] < 0)
                    {
                      refpartonB_PM = fabs(refparton_flavorForB[j])+8;
                    }
		}
	      else if(fabs(refparton_flavorForB[j]) == 21)
		{
		  refpartonB = 7;
		  refpartonB_PM = 7;
		}
	      else
		{
		  refpartonB = 0; 
		  refpartonB_PM = 0;
		}
	      if(refparton == -99 || refparton_PM == -99) continue;
	      if(refpartonB == -99 || refpartonB_PM == -99) continue;

	      // matched, unmatched, and total jet pT distributions

	      double Jet_CorrpT_refparton_ctbin[3] = {jet_pt_corr, (double)refparton, (double)ctbin};
	      double Jet_CorrpT_refpartonB_ctbin[3] = {jet_pt_corr, (double)refpartonB, (double)ctbin};
	      
	      if(jteta[j] > jet_eta_min_cut && jteta[j] < jet_eta_max_cut)
                {
                  if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
                    {
		      if((int)refpt[j] < 0) // unmatched
			{
			  hJet_UCorrpT_refparton_ctbin_W->Fill(Jet_CorrpT_refparton_ctbin, ptHatw);
                          hJet_UCorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, ptHatw);
			}
		      else if((int)refpt[j] >= 0) //matched
			{
			  hJet_MCorrpT_refparton_ctbin_W->Fill(Jet_CorrpT_refparton_ctbin, ptHatw);
                          hJet_MCorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, ptHatw);
			}
		      //total 
		      hJet_CorrpT_refparton_ctbin_W->Fill(Jet_CorrpT_refparton_ctbin, ptHatw);
		      hJet_CorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, ptHatw);
		    }
		}

	      if(refpt[j] <= 0.) continue; // discrad the unmacthed jet (refpt == -999.)
	      if(jet_pt_corr  <= 0.) continue;
	      if(rawpt[j] <= 0.) continue;

	      hJet_refpT_ctbin_pTCut0_noW->Fill(Jet_refpT_ctbin);
	      hJet_refpT_ctbin_pTCut0_noW->Fill(Jet_refpT_ctbin, ptHatw);
	      
	      //fill reco matched vector for correlation

	      TVector3 Matched_FilteredJet_CorrpT;
	      Matched_FilteredJet_CorrpT.SetPtEtaPhi(jet_pt_corr, jteta[j], jtphi[j]);
	      Matched_FilteredJet_CorrpT_Vec.push_back(Matched_FilteredJet_CorrpT);
	      
	      TVector3 Matched_FilteredWTAJet_CorrpT;
	      Matched_FilteredWTAJet_CorrpT.SetPtEtaPhi(jet_pt_corr, WTAeta[j], WTAphi[j]);
	      Matched_FilteredWTAJet_CorrpT_Vec.push_back(Matched_FilteredWTAJet_CorrpT);
	      
	      Matched_refparton_Vec.push_back(refparton);
	      Matched_refpartonB_Vec.push_back(refpartonB);

	      Matched_refparton_PM_Vec.push_back(refparton_PM);
	      Matched_refpartonB_PM_Vec.push_back(refpartonB_PM);

	      if(jteta[j] > jet_eta_min_cut && jteta[j] < jet_eta_max_cut)
		{
		  if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
		    {
		      int refptbin = hJetpTBin->FindBin(refpt[j]) -1;

		      double refparton_refpT[3] = {double(refpartonB_PM), double(refparton_PM), double(refptbin)};
		      //hqqbar_Scan_Reco_noW->Fill(refpartonB_PM, refparton_PM);
		      //hqqbar_Scan_Reco_W->Fill(refpartonB_PM, refparton_PM, ptHatw);

		      hqqbar_Scan_Reco_noW->Fill(refparton_refpT);
		      hqqbar_Scan_Reco_W->Fill(refparton_refpT, ptHatw);
		    }
		}

	      if(is_JES_JER)
		{
		  double JES_ratio_rawpT_vs_refpT = rawpt[j]/refpt[j];
		  double JES_ratio_CorrpT_vs_refpT = jet_pt_corr/refpt[j];
		  
		  double JER_ratio_rawpT_vs_refpT = (rawpt[j] - refpt[j])/refpt[j];
		  double JER_ratio_CorrpT_vs_refpT = (jet_pt_corr - refpt[j])/refpt[j];
		  
		  //JES
		  double Jes_RawpT_refpt_refparton_ctbin[4] = {JES_ratio_rawpT_vs_refpT, refpt[j], (double)refparton, (double)ctbin};
		  double Jes_CorrpT_refpt_refparton_ctbin[4] = {JES_ratio_CorrpT_vs_refpT, refpt[j], (double)refparton, (double)ctbin};
		  
		  double Jes_RawpT_refpt_refpartonB_ctbin[4] = {JES_ratio_rawpT_vs_refpT, refpt[j], (double)refpartonB, (double)ctbin};
		  double Jes_CorrpT_refpt_refpartonB_ctbin[4] = {JES_ratio_CorrpT_vs_refpT, refpt[j], (double)refpartonB, (double)ctbin};
		  
		  //JER
		  double Jer_RawpT_refpt_refparton_ctbin[4] = {JER_ratio_rawpT_vs_refpT, refpt[j], (double)refparton, (double)ctbin};
		  double Jer_CorrpT_refpt_refparton_ctbin[4] = {JER_ratio_CorrpT_vs_refpT, refpt[j], (double)refparton, (double)ctbin};
		  
		  double Jer_RawpT_refpt_refpartonB_ctbin[4] = {JER_ratio_rawpT_vs_refpT, refpt[j], (double)refpartonB, (double)ctbin};
		  double Jer_CorrpT_refpt_refpartonB_ctbin[4] = {JER_ratio_CorrpT_vs_refpT, refpt[j], (double)refpartonB, (double)ctbin};
		  
		  if(jteta[j] > jet_eta_min_cut && jteta[j] < jet_eta_max_cut)
		    {
		      //JES
		      //w/o pthat weight
		      hJes_rawpT_refpT_pTCut0_refparton_ctbin_noW->Fill(Jes_RawpT_refpt_refparton_ctbin);
		      hJes_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Fill(Jes_CorrpT_refpt_refparton_ctbin);
		      
		      hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Fill(Jes_RawpT_refpt_refpartonB_ctbin);
		      hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Fill(Jes_RawpT_refpt_refpartonB_ctbin);
		      //w/ pthat weight 
		      hJes_rawpT_refpT_pTCut0_refparton_ctbin_W->Fill(Jes_RawpT_refpt_refparton_ctbin, ptHatw);
		      hJes_CorrpT_refpT_pTCut0_refparton_ctbin_W->Fill(Jes_CorrpT_refpt_refparton_ctbin, ptHatw);
		      
		      hJes_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Fill(Jes_RawpT_refpt_refpartonB_ctbin, ptHatw);
		      hJes_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Fill(Jes_RawpT_refpt_refpartonB_ctbin, ptHatw);
		      
		      //JER
		      //w/o pthat weight
		      hJer_rawpT_refpT_pTCut0_refparton_ctbin_noW->Fill(Jer_RawpT_refpt_refparton_ctbin);
		      hJer_CorrpT_refpT_pTCut0_refparton_ctbin_noW->Fill(Jer_CorrpT_refpt_refparton_ctbin);
		      
		      hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_noW->Fill(Jer_RawpT_refpt_refpartonB_ctbin);
		      hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_noW->Fill(Jer_RawpT_refpt_refpartonB_ctbin);
		      
		      //w/ pthat weight
		      hJer_rawpT_refpT_pTCut0_refparton_ctbin_W->Fill(Jer_RawpT_refpt_refparton_ctbin, ptHatw);
		      hJer_CorrpT_refpT_pTCut0_refparton_ctbin_W->Fill(Jer_CorrpT_refpt_refparton_ctbin, ptHatw);
		      
		      hJer_rawpT_refpT_pTCut0_refpartonB_ctbin_W->Fill(Jer_RawpT_refpt_refpartonB_ctbin, ptHatw);
		      hJer_CorrpT_refpT_pTCut0_refpartonB_ctbin_W->Fill(Jer_RawpT_refpt_refpartonB_ctbin, ptHatw);
		      
		      //if(jtpt[j] > jet_pt_min_cut && jtpt[j] < jet_pt_max_cut)
		      //if(jtpt[j] > jet_pt_min_cut && jtpt[j] < jet_pt_max_cut && refpt[j] > jet_pt_min_cut && refpt[j] < jet_pt_max_cut)
		      if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
			{
			  //JES
			  //w/o pthat weight
			  hJes_rawpT_refpT_pTCut_refparton_ctbin_noW->Fill(Jes_RawpT_refpt_refparton_ctbin);
			  hJes_CorrpT_refpT_pTCut_refparton_ctbin_noW->Fill(Jes_CorrpT_refpt_refparton_ctbin);
			  
			  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jes_RawpT_refpt_refpartonB_ctbin);
			  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jes_CorrpT_refpt_refpartonB_ctbin);
			  
			  //w/ pthat weight
			  hJes_rawpT_refpT_pTCut_refparton_ctbin_W->Fill(Jes_RawpT_refpt_refparton_ctbin, ptHatw);
			  hJes_CorrpT_refpT_pTCut_refparton_ctbin_W->Fill(Jes_CorrpT_refpt_refparton_ctbin, ptHatw);
			  
			  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jes_RawpT_refpt_refpartonB_ctbin, ptHatw);
			  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jes_CorrpT_refpt_refpartonB_ctbin, ptHatw);
			  
			  //JER
			  //w/o pthat weight 		      
			  hJer_rawpT_refpT_pTCut_refparton_ctbin_noW->Fill(Jer_RawpT_refpt_refparton_ctbin);
			  hJer_CorrpT_refpT_pTCut_refparton_ctbin_noW->Fill(Jer_CorrpT_refpt_refparton_ctbin);
			  
			  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jer_RawpT_refpt_refpartonB_ctbin);
			  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jer_CorrpT_refpt_refpartonB_ctbin);
			  
			  //w/ pthat weight
			  hJer_rawpT_refpT_pTCut_refparton_ctbin_W->Fill(Jer_RawpT_refpt_refparton_ctbin, ptHatw);
			  hJer_CorrpT_refpT_pTCut_refparton_ctbin_W->Fill(Jer_CorrpT_refpt_refparton_ctbin, ptHatw);
			  
			  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jer_RawpT_refpt_refpartonB_ctbin, ptHatw);
			  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jer_CorrpT_refpt_refpartonB_ctbin, ptHatw);
			} // if jes and jer
		    } // Jet pt cut within is_MC
		} // Jet eta cut within is_MC
	    } // is_MC mc gen condition
	} // Jet loop end

      // filling leading and sub-leading quantities
    
      if(nref > 1) // condition for at least 2 jets in an event to find leading and subleading jets                    
	{
	  hEvents->AddBinContent(8,1);
	  
	  if(ld_RawpT != -999. && ld_RawEta != -999. && ld_RawPhi != -999. && ld_RawWTAEta != -999. && ld_RawWTAPhi != -999. && ld_CorrpT != -999. && sld_RawpT != -999. && sld_RawEta != -999. && sld_RawPhi != -999. && sld_RawWTAEta != -999. && ld_RawWTAPhi != -999. && sld_CorrpT != -999.)
	    {
	      double xj_CorrpT = sld_CorrpT/ld_CorrpT;
              double aj_CorrpT = (ld_CorrpT - sld_CorrpT)/(ld_CorrpT+sld_CorrpT);

	      double DPhi = TVector2::Phi_mpi_pi(ld_CorrPhi - sld_CorrPhi);
              double WTADPhi = TVector2::Phi_mpi_pi(ld_CorrWTAPhi - sld_CorrWTAPhi);

	      double Jet_CorrldpT_Eta_Phi_ctbin[4]={ld_CorrpT, ld_CorrEta, ld_CorrPhi, (double)ctbin};
	      double Jet_CorrsldpT_ETa_Phi_ctbin[4]={sld_CorrpT, sld_CorrEta, sld_CorrPhi, (double)ctbin};

	      double Jet_RawldpT_Eta_Phi_ctbin[4]={ld_RawpT, ld_RawEta, ld_RawPhi, (double)ctbin};
	      double Jet_RawsldpT_ETa_Phi_ctbin[4]={sld_RawpT, sld_RawEta, sld_RawPhi, (double)ctbin};

	      if(ld_CorrpT > leading_pT_min && ld_CorrpT < leading_pT_max && sld_CorrpT > subleading_pT_min && sld_CorrpT < subleading_pT_max)
		{
		  if(TMath::Abs(ld_CorrEta) < jet_eta_max_cut && TMath::Abs(sld_CorrEta) < jet_eta_max_cut)
		    {
		      //if(TMath::Abs(DPhi) > leading_subleading_deltaphi_min)
			{
			  hEvents->AddBinContent(9,1);
			  // w/o weight
			  hld_RawpT_Eta_Phi_ctbin_noW->Fill(Jet_RawldpT_Eta_Phi_ctbin);
			  hsld_RawpT_Eta_Phi_ctbin_noW->Fill(Jet_RawsldpT_ETa_Phi_ctbin);

			  hld_CorrpT_Eta_Phi_ctbin_noW->Fill(Jet_CorrldpT_Eta_Phi_ctbin);
			  hsld_CorrpT_Eta_Phi_ctbin_noW->Fill(Jet_CorrsldpT_ETa_Phi_ctbin);

			  // w/ weight
			  hld_RawpT_Eta_Phi_ctbin_W->Fill(Jet_RawldpT_Eta_Phi_ctbin, ptHatw);
                          hsld_RawpT_Eta_Phi_ctbin_W->Fill(Jet_RawsldpT_ETa_Phi_ctbin, ptHatw);

                          hld_CorrpT_Eta_Phi_ctbin_W->Fill(Jet_CorrldpT_Eta_Phi_ctbin, ptHatw);
                          hsld_CorrpT_Eta_Phi_ctbin_W->Fill(Jet_CorrsldpT_ETa_Phi_ctbin, ptHatw);

			} // Dphi condition
		    } // leading and subleading jet eta condition
		} // leading and subleading jet pt condition 
	    } // -999. condition
	} // nref > 1 condition

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jet loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // gen jet loop
      if(is_MC)
        {
	  double ld_genpT=-999., ld_genEta=-999., ld_genPhi=-999.; // leading Raw jet quantities                                     
	  double sld_genpT=-999., sld_genEta=-999, sld_genPhi=-999; // subleading Raw jet quantities          
	  double ld_genWTApT=-999., ld_genWTAEta=-999., ld_genWTAPhi=-999.; // leading Raw jet quantities                         
	  double sld_genWTApT=-999., sld_genWTAEta=-999, sld_genWTAPhi=-999; // subleading Raw jet quantities          
	  
	  for (int ign = 0; ign < (int)ngen; ign++) //gen Jet loop start                                        
            {
	      find_leading_subleading_Jets(gen_jtpt[ign], gen_jteta[ign], gen_jtphi[ign], ld_genpT, ld_genEta, ld_genPhi, sld_genpT, sld_genEta, sld_genPhi); // standard axis
	      
              find_leading_subleading_Jets(gen_jtpt[ign], gen_WTAeta[ign], gen_WTAphi[ign], ld_genWTApT, ld_genWTAEta, ld_genWTAPhi, sld_genWTApT, sld_genWTAEta, sld_genWTAPhi); // WTA axis
	      
	      //fill the vectors for correlations

	      TVector3 FilteredJet_GenpT;
	      FilteredJet_GenpT.SetPtEtaPhi(gen_jtpt[ign], gen_jteta[ign], gen_jtphi[ign]);
	      FilteredJet_GenpT_Vec.push_back(FilteredJet_GenpT);
	      
	      TVector3 FilteredWTAJet_GenpT;
	      FilteredWTAJet_GenpT.SetPtEtaPhi(gen_jtpt[ign], gen_WTAeta[ign], gen_WTAphi[ign]);
	      FilteredWTAJet_GenpT_Vec.push_back(FilteredWTAJet_GenpT);

	      Genrefparton_Vec.push_back(1);
	      GenrefpartonB_Vec.push_back(1);	      

	      //int nrfgn = 0;
	      int index_ref = -1;
	      int index_gen = -1;
	      double min_Dpt = 999999999.;

	      //std::cout<<std::endl;
	      for (int irf = 0; irf < (int)nref; irf++) //ref Jet loop start to calculate the parton flavour
		{
		  //std::cout<<"Gen pT And Ref pT are: "<<gen_jtpt[ign]<<"  "<<refpt[irf]<<"  "<<ign<<"  "<<irf<<std::endl;
		  double dpt = fabs(fabs(gen_jtpt[ign]) - fabs(refpt[irf]));

		  //if(gen_jtpt[ign] == refpt[irf])
		  if(dpt < min_Dpt)
		    {
		      min_Dpt = dpt;
		      index_ref = irf;
		      index_gen = ign;
		      //nrfgn++;
		    }
		}
	      //std::cout<<std::endl;
	      /*
	      if(nrfgn > 1)
		{
		  std::cout<<"More than one gen Jet matched with ref Jet, Please check"<<std::endl;
		}
	      */
	      //if(nrfgn == 1)
		//{
		  if(index_ref >= 0 && index_gen >= 0)
		    {
		      //std::cout<<"Matched Gen pT And Ref pT are: "<<gen_jtpt[index_gen]<<"  "<<refpt[index_ref]<<"  "<<index_gen<<"  "<<index_ref<<std::endl;
		      int refparton_gn = -99, refparton_gn_PM = -99;
		      if(fabs(refparton_flavor[index_ref]) >= 1 && fabs(refparton_flavor[index_ref]) <= 6)
			{
			  refparton_gn = fabs(refparton_flavor[index_ref]);
			  
			  if(refparton_flavor[index_ref] > 0)
			    {
			      refparton_gn_PM = fabs(refparton_flavor[index_ref]);
			    }
			  else if(refparton_flavor[index_ref] < 0)
			    {
			      refparton_gn_PM = fabs(refparton_flavor[index_ref])+8;
			    }
			}
		      else if(fabs(refparton_flavor[index_ref]) == 21)
			{
			  refparton_gn = 7;
			  refparton_gn_PM = 7;
			}
		      else 
			{
			  refparton_gn = 0;
			  refparton_gn_PM = 0;
			}
		      
		      int refpartonB_gn = -99, refpartonB_gn_PM = -99;
		      if(fabs(refparton_flavorForB[index_ref]) >= 1 && fabs(refparton_flavorForB[index_ref]) <= 6)
			{
			  refpartonB_gn = fabs(refparton_flavorForB[index_ref]);

			  if(refparton_flavorForB[index_ref] > 0)
                            {
                              refpartonB_gn_PM = fabs(refparton_flavorForB[index_ref]);
                            }
                          else if(refparton_flavorForB[index_ref] < 0)
                            {
                              refpartonB_gn_PM = fabs(refparton_flavorForB[index_ref])+8;
                            }
			}
		      else if(fabs(refparton_flavorForB[index_ref]) == 21)
			{
			  refpartonB_gn = 7;
			  refpartonB_gn_PM = 7;
			}
		      else
			{
			  refpartonB_gn = 0;
			  refpartonB_gn_PM = 0;
			}

		      if(refparton_gn == -99 || refparton_gn_PM == -99) continue;
		      if(refpartonB_gn == -99 || refpartonB_gn_PM == -99) continue;

		      TVector3 Matched_FilteredJet_GenpT;
		      Matched_FilteredJet_GenpT.SetPtEtaPhi(gen_jtpt[index_gen], gen_jteta[index_gen], gen_jtphi[index_gen]);
		      Matched_FilteredJet_GenpT_Vec.push_back(Matched_FilteredJet_GenpT);

		      TVector3 Matched_FilteredWTAJet_GenpT;
		      Matched_FilteredWTAJet_GenpT.SetPtEtaPhi(gen_jtpt[index_gen], gen_WTAeta[index_gen], gen_WTAphi[index_gen]);
		      Matched_FilteredWTAJet_GenpT_Vec.push_back(Matched_FilteredWTAJet_GenpT);

		      Matched_Genrefparton_Vec.push_back(refparton_gn);
		      Matched_GenrefpartonB_Vec.push_back(refpartonB_gn);

		      Matched_Genrefparton_PM_Vec.push_back(refparton_gn_PM);
		      Matched_GenrefpartonB_PM_Vec.push_back(refpartonB_gn_PM);

		      if(gen_jteta[index_gen] > jet_eta_min_cut && gen_jteta[index_gen] < jet_eta_max_cut)
			{
			  if(gen_jtpt[index_gen] > jet_pt_min_cut && gen_jtpt[index_gen] < jet_pt_max_cut)
			    {
			      int genptbin = hJetpTBin->FindBin(gen_jtpt[index_gen]) -1;

			      double refparton_genpT[3] = {double(refpartonB_gn_PM), double(refparton_gn_PM), double(genptbin)};
			      hqqbar_Scan_Gen_noW->Fill(refparton_genpT);
			      hqqbar_Scan_Gen_W->Fill(refparton_genpT, genweight);

			      //hqqbar_Scan_Gen_noW->Fill(refpartonB_gn_PM, refparton_gn_PM);
			      //hqqbar_Scan_Gen_W->Fill(refpartonB_gn_PM, refparton_gn_PM, genweight);
			    }
			}
		    }
		  //}
	      
	      double Jet_GenpT_Eta_Phi_ctbin[4] = {gen_jtpt[ign], gen_jteta[ign], gen_jtphi[ign], (double)ctbin};
	      double Jet_GenpT_WTAEta_WTAPhi_ctbin[4] = {gen_jtpt[ign], gen_WTAeta[ign], gen_WTAphi[ign], (double)ctbin};
	      
	      if(gen_jteta[ign] > jet_eta_min_cut && gen_jteta[ign] < jet_eta_max_cut)
		{
		  //w/o pthat weight
		  hJet_GenpT_Eta_Phi_ctbin_nopTCut_noW->Fill(Jet_GenpT_Eta_Phi_ctbin);
		  hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_noW->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin);
		  
		  //w/o pthat weight
		  hJet_GenpT_Eta_Phi_ctbin_nopTCut_W->Fill(Jet_GenpT_Eta_Phi_ctbin, genweight);
		  hJet_GenpT_WTAEta_WTAPhi_ctbin_nopTCut_W->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin, genweight);
		  
		  if(gen_jtpt[ign] > jet_pt_min_cut && gen_jtpt[ign] < jet_pt_max_cut)
		    {
		      //w/o pthat weight
		      hJet_GenpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_GenpT_Eta_Phi_ctbin);
		      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin);
		      
		      //w/o pthat weight
		      hJet_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_GenpT_Eta_Phi_ctbin, genweight);
		      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin, genweight);
		      
		    } // pt condition
		} // eta condition
            }// gen jet loop end
	  
	  //filling leading and subleading gen jet quantities
	  if(ngen > 1) // condition for at least 2 jets in an event to find leading and subleading jets                    
	    {
	      if(ld_genpT != -999. && ld_genEta != -999. && ld_genPhi != -999. && ld_genWTAEta != -999. && ld_genWTAPhi != -999. && sld_genpT != -999. && sld_genEta != -999. && sld_genPhi != -999. && sld_genWTAEta != -999. && ld_genWTAPhi != -999.) 
		{
		  double xj_genpT = sld_genpT/ld_genpT;
		  double aj_genpT = (ld_genpT - sld_genpT)/(ld_genpT+sld_genpT);
		  
		  double DPhi_gen = TVector2::Phi_mpi_pi(ld_genPhi - sld_genPhi);
		  double WTADPhi_gen = TVector2::Phi_mpi_pi(ld_genWTAPhi - sld_genWTAPhi);
		  
		  double Jet_genldpT_Eta_Phi_ctbin[4]={ld_genpT, ld_genEta, ld_genPhi, (double)ctbin};
		  double Jet_gensldpT_ETa_Phi_ctbin[4]={sld_genpT, sld_genEta, sld_genPhi, (double)ctbin};
		  
		  if(ld_genpT > leading_pT_min && ld_genpT < leading_pT_max && sld_genpT > subleading_pT_min && sld_genpT < subleading_pT_max)
		    {
		      if(TMath::Abs(ld_genEta) < jet_eta_max_cut && TMath::Abs(sld_genEta) < jet_eta_max_cut)
			{
			  //if(TMath::Abs(DPhi) > leading_subleading_deltaphi_min)
			  {
			    // w/o weight
			    hld_genpT_Eta_Phi_ctbin_noW->Fill(Jet_genldpT_Eta_Phi_ctbin);
			    hsld_genpT_Eta_Phi_ctbin_noW->Fill(Jet_gensldpT_ETa_Phi_ctbin);
			    
			    // w/ weight
			    hld_genpT_Eta_Phi_ctbin_W->Fill(Jet_genldpT_Eta_Phi_ctbin, genweight);
			    hsld_genpT_Eta_Phi_ctbin_W->Fill(Jet_gensldpT_ETa_Phi_ctbin, genweight);

			  } // Dphi condition
			} // leading and subleading jet eta condition
		    } // leading and subleading jet pt condition 
		} // -999. condition
	    } // ngen > 1 condition
	} // is MC condition

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      //~~~~~~~~~~~~~~~~:start correlation between E scheme and WTA axes in MC reco or data:~~~~~~~~~~~~~~~~~~~~~~~~~

      DeltaR_corr_EWTA(FilteredJet_CorrpT_Vec, FilteredWTAJet_CorrpT_Vec, ctbin, refparton_Vec, refpartonB_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, ptHatw, hNJets_noW, hNJets_W, hDR_EWTA_noW, hDR_EWTA_W, hDEta_EWTA_noW, hDEta_EWTA_W, hDPhi_EWTA_noW, hDPhi_EWTA_W, hJetpTBin);
      
      // start the correlation between E scheme and WTA axes in MC gen    
      if(is_MC)
	{
	  // reco matched 
	  DeltaR_corr_EWTA(Matched_FilteredJet_CorrpT_Vec, Matched_FilteredWTAJet_CorrpT_Vec, ctbin, Matched_refparton_Vec, Matched_refpartonB_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, ptHatw, hNJets_Matched_noW, hNJets_Matched_W, hDR_EWTA_Matched_noW, hDR_EWTA_Matched_W, hDEta_EWTA_Matched_noW, hDEta_EWTA_Matched_W, hDPhi_EWTA_Matched_noW, hDPhi_EWTA_Matched_W, hJetpTBin);

	  // quark antiquark
	  DeltaR_corr_EWTA_W(Matched_FilteredJet_CorrpT_Vec, Matched_FilteredWTAJet_CorrpT_Vec, ctbin, Matched_refparton_PM_Vec, Matched_refpartonB_PM_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, ptHatw, hNJets_Matched_PM_W, hDR_EWTA_Matched_PM_W, hJetpTBin);
	  
	  // gen
	  DeltaR_corr_EWTA(FilteredJet_GenpT_Vec, FilteredWTAJet_GenpT_Vec, ctbin, Genrefparton_Vec, GenrefpartonB_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, genweight, hNGenJets_noW, hNGenJets_W, hDR_Gen_EWTA_noW, hDR_Gen_EWTA_W, hDEta_Gen_EWTA_noW, hDEta_Gen_EWTA_W, hDPhi_Gen_EWTA_noW, hDPhi_Gen_EWTA_W, hJetpTBin);
	  
	  //gen ref matched
	  DeltaR_corr_EWTA(Matched_FilteredJet_GenpT_Vec, Matched_FilteredWTAJet_GenpT_Vec, ctbin, Matched_Genrefparton_Vec, Matched_GenrefpartonB_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, genweight, hNGenJets_Matched_noW, hNGenJets_Matched_W, hDR_Gen_EWTA_Matched_noW, hDR_Gen_EWTA_Matched_W, hDEta_Gen_EWTA_Matched_noW, hDEta_Gen_EWTA_Matched_W, hDPhi_Gen_EWTA_Matched_noW, hDPhi_Gen_EWTA_Matched_W, hJetpTBin);

	  // quark antiquark
	  DeltaR_corr_EWTA_W(Matched_FilteredJet_GenpT_Vec, Matched_FilteredWTAJet_GenpT_Vec, ctbin, Matched_Genrefparton_PM_Vec, Matched_GenrefpartonB_PM_Vec, jet_pt_min_cut, jet_pt_max_cut, jet_eta_min_cut, jet_eta_max_cut, genweight, hNGenJets_Matched_PM_W, hDR_Gen_EWTA_Matched_PM_W, hJetpTBin);

	} // is_MC 
      
      //~~~~~~~~~~~~~~~~:end correlation between E scheme and WTA axes in MC reco or data:~~~~~~~~~~~~~~~~~~~~~~~~~

      //clear all the vectors filled in event wise

      FilteredJet_jtEta.clear();
      FilteredJet_jtPhi.clear();
      FilteredJet_WTAEta.clear();
      FilteredJet_WTAPhi.clear();

      // reco jets vector       
      FilteredJet_RawpT_Vec.clear();
      FilteredJet_CorrpT_Vec.clear();
      FilteredWTAJet_RawpT_Vec.clear();
      FilteredWTAJet_CorrpT_Vec.clear();

      refparton_Vec.clear();
      refpartonB_Vec.clear();

      //reco matched jets vector   
      Matched_FilteredJet_CorrpT_Vec.clear();
      Matched_FilteredWTAJet_CorrpT_Vec.clear();

      Matched_refparton_Vec.clear();
      Matched_refpartonB_Vec.clear();

      Matched_refparton_PM_Vec.clear();
      Matched_refpartonB_PM_Vec.clear();

      //gen jets vector   
      FilteredJet_GenpT_Vec.clear();
      FilteredWTAJet_GenpT_Vec.clear();

      Genrefparton_Vec.clear();
      GenrefpartonB_Vec.clear();

      //gen matched jets vector   
      Matched_FilteredJet_GenpT_Vec.clear();
      Matched_FilteredWTAJet_GenpT_Vec.clear();

      Matched_Genrefparton_Vec.clear();
      Matched_GenrefpartonB_Vec.clear();

      Matched_Genrefparton_PM_Vec.clear();
      Matched_GenrefpartonB_PM_Vec.clear();

    } // events loop end

  std::cout<<"Total nref is: "<<count<<std::endl;

  delete hea_tree;
  delete hlt_tree;
  delete ski_tree;
  delete jet_tree;

  //delete trk_tree;
  //delete gen_tree;

  
  std::string outfilename = Form("%s/%s_Outfile_NoRefpTCut_EventpTHatFilter_pTHat%1.1f_JetpT%1.1f_LdpT%1.1f_SldpT%1.1f_NoTrigger_%d",out_file.Data(), colliding_system.Data(), pthat_cut, jet_pt_min_cut, leading_pT_min, subleading_pT_min, itxtoutFile);

  std::replace(outfilename.begin(), outfilename.end(), '.', 'p'); // replace . to p
  std::replace(outfilename.begin(), outfilename.end(), '-', 'N'); // replace - to N for negative

  TFile* fout = new TFile(Form("%s.root", outfilename.c_str()), "recreate");	 

  //TFile* fout = new TFile(Form("%s/%s_Outfile_NoRefpTCut_pTHat%1.1f_JetpT%1.1f_LdpT%1.1f_SldpT%1.1f_NoTrigger_%d.root",out_file.Data(), colliding_system.Data(), pthat_cut, jet_pt_min_cut, leading_pT_min, subleading_pT_min, itxtoutFile), "recreate");
  fout->mkdir("Event_Hist");
  fout->cd("Event_Hist");
  Write_Event_hist(is_MC);

  fout->mkdir("Jet_QA_Hist");
  fout->cd("Jet_QA_Hist");
  Write_Jet_QA_hist(is_MC);

  if(is_MC)
    {
      fout->mkdir("JES_JER_Hist");
      fout->cd("JES_JER_Hist");
      Write_JES_JER_hist(is_JES_JER);
    }

  fout->mkdir("Corr_Hist");
  fout->cd("Corr_Hist");
  Write_Corr_hist(is_MC);
  
  fout->Write();
  fout->Close();
}

// main program
int main(int argc, char **argv)
{
  using namespace std;

  TString inputfile;
  int itxtout;
  TString outfile;
  TString coll_sys;

  if(argc == 1)
    {
      std::cout<<"You did not pass any argument to the code other than the program name"<<std::endl;
      std::cout<<"You need to pass 4 arguments including the program name"<<std::endl;
    }

  if(argc >=1 && argc <=4)
    {
      std::cout<<"Only "<<argc<<" arguments you have given including the program name"<<std::endl;
      std::cout<<"You need to pass 4 arguments including the program name"<<std::endl;
    }

  if(argc == 5)
    {
      std::cout<<std::endl;
      std::cout<<"You have given 5 arguments including the program name;  Your program will run"<<std::endl;
      std::cout<<std::endl;

      inputfile = argv[1];
      itxtout = atoi(argv[2]);
      outfile = argv[3];
      coll_sys = argv[4];

      Tree_Analyzer_MAOD(inputfile, itxtout, outfile, coll_sys);
    }
  return 0;
}
