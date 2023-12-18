#include "call_libraries.h"

bool pTHatFilter(double jetPt, double pthat)
{
  bool result = false;
  
  if(jetPt < 1.702*pthat + 9.701) result = true;
  
  return result;
}

void find_leading_subleading_Jets(double pt, double eta, double phi, double &leadpt, double &leadeta, double &leadphi, double &sublpt, double &subleta, double &sublphi)
{
  if( pt > leadpt )
    {
      sublpt = leadpt;
      leadpt = pt;
      leadeta = eta;
      leadphi = phi;
    }
  else if( pt > sublpt)
    {
      sublpt = pt;
      subleta = eta;
      sublphi = phi;
    }
}

void find_leading_Jets(double pt, double eta, double phi, int index, double refpt, int falvorB, double &leadpt, double &leadeta, double &leadphi, int &leadindex, double &leadrefpt, int &leadfalvorB)
{
  if( pt > leadpt )
    {
      leadpt = pt;
      leadeta = eta;
      leadphi = phi;
      leadindex = index;
      leadrefpt = refpt;
      leadfalvorB = falvorB;
    }
}

void find_leading_Tracks(double pt, double eta, double phi, int index, double &leadpt, double &leadeta, double &leadphi, int &leadindex)
{
  if( pt > leadpt )
    {
      leadpt = pt;
      leadeta = eta;
      leadphi = phi;
      leadindex = index;
    }
}

void DeltaR_corr_EWTA(std::vector<TVector3> Tvec1, std::vector<TVector3> Tvec2, int centBin, std::vector<int> refparton, std::vector<int> refpartonB, double pTmin, double pTmax, double Etamin, double Etamax, double Evtweight, THnSparseD* sNjet_noW, THnSparseD* sNjet_W, THnSparseD* hDR_noW, THnSparseD* hDR_W, THnSparseD* hDEta_noW, THnSparseD* hDEta_W, THnSparseD* hDPhi_noW, THnSparseD* hDPhi_W, TH1D* hpTBin, TH1D* hptdist, TF1* fpt, bool ismc, bool is_weight, bool isptw)
{

  //std::cout<<"size of :"<<Tvec1.size()<<"  "<<Tvec2.size()<<"  "<<refparton.size()<<std::endl;

  if(Tvec1.size() != Tvec2.size())
    {
      std::cout<<"~~~~~~~~:Two vectors sizes are not same (something wrong), Please check:~~~~~~~~~"<<std::endl;
    }
  if(Tvec1.size() != refparton.size())
    {
      std::cout<<"~~~~~~~~:refparton vec size are not same as Tvec1(something wrong), Please check:~~~~~~~~~"<<std::endl;
    }

  const int xndim = 4;

  for(int i1 = 0; i1 < (int)Tvec1.size(); i1++)
    {
      TVector3 vec1 = Tvec1[i1];
      double pT1 = vec1.Pt();
      double Eta1 = vec1.Eta();
      double Phi1 = vec1.Phi();
      int refpar = refparton[i1];
      int refparB = refpartonB[i1];

      int ptbin = hpTBin->FindBin(pT1) -1;

      double corrpt_w = Evtweight;      

      if(ismc && is_weight && isptw)
	{
	  corrpt_w = Evtweight*(1./fpt->Eval(pT1));
	}

      if((pT1 > pTmin && pT1 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  hptdist->Fill(pT1, corrpt_w);

	  double xaxis[xndim] = {1, (double)ptbin, (double)refparB, (double)centBin};
	  sNjet_noW->Fill(xaxis);
	  sNjet_W->Fill(xaxis, corrpt_w);
	}
      
      TVector3 vec2 = Tvec2[i1];
      double pT2 = vec2.Pt();
      double Eta2 = vec2.Eta();
      double Phi2 = vec2.Phi();
      
      if((pT1 > pTmin && pT1 < pTmax) && (pT2 > pTmin && pT2 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  double DEta = Eta1 - Eta2;
	  double DPhi = TVector2::Phi_mpi_pi(Phi1 - Phi2);
	  double DR = TMath::Sqrt(pow(DEta,2) + pow(DPhi,2));
	  
	  double DR_axis[xndim] = {DR, (double)ptbin, (double)refparB, (double)centBin};
	  hDR_noW->Fill(DR_axis);
	  hDR_W->Fill(DR_axis, corrpt_w);
	  
	  double DEta_axis[xndim] = {DEta, (double)ptbin, (double)refparB, (double)centBin};
	  hDEta_noW->Fill(DEta_axis);
	  hDEta_W->Fill(DEta_axis, corrpt_w);
	  
	  double DPhi_axis[xndim] = {DPhi, (double)ptbin, (double)refparB, (double)centBin};
	  hDPhi_noW->Fill(DPhi_axis);
	  hDPhi_W->Fill(DPhi_axis, corrpt_w);
	}
    }
}

void DeltaR_corr_EWTA_W(std::vector<TVector3> Tvec1, std::vector<TVector3> Tvec2, int centBin, std::vector<int> refparton, std::vector<int> refpartonB, double pTmin, double pTmax, double Etamin, double Etamax, double Evtweight, THnSparseD* sNjet_W, THnSparseD* hDR_W, TH1D* hpTBin, TH1D* hptdist, TF1* fpt, bool ismc, bool is_weight, bool isptw)
{
  //std::cout<<"size of :"<<Tvec1.size()<<"  "<<Tvec2.size()<<"  "<<refparton.size()<<std::endl;

  if(Tvec1.size() != Tvec2.size())
    {
      std::cout<<"~~~~~~~~:Two vectors sizes are not same (something wrong), Please check:~~~~~~~~~"<<std::endl;
    }
  if(Tvec1.size() != refparton.size())
    {
      std::cout<<"~~~~~~~~:refparton vec size are not same as Tvec1(something wrong), Please check:~~~~~~~~~"<<std::endl;
    }

  const int xndim = 4;

  for(int i1 = 0; i1 < (int)Tvec1.size(); i1++)
    {
      TVector3 vec1 = Tvec1[i1];
      double pT1 = vec1.Pt();
      double Eta1 = vec1.Eta();
      double Phi1 = vec1.Phi();
      int refpar = refparton[i1];
      int refparB = refpartonB[i1];

      int ptbin = hpTBin->FindBin(pT1) -1;
      
      double corrpt_w = Evtweight;

      if(ismc && is_weight && isptw)
        {
          corrpt_w = Evtweight*(1./fpt->Eval(pT1));
	}

      if((pT1 > pTmin && pT1 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  hptdist->Fill(pT1, corrpt_w);

	  double xaxis[xndim] = {1, (double)ptbin, (double)refparB, (double)centBin};
	  sNjet_W->Fill(xaxis, corrpt_w);
	}
      
      TVector3 vec2 = Tvec2[i1];
      double pT2 = vec2.Pt();
      double Eta2 = vec2.Eta();
      double Phi2 = vec2.Phi();
      
      if((pT1 > pTmin && pT1 < pTmax) && (pT2 > pTmin && pT2 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  double DEta = Eta1 - Eta2;
	  double DPhi = TVector2::Phi_mpi_pi(Phi1 - Phi2);
	  double DR = TMath::Sqrt(pow(DEta,2) + pow(DPhi,2));
	  
	  double DR_axis[xndim] = {DR, (double)ptbin, (double)refparB, (double)centBin};
	  hDR_W->Fill(DR_axis, corrpt_w);
	}
    }
}

void DeltaR_corr_EWTA_WTACorrectionOnGen(std::vector<TVector3> Tvec1, std::vector<TVector3> Tvec2, int centBin, std::vector<int> refparton, std::vector<int> refpartonB, double pTmin, double pTmax, double Etamin, double Etamax, double Evtweight, THnSparseD* sNjet_noW, THnSparseD* sNjet_W, THnSparseD* hDR_noW, THnSparseD* hDR_W, THnSparseD* hDEta_noW, THnSparseD* hDEta_W, THnSparseD* hDPhi_noW, THnSparseD* hDPhi_W, TH1D* hpTBin, TH1D* hptdist, TF1* fpt, bool ismc, bool is_weight, bool isptw, bool is_EWTA_Corr, TH1D* hDR_EWTA_W[NCentbin][jtpT_nbins])
{

  //std::cout<<"size of :"<<Tvec1.size()<<"  "<<Tvec2.size()<<"  "<<refparton.size()<<std::endl;

  if(Tvec1.size() != Tvec2.size())
    {
      std::cout<<"~~~~~~~~:Two vectors sizes are not same (something wrong), Please check:~~~~~~~~~"<<std::endl;
    }
  if(Tvec1.size() != refparton.size())
    {
      std::cout<<"~~~~~~~~:refparton vec size are not same as Tvec1(something wrong), Please check:~~~~~~~~~"<<std::endl;
    }

  const int xndim = 4;

  for(int i1 = 0; i1 < (int)Tvec1.size(); i1++)
    {
      TVector3 vec1 = Tvec1[i1];
      double pT1 = vec1.Pt();
      double Eta1 = vec1.Eta();
      double Phi1 = vec1.Phi();
      int refpar = refparton[i1];
      int refparB = refpartonB[i1];

      int ptbin = hpTBin->FindBin(pT1) -1;

      double corrpt_w = Evtweight;      

      if(ismc && is_weight && isptw)
	{
	  corrpt_w = Evtweight*(1./fpt->Eval(pT1));
	}

      if((pT1 > pTmin && pT1 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  hptdist->Fill(pT1, corrpt_w);

	  double xaxis[xndim] = {1, (double)ptbin, (double)refparB, (double)centBin};
	  sNjet_noW->Fill(xaxis);
	  sNjet_W->Fill(xaxis, corrpt_w);
	}
      
      TVector3 vec2 = Tvec2[i1];
      double pT2 = vec2.Pt();
      double Eta2 = vec2.Eta();
      double Phi2 = vec2.Phi();
      
      if((pT1 > pTmin && pT1 < pTmax) && (pT2 > pTmin && pT2 < pTmax) && (Eta1 > Etamin && Eta1 < Etamax))
	{
	  double DEta = Eta1 - Eta2;
	  double DPhi = TVector2::Phi_mpi_pi(Phi1 - Phi2);
	  double DR = TMath::Sqrt(pow(DEta,2) + pow(DPhi,2));
	  //std::cout<<"DR before: "<<DR<<"  "<<pT1<<"  "<< ptbin<<"  "<<centBin<<std::endl;
	  if(ismc && is_EWTA_Corr)
	    {
	      gRandom->SetSeed(0);
	      gRandom = new TRandom3(0);

	      //DR = DR + (hDR_EWTA_W[centBin][ptbin]->GetRandom());
	      DR = DR*(hDR_EWTA_W[centBin][ptbin]->GetRandom());
	      //while(DR < 0 ) DR = DR + (hDR_EWTA_W[centBin][ptbin]->GetRandom());
	      //std::cout<<"DR after: "<<DR<<"  "<<pT1<<"  "<< ptbin<<"  "<<centBin<<"  "<<hDR_EWTA_W[centBin][ptbin]->GetRandom()<<std::endl;
	      //std::cout<<"DR after: "<<DR<<"  "<<std::endl;
	    }


	  double DR_axis[xndim] = {DR, (double)ptbin, (double)refparB, (double)centBin};
	  hDR_noW->Fill(DR_axis);
	  hDR_W->Fill(DR_axis, corrpt_w);
	  
	  double DEta_axis[xndim] = {DEta, (double)ptbin, (double)refparB, (double)centBin};
	  hDEta_noW->Fill(DEta_axis);
	  hDEta_W->Fill(DEta_axis, corrpt_w);
	  
	  double DPhi_axis[xndim] = {DPhi, (double)ptbin, (double)refparB, (double)centBin};
	  hDPhi_noW->Fill(DPhi_axis);
	  hDPhi_W->Fill(DPhi_axis, corrpt_w);
	}
    }
}

void DeltaR_corr_Gen_Reco(std::vector<TVector3> Tvec1_gn, std::vector<TVector3> Tvec2_gn, std::vector<TVector3> Tvec1_rc, std::vector<TVector3> Tvec2_rc, std::vector<int> refpartonB_gn, std::vector<int> refpartonB_rc, int centBin, double pTmin, double pTmax, double Etamin, double Etamax, double Evtweight, THnSparseD* hDR_gn_rc, THnSparseD* hDR_gn_rc_ratio, THnSparseD* hR_E_gn_rc, THnSparseD* hR_WTA_gn_rc, THnSparseD* hEta_E_gn_rc, THnSparseD* hPhi_E_gn_rc, THnSparseD* hEta_WTA_gn_rc, THnSparseD* hPhi_WTA_gn_rc, TH1D* hpTBin)
{
  if(Tvec1_gn.size() != Tvec2_gn.size())
    {
      std::cout<<"~~~~~~~~:Two Gen vectors sizes are not same (something wrong), Please check:~~~~~~~~~"<<std::endl;
    }
  if(Tvec1_rc.size() != Tvec2_rc.size())
    {
      std::cout<<"~~~~~~~~:Two Reco vectors size are not same (something wrong), Please check:~~~~~~~~~"<<std::endl;
    }

  for(int i1 = 0; i1 < (int)Tvec1_gn.size(); i1++)
    {
      TVector3 vec1_gn = Tvec1_gn[i1];
      double pT1_gn = vec1_gn.Pt();
      double Eta1_gn = vec1_gn.Eta();
      double Phi1_gn = vec1_gn.Phi();
      int refparB_gn = refpartonB_gn[i1];

      int ptbin_gn = hpTBin->FindBin(pT1_gn) -1;

      TVector3 vec2_gn = Tvec2_gn[i1];
      double pT2_gn = vec2_gn.Pt();
      double Eta2_gn = vec2_gn.Eta();
      double Phi2_gn = vec2_gn.Phi();
      
      for(int i1 = 0; i1 < (int)Tvec1_rc.size(); i1++)
	{
	  TVector3 vec1_rc = Tvec1_rc[i1];
	  double pT1_rc = vec1_rc.Pt();
	  double Eta1_rc = vec1_rc.Eta();
	  double Phi1_rc = vec1_rc.Phi();
	  int refparB_rc = refpartonB_rc[i1];
	  
	  int ptbin_rc = hpTBin->FindBin(pT1_rc) -1;

	  TVector3 vec2_rc = Tvec2_rc[i1];
	  double pT2_rc = vec2_rc.Pt();
	  double Eta2_rc = vec2_rc.Eta();
	  double Phi2_rc = vec2_rc.Phi();


	  if((pT1_gn > pTmin && pT1_gn < pTmax) && (pT2_gn > pTmin && pT2_gn < pTmax) && (Eta1_gn > Etamin && Eta1_gn < Etamax))
	    {
	      if((pT1_rc > pTmin && pT1_rc < pTmax) && (pT2_rc > pTmin && pT2_rc < pTmax) && (Eta1_rc > Etamin && Eta1_rc < Etamax))
	      if((Eta1_rc > Etamin && Eta1_rc < Etamax))
		{
		  double DEta_gn = Eta1_gn - Eta2_gn;
		  double DPhi_gn = TVector2::Phi_mpi_pi(Phi1_gn - Phi2_gn);
		  double DR_gn = TMath::Sqrt(pow(DEta_gn,2) + pow(DPhi_gn,2));
		  
		  double DEta_rc = Eta1_rc - Eta2_rc;
		  double DPhi_rc = TVector2::Phi_mpi_pi(Phi1_rc - Phi2_rc);
		  double DR_rc = TMath::Sqrt(pow(DEta_rc,2) + pow(DPhi_rc,2));
		  
		  double DR_axis[5] = {DR_gn, DR_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};

		  double Gen_Reco_Ratio = DR_rc/DR_gn;
		  double DR_axis_Ratio[4] = {Gen_Reco_Ratio, (double)ptbin_gn, (double)refparB_gn, (double)centBin};

		  hDR_gn_rc->Fill(DR_axis, Evtweight);
		  hDR_gn_rc_ratio->Fill(DR_axis_Ratio, Evtweight);

		  double R_E_gn = TMath::Sqrt(pow(Eta1_gn,2)+pow(Phi1_gn,2));
		  double R_WTA_gn = TMath::Sqrt(pow(Eta2_gn,2)+pow(Phi2_gn,2));
		  
		  double R_E_rc = TMath::Sqrt(pow(Eta1_rc,2)+pow(Phi1_rc,2));
		  double R_WTA_rc = TMath::Sqrt(pow(Eta2_rc,2)+pow(Phi2_rc,2));
		  
		  double R_E_axis[5] = {R_E_gn, R_E_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};
		  double R_WTA_axis[5] = {R_WTA_gn, R_WTA_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};
		  hR_E_gn_rc->Fill(R_E_axis, Evtweight);
		  hR_WTA_gn_rc->Fill(R_WTA_axis, Evtweight);

		  double Eta_E_axis[5] = {Eta1_gn, Eta1_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};
		  double Phi_E_axis[5] = {Phi1_gn, Phi1_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};

		  double Eta_WTA_axis[5] = {Eta2_gn, Eta2_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};
		  double Phi_WTA_axis[5] = {Phi2_gn, Phi2_rc, (double)ptbin_gn, (double)refparB_gn, (double)centBin};
		  
		  hEta_E_gn_rc->Fill(Eta_E_axis, Evtweight);
		  hPhi_E_gn_rc->Fill(Phi_E_axis, Evtweight);

		  hEta_WTA_gn_rc->Fill(Eta_WTA_axis, Evtweight);
		  hPhi_WTA_gn_rc->Fill(Phi_WTA_axis, Evtweight);
		}
	    }
	}
    }
}
