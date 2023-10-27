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
