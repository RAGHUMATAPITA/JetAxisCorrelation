#include "call_libraries.h"  // call libraries from ROOT and C++

//vectors for checking
std::vector<double> FilteredJet_jtEta;
std::vector<double> FilteredJet_jtPhi;
std::vector<double> FilteredJet_WTAEta;
std::vector<double> FilteredJet_WTAPhi;

// reco jets vector
std::vector<TVector3> FilteredJet_CorrpT_Vec;
std::vector<TVector3> FilteredJet_RawpT_Vec;
std::vector<TVector3> FilteredWTAJet_CorrpT_Vec;
std::vector<TVector3> FilteredWTAJet_RawpT_Vec;
std::vector<int> refparton_Vec;
std::vector<int> refpartonB_Vec;

//reco matched jets vector
std::vector<TVector3> Matched_FilteredJet_CorrpT_Vec;
std::vector<TVector3> Matched_FilteredWTAJet_CorrpT_Vec;
std::vector<int> Matched_refparton_Vec;
std::vector<int> Matched_refpartonB_Vec;
std::vector<int> Matched_refparton_PM_Vec;
std::vector<int> Matched_refpartonB_PM_Vec;

//gen jets vector
std::vector<TVector3> FilteredJet_GenpT_Vec;
std::vector<TVector3> FilteredWTAJet_GenpT_Vec;
std::vector<int> Genrefparton_Vec;
std::vector<int> GenrefpartonB_Vec;

//gen matched jets vector
std::vector<TVector3> Matched_FilteredJet_GenpT_Vec;
std::vector<TVector3> Matched_FilteredWTAJet_GenpT_Vec;
std::vector<int> Matched_Genrefparton_Vec;
std::vector<int> Matched_GenrefpartonB_Vec;
std::vector<int> Matched_Genrefparton_PM_Vec;
std::vector<int> Matched_GenrefpartonB_PM_Vec;
