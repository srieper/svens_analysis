#ifndef MDSM_LFV_HMuE_CL_H
#define MDSM_LFV_HMuE_CL_H

//#ifndef MDSM_MC_CL_H
//#include "MDSM_Mc_CL.h"
//#endif
#include "MDSM_Structs.h"
#include "MDSM_Mc_CL.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include <vector>
#include <iostream>
#include <math.h>
#include "TCanvas.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "MDSM_master_HCL.h"
#include <memory>
#include <map>
#include <string>
#include <algorithm>
#include "TTree.h"
#include "tau_decay.h"
#include "Sven.h"
#include "../../../HHKinFit/interface/HHKinFitSingleHMaster.h"

using std::string;

class MDSM_LFV_HMuE_CL : public TObject{
 public:
  MDSM_LFV_HMuE_CL();
  //~MDSM_LFV_HMuE_CL();

  //functions to be called
  void analyze(MDSM_master_HCL* event,MDSM_Mc_CL* mcs);
  void write();

  void define_histos();
  void add_histos_to_list();

  bool blind;

  //some constants
  double pi;
  double mass_higgs;
  double weight;
 
  TObjArray Hlist;
 
  //outfile
  TFile* outfile;

 private:
  MDSM_master_HCL* event;

  //some standard histos
  std::map<std::string ,TH2F*> h2_Easymm;
  std::map<std::string ,TH2F*> h2_z_zbar;
  std::map<std::string ,TH1F*> h_Easymm;
  std::map<std::string ,TH1F*> h_;
  std::map<std::string ,TH1F*> h_bg_rejec;
  std::map<std::string, TTree*> trees;

  const std::map<std::string, Sven::DECAY_MODES> decay_modes = {{"rho", Sven::DECAY_MODES::RHO}, {"e", Sven::DECAY_MODES::E}, {"mu", Sven::DECAY_MODES::MU}, {"pi", Sven::DECAY_MODES::PI}, {"a1", Sven::DECAY_MODES::A1}, {"had", Sven::DECAY_MODES::HAD}};
  const std::map<int,int> bmass = {{25, 125}, {23,91}};

  HHKinFitSingleHMaster *kinFit, *kinFit_gen;

  TLorentzVector (HHKinFitSingleHMaster::*fit_part)(Int_t) const;
  TLorentzVector (HHKinFitSingleHMaster::*fit_antipart)(Int_t) const;
  TLorentzVector (HHKinFitSingleHMaster::*GetFit[2])(Int_t) const {0,0};

  void add_map_to_list(std::map<std::string,TH2F*>);
  void add_map_to_list(std::map<std::string,TH1F*>);

  void norm_histos();

  void fill_histos_decay_mode(const std::vector<tau_decay>& gen_taus, const ROOT::Math::XYZTVector* boson_vec, Int_t boson_int=0, const HHKinFitSingleHMaster* fitted_tau=0, const HHKinFitSingleHMaster* fitted_tau_gen=0);
  template<class V1, class V2, class V3, class V4>
  void fill_histos_calculate_values(const V1& vis1, const V2& vis2, const V3& tau1, const V4& tau2, std::map<std::string,double>& values, std::map<std::string,bool>& bool_values, bool vis_only, const std::string& name1, const std::string& name2, const std::string& addition, bool keep_particle_order);
  void fill_histos_match(const std::string& name, const std::map<std::string,double>& values, const std::map<std::string,bool>& bool_values, bool Easymm=true);
  void fill_histos_cut(const std::string& name, const std::map<std::string,double>& values, const std::map<std::string,bool>& bool_values, bool Fill_Easymm);
  void fill_histos_call_create_and_fill(const std::string& name, const std::map<std::string,double>& values, bool Fill_Easymm);
  void fill_bg_rejec_histos(const std::string& name, const std::map<std::string,double>& values);

  void pt_lowest_highest_histos(const std::vector<tau_decay>&);

  struct branch_data {
    Double_t weight, chi2_asymm, chi2H, chi2Z, convergenceH, convergenceZ, pt_tau_low, pt_tau_high;
    Int_t n_reco_part_n, n_reco_part_p, best_fit_mass, n_jets, boson_int;
    Bool_t match, reco, kin_fit;
  };
};

namespace Sven {
  const TH_data z_zbar("h2_z_zbar_", "x=E_{A^{-}}/E_{#tau^{-}}, y=E_{B^{+}}/E_{#tau^{+}}", 40, 0, 2);
  const TH_data Easymm("h2_Easymm_", "(E_{#tau^{-}} - E_{#tau^{+}})/(E_{#tau^{-}} + E_{#tau^{-}})", 24, 1.2);
  const TH_data bg_rejec("h_event_rejection_", "event rejection", 41, 0, 2.05);
  const TH_data zzbar_sum("h_", "z+zbar", 41, 0, 2.05);
  const TH_data th_data_pt_vis("h_", "pt vis", 200, 0, 200);
}

#endif
