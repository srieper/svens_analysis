#include "MDSM_LFV_HMuE_CL.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <memory>
#include "TClonesArray.h"

//#define DEBUG_H_ditau

using std::cout;
using std::endl;
using std::cerr;
using std::to_string;
using std::pair;
using std::map;
using std::vector;
using std::string;
using std::out_of_range;
using Sven::asymm;
using Sven::create_and_fill_histos;
using Sven::createTH1F;
using Sven::createTH2F;
using Sven::to_TLorentzVector;
using Sven::DECAY_MODES;
using Sven::operator<<;
using Sven::th_data_pt_vis;
using Sven::z_zbar;
using Sven::Easymm;
using Sven::bg_rejec;
using Sven::zzbar_sum;

//constructor
MDSM_LFV_HMuE_CL::MDSM_LFV_HMuE_CL() :outfile(0), event(0), kinFit(0), kinFit_gen(0), fit_part(0), fit_antipart(0) {
  blind=false;
  pi=3.14159;
  mass_higgs=125; //[GeV]

  Hlist=  TObjArray(0);

  //define histo in different function
  define_histos();
  //add_histos_to_list();
}

//MDSM_LFV_HMuE_CL::~MDSM_LFV_HMuE_CL() {
//}


void MDSM_LFV_HMuE_CL::write() {
  static bool histos_added_to_list=false;
  if(!histos_added_to_list) {
    add_histos_to_list();
    histos_added_to_list=true;
  }

  outfile->cd();
  outfile->mkdir("LFV_HMuE_histos");
  outfile->cd("LFV_HMuE_histos");
  norm_histos();
  Hlist.Write();
  outfile->cd();

  outfile->mkdir("H_ditau_trees");
  outfile->cd("H_ditau_trees");
  for(map<string,TTree*>::const_iterator iter=trees.begin(); iter!=trees.end(); ++iter) {
    iter->second->Write();
  }

  outfile->cd();

  return;
}



//do analysis
void MDSM_LFV_HMuE_CL::analyze(MDSM_master_HCL* akt_event,MDSM_Mc_CL* mcs){
#ifdef DEBUG_H_ditau
cout << "analyze called\n";
#endif
  event = akt_event;
  weight=event->weight;

  if(kinFit) {
    delete kinFit;
    kinFit=0;
  }
  if(kinFit_gen) {
    delete kinFit_gen;
    kinFit_gen=0;
  }

  //MC-tau's + daughters in gen_taus einlesen
  vector<tau_decay> gen_taus;
  const ROOT::Math::XYZTVector* boson_vec=0;
  int boson_int=0;
  for (int i=0; i<mcs->get_vector_size(); i++) {
    mcs->set_akt_values(i);

    if(mcs->status.akt_val!=3) continue;
    if(fabs(mcs->pdgid.akt_val)==15){
      gen_taus.push_back(tau_decay(mcs, i, event));
    }
    else if(mcs->pdgid.akt_val==25 || mcs->pdgid.akt_val==23) {
      if(boson_int) cout << "Zweites Higgs oder Z gefunden\n";
      boson_vec = mcs->get_vec_val(i);
      boson_int = mcs->pdgid.akt_val;
    }
  }

  //Zuordnung von tau_gen zu tau_reco
  if(event->t_tau.size()==2 && gen_taus.size()==2 && (gen_taus[0].is_particle()==event->t_tau_charge[0] || gen_taus[1].is_particle()==event->t_tau_charge[1]))
    iter_swap(gen_taus.begin(), gen_taus.begin()+1);
  if(gen_taus.size()==2) {
    fit_part     = gen_taus[0].is_particle()?&HHKinFitSingleHMaster::getTau1Fitted:&HHKinFitSingleHMaster::getTau2Fitted;
    fit_antipart = gen_taus[1].is_particle()?&HHKinFitSingleHMaster::getTau1Fitted:&HHKinFitSingleHMaster::getTau2Fitted;
    GetFit[0] = &HHKinFitSingleHMaster::getTau1Fitted;
    GetFit[1] = &HHKinFitSingleHMaster::getTau2Fitted;
  }

  //Fit
  if(gen_taus.size()==2 && gen_taus[0].reco_vec() && gen_taus[1].reco_vec()) {

    TLorentzVector met_vec      = to_TLorentzVector(akt_event->MVA_MET_TT);
    TLorentzVector tau1_vis     = to_TLorentzVector(*gen_taus[0].reco_vec());
    TLorentzVector tau2_vis     = to_TLorentzVector(*gen_taus[1].reco_vec());
    TLorentzVector tau1_gen_vis = to_TLorentzVector(*gen_taus[0].vis_vec());
    TLorentzVector tau2_gen_vis = to_TLorentzVector(*gen_taus[1].vis_vec());

    kinFit     = new HHKinFitSingleHMaster(&tau1_vis, &tau2_vis);
    kinFit_gen = new HHKinFitSingleHMaster(&tau1_gen_vis, &tau2_gen_vis);

    TMatrixD met_cov = TMatrixD(2,2);
    met_cov[0][0]    = akt_event->MVA_MET_TT_cov_0_0;
    met_cov[1][0]    = akt_event->MVA_MET_TT_cov_1_0;
    met_cov[0][1]    = akt_event->MVA_MET_TT_cov_0_1;
    met_cov[1][1]    = akt_event->MVA_MET_TT_cov_1_1;

    std::vector<Int_t> hypo_mh = {bmass.at(25), bmass.at(23)};//, 86, 96, 105, 120, 130}; //not more than 20 entries!

    kinFit->setAdvancedBalance(&met_vec, met_cov);
    kinFit->addMhHypothesis(hypo_mh);
    kinFit->doFullFit();

    try {
      Double_t fit_prob = kinFit->getFitProb(bmass.at(boson_int));
      Double_t chi2     = kinFit->getChi2(bmass.at(boson_int));
      h_["fit_prob"]->Fill(fit_prob);
      h_["chi2_dist"]->Fill(chi2);

      fit_prob = kinFit->getFitProb(kinFit->getBestHypoFullFit());
      chi2     = kinFit->getBestChi2FullFit();
      h_["fit_prob_best"]->Fill(fit_prob);
      h_["chi2_dist_best"]->Fill(chi2);
      //if(int bh = kinFit->getBestHypoFullFit() < 0)   //??? test/debug?
        //cout << "BestHypoFullFit(): " << bh << endl;

      h_["chi2_asymm"]->Fill(Sven::asymm(kinFit->getChi2(bmass.at(25)), kinFit->getChi2(bmass.at(23))));
      //for more than 2 masses
      //int bh_hz = kinFit->getChi2(bmass.at(25))<kinFit->getChi2(bmass.at(23))?bmass.at(25):bmass.at(23);
      //h2_Easymm["bestmass_hzbest"]->Fill(bh_hz, kinFit->getBestHypoFullFit());

      h_["best_fit_mass"]->Fill(kinFit->getBestHypoFullFit());
      h_["convergence_higgs"]->Fill(kinFit->GetConvergence(bmass.at(25)));
      h_["convergence_z"]->Fill(kinFit->GetConvergence(bmass.at(23)));
      h2_Easymm["best_fit_mass_ptH"]->Fill(kinFit->getBestHypoFullFit(), boson_vec->pt());
      Double_t dphi = Sven::DeltaPhi(*gen_taus[0].reco_vec(), *gen_taus[1].reco_vec());
      h2_Easymm["best_fit_mass_phi"]->Fill(kinFit->getBestHypoFullFit(), fabs(dphi));
      h2_Easymm["chi2H_chi2Z"]->Fill(kinFit->getChi2(bmass.at(25)), kinFit->getChi2(bmass.at(23)));
    } catch(out_of_range) {}

    kinFit_gen->setAdvancedBalance(&met_vec, met_cov);
    kinFit_gen->addMhHypothesis(hypo_mh);
    kinFit_gen->doFullFit();

    try {
      Double_t fit_prob = kinFit_gen->getFitProb(bmass.at(boson_int));
      Double_t chi2     = kinFit_gen->getChi2(bmass.at(boson_int));
      h_["chi2_dist_genfit"]->Fill(chi2);
      h_["fit_prob_genfit"]->Fill(fit_prob);
    } catch(out_of_range) {}
  }
#ifdef DEBUG_H_ditau
cout << "after fit\n";
#endif

  //Histos
  if(gen_taus.size()==2)
  {
    try {
      fill_histos_decay_mode(gen_taus, boson_vec, boson_int, kinFit, kinFit_gen);
    } catch(ErrorSven::missing_fit) {
#ifdef DEBUG_H_ditau
cout << "missing fit caught in analyze\n";
#endif
    } catch(out_of_range oor) {
      cerr << "Out of range in 'analyse'-1: " << oor.what() << "\n";
    } catch(ErrorSven::somethingwenthorriblywrong err) {
      cerr << "Something went terribly wrong (in function \"" << err.func << "\")\n";
    }
#ifdef DEBUG_H_ditau
cout << "after majority of histos\n";
#endif

    pt_lowest_highest_histos(gen_taus);

    //Betrachte nur Zerfaelle in Pionen
    if(gen_taus[0].n_K() || gen_taus[1].n_K()) {
#ifdef DEBUG_H_ditau
cout << "Kaon end analyze\n";
#endif
      h_["K"]->Fill(1);
      return;
    }
    if(gen_taus[0].n_other() || gen_taus[1].n_other()) {
      cout << "unbekanntes Teilchen gefunden\n";
      return;
    }
    if(gen_taus[0].n_mu() || gen_taus[1].n_mu()) 
#ifdef DEBUG_H_ditau
cout << "muon found, end analyze\n";
#endif
      return;

#ifdef DEBUG_H_ditau
cout << "bofore majority of individual histos\n";
#endif
    int part = int(gen_taus[1].is_particle());
    int antipart = int(gen_taus[0].is_particle());
    if(part == antipart) cout << "part == antipart == " << part << " => " << (part?("2 tau-"):("2 tau+")) << endl;
    double z_gen = gen_taus[part].vis_vec()->E()/gen_taus[part].mother_vec()->E();
    double zbar_gen = gen_taus[antipart].vis_vec()->E()/gen_taus[antipart].mother_vec()->E();

    try {
      //2x tau->rho
      //switch(decay_mode.at(gen_taus[part].decay_mode()) + decay_mode.at(gan_taus[antipart].decay_mode())) ...
      if(gen_taus[part].decay_mode()=="rho" && gen_taus[antipart].decay_mode()=="rho") {
        double Easymm_rhop_gen = fabs(asymm(gen_taus[antipart].get_particle_vec(111).E(),gen_taus[antipart].get_particle_vec(211).E()));
        double Easymm_rhon_gen = fabs(asymm(gen_taus[part].get_particle_vec(111).E(),gen_taus[part].get_particle_vec(-211).E()));
        h2_Easymm["asymm_2rho_gen_all"]->Fill(Easymm_rhon_gen, Easymm_rhop_gen);
        if(z_gen>.4 && zbar_gen>.4)
          h2_Easymm["asymm_2rho_gen_all_z>.4"]->Fill(Easymm_rhon_gen, Easymm_rhop_gen);
      }
      //2x tau->a1
      else if(gen_taus[part].decay_mode()=="a1" && gen_taus[antipart].decay_mode()=="a1") {
        if(!gen_taus[part].n_pi0() && !gen_taus[antipart].n_pi0()) {
          try {
            pair<double,double> taun_fracs = gen_taus[part].get_Efrac_a1to3pi();
            pair<double,double> taup_fracs = gen_taus[antipart].get_Efrac_a1to3pi();
            h2_Easymm["asymm_2a1to3pi_gen_all"]->Fill(fabs(Sven::asymm(taun_fracs.first, taun_fracs.second)), fabs(Sven::asymm(taup_fracs.first, taup_fracs.second)));
          } catch(tau_decay::wrong_particle) {
            cerr << "Falsches Teilchen: kein a1->3pi!" << endl;
          }
        }
        try {
          pair<double,double> taun_fracs = gen_taus[part].get_Efrac_a1();
          pair<double,double> taup_fracs = gen_taus[antipart].get_Efrac_a1();
          h2_Easymm["asymm_2a1_gen_all"]->Fill(fabs(Sven::asymm(taun_fracs.first, taun_fracs.second)), fabs(Sven::asymm(taup_fracs.first, taup_fracs.second)));
        } catch(tau_decay::wrong_particle) {
          cerr << "Falsches Teilchen: kein a1!" << endl;
        }
      }
    } catch(tau_decay::unknown_decay_mode) {
    }
  }
  else {
    cout << "gen_taus.size() == " << gen_taus.size() << " != 2" << endl;
    return;
  }


  if(event->t_tau.size()==2 && gen_taus[0].reco_vec() && gen_taus[1].reco_vec()){
  

    for(int i=0; i<2; i++) {
      h_["deltaR_gen"]->Fill(gen_taus[i].get_deltaR_gen());
      h_["deltaR_reco"]->Fill(gen_taus[i].get_deltaR());
    }

    //Nr von Teilchen und Antiteilchen
    int part = int(gen_taus[1].is_particle());
    int antipart = int(gen_taus[0].is_particle());

    //matching
    if(gen_taus[0].get_deltaR()<.1 && gen_taus[1].get_deltaR()<.1 && kinFit) {
      try {
        h_["fit_tau1_res"]->Fill(Sven::relValue2(gen_taus[part].mother_vec()->E(), (kinFit->*fit_part)(bmass.at(boson_int)).E()).first);
        h_["fit_tau2_res"]->Fill(Sven::relValue2(gen_taus[antipart].mother_vec()->E(), (kinFit->*fit_antipart)(bmass.at(boson_int)).E()).first);
        h_["fit_tau1_res_best"]->Fill(Sven::relValue2(gen_taus[part].mother_vec()->E(), (kinFit->*fit_part)(-1).E()).first);
        h_["fit_tau2_res_best"]->Fill(Sven::relValue2(gen_taus[antipart].mother_vec()->E(), (kinFit->*fit_antipart)(-1).E()).first);
      } catch(out_of_range oor) {//why? bmass?
      }

      double z_gen = gen_taus[part].vis_vec()->E()/gen_taus[part].mother_vec()->E();
      double zbar_gen = gen_taus[antipart].vis_vec()->E()/gen_taus[antipart].mother_vec()->E();

     try {
      //2x tau->rho
      if(gen_taus[part].decay_mode()=="rho" && gen_taus[antipart].decay_mode()=="rho") {
        double Easymm_rhop_gen = fabs(asymm(gen_taus[antipart].get_particle_vec(111).E(),gen_taus[antipart].get_particle_vec(211).E()));
        double Easymm_rhon_gen = fabs(asymm(gen_taus[part].get_particle_vec(111).E(),gen_taus[part].get_particle_vec(-211).E()));
        h2_Easymm["asymm_2rho_gen_match"]->Fill(Easymm_rhon_gen, Easymm_rhop_gen);
        if(z_gen>.4 && zbar_gen>.4)
          h2_Easymm["asymm_2rho_gen_match_z>.4"]->Fill(Easymm_rhon_gen, Easymm_rhop_gen);
      }
      //2x tau->a1
      else if(gen_taus[part].decay_mode()=="a1" && gen_taus[antipart].decay_mode()=="a1") {
        if(!gen_taus[part].n_pi0() && !gen_taus[antipart].n_pi0()) {
          try {
            pair<double,double> taun_fracs = gen_taus[part].get_Efrac_a1to3pi();
            pair<double,double> taup_fracs = gen_taus[antipart].get_Efrac_a1to3pi();
            h2_Easymm["asymm_2a1to3pi_gen_match"]->Fill(fabs(Sven::asymm(taun_fracs.first, taun_fracs.second)), fabs(Sven::asymm(taup_fracs.first, taup_fracs.second)));
          } catch(tau_decay::wrong_particle) {
            cerr << "Falsches Teilchen: kein a1->3pi!" << endl;
          }
        }
        try {
          pair<double,double> taun_fracs = gen_taus[part].get_Efrac_a1();
          pair<double,double> taup_fracs = gen_taus[antipart].get_Efrac_a1();
          h2_Easymm["asymm_2a1_gen_match"]->Fill(fabs(Sven::asymm(taun_fracs.first, taun_fracs.second)), fabs(Sven::asymm(taup_fracs.first, taup_fracs.second)));
        } catch(tau_decay::wrong_particle) {
          cerr << "Falsches Teilchen: kein a1!" << endl;
        }
      }
     } catch(tau_decay::unknown_decay_mode) {
     }
    }
  }
#ifdef DEBUG_H_ditau
cout << "end of analyze\n";
#endif
  
 return;
}





void MDSM_LFV_HMuE_CL::define_histos(){

  using Sven::TH_data;
 
  ////////////////////////////
  //now actual histos

  {
    TH_data Easymmasymm("h2_Easymm", "Energieasymmetrien der Zerfallsprodukte gegeneinander", 24, 0, 1.2);
    vector<string> keys2_aa = {"asymm_2rho_gen_all", "asymm_2rho_gen_all_z>.4", "asymm_2rho_gen_match", "asymm_2rho_gen_match_z>.4", "asymm_2a1_gen_all", "asymm_2a1_gen_match", "asymm_2a1to3pi_gen_all", "asymm_2a1to3pi_gen_match"};
    for(vector<string>::const_iterator iter=keys2_aa.begin(); iter!=keys2_aa.end(); iter++) {
      createTH2F(&h2_Easymm[*iter], Easymmasymm, *iter);
    }
  }


  h_["K"] = new TH1F("h_K", "K", 3, 0, 2);
  h_["deltaR_gen"] = new TH1F("h_deltaR_gen", "deltaR_gen", 56, 0, 3.5);
  h_["deltaR_reco"] = new TH1F("h_deltaR_reco", "deltaR_reco", 56, 0, 3.5);
  h_["fit_tau1_res"] = new TH1F("h_hit_tau1_res", "tau1 resolution (after fit)", 24, -1.2, 1.2);
  h_["fit_tau2_res"] = new TH1F("h_hit_tau2_res", "tau2 resolution (after fit)", 24, -1.2, 1.2);
  h_["fit_tau1_res_best"] = new TH1F("h_hit_tau1_res_best", "tau1 resolution (after fit, best fit)", 24, -1.2, 1.2);
  h_["fit_tau2_res_best"] = new TH1F("h_hit_tau2_res_best", "tau2 resolution (after fit, best fit)", 24, -1.2, 1.2);
  h_["fit_prob"] = new TH1F("h_fit_prob", "fit probability distribution", 42,0, 1.05);
  h_["fit_prob_best"] = new TH1F("h_fit_prob_best", "fit probability distribution (best fit)", 42,0, 1.05);
  h_["fit_prob_genfit"] = new TH1F("h_fit_prob_genfit", "fit probability distribution", 42,0, 1.05);
  h_["chi2_dist"] = new TH1F("h_chi2_dist", "#chi^{2} distribution", 30, 0, 30);
  h_["chi2_dist_best"] = new TH1F("h_chi2_dist_best", "#chi^{2} distribution (best fit)", 30, 0, 30);
  h_["chi2_dist_genfit"] = new TH1F("h_chi2_dist_genfit", "#chi^{2} distribution", 30, 0, 30);
  h_["best_fit_mass"] = new TH1F("h_best_fit_mass", "Massenhypothese des besten fits", 55, 80, 135);
  h_["chi2_asymm"] = new TH1F("chi2_asymm", "#chi^{2}_{asym}", 22, -1.1, 1.1);
  h_["convergence_higgs"] = new TH1F("convergence_higgs", "", 12, -2, 10);
  h_["convergence_z"] = new TH1F("convergence_z", "", 12, -2, 10);
  h2_Easymm["best_fit_mass_ptH"] = new TH2F("h2_best_fit_mass_ptH", "Massenhypothese des besten fits gegen pt_{Higgs}", 10, 80, 130, 600, 0, 300);
  h2_Easymm["best_fit_mass_phi"] = new TH2F("h2_best_fit_mass_phi", "Massenhypothese des besten fits gegen #phi", 10, 80, 130, 20, 0, 4);
  h2_Easymm["chi2H_chi2Z"] = new TH2F("h2_chi2H_chi2Z", "#chi^{2}_{Higgs} gegen #chi^{2}_{Z}", 100, 0, 50, 100, 0, 50);
  h2_Easymm["bestmass_hzbest"] = new TH2F("bestmass_hzbest", "Mass best fit vs. Mass best fit Higgs/Z", 50, 85, 135, 50, 85, 135);

  return;
}



void MDSM_LFV_HMuE_CL::add_histos_to_list(){
  ///how to fill histos
 
  add_map_to_list(h2_Easymm);
  add_map_to_list(h2_z_zbar);
  add_map_to_list(h_Easymm);
  add_map_to_list(h_bg_rejec);
  add_map_to_list(h_);

  //setdirectory to 0for all histos
  for(int n=0; n< Hlist.GetEntries(); n++){
    TH1F *hd=(TH1F*)(Hlist[n]);
    hd->SetDirectory(0);
    //hd->Sumw2();
  }

  return;
}

//Fill histos
void MDSM_LFV_HMuE_CL::fill_histos_decay_mode(const std::vector<tau_decay>& gen_taus, const ROOT::Math::XYZTVector* boson_vec, Int_t boson_int, const HHKinFitSingleHMaster* fitted_tau, const HHKinFitSingleHMaster* fitted_tau_gen) {
#ifdef DEBUG_H_ditau
cout << "fill_histos_decay_mode called\n";
#endif

  int part = int(gen_taus.at(1).is_particle());
  int antipart = int(gen_taus.at(0).is_particle());
  if(part == antipart) cout << "part == antipart == " << part << " => " << (part?("2 tau-"):("2 tau+")) << endl;
  bool known_decaymode=true;

  string decay_mode_taun, decay_mode_taup;
  try {
    decay_mode_taun = gen_taus.at(part).decay_mode();
    decay_mode_taup = gen_taus.at(antipart).decay_mode();
  }
  catch(tau_decay::unknown_decay_mode) {
    known_decaymode=false;
      decay_mode_taun = decay_mode_taup = "unknown";
  }
  if(decay_mode_taun=="had" || decay_mode_taup=="had")
    known_decaymode=false;
  string name = decay_mode_taun + string("m_") + decay_mode_taup  + string("p");
#ifdef DEBUG_H_ditau
cout << "decay mode: " << name << endl;
#endif
  map<string,double> values;
  map<string,bool> bool_values;
  values["weight"] = event->weight;
  values["chi2_asymm"] = fitted_tau?Sven::asymm(fitted_tau->getChi2(bmass.at(25)), fitted_tau->getChi2(bmass.at(23))):1000;
  bool_values["match"] = (gen_taus.at(0).reco_vec() && gen_taus.at(1).reco_vec() && gen_taus.at(0).get_deltaR()<.1 && gen_taus.at(1).get_deltaR()<.1);

  if(known_decaymode) {
    string name2 = "";
    bool name2_keep_particle_order=false;
    if(decay_modes.at(decay_mode_taun)<decay_modes.at(decay_mode_taup)) {
      name2 = decay_mode_taun + string("_") + decay_mode_taup;
      name2_keep_particle_order=true;
    }
    else if(decay_modes.at(decay_mode_taun)>decay_modes.at(decay_mode_taup)) {
      name2 = decay_mode_taup + string("_") + decay_mode_taun;
      name2_keep_particle_order=false;
    }

    bool_values["cut"] = (boson_vec->pt()<20 && fabs(gen_taus.at(0).mother_vec()->Eta())<.5 && fabs(gen_taus.at(1).mother_vec()->Eta())<.5);
    bool_values["1 jet"] = event->n_jets==1;
    bool_values["boosted"] = boson_vec->pt()>100;
    bool_values["highly boosted"] = boson_vec->pt()>170;
    bool_values["full_had"] = (gen_taus.at(part).is_hadronic() && gen_taus.at(antipart).is_hadronic()); //!true if !full_pion
    bool_values["full_pion"] = bool_values["full_had"];
    for(int i=0; i<2 && bool_values.at("full_pion"); ++i) {
      try {
        switch(decay_modes.at(gen_taus.at(i).decay_mode())) {
          case DECAY_MODES::RHO:
          case DECAY_MODES::PI:
          case DECAY_MODES::A1:
            break;
          default:
            bool_values["full_pion"] = false;
            break;
        }
      } catch(tau_decay::unknown_decay_mode) {//should never happen, only known-_known+ here
        bool_values["full_pion"] = false;
      } catch(out_of_range oor) {
        cerr << "Out of range 'fill_histos_decay_mode'-0\n";
      }
    }

    //use generator
    try {
      fill_histos_calculate_values(*gen_taus.at(part).vis_vec(), *gen_taus.at(antipart).vis_vec(), *gen_taus.at(part).mother_vec(), *gen_taus.at(antipart).mother_vec(), values, bool_values, true, name, name2, "_gen", name2_keep_particle_order);
    } catch(out_of_range oor) {
      cerr << "Out of range 'fill_histos_decay_mode'-1\n";
    }

    //use generator + KinFit
    //if(fitted_tau_gen && fitted_tau_gen->size()>0 && fitted_tau_gen->at(0).size()==2)
      //fill_histos_calculate_values(*gen_taus.at(part).vis_vec(), *gen_taus.at(antipart).vis_vec(), fitted_tau_gen->at(0)[part], fitted_tau_gen->at(0)[antipart], values, bool_values, false, name, name2, "_gen_kin", name2_keep_particle_order);

    //use generator + KinFit(best) -- fix before use
    //if(fitted_tau_gen_best && fitted_tau_gen_best->size()==2)
      //fill_histos_calculate_values(*gen_taus.at(part).vis_vec(), *gen_taus.at(antipart).vis_vec(), (*fitted_tau_gen_best)[part], (*fitted_tau_gen_best)[antipart], values, bool_values, false, name, name2, "_gen_kin_best", name2_keep_particle_order);

    if(gen_taus.at(part).reco_vec() && gen_taus.at(antipart).reco_vec()) {
      bool_values["boosted"] = (*gen_taus.at(part).reco_vec() + *gen_taus.at(antipart).reco_vec() + event->MVA_MET_TT).pt()>100;
      bool_values["highly boosted"] = (*gen_taus.at(part).reco_vec() + *gen_taus.at(antipart).reco_vec() + event->MVA_MET_TT).pt()>170;

      //use reco
      if(bool_values.at("match"))
        fill_histos_calculate_values(*gen_taus.at(part).reco_vec(), *gen_taus.at(antipart).reco_vec(), *gen_taus.at(part).mother_vec(), *gen_taus.at(antipart).mother_vec(), values, bool_values, true, name, name2, "_reco_match", name2_keep_particle_order);

      //use reco + KinFit
      if(fitted_tau) {
        if(fitted_tau && boson_vec && gen_taus.at(part).reco_vec() && gen_taus.at(antipart).reco_vec() && (fitted_tau->GetConvergence(bmass.at(boson_int))==1 || fitted_tau->GetConvergence(bmass.at(boson_int))==2)) //why boson_vec?
          fill_histos_calculate_values(*gen_taus.at(part).reco_vec(), *gen_taus.at(antipart).reco_vec(), (fitted_tau->*fit_part)(bmass.at(boson_int)), (fitted_tau->*fit_antipart)(bmass.at(boson_int)), values, bool_values, false, name, name2, "_reco_kin", name2_keep_particle_order);

        //use reco + KinFit(best)
        if(fitted_tau && gen_taus.at(part).reco_vec() && gen_taus.at(antipart).reco_vec() && (fitted_tau->GetConvergence(-1)==1 || fitted_tau->GetConvergence(-1)==2)) {
          try {
            fill_histos_calculate_values(*gen_taus.at(part).reco_vec(), *gen_taus.at(antipart).reco_vec(), (fitted_tau->*fit_part)(-1), (fitted_tau->*fit_antipart)(-1), values, bool_values, false, name, name2, "_reco_kin_best", name2_keep_particle_order);
            fill_histos_calculate_values(*gen_taus.at(part).reco_vec(), *gen_taus.at(antipart).reco_vec(), (fitted_tau->*fit_part)(-1), (fitted_tau->*fit_antipart)(-1), values, bool_values, false, name, name2, "_reco_kin_best"+to_string(fitted_tau->getBestHypoFullFit()), name2_keep_particle_order);
            if(values.at("chi2_asymm")<-.7)
              fill_histos_calculate_values(*gen_taus.at(part).reco_vec(), *gen_taus.at(antipart).reco_vec(), (fitted_tau->*fit_part)(-1), (fitted_tau->*fit_antipart)(-1), values, bool_values, false, name, name2, "_reco_kin_best_chi2a_lt.7", name2_keep_particle_order);
          } catch(out_of_range) {
            cerr << "Out of range 'fill_histos_decay_mode'-4\n";
          }
        }
      }
    }
  }
  else if(decay_mode_taun=="unknown" || decay_mode_taup=="unknown")
    name="unknown";
#ifdef DEBUG_H_ditau
cout << "right befor trees are fillde\n";
#endif
  if(fitted_tau && fitted_tau->getBestHypoFullFit()==-1)
    fitted_tau=0;
  if(fitted_tau) name+="_fit";
  else name+="_nofit";
  vector<int> npart = {part, antipart};
  const static ROOT::Math::XYZTVector *p_tau_gen[2], *p_tau_reco[2], *p_tau_vis[2], *p_MET, *p_higgs;
  static TLorentzVector fit_H[2], fit_best[2], fit_Z[2];
  static branch_data data;
  //int Z_mass = boson_int==25?91:125;
  for(int i=0; i<2; ++i) {
    p_tau_gen[i]  = gen_taus.at(npart[i]).mother_vec();
    p_tau_reco[i] = gen_taus.at(npart[i]).reco_vec();
    p_tau_vis[i]  = gen_taus.at(npart[i]).vis_vec();
    //fit_H[i]    = fitted_tau?(fitted_tau->*GetFit[npart[i]])(bmass.at(boson_int)):TLorentzVector(0,0,0,0);
    fit_H[i]    = fitted_tau?(fitted_tau->*GetFit[npart[i]])(bmass.at(25)):TLorentzVector(0,0,0,0);
    fit_best[i] = fitted_tau?(fitted_tau->*GetFit[npart[i]])(-1):TLorentzVector(0,0,0,0);
    //fit_Z[i]    = fitted_tau?(fitted_tau->*GetFit[npart[i]])(Z_mass):TLorentzVector(0,0,0,0);
    fit_Z[i]    = fitted_tau?(fitted_tau->*GetFit[npart[i]])(bmass.at(23)):TLorentzVector(0,0,0,0);
  }
  p_MET = &event->MVA_MET_TT;
  p_higgs = boson_vec;

  try {
    data.weight       = event->weight;
    data.chi2_asymm   = fitted_tau?values.at("chi2_asymm"):1000;
    //data.chi2H        = fitted_tau?fitted_tau->getChi2(bmass.at(boson_int)):-1000;
    //data.chi2Z        = fitted_tau?fitted_tau->getChi2(Z_mass):-1000;
    //data.convergenceH = fitted_tau?fitted_tau->GetConvergence(bmass.at(boson_int)):1000;
    //data.convergenceZ = fitted_tau?fitted_tau->GetConvergence(Z_mass):1000;
    data.chi2H        = fitted_tau?fitted_tau->getChi2(bmass.at(25)):-1000;
    data.chi2Z        = fitted_tau?fitted_tau->getChi2(bmass.at(23)):-1000;
    data.convergenceH = fitted_tau?fitted_tau->GetConvergence(bmass.at(25)):1000;
    data.convergenceZ = fitted_tau?fitted_tau->GetConvergence(bmass.at(23)):1000;
    data.pt_tau_low  = gen_taus.at(0).vis_vec()->pt();
    data.pt_tau_high = gen_taus.at(1).vis_vec()->pt();
    if(data.pt_tau_low>data.pt_tau_high) {
      data.pt_tau_low  = data.pt_tau_high;
      data.pt_tau_high = gen_taus.at(0).vis_vec()->pt();
    }
    data.n_reco_part_n = gen_taus.at(part).n_reco_particles();
    data.n_reco_part_p = gen_taus.at(antipart).n_reco_particles();
    data.best_fit_mass = fitted_tau?fitted_tau->getBestHypoFullFit():0;
    data.n_jets        = event->n_jets;
    data.boson_int     = boson_int;
    data.match         = bool_values.at("match");
    data.reco          = (gen_taus.at(part).reco_vec() && gen_taus.at(antipart).reco_vec());
    data.kin_fit       = fitted_tau;
  } catch(out_of_range) {
    cerr << "Out of range 'fill_histos_decay_mode'-7\n";
  }

  try {
    trees.at(name)->Fill();
  } catch(out_of_range) {
    trees[name] = new TTree(name.c_str(), name.c_str());
    trees[name]->Branch("taun_gen", &p_tau_gen[0]);
    trees[name]->Branch("taup_gen", &p_tau_gen[1]);
    trees[name]->Branch("taun_reco", &p_tau_reco[0]);
    trees[name]->Branch("taup_reco", &p_tau_reco[1]);
    trees[name]->Branch("taun_gen_vis", &p_tau_vis[0]);
    trees[name]->Branch("taup_gen_vis", &p_tau_vis[1]);
    trees[name]->Branch("MET", &p_MET);
    trees[name]->Branch("Boson", &p_higgs);
    trees[name]->Branch("taun_kin_H", &fit_H[0]);
    trees[name]->Branch("taup_kin_H", &fit_H[1]);
    trees[name]->Branch("taun_kin_best", &fit_best[0]);
    trees[name]->Branch("taup_kin_best", &fit_best[1]);
    trees[name]->Branch("taun_kin_Z", &fit_Z[0]);
    trees[name]->Branch("taup_kin_Z", &fit_Z[1]);
    trees[name]->Branch("scalars", &data.weight, "weight/D:chi2_asymm:chi2H:chi2Z:convergenceH:convergenceZ:pt_tau_low:pt_tau_high:n_reco_part_n/I:n_reco_part_p:best_fit_mass:n_jets:boson_int:match/O:reco:kin_fit");
    trees.at(name)->Fill();
  }
#ifdef DEBUG_H_ditau
cout << "end of fill_histos_decay_mode\n";
#endif
}

void MDSM_LFV_HMuE_CL::fill_histos_match(const string& name, const map<string,double>& values, const map<string,bool>& bool_values, bool Fill_Easymm) {
  try {
    fill_histos_cut(name+string("_all"), values, bool_values, Fill_Easymm);
    if(bool_values.at("match"))
      fill_histos_cut(name+string("_match"), values, bool_values, Fill_Easymm);
  } catch(std::out_of_range) {
    throw ErrorSven::somethingwenthorriblywrong("fill_histos_match");
  }
}

template<class V1, class V2, class V3, class V4>
void MDSM_LFV_HMuE_CL::fill_histos_calculate_values(const V1& vis1, const V2& vis2, const V3& tau1, const V4& tau2, std::map<std::string,double>& values, map<string,bool>& bool_values, bool vis_only, const std::string& name1, const std::string& name2, const std::string& addition, bool keep_particle_order) {
  if(!vis_only)
    values["Easymm"] = Sven::asymm(vis1.E(), vis2.E());
  else
    values["Easymm"] = 1000;
  try {
    values["z"]       = vis1.E()/tau1.E();
    values["zbar"]    = vis2.E()/tau2.E();
    values["z_max"]   = Double_t(10)/9+1e-9;
    values["zasymm"]  = Sven::asymm(values.at("z"), values.at("zbar"));
    values["pt_low"]  = vis1.pt()<vis2.pt()?vis1.pt():vis2.pt();
    values["pt_high"] = vis1.pt()<vis2.pt()?vis2.pt():vis1.pt();
    bool_values["eta23"]         = (fabs(vis1.Eta())<2.3 && fabs(vis2.Eta())<2.3);
    bool_values["eta21"]         = (fabs(vis1.Eta())<2.1 && fabs(vis2.Eta())<2.1);
    bool_values["reco_cuts"]     = (bool_values.at("eta23") && vis1.pt()>20 && vis2.pt()>20);
    bool_values["cut_1x_pt20"]   = (vis1.pt()>20 || vis2.pt()>20);
    //bool_values["cut_1x_pt30"]   = (vis1.pt()>30 || vis2.pt()>30);
    //bool_values["cut_1x_pt40"]   = (vis1.pt()>40 || vis2.pt()>40);
    //bool_values["cut_1x_pt45"]   = (vis1.pt()>45 || vis2.pt()>45);
    bool_values["cut_2x_pt45"]   = (vis1.pt()>45 && vis2.pt()>45);
    bool_values["baseline thth"] = bool_values.at("eta21") && bool_values.at("cut_2x_pt45");
  } catch(out_of_range oor) {
    cerr << "out of range in 'fill_histos_calculate_values'-1\n";
  }

  /*if(values.at("z")>z_max || values.at("zbar")>z_max) {
    static int count = 0;
    ++count;
    if(!count%25)
      cout << "count == " << count << endl;
    cout << name1+addition << endl;
    cout << "z == " << values.at("z") << "  zbar == " << values.at("zbar") << "   z_max == " << z_max << endl;
    cout << "vis1:  " << vis1 << endl;
    cout << "tau1:  " << tau1 << endl;
    cout << "vis2:  " << vis2 << endl;
    cout << "tau2:  " << tau2 << endl << endl;
  }*/
  if(values.at("z")>values.at("z_max") || values.at("zbar")>values.at("z_max"))
    return;

  //for(int i=20; i<=45; i+=5)
    //bool_values["cut_1x_pt" + to_string(i)] = (vis1.pt()>i || vis2.pt()>i);
  for(int i=10; i<=30; i+=2) {
    string cut = ("cut_2x_ptmin_" + to_string(i));
    bool_values[cut] = (vis1.pt()>i && vis2.pt()>i);
           cut = ("cut_ptmin_1x30_1x" + to_string(i));
    bool_values[cut] = ((vis1.pt()>30 && vis2.pt()>i) || (vis2.pt()>30 && vis1.pt()>i));
  }

  fill_histos_match(name1+addition, values, bool_values, !vis_only);
  //if(bool_values.at("full_had"))  //function not called if full_had and not full_pion
    //fill_histos_match(string("full_had")+addition, values, bool_values, !vis_only);
  if(bool_values.at("full_pion"))
    fill_histos_match(string("full_pion")+addition, values, bool_values, !vis_only);
  //fill_histos_match(name1+addition+string("_")+to_string(values.at("mass_h")), values, bool_values, !vis_only);

  if(keep_particle_order) {
    fill_histos_match(name2+addition, values, bool_values, !vis_only);
    //fill_histos_match(name2+addition+string("_")+to_string(values.at("mass_h")), values, bool_values, !vis_only);
  }
  else if(name2!="") {
    if(!vis_only)
      values["Easymm"] = Sven::asymm(vis2.E(),vis1.E());
    values["z"] = vis2.E()/tau2.E();
    values["zbar"] = vis1.E()/tau1.E();
    values["zasymm"] = Sven::asymm(values.at("z"), values.at("zbar"));
    fill_histos_match(name2+addition, values, bool_values, !vis_only);
    //fill_histos_match(name2+addition+string("_")+to_string(values.at("mass_h")), values, bool_values, !vis_only);
  }
}

void MDSM_LFV_HMuE_CL::fill_histos_cut(const string& name, const map<string,double>& values, const map<string,bool>& bool_values, bool Fill_Easymm) {
  try {
    fill_histos_call_create_and_fill(name, values, Fill_Easymm);
    if(bool_values.at("cut"))
      fill_histos_call_create_and_fill(name+string("_cut"), values, Fill_Easymm);
    if(bool_values.at("reco_cuts"))
      fill_histos_call_create_and_fill(name+string("_reco_cuts"), values, Fill_Easymm);
    if(bool_values.at("cut_1x_pt20"))
      fill_histos_call_create_and_fill(name+string("_cut_1x_pt20"), values, Fill_Easymm);
    if(bool_values.at("eta23"))
      for(int i=10; i<=30; i+=2) {
        string i_str="cut_2x_ptmin_"+to_string(i);
        //if(bool_values.at(i_str))
          //fill_histos_call_create_and_fill(name+string("_")+i_str, values, Fill_Easymm);
        i_str=("cut_ptmin_1x30_1x" + to_string(i));
        if(bool_values.at(i_str))
          fill_histos_call_create_and_fill(name+string("_")+i_str, values, Fill_Easymm);
      }
    if(bool_values.at("baseline thth"))
      fill_histos_call_create_and_fill(name+string("_baseline"), values, Fill_Easymm);
    if(bool_values.at("1 jet") && bool_values.at("boosted")) {
      fill_histos_call_create_and_fill(name+string("_1jet_boosted"), values, Fill_Easymm);
      if(bool_values.at("baseline thth"))
        fill_histos_call_create_and_fill(name+string("_baseline_1jet_boosted"), values, Fill_Easymm);
      if(bool_values.at("highly boosted")) {
        fill_histos_call_create_and_fill(name+string("_1jet_highly_boosted"), values, Fill_Easymm);
        if(bool_values.at("baseline thth"))
          fill_histos_call_create_and_fill(name+string("_baseline_1jet_highly_boosted"), values, Fill_Easymm);
      }
    }
    if(bool_values.at("reco_cuts"))
      for(int i=20; i<=45 && values.at("pt_high")>i; i+=5)
        fill_histos_call_create_and_fill(name+string("_cut_pt_20_" + to_string(i)), values, Fill_Easymm);
  } catch(std::out_of_range) {
    throw ErrorSven::somethingwenthorriblywrong("fill_histos_cut");
  }
}

void MDSM_LFV_HMuE_CL::fill_histos_call_create_and_fill(const string& name, const map<string,double>& values, bool Fill_Easymm) {
  try {
    double w = values.at("weight");
    create_and_fill_histos(h2_z_zbar, name, values.at("z"), values.at("zbar"), z_zbar, w);
    create_and_fill_histos(h2_z_zbar, string("pt_low_high_")+name, values.at("pt_low"), values.at("pt_high"), th_data_pt_vis, w);
    //create_and_fill_histos(h_Easymm, string("zasymm_")+name, values.at("zasymm"), Easymm, w);
    create_and_fill_histos(h_Easymm, string("z+zbar_")+name, values.at("z")+values.at("zbar"), zzbar_sum, w);
    create_and_fill_histos(h_Easymm, string("z+zbar2_")+name, fabs(values.at("z")+values.at("zbar")-1), zzbar_sum, w);
    fill_bg_rejec_histos(name, values);
    //if(Fill_Easymm)
      //create_and_fill_histos(h_Easymm, name, values.at("Easymm"), Easymm, w);
  } catch(std::out_of_range) {
    throw ErrorSven::somethingwenthorriblywrong("fill_histos_call_create_and_fill");
  }
}

void MDSM_LFV_HMuE_CL::fill_bg_rejec_histos(const string& name, const map<string,double>& values) {
  try {
    double w = values.at("weight");
    Double_t sum_zzbar = values.at("z")+values.at("zbar");
    for(int i=int(sum_zzbar*20); i>=0; --i)
      create_and_fill_histos(h_bg_rejec, name, i*.05, bg_rejec, w);
    for(int i=int(fabs(sum_zzbar-1)*20); i>=0; --i)
      create_and_fill_histos(h_bg_rejec, string("2_")+name, i*0.05, bg_rejec, w);
  } catch(std::out_of_range) {
    throw ErrorSven::somethingwenthorriblywrong("fill_bg_rejec_histos");
  }
}

void MDSM_LFV_HMuE_CL::pt_lowest_highest_histos(const std::vector<tau_decay>& gen_taus) {
#ifdef DEBUG_H_ditau
cout << "'pt_lowest_highest_histos' called\n";
#endif
  double pt_gen_lowest, pt_gen_highest;
  if((pt_gen_lowest=gen_taus.at(0).vis_vec()->pt())>(pt_gen_highest=gen_taus.at(1).vis_vec()->pt())) {
    pt_gen_lowest=gen_taus.at(1).vis_vec()->pt();
    pt_gen_highest=gen_taus.at(0).vis_vec()->pt();
  }
  create_and_fill_histos(h_, "pt_gen_lowest", pt_gen_lowest, th_data_pt_vis);
  create_and_fill_histos(h_, "pt_gen_highest", pt_gen_highest, th_data_pt_vis);
  if(gen_taus.at(0).is_hadronic() || gen_taus.at(1).is_hadronic()) {
    create_and_fill_histos(h_, "pt_gen_lowest_semi_had", pt_gen_lowest, th_data_pt_vis);
    create_and_fill_histos(h_, "pt_gen_highest_semi_had", pt_gen_highest, th_data_pt_vis);
    if(gen_taus.at(0).is_hadronic() && gen_taus.at(1).is_hadronic()) {
      create_and_fill_histos(h_, "pt_gen_lowest_full_had", pt_gen_lowest, th_data_pt_vis);
      create_and_fill_histos(h_, "pt_gen_highest_full_had", pt_gen_highest, th_data_pt_vis);
    }
    bool pionic[] = {false, false};
    for(int i=0; i<2; ++i) {
      try {
        switch (decay_modes.at(gen_taus.at(i).decay_mode())) {
          case DECAY_MODES::RHO:
          case DECAY_MODES::A1:
          case DECAY_MODES::PI:
            pionic[i]=true;
            break;
        }
      } catch(tau_decay::unknown_decay_mode) {
      }
    }
    if(pionic[0] && pionic[1]) {
      create_and_fill_histos(h_, "pt_gen_lowest_full_pion", pt_gen_lowest, th_data_pt_vis);
      create_and_fill_histos(h_, "pt_gen_highest_full_pion", pt_gen_highest, th_data_pt_vis);
    }
    else if(pionic[0] || pionic[1]) { //inconsistent!!
      create_and_fill_histos(h_, "pt_gen_lowest_semi_pion", pt_gen_lowest, th_data_pt_vis);
      create_and_fill_histos(h_, "pt_gen_highest_semi_pion", pt_gen_highest, th_data_pt_vis);
    }
  }
#ifdef DEBUG_H_ditau
cout << "end of 'pt_lowest_highest_histos'\n";
#endif
}

void MDSM_LFV_HMuE_CL::add_map_to_list(map<std::string,TH2F*> THmap) {
  for(map<std::string,TH2F*>::const_iterator iter=THmap.begin(); iter!=THmap.end(); iter++) {
    Hlist.Add(iter->second);
  }
}

void MDSM_LFV_HMuE_CL::add_map_to_list(map<std::string,TH1F*> THmap) {
  for(map<std::string,TH1F*>::const_iterator iter=THmap.begin(); iter!=THmap.end(); iter++) {
    Hlist.Add(iter->second);
  }
}

void MDSM_LFV_HMuE_CL::norm_histos() {
  for(map<std::string,TH2F*>::const_iterator iter=h2_Easymm.begin(); iter!=h2_Easymm.end(); iter++) {
    Double_t norm = 1;
    Double_t scale = norm/(iter->second->Integral());
    iter->second->Scale(scale);
  }
  for(map<std::string,TH2F*>::const_iterator iter=h2_z_zbar.begin(); iter!=h2_z_zbar.end(); iter++) {
    Double_t norm = 1;
    Double_t scale = norm/(iter->second->Integral());
    iter->second->Scale(scale);
  }
  for(map<std::string,TH1F*>::const_iterator iter=h_Easymm.begin(); iter!=h_Easymm.end(); iter++) {
    Double_t norm = 1;
    Double_t scale = norm/(iter->second->Integral());
    iter->second->Scale(scale);
  }
  /*for(map<string,TH1F*>::const_iterator iter=h_.begin(); iter!=h_.end(); iter++) {
    Double_t norm = 1;
    Double_t scale = norm/(iter->second->Integral());
    iter->second->Scale(scale);
  }*/
  /*for(map<string,TH1F*>::const_iterator iter=h_bg_rejec.begin(); iter!=h_bg_rejec.end(); ++iter) {
    Double_t norm = 1;
    Double_t scale = norm/(iter->second->Integral());
    iter->second->Scale(scale);
  }*/
}
