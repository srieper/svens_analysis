#ifndef SVEN__STUFF
#define SVEN__STUFF

#include <utility>
#include <string>
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include <map>
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include <ostream>

using std::string;

namespace Sven {
  const double pi=3.14159265358979323846264338328;

  inline double asymm(double x1, double x2) { return ((x1-x2)/(x1+x2)); }
  inline std::pair<double,double> relValue2(double x1, double x2) { return std::make_pair((x1-x2)/x1 ,(x1-x2)/x2); }

  template <class VectorO, class VectorI>
  inline void convert_vector(VectorO& out, const VectorI& in) { out.SetXYZT(in.x(), in.y(), in.z(), in.t()); }
  inline void XYZT_to_TLorentzVector(TLorentzVector* out, const ROOT::Math::XYZTVector* in) { out->SetXYZT(in->x(), in->y(), in->z(), in->t()); }
  inline TLorentzVector to_TLorentzVector(const ROOT::Math::XYZTVector& in) { return TLorentzVector(in.x(), in.y(), in.z(), in.t()); }

  struct TH_data {
    std::string name, title;
    double xlow, xup, ylow, yup;
    int nbinsx, nbinsy;
    TH_data() :name(""), title(""), xlow(0), xup(0), ylow(0), yup(0), nbinsx(0), nbinsy(0) {}
    TH_data(std::string nn, std::string tt, int nx, double xl, double xu, int ny, double yl, double yu)
      :name(nn), title(tt), xlow(xl), xup(xu), ylow(yl), yup(yu), nbinsx(nx), nbinsy(ny) {}
    TH_data(std::string nn, std::string tt, int nx, double xl, double xu)
      :name(nn), title(tt), xlow(xl), xup(xu), ylow(xl), yup(xu), nbinsx(nx), nbinsy(nx) {}
    TH_data(std::string nn, std::string tt, int nx, double xu)
      :name(nn), title(tt), xlow(-xu), xup(xu), ylow(-xu), yup(xu), nbinsx(nx), nbinsy(nx) {}
  };

  enum DECAY_MODES {RHO=1, E=2, MU=4, PI=8, A1=16, HAD=32};

  inline void createTH2F(TH2F** histo, const Sven::TH_data& thdata, const std::string& name2="") {
    *histo = new TH2F((thdata.name+name2).c_str(), thdata.title.c_str(), thdata.nbinsx, thdata.xlow, thdata.xup, thdata.nbinsy, thdata.ylow, thdata.yup); }

  inline void createTH1F(TH1F** histo, const Sven::TH_data& thdata, const std::string& name2="") {
    *histo = new TH1F((thdata.name+name2).c_str(), thdata.title.c_str(), thdata.nbinsx, thdata.xlow, thdata.xup); }

  inline void create_and_fill_histos(std::map<std::string, TH1F*>& THmap, const std::string& name, double x, const Sven::TH_data& thdata, double w=1) {
    try {
      THmap.at(name)->Fill(x,w);  
    } catch(std::out_of_range) {
      createTH1F(&THmap[name], thdata, name);
      THmap.at(name)->Fill(x,w);
    }
  }

  inline void create_and_fill_histos(std::map<std::string, TH2F*>& THmap, const std::string& name, double x, double y, const Sven::TH_data& thdata, double w=1) {
    try {
      THmap.at(name)->Fill(x,y,w);
    } catch(std::out_of_range) {
      createTH2F(&THmap[name], thdata, name);
      THmap.at(name)->Fill(x,y,w);
    }
  }

  std::ostream& operator<<(std::ostream&, TLorentzVector);
  std::ostream& operator<<(std::ostream&, ROOT::Math::XYZTVector);

  template <class Vector1, class Vector2>
  inline double DeltaPhi( const Vector1 & v1, const Vector2 & v2) {
    double dphi = v2.Phi() - v1.Phi();
    if ( dphi > pi ) {
      dphi -= 2.0*pi;
    } else if ( dphi <= -pi ) {
      dphi += 2.0*pi;
    }
    return dphi;
  }

  template <class Vector1, class Vector2>
  inline double DeltaR2( const Vector1 & v1, const Vector2 & v2) {
    double dphi = DeltaPhi(v1,v2);
    double deta = v2.Eta() - v1.Eta();
    return dphi*dphi + deta*deta;
  }

  template <class Vector1, class Vector2>
  inline double DeltaR( const Vector1 & v1, const Vector2 & v2) {
    return std::sqrt( DeltaR2(v1,v2) );
  }
}

namespace ErrorSven {
  struct missing_fit{};
  struct boson_error{
    int id;
    boson_error(int i) :id(i) {}
  };
  struct gen_error{};
  struct somethingwenthorriblywrong{
    const char* func;
    somethingwenthorriblywrong() :func("") {}
    somethingwenthorriblywrong(const char* f) :func(f) {}
  };
}

#endif
