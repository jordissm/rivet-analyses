// EHIJING_SMASH_2026_DNDPtDZ.cc
//
// Rivet analysis to produce (1/Nevt) dN/(dpT dz_h) vs pT in fixed z_h slices
// from SMASH final states + per-event DIS metadata.
//
// Species: pi+, pi-, K+, K-
// z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)
//
// Conventions:
//  - z_h computed invariantly in LAB: z_h = (P·p_h)/(P·q)
//  - pT computed w.r.t. q direction in chosen frame: LAB|TRF|BREIT
//
// Config via env vars (mirrors your existing approach):
//   RIVET_METAFILE         (required)
//   RIVET_FRAME            (default TRF)
//   RIVET_VETO_SPECTATORS  (default 1)
//   RIVET_SPECTATOR_PMAX   (default 0.30 GeV)
//   RIVET_SPECTATOR_TKMAX  (default -1 disabled)
//
// Optional pT histogram axis config:
//   RIVET_PT_MIN   (default 0.0 GeV)
//   RIVET_PT_MAX   (default 2.0 GeV)
//   RIVET_PT_NBINS (default 40)
//
// Output histograms (Histo1D):
//   /EHIJING_SMASH_DNDPtDZ/dN_dpTdz_pip_z02_03
//   /EHIJING_SMASH_DNDPtDZ/dN_dpTdz_pip_z03_04
//   /EHIJING_SMASH_DNDPtDZ/dN_dpTdz_pip_z04_06
//   /EHIJING_SMASH_DNDPtDZ/dN_dpTdz_pip_z06_08
//   ... similarly for pim, kp, km
//
// Normalization in finalize:
//  - divide by Nev = sumW()
//  - convert to density in pT and z: scale by 1/(dpT * dz)

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Math/LorentzTrans.hh"

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <array>

namespace Rivet {

  struct MetaDIS {
    int event = -1;
    int Z = 0, A = 0;
    double xB = NAN, Q2 = NAN, y = NAN, nu = NAN;
    FourMomentum P;
    FourMomentum q;
  };

  enum class FrameChoice { LAB, TRF, BREIT };

  static inline FrameChoice parseFrame(const std::string& s) {
    if (s == "LAB")   return FrameChoice::LAB;
    if (s == "TRF")   return FrameChoice::TRF;
    if (s == "BREIT") return FrameChoice::BREIT;
    throw UserError("Unknown FRAME='" + s + "'. Use LAB, TRF, or BREIT.");
  }

  static inline std::string frameName(FrameChoice f) {
    if (f == FrameChoice::LAB) return "LAB";
    if (f == FrameChoice::TRF) return "TRF";
    return "BREIT";
  }

  static inline std::string getenv_str(const char* key) {
    const char* v = std::getenv(key);
    return v ? std::string(v) : std::string();
  }

  static inline bool getenv_double(const char* key, double& out) {
    const std::string s = getenv_str(key);
    if (s.empty()) return false;
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || !std::isfinite(v)) return false;
    out = v;
    return true;
  }

  static inline bool getenv_int(const char* key, int& out) {
    const std::string s = getenv_str(key);
    if (s.empty()) return false;
    char* end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str()) return false;
    out = (int)v;
    return true;
  }

  static inline bool getenv_bool_default(const char* key, bool defval) {
    const std::string s = getenv_str(key);
    if (s.empty()) return defval;
    if (s == "1" || s == "true" || s == "TRUE" || s == "yes" || s == "YES") return true;
    if (s == "0" || s == "false" || s == "FALSE" || s == "no"  || s == "NO")  return false;
    return defval;
  }

  static std::vector<std::string> readConcatenatedJsonObjects(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw UserError("Cannot open METAFILE='" + path + "'");

    std::vector<std::string> objs;
    std::string line, buf;
    int depth = 0;
    bool inStr = false;
    bool escape = false;
    bool started = false;

    while (std::getline(fin, line)) {
      for (char c : line) {
        if (!started) {
          if (std::isspace(static_cast<unsigned char>(c))) continue;
          if (c == '{') {
            started = true;
            depth = 1;
            buf.clear();
            buf.push_back(c);
            inStr = false; escape = false;
          }
          continue;
        } else {
          buf.push_back(c);

          if (escape) { escape = false; continue; }
          if (c == '\\') { if (inStr) escape = true; continue; }
          if (c == '"') { inStr = !inStr; continue; }
          if (inStr) continue;

          if (c == '{') depth++;
          else if (c == '}') depth--;

          if (depth == 0) {
            objs.push_back(buf);
            buf.clear();
            started = false;
          }
        }
      }
      if (started) buf.push_back('\n');
    }

    if (started || depth != 0) {
      throw UserError("METAFILE appears truncated/unbalanced braces: '" + path + "'");
    }
    return objs;
  }

  static inline FourMomentum readVec4_pxpy_pz_E(const boost::property_tree::ptree& arr) {
    std::vector<double> v;
    v.reserve(4);
    for (const auto& kv : arr) v.push_back(kv.second.get_value<double>());
    if (v.size() != 4) throw UserError("Expected 4-vector array of length 4");
    return FourMomentum(v[3], v[0], v[1], v[2]);
  }

  static std::unordered_map<int, MetaDIS> loadMeta(const std::string& path) {
    std::unordered_map<int, MetaDIS> out;

    const auto objs = readConcatenatedJsonObjects(path);
    for (const std::string& js : objs) {
      std::stringstream ss(js);
      boost::property_tree::ptree pt;
      boost::property_tree::read_json(ss, pt);

      MetaDIS m;
      m.event = pt.get<int>("event");
      m.Z     = pt.get<int>("Z");
      m.A     = pt.get<int>("A");
      m.xB    = pt.get<double>("xB");
      m.Q2    = pt.get<double>("Q2");
      m.y     = pt.get<double>("y");
      m.nu    = pt.get<double>("nu");

      m.P = readVec4_pxpy_pz_E(pt.get_child("P4"));
      m.q = readVec4_pxpy_pz_E(pt.get_child("q4"));

      out.emplace(m.event, m);
    }

    if (out.empty()) throw UserError("METAFILE had zero parsed events: '" + path + "'");
    return out;
  }

  static inline void boostBy(FourMomentum& p, const Vector3& beta) {
    const LorentzTransform tr = LorentzTransform::mkFrameTransformFromBeta(beta);
    p = tr.transform(p);
  }

  static inline double pT2_wrt_q(const FourMomentum& ph, const FourMomentum& q) {
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return 0.0;

    const Vector3 nhat = qv / qmag;
    const Vector3 pv = ph.p3();
    const double ppar = dot(pv, nhat);
    const Vector3 pT = pv - ppar * nhat;
    return pT.mod2();
  }

  static inline void toFrame(FrameChoice frame, FourMomentum& ph, FourMomentum& q, FourMomentum& P) {
    if (frame == FrameChoice::LAB) return;

    // TRF: boost so target is at rest
    const double PE = P.E();
    if (PE <= 0) return;
    const Vector3 betaP = P.p3() / PE;
    const Vector3 bTRF  = -betaP;

    boostBy(ph, bTRF);
    boostBy(q,  bTRF);
    boostBy(P,  bTRF);

    if (frame == FrameChoice::TRF) return;

    // Breit: boost along q so that q^0 = 0
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return;

    const double beta = q.E() / qmag;
    if (!std::isfinite(beta) || std::abs(beta) >= 1.0) return;

    const Vector3 nhat = qv / qmag;
    const Vector3 bBreit = beta * nhat;

    boostBy(ph, bBreit);
    boostBy(q,  bBreit);
    boostBy(P,  bBreit);
  }

  static inline void toGammaNCM(FourMomentum& ph, FourMomentum& q, FourMomentum& P) {
    const FourMomentum W = P + q;
    const double WE = W.E();
    if (WE <= 0) return;
    const Vector3 betaW = W.p3() / WE;
    const Vector3 bCM   = -betaW;
    boostBy(ph, bCM);
    boostBy(q,  bCM);
    boostBy(P,  bCM);
  }

  static inline double pL_wrt_q(const FourMomentum& ph, const FourMomentum& q) {
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return 0.0;
    const Vector3 nhat = qv / qmag;
    return dot(ph.p3(), nhat);
  }

  class EHIJING_SMASH_DNDPtDZ : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(EHIJING_SMASH_DNDPtDZ);

    void init() override {
      declare(FinalState(), "FS");

      // METAFILE
      _metafile = getenv_str("RIVET_METAFILE");
      if (_metafile.empty()) _metafile = getOption<string>("METAFILE", "");
      if (_metafile.empty()) {
        throw UserError("You must provide METAFILE. Set env var RIVET_METAFILE=/path/run.meta.json");
      }

      // FRAME
      std::string frameStr = getenv_str("RIVET_FRAME");
      if (frameStr.empty()) frameStr = getOption<string>("FRAME", "TRF");
      _frame = parseFrame(frameStr);

      // Spectator veto config
      _vetoSpectators = getenv_bool_default("RIVET_VETO_SPECTATORS", true);
      _spectatorPmax = 0.30;
      (void) getenv_double("RIVET_SPECTATOR_PMAX", _spectatorPmax);
      _spectatorTkmax = -1.0;
      (void) getenv_double("RIVET_SPECTATOR_TKMAX", _spectatorTkmax);
      if (_spectatorPmax < 0) throw UserError("RIVET_SPECTATOR_PMAX must be >= 0");

      // pT axis config
      _ptMin = 0.0; _ptMax = 2.0; _ptNBins = 40;
      (void) getenv_double("RIVET_PT_MIN", _ptMin);
      (void) getenv_double("RIVET_PT_MAX", _ptMax);
      (void) getenv_int("RIVET_PT_NBINS", _ptNBins);
      if (!std::isfinite(_ptMin) || !std::isfinite(_ptMax) || _ptMax <= _ptMin) {
        throw UserError("Invalid pT axis: require PT_MAX > PT_MIN.");
      }
      if (_ptNBins <= 0) throw UserError("RIVET_PT_NBINS must be > 0");

      _meta = loadMeta(_metafile);

      // Book histograms: [species][zbin]
      for (size_t is = 0; is < NSPEC; ++is) {
        for (size_t iz = 0; iz < NZ; ++iz) {
          const std::string name = "dN_dptdz_" + std::string(specTag(is)) + "_" + std::string(zTag(iz));
          book(_h[is][iz], name, _ptNBins, _ptMin, _ptMax);
        }
      }

      MSG_INFO("Loaded " << _meta.size() << " metadata entries from " << _metafile);
      MSG_INFO("Frame choice (for pT): " << frameName(_frame));
      MSG_INFO("pT axis: nbins=" << _ptNBins << " range=[" << _ptMin << "," << _ptMax << "] GeV");
      MSG_INFO("z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)");
      if (_vetoSpectators) {
        MSG_INFO("Spectator veto ENABLED (gamma*-N CM): veto nucleons with pL<0 and "
                 << "[ |p|<" << _spectatorPmax << " GeV"
                 << (_spectatorTkmax >= 0 ? (std::string(" OR (E-m)<") + std::to_string(_spectatorTkmax) + " GeV") : std::string(""))
                 << " ]");
      } else {
        MSG_INFO("Spectator veto DISABLED.");
      }
    }

    void analyze(const Event& event) override {
      _nEventsSeen++;

      const FinalState& fs = apply<FinalState>(event, "FS");
      const auto* ge = event.genEvent();
      const int evnum = ge ? ge->event_number() : -1;

      auto it = _meta.find(evnum);
      if (it == _meta.end()) {
        vetoEvent;
      }
      _nEventsWithMeta++;

      const MetaDIS& m = it->second;

      const double Pdotq_lab = m.P.dot(m.q);
      if (!std::isfinite(Pdotq_lab) || Pdotq_lab == 0) {
        _nBadPdotq++;
        vetoEvent;
      }

      for (const Particle& p : fs.particles()) {
        const int pid = p.pid();
        if (PID::isLepton(pid)) continue;
        if (pid == 22) continue;

        const int is = speciesIndex(pid);
        if (is < 0) continue;

        // z_h (LAB invariant)
        const double zh = (m.P.dot(p.momentum())) / Pdotq_lab;
        if (!std::isfinite(zh)) continue;

        const int iz = zbinIndex(zh);
        if (iz < 0) continue;

        // Spectator veto
        if (_vetoSpectators && isSpectatorNucleon_(pid, p.momentum(), m.P, m.q)) {
          _nSpectatorVeto++;
          continue;
        }

        // pT w.r.t q in requested frame
        FourMomentum ph = p.momentum();
        FourMomentum q  = m.q;
        FourMomentum P  = m.P;
        toFrame(_frame, ph, q, P);

        const double pt2 = pT2_wrt_q(ph, q);
        if (!std::isfinite(pt2) || pt2 < 0) continue;
        const double pt = std::sqrt(pt2);

        _h[is][iz]->fill(pt, 1.0);
        _nFilled++;
      }
    }

    void finalize() override {
      MSG_INFO("==== Debug summary ====");
      MSG_INFO("Events seen:          " << _nEventsSeen);
      MSG_INFO("Events w/ metadata:   " << _nEventsWithMeta);
      MSG_INFO("Bad/zero P·q:         " << _nBadPdotq);
      if (_vetoSpectators) MSG_INFO("Spectator vetoed:     " << _nSpectatorVeto);
      MSG_INFO("Filled entries:       " << _nFilled);
      MSG_INFO("=======================");

      // Per-event normalization
      const double Nev = (sumW() > 0) ? sumW() : 1.0;
      const double perEvent = 1.0 / Nev;

      // Convert to density: (1/Nev) dN/(dpT dz)
      for (size_t is = 0; is < NSPEC; ++is) {
        for (size_t iz = 0; iz < NZ; ++iz) {
          if (!_h[is][iz]) continue;

          scale(_h[is][iz], perEvent);

          const double dz = (zEdges[iz+1] - zEdges[iz]);
          if (!(dz > 0)) continue;

          for (auto& bin : _h[is][iz]->bins()) {
            const double dpt = bin.xWidth();
            if (dpt > 0) bin.scaleW(1.0 / (dpt * dz));
          }
        }
      }
    }

  private:
    // Species mapping
    static constexpr size_t NSPEC = 4;
    static constexpr size_t NZ = 4;

    // z edges: [0.2, 0.3, 0.4, 0.6, 0.8]
    static constexpr double zEdges[NZ+1] = {0.2, 0.3, 0.4, 0.6, 0.8};

    static inline const char* specTag(size_t i) {
      switch (i) {
        case 0: return "pip";
        case 1: return "pim";
        case 2: return "kp";
        case 3: return "km";
        default: return "unk";
      }
    }

    static inline const char* zTag(size_t i) {
      switch (i) {
        case 0: return "zh0p2_0p3";
        case 1: return "zh0p3_0p4";
        case 2: return "zh0p4_0p6";
        case 3: return "zh0p6_0p8";
        default: return "zhX";
      }
    }

    static inline int speciesIndex(int pid) {
      switch (pid) {
        case  211: return 0; // pi+
        case -211: return 1; // pi-
        case  321: return 2; // K+
        case -321: return 3; // K-
        default:   return -1;
      }
    }

    static inline int zbinIndex(double zh) {
      // open intervals (0.2,0.3), etc. as you requested
      for (size_t i = 0; i < NZ; ++i) {
        if (zh > zEdges[i] && zh < zEdges[i+1]) return (int)i;
      }
      return -1;
    }

    bool isSpectatorNucleon_(int pid,
                             const FourMomentum& ph_lab,
                             const FourMomentum& P_lab,
                             const FourMomentum& q_lab) {
      const int apid = std::abs(pid);
      if (!(apid == 2212 || apid == 2112)) return false;

      FourMomentum ph = ph_lab;
      FourMomentum P  = P_lab;
      FourMomentum q  = q_lab;
      toGammaNCM(ph, q, P);

      const double pL = pL_wrt_q(ph, q);
      const double pmag = ph.p3().mod();
      const double m = ph.mass();
      const double Tk = ph.E() - m;

      if (!(pL < 0)) return false;

      bool slow = false;
      if (pmag < _spectatorPmax) slow = true;
      if (_spectatorTkmax >= 0.0 && Tk < _spectatorTkmax) slow = true;

      return slow;
    }

  private:
    // Config
    std::string _metafile;
    FrameChoice _frame = FrameChoice::TRF;

    bool _vetoSpectators = false;
    double _spectatorPmax = 0.30;
    double _spectatorTkmax = -1.0;

    double _ptMin = 0.0, _ptMax = 2.0;
    int _ptNBins = 40;

    // Data
    std::unordered_map<int, MetaDIS> _meta;

    // Histos
    std::array<std::array<Histo1DPtr, NZ>, NSPEC> _h;

    // Debug
    size_t _nEventsSeen = 0;
    size_t _nEventsWithMeta = 0;
    size_t _nBadPdotq = 0;
    size_t _nSpectatorVeto = 0;
    size_t _nFilled = 0;
  };

  constexpr double EHIJING_SMASH_DNDPtDZ::zEdges[NZ+1];

  RIVET_DECLARE_PLUGIN(EHIJING_SMASH_DNDPtDZ);

} // namespace Rivet
