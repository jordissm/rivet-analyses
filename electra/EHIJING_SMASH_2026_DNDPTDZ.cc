// EHIJING_SMASH_2026_DNDPtDZ.cc
//
// Rivet analysis to produce (1/Nevt) dN/(dpT dz_h) vs pT in fixed z_h slices
// from SMASH final states + per-event DIS metadata.
//
// Species: pi+, pi-, K+, K-
// z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)
//
// Config via env vars (mirrors your existing approach):
//   RIVET_METAFILE         (required)
//   RIVET_FRAME            (default BREIT)
//
// Optional pT histogram axis config:
//   RIVET_PT_MIN   (default 0.0 GeV)
//   RIVET_PT_MAX   (default 1.1 GeV)
//   RIVET_PT_NBINS (default 20)
//
// Output histograms (Histo1D):
//   /EHIJING_SMASH_DNDPTDZ/dN_dpTdz_pip_z02_03
//   /EHIJING_SMASH_DNDPTDZ/dN_dpTdz_pip_z03_04
//   /EHIJING_SMASH_DNDPTDZ/dN_dpTdz_pip_z04_06
//   /EHIJING_SMASH_DNDPTDZ/dN_dpTdz_pip_z06_08
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

  static inline FourMomentum readVec4_px_py_pz_E(const boost::property_tree::ptree& arr) {
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

      MetaDIS eventKinematics;
      eventKinematics.event = pt.get<int>("event");
      eventKinematics.Z     = pt.get<int>("Z");
      eventKinematics.A     = pt.get<int>("A");
      eventKinematics.xB    = pt.get<double>("xB");
      eventKinematics.Q2    = pt.get<double>("Q2");
      eventKinematics.y     = pt.get<double>("y");
      eventKinematics.nu    = pt.get<double>("nu");

      eventKinematics.P = readVec4_px_py_pz_E(pt.get_child("P4"));
      eventKinematics.q = readVec4_px_py_pz_E(pt.get_child("q4"));

      out.emplace(eventKinematics.event, eventKinematics);
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

  static inline double pL_wrt_q(const FourMomentum& ph, const FourMomentum& q) {
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return 0.0;
    const Vector3 nhat = qv / qmag;
    return dot(ph.p3(), nhat);
  }

  class EHIJING_SMASH_DNDPTDZ : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(EHIJING_SMASH_DNDPTDZ);

    void init() override {
      declare(FinalState(), "FS");

      // METAFILE
      _metafile = getenv_str("RIVET_METAFILE");
      if (_metafile.empty()) _metafile = getOption<string>("METAFILE", "");
      if (_metafile.empty()) {
        throw UserError("You must provide METAFILE. Set env var RIVET_METAFILE=/path/DISKinematics.meta.json");
      }

      // Set default target reference frame
      std::string frameStr = getenv_str("RIVET_FRAME");
      if (frameStr.empty()) frameStr = getOption<string>("FRAME", "BREIT");
      _frame = parseFrame(frameStr);

      // pT axis config
      _ptMin = 0.0; _ptMax = 1.1; _ptNBins = 20;
      (void) getenv_double("RIVET_PT_MIN", _ptMin);
      (void) getenv_double("RIVET_PT_MAX", _ptMax);
      (void) getenv_int("RIVET_PT_NBINS", _ptNBins);
      if (!std::isfinite(_ptMin) || !std::isfinite(_ptMax) || _ptMax <= _ptMin) {
        throw UserError("Invalid pT axis: require PT_MAX > PT_MIN.");
      }
      if (_ptNBins <= 0) throw UserError("RIVET_PT_NBINS must be > 0");

      // Load kinematics metadata file
      _meta = loadMeta(_metafile);

      // Book histograms: [species][zbin]
      for (size_t is = 0; is < NSPEC; ++is) {
        for (size_t iz = 0; iz < NZ; ++iz) {
          const std::string name = "dN_dptdz_" + std::string(specTag(is)) + "_" + std::string(zTag(iz));
          book(_h[is][iz], name, _ptNBins, _ptMin, _ptMax);
        }
      }

      // Print config information
      MSG_INFO("Loaded " << _meta.size() << " metadata entries from " << _metafile);
      MSG_INFO("Frame choice (for pT): " << frameName(_frame));
      MSG_INFO("pT axis: nbins=" << _ptNBins << " range=[" << _ptMin << "," << _ptMax << "] GeV");
      MSG_INFO("z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)");
    }

    void analyze(const Event& event) override {

      // Counter for total number of events found in the input directory.
      _nEventsSeen++;

      const FinalState& fs = apply<FinalState>(event, "FS");
      const auto* ge = event.genEvent();
      const int evnum = ge ? ge->event_number() : -1;

      // Find metadata corresponding to the event number.
      // Veto event if metadata is not found.
      auto it = _meta.find(evnum);
      if (it == _meta.end()) {
        _nEventsVetoed++;
        vetoEvent;
      }
      _nEventsWithMeta++;

      // Load event DIS kinematics.
      const MetaDIS& eventKinematics = it->second;

      // Compute P·q from DIS kinematics metadata.
      //    P: four-momentum of the struck nucleon
      //    q: four-momentum of the exchanged virtual photon
      const double Pdotq = eventKinematics.P.dot(eventKinematics.q);
      // Veto events with problematic P·q.
      if (!std::isfinite(Pdotq) || Pdotq == 0) {
        _nBadPdotq++;
        _nEventsVetoed++;
        vetoEvent;
      }

      // Loop over all particles in the event.
      for (const Particle& particle : fs.particles()) {

        // Ignore leptons and photons.
        const int pid = particle.pid();
        if (PID::isLepton(pid)) continue;
        if (pid == 22) continue;

        // Set internal index for particles of interest.
        const int is = speciesIndex(pid);
        // Ignore particles that are not of interest.
        if (is < 0) continue;

        // Find observed hadron four-momentum ph.
        FourMomentum ph = particle.momentum();

        // Compute momentum fraction zh from DIS kinematics metadata.
        // zh = P·p_h / P·q
        //    P: Four-momentum of the struck nucleon
        //    q: Four-momentum of the exchanged virtual photon
        //    ph: Four-momentum of the observed hadron
        const double zh = (eventKinematics.P.dot(ph)) / Pdotq;
        // Compute momentum fraction from Eh/ν.
        // zh = Eh / ν (in TARGET REST FRAME)
        //    Eh: Energy of the observed hadron
        //    ν: Energy of the exchanged virtual photon
        const double zh_trf = ph.E() / eventKinematics.q.E();
        // Compare both ways of computing zh.
        MSG_INFO("[DEBUG] Invariant zh: " << zh << "; TRF zh: " << zh_trf);
        // Ignore particles with bad zh.
        if (!std::isfinite(zh)) continue;

        // Find bin corresponding to the observed zh.
        const int iz = zbinIndex(zh);
        if (iz < 0) continue;

        // Boost four-momenta to the frame of reference of interest.
        FourMomentum q  = eventKinematics.q;
        FourMomentum P  = eventKinematics.P;
        toFrame(_frame, ph, q, P);
        
        // Compute transverse-momentum pT w.r.t q in the requested frame.
        const double pT2 = pT2_wrt_q(ph, q);
        // Ignore hadrons with bad pT.
        if (!std::isfinite(pT2) || pT2 < 0) continue;
        const double pT = std::sqrt(pT2);

        // Add hadron to the histogram according to its
        // species and momentum fraction.
        _h[is][iz]->fill(pT, 1.0);
        _nFilled++;
      }
    }

    void finalize() override {
      MSG_INFO("======= Summary =======");
      MSG_INFO("Events seen:          " << _nEventsSeen);
      MSG_INFO("Events w/ metadata:   " << _nEventsWithMeta);
      MSG_INFO("Bad/zero P·q:         " << _nBadPdotq);
      MSG_INFO("Events vetoed:     " << _nEventsVetoed);
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
            const double dpT = bin.xWidth();
            if (dpT > 0) bin.scaleW(1.0 / (dpT * dz));
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
      for (size_t i = 0; i < NZ; ++i) {
        if (zh > zEdges[i] && zh <= zEdges[i+1]) return (int)i;
      }
      return -1;
    }


  private:
    // Config
    std::string _metafile;
    FrameChoice _frame = FrameChoice::TRF;

    double _ptMin = 0.0, _ptMax = 1.1;
    int _ptNBins = 20;

    // Data
    std::unordered_map<int, MetaDIS> _meta;

    // Histos
    std::array<std::array<Histo1DPtr, NZ>, NSPEC> _h;

    // Debug
    size_t _nEventsSeen = 0;
    size_t _nEventsWithMeta = 0;
    size_t _nBadPdotq = 0;
    size_t _nEventsVetoed = 0;
    size_t _nFilled = 0;
  };

  constexpr double EHIJING_SMASH_DNDPTDZ::zEdges[NZ+1];

  RIVET_DECLARE_PLUGIN(EHIJING_SMASH_DNDPTDZ);

} // namespace Rivet
