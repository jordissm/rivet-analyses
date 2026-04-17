// EHIJING_SMASH_DNDZ.cc
//
// Rivet analysis to produce dN/dz_h (optionally restricted to a pT^2 window)
// from SMASH final states + per-event DIS metadata.
//
// IMPORTANT (your requested behavior):
//   - Infer the *original eHIJING event number* from the *input HepMC file path*,
//     by extracting the token `evt_XXXXXX` (6 digits) from that path.
//   - Use that inferred event id as the lookup key into the metadata map for ALL
//     SMASH replicas in that folder/file.
//   - Still supports Option A attribute-based keying if present ("ehijing_event"/"orig_event"),
//     but path-based inference is the required fallback when replicas have renumbered events.
//
// How to provide the path to the Rivet plugin:
//   - Set env var: RIVET_HEPMC_PATH=/full/path/to/.../evt_000199/.../SMASH_HepMC_....asciiv3
//     OR pass analysis option: -a EHIJING_SMASH_DNDZ:HEPMC_PATH=/full/path/...
//
// Why this is necessary:
//   - Rivet plugins generally do NOT get the input filename from the Event object,
//     so we must receive it via env/option.
//
// Other notes unchanged below.

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Math/LorentzTrans.hh"

// HepMC3 ancestry access
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// HepMC3 attributes (Option A)
#include "HepMC3/Attribute.h"

// Boost property_tree
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// Std
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <regex>

namespace Rivet {

  // =========================
  // HepMC3 parent utilities
  // =========================

  // Return true if particle p has an IMMEDIATE parent whose PDG id is in `parentPids`.
  // i.e. checks only incoming particles to p.production_vertex().
  static inline bool hasDirectParentPid_HepMC3(const Particle& p,
                                              const std::unordered_set<int>& parentPids) {
    const auto gp = p.genParticle();
    if (!gp) return false;

    const auto vtx = gp->production_vertex();
    if (!vtx) return false;

    for (const auto& parentPtr : vtx->particles_in()) {
      const HepMC3::GenParticle* parent = parentPtr.get();
      if (!parent) continue;
      if (parentPids.count(parent->pid())) return true;
    }
    return false;
  }

  static inline bool isDirectFromVectorMesonOrKstar_HepMC3(const Particle& p) {
    // rho, omega, phi, K*
    static const std::unordered_set<int> VMK = {
      113, 213, -213, 223, 333
      // 313, -313, 323, -323
    };
    return hasDirectParentPid_HepMC3(p, VMK);
  }

  // =========================
  // eHIJING event id inference
  // =========================

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

  static inline bool getenv_bool_default(const char* key, bool defval) {
    const std::string s = getenv_str(key);
    if (s.empty()) return defval;
    if (s == "1" || s == "true" || s == "TRUE" || s == "yes" || s == "YES") return true;
    if (s == "0" || s == "false" || s == "FALSE" || s == "no"  || s == "NO")  return false;
    return defval;
  }

  /// Extract evt_XXXXXX (6 digits) from a path, return integer or -1 if not found.
  static inline int extractEvtFromPath(const std::string& path) {
    // Match ".../evt_000199/..." anywhere in the string
    static const std::regex re(R"(evt_(\d{6}))");
    std::smatch m;
    if (!std::regex_search(path, m, re)) return -1;
    try {
      return std::stoi(m[1].str());
    } catch (...) {
      return -1;
    }
  }

  /// Option A attribute-based key (if SMASH writes it).
  static inline int metaKeyFromEventAttribute(const HepMC3::GenEvent* ge) {
    if (!ge) return -1;
    if (auto a = ge->attribute<HepMC3::IntAttribute>("ehijing_event")) return a->value();
    if (auto a = ge->attribute<HepMC3::IntAttribute>("orig_event"))    return a->value();
    return -1;
  }

  // =========================
  // Analysis code
  // =========================

  struct MetaDIS {
    int event = -1;
    int Z = 0, A = 0;
    double xB = NAN, Q2 = NAN, y = NAN, nu = NAN;
    FourMomentum P; // target 4-mom (E,px,py,pz)
    FourMomentum q; // photon 4-mom (E,px,py,pz)
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

  /// Read concatenated JSON objects from file (no surrounding array).
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
    // Expect [px, py, pz, E]
    std::vector<double> v;
    v.reserve(4);
    for (const auto& kv : arr) v.push_back(kv.second.get_value<double>());
    if (v.size() != 4) throw UserError("Expected 4-vector array of length 4");
    return FourMomentum(v[3], v[0], v[1], v[2]); // (E,px,py,pz)
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

  /// Boost p by velocity-vector `beta` (|beta|<1) using Rivet LorentzTransform.
  static inline void boostBy(FourMomentum& p, const Vector3& beta) {
    const LorentzTransform tr = LorentzTransform::mkFrameTransformFromBeta(beta);
    p = tr.transform(p);
  }

  /// Compute pT^2 w.r.t. photon direction in the *current* frame of q and ph
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

  /// Apply frame choice boosts to (ph, q, P) in-place.
  static inline void toFrame(FrameChoice frame, FourMomentum& ph, FourMomentum& q, FourMomentum& P) {
    if (frame == FrameChoice::LAB) return;

    // 1) TRF: boost so that target P is at rest
    const double PE = P.E();
    if (PE <= 0) return;

    const Vector3 betaP = P.p3() / PE;
    const Vector3 bTRF  = -betaP;

    boostBy(ph, bTRF);
    boostBy(q,  bTRF);
    boostBy(P,  bTRF);

    if (frame == FrameChoice::TRF) return;

    // 2) Breit: boost along q direction so that q^0 = 0
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return;

    const double beta = q.E() / qmag;
    if (!std::isfinite(beta) || std::abs(beta) >= 1.0) return;

    const Vector3 nhat = qv / qmag;
    const Vector3 bBreit = beta * nhat;

    // Sign convention kept as in your original code
    boostBy(ph, bBreit);
    boostBy(q,  bBreit);
    boostBy(P,  bBreit);
  }

  /// Boost (P,q,ph) into the gamma*-N CM frame: rest frame of W = P + q
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

  /// Helper: return longitudinal momentum along q direction in current frame
  static inline double pL_wrt_q(const FourMomentum& ph, const FourMomentum& q) {
    const Vector3 qv = q.p3();
    const double qmag = qv.mod();
    if (qmag == 0) return 0.0;
    const Vector3 nhat = qv / qmag;
    return dot(ph.p3(), nhat);
  }


  class EHIJING_SMASH_2026_DNDZ : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(EHIJING_SMASH_2026_DNDZ);

    void init() override {
      declare(FinalState(), "FS");

      // METAFILE
      _metafile = getenv_str("RIVET_METAFILE");
      if (_metafile.empty()) _metafile = getOption<string>("METAFILE", "");
      if (_metafile.empty()) {
        throw UserError("You must provide METAFILE. Set env var RIVET_METAFILE=/path/run.meta.json");
      }

      // INPUT HepMC path (needed to infer evt_XXXXXX)
      _hepmcPath = getenv_str("RIVET_HEPMC_PATH");
      if (_hepmcPath.empty()) _hepmcPath = getOption<string>("HEPMC_PATH", "");
      _evtFromPath = extractEvtFromPath(_hepmcPath);

      if (_evtFromPath < 0) {
        throw UserError(
          "Could not infer eHIJING event id from input path.\n"
          "Expected the HepMC file path to contain 'evt_XXXXXX' (6 digits).\n"
          "Provide it via env var RIVET_HEPMC_PATH=/.../evt_000199/.../file.hepmc\n"
          "or analysis option :HEPMC_PATH=/.../evt_000199/.../file.hepmc"
        );
      }

      // FRAME (for pT^2)
      std::string frameStr = getenv_str("RIVET_FRAME");
      if (frameStr.empty()) frameStr = getOption<string>("FRAME", "TRF");
      _frame = parseFrame(frameStr);

      // Optional pT^2 window (GeV^2)
      _usePt2Min = getenv_double("RIVET_PT2_MIN", _pt2Min);
      _usePt2Max = getenv_double("RIVET_PT2_MAX", _pt2Max);
      if (_usePt2Min && _usePt2Max && _pt2Max < _pt2Min) {
        throw UserError("RIVET_PT2_MAX < RIVET_PT2_MIN (swap them or fix values).");
      }

      // Spectator veto config
      _vetoSpectators = getenv_bool_default("RIVET_VETO_SPECTATORS", false);
      _spectatorPmax = 0.30;   // GeV
      (void) getenv_double("RIVET_SPECTATOR_PMAX", _spectatorPmax);
      _spectatorTkmax = -1.0;  // GeV (disabled)
      (void) getenv_double("RIVET_SPECTATOR_TKMAX", _spectatorTkmax);

      if (_spectatorPmax < 0) throw UserError("RIVET_SPECTATOR_PMAX must be >= 0");

      _meta = loadMeta(_metafile);

      // Output histograms: dN/dz_h per species (density in finalize)
      book(_h_dN_dzh_pip,  "dN_dzh_pip",  10, 0.1, 1.1);
      book(_h_dN_dzh_pim,  "dN_dzh_pim",  10, 0.1, 1.1);
      book(_h_dN_dzh_kp,   "dN_dzh_kp",   10, 0.1, 1.1);
      book(_h_dN_dzh_km,   "dN_dzh_km",   10, 0.1, 1.1);
      book(_h_dN_dzh_p,    "dN_dzh_p",    10, 0.1, 1.1);
      book(_h_dN_dzh_pbar, "dN_dzh_pbar", 10, 0.1, 1.1);

      MSG_INFO("Loaded " << _meta.size() << " metadata entries from " << _metafile);
      MSG_INFO("Inferred eHIJING event id from HepMC path: evt_" << std::setw(6) << std::setfill('0') << _evtFromPath);
      MSG_INFO("HepMC path used for inference: " << _hepmcPath);

      MSG_INFO("Frame choice (for pT^2): " << frameName(_frame));
      if (_usePt2Min || _usePt2Max) {
        MSG_INFO("pT^2 window enabled: "
                 << (_usePt2Min ? std::to_string(_pt2Min) : std::string("-inf"))
                 << " <= pT^2 < "
                 << (_usePt2Max ? std::to_string(_pt2Max) : std::string("+inf"))
                 << "  (GeV^2)");
      } else {
        MSG_INFO("pT^2 window disabled (integrating over all pT^2).");
      }

      MSG_INFO("Applying DIS cuts ALWAYS: Q2>1.0 GeV^2, W2>10.0 GeV^2, 0.1<y<0.85, 0.023<xB<0.60");
      MSG_INFO("Applying hadron momentum cut ALWAYS (LAB): 2.0<|p_h|<15.0 GeV");
      MSG_INFO("Applying VM subtraction ONLY for pi/K and ONLY for immediate parent in {rho,omega,phi,K*}");

      MSG_INFO("Metadata keying policy:");
      MSG_INFO("  (1) If HepMC3 event has IntAttribute 'ehijing_event' or 'orig_event', use it");
      MSG_INFO("  (2) Else use evt_XXXXXX inferred from HepMC file path (this is your required behavior)");
      MSG_INFO("  (3) Else fallback to HepMC event_number() (should rarely happen now)");

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
      const int raw_evnum = ge ? ge->event_number() : -1;

      // (1) Option A attribute if present, else (2) evt_XXXXXX from path, else (3) raw event_number
      int meta_key = metaKeyFromEventAttribute(ge);
      if (meta_key < 0) meta_key = _evtFromPath;
      if (meta_key < 0) meta_key = raw_evnum;

      auto it = _meta.find(meta_key);
      if (it == _meta.end()) {
        MSG_WARNING("No metadata for meta_key=" << meta_key
                    << " (raw event_number=" << raw_evnum
                    << ", evtFromPath=" << _evtFromPath << ") -> veto event");
        _nEventsNoMeta++;
        vetoEvent;
      }
      _nEventsWithMeta++;

      const MetaDIS& m = it->second;

      // Require finite DIS metadata
      if (!std::isfinite(m.Q2) || !std::isfinite(m.y) || !std::isfinite(m.xB)) {
        _nEventsVetoDIS++;
        vetoEvent;
      }

      // Compute P·q and W^2
      const double Pdotq_lab = m.P.dot(m.q);
      if (!std::isfinite(Pdotq_lab) || Pdotq_lab == 0) {
        _nBadPdotq++;
        vetoEvent;
      }

      const double W2 = (m.P + m.q).mass2();

      // Event-level DIS cuts (HERMES-like)
      if (m.Q2 <= 1.0) { _nEventsVetoDIS++; vetoEvent; }
      if (!std::isfinite(W2) || W2 <= 10.0) { _nEventsVetoDIS++; vetoEvent; }
      if (!(m.y > 0.10 && m.y < 0.85)) { _nEventsVetoDIS++; vetoEvent; }
      if (!(m.xB > 0.023 && m.xB < 0.60)) { _nEventsVetoDIS++; vetoEvent; }
      _nEventsAfterDIS++;

      static bool printed = false;

      for (const Particle& p : fs.particles()) {
        _nPartsTotal++;

        const int pid = p.pid();
        if (PID::isLepton(pid)) continue;
        if (pid == 22) continue; // gamma
        _nAfterPid++;

        // Select only requested hadron species
        Histo1DPtr* h = histForPid_(pid);
        if (h == nullptr) continue;
        _nAfterSpecies++;

        // Hadron momentum cut (LAB)
        const double Ph_lab = p.momentum().p3().mod();
        if (!std::isfinite(Ph_lab) || Ph_lab < 0.0 || Ph_lab > 15.0) {
          _nHadVetoPhMom++;
          continue;
        }

        // Vector-meson subtraction: ONLY for pions/kaons, and ONLY if immediate parent is VM/K*
        // {
        //   const int apid = std::abs(pid);
        //   const bool isPiK = (apid == 211 || apid == 321);
        //   if (isPiK && isDirectFromVectorMesonOrKstar_HepMC3(p)) {
        //     _nVetoVM++;
        //     continue;
        //   }
        // }

        // Invariant z_h (LAB)
        const double zh = (m.P.dot(p.momentum())) / Pdotq_lab;
        if (std::isfinite(zh)) {
          _nZhFinite++;
          _maxZhSeen = std::max(_maxZhSeen, zh);
        }
        if (!std::isfinite(zh) || zh <= 0.0 || zh >= 1.0) continue;
        _nZhInRange++;

        // One-time frame diagnostics
        if (!printed) {
          printFrameDiagnostics_(m, p.momentum(), _frame);
          printed = true;
        }

        // Spectator nucleon veto (gamma*-N CM)
        if (_vetoSpectators && isSpectatorNucleon_(pid, p.momentum(), m.P, m.q)) {
          _nSpectatorVeto++;
          continue;
        }

        // Optional pT^2 window (computed in chosen frame)
        if (_usePt2Min || _usePt2Max) {
          FourMomentum ph = p.momentum();
          FourMomentum q  = m.q;
          FourMomentum P  = m.P;

          toFrame(_frame, ph, q, P);

          const double pt2 = pT2_wrt_q(ph, q);
          if (std::isfinite(pt2)) {
            _nPt2Finite++;
            _maxPt2Seen = std::max(_maxPt2Seen, pt2);
          }
          if (!std::isfinite(pt2) || pt2 < 0.0) continue;

          if (_usePt2Min && pt2 < _pt2Min) continue;
          if (_usePt2Max && pt2 >= _pt2Max) continue;

          _nPt2InWindow++;
        }

        // Fill dN/dz_h per species
        (*h)->fill(zh, 1.0);
        _nFilledTotal++;

        bumpSpeciesCounter_(pid);
      }
    }


    void finalize() override {
      MSG_INFO("==== Debug summary ====");
      MSG_INFO("Events seen:            " << _nEventsSeen);
      MSG_INFO("Events w/ metadata:     " << _nEventsWithMeta);
      MSG_INFO("Events w/o metadata:    " << _nEventsNoMeta);
      MSG_INFO("Bad/zero P·q:           " << _nBadPdotq);
      MSG_INFO("Events vetoed by DIS:   " << _nEventsVetoDIS);
      MSG_INFO("Events after DIS cuts:  " << _nEventsAfterDIS);
      MSG_INFO("Particles total:        " << _nPartsTotal);
      MSG_INFO("After PID cuts:         " << _nAfterPid);
      MSG_INFO("After species filter:   " << _nAfterSpecies);
      MSG_INFO("Hadrons vetoed by Ph:   " << _nHadVetoPhMom);
      MSG_INFO("z_h finite:             " << _nZhFinite);
      MSG_INFO("0<z_h<1:                " << _nZhInRange);
      if (_vetoSpectators) {
        MSG_INFO("Spectator vetoed:       " << _nSpectatorVeto);
      }
      if (_usePt2Min || _usePt2Max) {
        MSG_INFO("pT2 finite:             " << _nPt2Finite);
        MSG_INFO("pT2 in window:          " << _nPt2InWindow);
        MSG_INFO("Max pT2 seen:           " << _maxPt2Seen);
      }
      MSG_INFO("Filled (total):         " << _nFilledTotal);
      MSG_INFO("Filled pi+:             " << _nFilled_pip);
      MSG_INFO("Filled pi-:             " << _nFilled_pim);
      MSG_INFO("Filled K+:              " << _nFilled_kp);
      MSG_INFO("Filled K-:              " << _nFilled_km);
      MSG_INFO("Filled p:               " << _nFilled_p);
      MSG_INFO("Filled pbar:            " << _nFilled_pbar);
      MSG_INFO("Max z_h seen:           " << _maxZhSeen);
      MSG_INFO("=======================");

      // Normalize: per Rivet-kept event (replicas count as separate events here)
      const double norm = (sumW() > 0) ? 1.0 / sumW() : 1.0;

      normalizeToDensity_(_h_dN_dzh_pip,  norm);
      normalizeToDensity_(_h_dN_dzh_pim,  norm);
      normalizeToDensity_(_h_dN_dzh_kp,   norm);
      normalizeToDensity_(_h_dN_dzh_km,   norm);
      normalizeToDensity_(_h_dN_dzh_p,    norm);
      normalizeToDensity_(_h_dN_dzh_pbar, norm);
    }

  private:
    Histo1DPtr* histForPid_(int pid) {
      switch (pid) {
        case  211:  return &_h_dN_dzh_pip;   // pi+
        case -211:  return &_h_dN_dzh_pim;   // pi-
        case  321:  return &_h_dN_dzh_kp;    // K+
        case -321:  return &_h_dN_dzh_km;    // K-
        case  2212: return &_h_dN_dzh_p;     // p
        case -2212: return &_h_dN_dzh_pbar;  // anti-p
        default:    return nullptr;
      }
    }

    void bumpSpeciesCounter_(int pid) {
      switch (pid) {
        case  211:  _nFilled_pip++;  break;
        case -211:  _nFilled_pim++;  break;
        case  321:  _nFilled_kp++;   break;
        case -321:  _nFilled_km++;   break;
        case  2212: _nFilled_p++;    break;
        case -2212: _nFilled_pbar++; break;
        default: break;
      }
    }

    void normalizeToDensity_(Histo1DPtr& h, double norm) {
      if (!h) return;
      scale(h, norm);
      for (auto& bin : h->bins()) {
        const double w = bin.xWidth();
        if (w > 0) bin.scaleW(1.0 / w);
      }
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

    void printFrameDiagnostics_(const MetaDIS& m,
                               const FourMomentum& ph_lab,
                               FrameChoice frameRequested) {
      FourMomentum P = m.P;
      FourMomentum q = m.q;
      FourMomentum ph = ph_lab;

      auto p3mag = [](const FourMomentum& x){ return x.p3().mod(); };
      auto betaVec = [](const FourMomentum& x){
        const double E = x.E();
        return (E != 0.0) ? (x.p3() / E) : Vector3(0,0,0);
      };

      const double q2_lab    = q.mass2();
      const double ph2_lab   = ph.mass2();
      const double P2_lab    = P.mass2();
      const double Pdotq_lab = P.dot(q);

      MSG_INFO("========== Frame diagnostics (requested=" << frameName(frameRequested) << ") ==========");

      auto dump = [&](const std::string& tag,
                      const FourMomentum& Pnow,
                      const FourMomentum& qnow,
                      const FourMomentum& phnow) {

        const double q2  = qnow.mass2();
        const double ph2 = phnow.mass2();
        const double P2  = Pnow.mass2();
        const double Pdotq = Pnow.dot(qnow);

        const Vector3 betaP = betaVec(Pnow);
        const double betaPmag = betaP.mod();

        MSG_INFO("[" << tag << "]");
        MSG_INFO("  P:  E=" << Pnow.E()  << " |p|=" << p3mag(Pnow)
                 << "  m2=" << P2
                 << "  beta=" << betaPmag);
        MSG_INFO("  q:  E=" << qnow.E()  << " |q|=" << p3mag(qnow)
                 << "  q2=" << q2);
        MSG_INFO("  ph: E=" << phnow.E() << " |p|=" << p3mag(phnow)
                 << "  m2=" << ph2);

        MSG_INFO("  Invariant deltas: "
                 << "d(q2)=" << (q2 - q2_lab)
                 << "  d(P2)=" << (P2 - P2_lab)
                 << "  d(ph2)=" << (ph2 - ph2_lab)
                 << "  d(P·q)=" << (Pdotq - Pdotq_lab));

        MSG_INFO("  TRF check: |P.p3|=" << p3mag(Pnow)
                 << "  (should be ~0 in TRF)");

        const double Q = (q2 < 0) ? std::sqrt(-q2) : NAN;
        MSG_INFO("  Breit check: q.E=" << qnow.E()
                 << "  Q=sqrt(-q2)=" << Q
                 << "  |q|=" << p3mag(qnow)
                 << "  (in BREIT: q.E~0 and |q|~Q)");

        const double qmag = p3mag(qnow);
        const double betaBreit = (qmag > 0) ? (qnow.E() / qmag) : NAN;
        MSG_INFO("  BetaBreit candidate (q.E/|q|)=" << betaBreit
                 << "  (DIS should satisfy |beta|<1)");
      };

      dump("LAB", P, q, ph);

      {
        if (P.E() != 0.0) {
          const Vector3 bTRF = -(P.p3() / P.E());
          boostBy(P,  bTRF);
          boostBy(q,  bTRF);
          boostBy(ph, bTRF);
          dump("TRF", P, q, ph);
        } else {
          MSG_INFO("[TRF] skipped: P.E=0");
        }
      }

      {
        const Vector3 qv = q.p3();
        const double qmag = qv.mod();
        if (qmag > 0) {
          const Vector3 nhat = qv / qmag;
          const double beta = q.E() / qmag;
          const Vector3 bBreit = beta * nhat;

          boostBy(P,  bBreit);
          boostBy(q,  bBreit);
          boostBy(ph, bBreit);
          dump("BREIT", P, q, ph);
        } else {
          MSG_INFO("[BREIT] skipped: |q|=0 after TRF");
        }
      }

      MSG_INFO("=====================================================");
    }

  private:
    // Config
    std::string _metafile;
    FrameChoice _frame = FrameChoice::TRF;

    // Path-based inference
    std::string _hepmcPath;
    int _evtFromPath = -1;

    bool _usePt2Min = false;
    bool _usePt2Max = false;
    double _pt2Min = 0.0;
    double _pt2Max = 0.0;

    // Spectator veto config
    bool _vetoSpectators = true;
    double _spectatorPmax = 0.30;  // GeV
    double _spectatorTkmax = -1.0; // GeV (disabled if <0)

    // Data
    std::unordered_map<int, MetaDIS> _meta;

    // Species histograms
    Histo1DPtr _h_dN_dzh_pip;
    Histo1DPtr _h_dN_dzh_pim;
    Histo1DPtr _h_dN_dzh_kp;
    Histo1DPtr _h_dN_dzh_km;
    Histo1DPtr _h_dN_dzh_p;
    Histo1DPtr _h_dN_dzh_pbar;

    // Debug counters
    size_t _nEventsSeen = 0;
    size_t _nEventsWithMeta = 0;
    size_t _nEventsNoMeta = 0;
    size_t _nBadPdotq = 0;

    size_t _nEventsAfterDIS = 0;
    size_t _nEventsVetoDIS  = 0;

    size_t _nPartsTotal = 0;
    size_t _nAfterPid = 0;
    size_t _nAfterSpecies = 0;

    size_t _nHadVetoPhMom = 0;
    size_t _nVetoVM = 0;

    size_t _nZhFinite = 0;
    size_t _nZhInRange = 0;

    size_t _nPt2Finite = 0;
    size_t _nPt2InWindow = 0;

    size_t _nSpectatorVeto = 0;

    size_t _nFilledTotal = 0;
    size_t _nFilled_pip  = 0;
    size_t _nFilled_pim  = 0;
    size_t _nFilled_kp   = 0;
    size_t _nFilled_km   = 0;
    size_t _nFilled_p    = 0;
    size_t _nFilled_pbar = 0;

    double _maxZhSeen = 0.0;
    double _maxPt2Seen = 0.0;
  };

  RIVET_DECLARE_PLUGIN(EHIJING_SMASH_2026_DNDZ);

} // namespace Rivet
