// EHIJING_SMASH_DNDXBDZ.cc
//
// Rivet analysis to produce (1/Nevt) dN/(dx_B dz_h) vs x_B in fixed z_h slices
// from SMASH final states + per-event DIS metadata.
//
// Species: pi+, pi-, K+, K-
// z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)
//
// Config via env vars:
//   RIVET_METAFILE         (required)
//
// Output histograms (Histo1D):
//   /EHIJING_SMASH_DNDXBDZ/dN_dxbdz_pip_zh0p2_0p3
//   /EHIJING_SMASH_DNDXBDZ/dN_dxbdz_pip_zh0p3_0p4
//   /EHIJING_SMASH_DNDXBDZ/dN_dxbdz_pip_zh0p4_0p6
//   /EHIJING_SMASH_DNDXBDZ/dN_dxbdz_pip_zh0p6_0p8
//   ... similarly for pim, kp, km
//
// Normalization in finalize:
//  - divide by Nev = sumW()
//  - convert to density in x_B and z_h: scale by 1/(dx_B * dz)

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
        out = static_cast<int>(v);
        return true;
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
                    inStr = false;
                    escape = false;
                }
                continue;
                } else {
                    buf.push_back(c);

                    if (escape) {
                        escape = false;
                        continue;
                    }
                    if (c == '\\') {
                        if (inStr) escape = true;
                        continue;
                    }
                    if (c == '"') {
                        inStr = !inStr;
                        continue;
                    }
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

    class EHIJING_SMASH_DNDXBDZ : public Analysis {

    public:

        RIVET_DEFAULT_ANALYSIS_CTOR(EHIJING_SMASH_DNDXBDZ);

        void init() override {
            declare(FinalState(), "FS");

            // METAFILE
            _metafile = getenv_str("RIVET_METAFILE");
            if (_metafile.empty()) _metafile = getOption<string>("METAFILE", "");
            if (_metafile.empty()) {
                throw UserError("You must provide METAFILE. Set env var RIVET_METAFILE=/path/DISKinematics.meta.jsonl");
            }

            // x_B axis config
            _xbEdges = {0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6};
            if (_xbEdges.size() < 2) {
                throw UserError("x_B bin edge list must contain at least two entries.");
            }
            for (size_t i = 1; i < _xbEdges.size(); ++i) {
                if (!std::isfinite(_xbEdges[i-1]) || !std::isfinite(_xbEdges[i]) || _xbEdges[i] <= _xbEdges[i-1]) {
                    throw UserError("x_B bin edges must be strictly increasing and finite.");
                }
            }

            // Load kinematics metadata file
            _meta = loadMeta(_metafile);

            // Book histograms: [species][zbin]
            for (size_t is = 0; is < NSPEC; ++is) {
                for (size_t iz = 0; iz < NZ; ++iz) {
                    const std::string name = "dN_dxbdz_" + std::string(specTag(is)) + "_" + std::string(zTag(iz));
                    book(_h[is][iz], name, _xbEdges);
                }
            }
            book(_hDIS, "DIS_xB", _xbEdges);

            MSG_INFO("Loaded " << _meta.size() << " metadata entries from " << _metafile);
            MSG_INFO("x_B axis edges: 0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6");
            MSG_INFO("z_h slices: (0.2,0.3), (0.3,0.4), (0.4,0.6), (0.6,0.8)");
        } // init()

        void analyze(const Event& event) override {
            _nEventsSeen++;

            const FinalState& fs = apply<FinalState>(event, "FS");
            const auto* ge = event.genEvent();
            const int evnum = ge ? ge->event_number() : -1;

            // Find metadata corresponding to the event number.
            auto it = _meta.find(evnum);
            if (it == _meta.end()) {
                _nEventsVetoed++;
                vetoEvent;
            }
            _nEventsWithMeta++;

            const MetaDIS& eventKinematics = it->second;

            // Event-level x_B
            const double xB = eventKinematics.xB;
            if (!std::isfinite(xB)) {
                _nEventsVetoed++;
                vetoEvent;
            }

            // Compute P·q from DIS kinematics metadata.
            const double Pdotq = eventKinematics.P.dot(eventKinematics.q);

            _hDIS->fill(xB, 1.0);

            for (const Particle& particle : fs.particles()) {
                const int pid = particle.pid();

                // Ignore leptons and photons
                if (PID::isLepton(pid)) {
                    continue;
                }
                if (pid == 22) {
                    continue;
                }

                // Keep only requested species
                const int is = speciesIndex(pid);
                if (is < 0) {
                    continue;
                }

                FourMomentum ph = particle.momentum();

                // z_h = P·p_h / P·q
                const double zh = (eventKinematics.P.dot(ph)) / Pdotq;
                if (!std::isfinite(zh)) {
                    continue;
                }

                // Optional cross-check in TRF
                FourMomentum ph_trf = ph;
                FourMomentum q_trf  = eventKinematics.q;
                FourMomentum P_trf  = eventKinematics.P;
                toFrame(FrameChoice::TRF, ph_trf, q_trf, P_trf);
                const double zh_trf = ph_trf.E() / q_trf.E();
                if (std::abs(zh - zh_trf) > 1e-6) {
                    MSG_WARNING("Hadron momentum fraction from P·p_h / P·q: "
                        << zh << " does not match calculation from E_h / nu: " << zh_trf);
                }

                // z_h slice selection
                const int iz = zbinIndex(zh);
                if (iz < 0) {
                    continue;
                }

                // Same hadron momentum cut as in your current analysis
                const double ph_abs = ph.p3().mod();
                if (!std::isfinite(ph_abs) || ph_abs < 2.0 || ph_abs > 15.0) {
                    continue;
                }

                // Fill x_B histogram for this z_h slice and species
                _h[is][iz]->fill(xB, 1.0);
                _nFilled++;
            }
        } // analyze()

        void finalize() override {
            MSG_INFO("============ Summary ============");
            MSG_INFO("Events seen:          " << _nEventsSeen);
            MSG_INFO("Events w/ metadata:   " << _nEventsWithMeta);
            MSG_INFO("Events vetoed:        " << _nEventsVetoed);
            MSG_INFO("Filled entries:       " << _nFilled);
            MSG_INFO("=================================");

            // Per-event normalization
            const double Nev = (sumW() > 0) ? sumW() : 1.0;
            MSG_INFO("Normalizing by Nev = " << Nev);

            // Convert to density: (1/Nev) dN/(dx_B dz)
            for (size_t is = 0; is < NSPEC; ++is) {
                for (size_t iz = 0; iz < NZ; ++iz) {
                    if (!_h[is][iz]) continue;

                    const double dz = (zEdges[iz+1] - zEdges[iz]);
                    if (!(dz > 0)) continue;

                    for (size_t ib = 0; ib < _h[is][iz]->numBins(); ++ib) {
                        const double nDIS = _hDIS->bin(ib).sumW();
                        if (nDIS > 0.0) {
                            _h[is][iz]->bin(ib).scaleW(1.0 / (nDIS * dz));
                        } else {
                            _h[is][iz]->bin(ib).scaleW(0.0);
                        }
                    }
                }
            }
        } // finalize()

    private:

        static constexpr size_t NSPEC = 4;
        static constexpr size_t NZ = 4;

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
                if (zh > zEdges[i] && zh < zEdges[i+1]) return static_cast<int>(i);
            }
            return -1;
        }

    private:
        std::string _metafile;

        std::vector<double> _xbEdges = {0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6};

        std::unordered_map<int, MetaDIS> _meta;
        std::array<std::array<Histo1DPtr, NZ>, NSPEC> _h;
        Histo1DPtr _hDIS;

        size_t _nEventsSeen = 0;
        size_t _nEventsWithMeta = 0;
        size_t _nEventsVetoed = 0;
        size_t _nFilled = 0;
    };

    constexpr double EHIJING_SMASH_DNDXBDZ::zEdges[NZ+1];

    RIVET_DECLARE_PLUGIN(EHIJING_SMASH_DNDXBDZ);

    } // namespace Rivet