// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

  /// @brief pT spectra for π±, K±, p/pbar with |eta|<0.5; mean pT and total charged (π±+K±+p/pbar)
  class SMASH_2023_I2693474 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(SMASH_2023_I2693474);

    void init() {
      // Final-state particles in |eta| < 0.5
      // const FinalState fs(Cuts::abseta < 5.0);
      const FinalState fs;
      declare(fs, "FS");

      // pT spectra histograms in (0, 20) GeV with e.g. 40 bins
      book(_h_pt_pion["dNdpT_pion"], "dNdpT_pion",   50, 0.0, 5.0);
      book(_h_pt_kaon["dNdpT_kaon"], "dNdpT_kaon",   50, 0.0, 5.0);
      book(_h_pt_proton["dNdpT_proton"], "dNdpT_proton", 50, 0.0, 5.0);

      // Counters for species and the combined charged total (π± + K± + p/pbar)
      book(_N_pi, "N_pi");
      book(_N_K, "N_K");
      book(_N_p, "N_p"); // here: p and pbar only
      book(_N_total, "N_total"); // sum of the 3 above
    }

    void analyze(const Event& event) {
      const Particles& parts = apply<FinalState>(event, "FS").particles();

      for (const Particle& p : parts) {
        const int apid = std::abs(p.pid());
        const double pt = p.pT()/GeV;

        // π±
        if (apid == 211) {
          _h_pt_pion["dNdpT_pion"]->fill(p.pT()/GeV);
          _N_pi->fill();
          _N_total->fill();
          continue;
        }

        // K±
        if (apid == 321) {
          _h_pt_kaon["dNdpT_kaon"]->fill(p.pT()/GeV);
          _N_K->fill();
          _N_total->fill();
          continue;
        }

        // p / pbar (charged nucleons in this context)
        if (apid == 2212) {
          _h_pt_proton["dNdpT_proton"]->fill(p.pT()/GeV);
          _N_p->fill();
          _N_total->fill();
          continue;
        }
      }
    }

    void finalize() {
      // Normalize spectra to per-event yields (dN/dpT)
      const double invNev = 1.0 / sumOfWeights();
      if (sumOfWeights() > 0.0) {
        _h_pt_pion["dNdpT_pion"]->scaleW(invNev);
        _h_pt_kaon["dNdpT_kaon"]->scaleW(invNev);
        _h_pt_proton["dNdpT_proton"]->scaleW(invNev);
      }

    }

    /// @name Histograms and Counters
    /// @{
    map<string, Histo1DPtr> _h_pt_pion, _h_pt_kaon, _h_pt_proton;
    CounterPtr _N_pi, _N_K, _N_p, _N_total;
    /// @}

  };

  DECLARE_RIVET_PLUGIN(SMASH_2023_I2693474);

}
