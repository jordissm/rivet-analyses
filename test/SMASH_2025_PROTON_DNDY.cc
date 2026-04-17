// -*- C++ -*-
// SMASH_2025_PROTON_DNDY
// Rivet analysis: dN/dy of final-state protons as a function of rapidity

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  /// dN/dy of protons
  class SMASH_2025_PROTON_DNDY : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SMASH_2025_PROTON_DNDY);

    /// Book projections and histograms
    void init() override {
      // Final-state particles (stable)
      const FinalState fs;
      declare(fs, "FS");

      // Histogram for dN/dy of protons
      const int    N_BINS = 40;
      const double Y_MIN  = -4.12;
      const double Y_MAX  =  3.88;

      // Rivet 3 style: use book(_h, ...)
      book(_h_dN_dy_p, "dN_dy_p", N_BINS, Y_MIN, Y_MAX);
    }

    /// Per-event analysis
    void analyze(const Event& event) override {
      const FinalState& fs = apply<FinalState>(event, "FS");

      // Final-state **protons only** (PDG ID = 2212)
      for (const Particle& p : fs.particles(Cuts::pid == 2212)) {
        const double y = p.rapidity();   // true rapidity
        _h_dN_dy_p->fill(y);            // Rivet 3: no explicit weight
      }
    }

    /// Normalize to (1 / N_events) dN/dy
    void finalize() override {
      if (sumOfWeights() > 0.0) {
        scale(_h_dN_dy_p, 1.0 / sumOfWeights());
      }
    }

  private:
    Histo1DPtr _h_dN_dy_p;
  };

  // Plugin hook
  RIVET_DECLARE_PLUGIN(SMASH_2025_PROTON_DNDY);

}
