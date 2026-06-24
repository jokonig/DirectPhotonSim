
namespace pmh {
  const double hbarc   = 0.1973269804;   // GeV*fm 
 
  // Fireball
  const double T0      = 0.300;          // GeV   Initial temperature
  const double Tc      = 0.155;          // GeV   QGP hadron gas transition
  const double Tfo     = 0.120;          // GeV   Freeze-out
  const double tau0    = 0.4;            // fm/c  
  const double R0      = 3.0;            // fm    initial radius
  const double vT      = 0.5;            // c     expansion rate
  const double n       = 4.0;            // power of temperature dependence
 
}

double PhotonProb(double T, double pT){
    return pow(T, pmh::n) * exp(-pT / T);
}


void SimpleSim(){
    cout << "PhotonProb(0.3, 0.2) = " << PhotonProb(0.3, 0.2) << endl;
    cout << "PhotonProb(0.1, 0.2) = " << PhotonProb(0.1, 0.2) << endl;

    TF1* fPhotonProb = new TF1("PhotonProb", "[2]*pow([0], [1]) * exp(-x / [0])", 0, 10); // [2] is proportional to the volume! 
    fPhotonProb->SetParameters(pmh::T0, pmh::n, 100.);

    TF1 *fPoisson = new TF1("fPoisson", "TMath::Poisson(x,[0])", 0, 20);
    fPoisson->SetParameter(0, 5.0);   // mean λ = 5

    const double timestep = 0.01; // fm/c
    const double tmax     = 1.0; // fm/c
    const int    nsteps   = (int)(tmax / timestep);

    TH2F* hPtVsTime = new TH2F("hPtVsTime", "hPtVsTime", 100, 0., 10., nsteps, 0., tmax);
    TH2F* hPtVsTemp = new TH2F("hPtVsTemp", "hPtVsTemp", 100, 0., 10., nsteps, 0., pmh::T0);

    for(int iEv = 0; iEv < 100; ++iEv){
        double V0 = -1;
        for(int i = 0; i < nsteps; ++i){
            const double RT = pmh::R0 * (1 + timestep * i * pmh::vT);
            const double RL = RT; // timestep; // one could also assume expansion in longitudinal direction = c?? 

            const double Volume = M_PI*RT*RT*RL;   // of course its a cylinder
            if(V0 < 0) {
                V0 = Volume;
            }

            const double Temperature = pmh::T0 * (V0/Volume);
            cout << "Volume " << Volume << "  Temperature " << Temperature << endl;

            fPhotonProb->SetParameter(0, Temperature);
            fPhotonProb->SetParameter(2, 100*Volume); // Faktor 100 is for more statistics

            fPoisson->SetParameter(0, fPhotonProb->Integral(0, 10));
            const int nPhotons = static_cast<int>(fPoisson->GetRandom()); // this needs to have some randomness and transformed to int
            // cout << "nPhotons " << nPhotons << endl;

            for(int ip = 0; ip < nPhotons; ++ip){
                const double pTPhoton = fPhotonProb->GetRandom();
                cout << "Wow, a photon (nPhotons = " << nPhotons << ")"  << endl;
                hPtVsTime->Fill(pTPhoton, i*timestep);
                hPtVsTemp->Fill(pTPhoton, Temperature);
            }

            // cout << "fPhotonProb " << fPhotonProb->Integral(0, 10) << endl;
        }
    }

    TFile *fout = TFile::Open("pmh.root", "Recreate");
    fout->cd();
    hPtVsTime->Write();
    hPtVsTemp->Write();
    fout->Close();
}