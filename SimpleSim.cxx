#include <TRandom.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TH2F.h>
#include <TF3.h>
#include <TDirectory.h>
#include <TString.h>
#include <cmath>
#include <vector>
#include <iostream>

namespace pmh {
  const double hbarc   = 0.1973269804;   // GeV*fm

  // Fireball
  const double T0      = 0.300;          // GeV
  const double Tc      = 0.155;          // GeV
  const double Tfo     = 0.120;          // GeV
  const double tau0    = 0.4;            // fm/c
  const double R0      = 3.0;            // fm
  const double vT      = 0.5;            // c
  const double n       = 4.0;            // power of temperature dependence, really ok?

  const double lambdaInj = 0.5;          

  const double norm    = 700.0;
}

double FitBE(double* x, double* par){
  const double q = x[0];
  return 1.0 + par[0]*exp(-par[1]*par[1]*q*q/(pmh::hbarc*pmh::hbarc));
}

inline double SamplePhotonMomentum(double T){
  double u = gRandom->Uniform() * gRandom->Uniform() * gRandom->Uniform();
  return -T * log(u + 1e-300);   // guard the (measure-zero) u=0 case
}


TH1D* SlicedC(TH3F* sig, TH3F* ref, char proj, int nCut, const char* name){
  sig->GetXaxis()->SetRange(); sig->GetYaxis()->SetRange(); sig->GetZaxis()->SetRange();
  ref->GetXaxis()->SetRange(); ref->GetYaxis()->SetRange(); ref->GetZaxis()->SetRange();
  TH1D *s=nullptr,*r=nullptr;
  if(proj=='x'){
    sig->GetYaxis()->SetRange(1,nCut); sig->GetZaxis()->SetRange(1,nCut);
    ref->GetYaxis()->SetRange(1,nCut); ref->GetZaxis()->SetRange(1,nCut);
    s=sig->ProjectionX(Form("%s_sig",name)); r=ref->ProjectionX(Form("%s_ref",name));
  } else if(proj=='y'){
    sig->GetXaxis()->SetRange(1,nCut); sig->GetZaxis()->SetRange(1,nCut);
    ref->GetXaxis()->SetRange(1,nCut); ref->GetZaxis()->SetRange(1,nCut);
    s=sig->ProjectionY(Form("%s_sig",name)); r=ref->ProjectionY(Form("%s_ref",name));
  } else {
    sig->GetXaxis()->SetRange(1,nCut); sig->GetYaxis()->SetRange(1,nCut);
    ref->GetXaxis()->SetRange(1,nCut); ref->GetYaxis()->SetRange(1,nCut);
    s=sig->ProjectionZ(Form("%s_sig",name)); r=ref->ProjectionZ(Form("%s_ref",name));
  }
  // reset so later projections start clean
  sig->GetXaxis()->SetRange(); sig->GetYaxis()->SetRange(); sig->GetZaxis()->SetRange();
  ref->GetXaxis()->SetRange(); ref->GetYaxis()->SetRange(); ref->GetZaxis()->SetRange();
  TH1D* c=(TH1D*)s->Clone(name); c->Divide(r);
  return c;
}

void printProgress(double_t progress)
{
  int barWidth = 50;
  cout << "         [";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)       cout << "=";
    else if (i == pos) cout << ">";
    else               cout << " ";
  }
  cout << "] " << int32_t(progress * 100.0) << " %\r" << flush;
  if (progress >= 1.) cout << endl;
}

void SimpleSim(){

  const double timestep = 0.05; // fm/c
  const double tmax     = 15.0; // fm/c
  const int    nsteps   = (int)(tmax / timestep);
  const int    nEvents  = 1e3; 

  const int    nQAevents = 500; 

  //******************************************************
  // Histograms
  //******************************************************


  TH2F* hPtVsTime = new TH2F("hPtVsTime",";p_{T} (GeV);t-#tau_{0} (fm/c)",100,0.,10.,nsteps,0.,tmax);
  TH2F* hPtVsTemp = new TH2F("hPtVsTemp",";p_{T} (GeV);T (GeV)",100,0.,10.,nsteps,0.,pmh::T0);
  TH1F* hP   = new TH1F("hP",  ";p (GeV);Counts",100,0,10);
  TH1F* hPt  = new TH1F("hPt", ";p_{T} (GeV);Counts",100,0,10);
  TH1F* hEta = new TH1F("hEta",";#eta;Counts",100,-5,5);
  TH1F* hPhi = new TH1F("hPhi",";#varphi;Counts",64,0,TMath::TwoPi());

  TH1F* hTvsTime = new TH1F("hTvsTime",";t-#tau_{0} (fm/c);T (GeV)",nsteps,0.,tmax);
  TH1F* hVvsTime = new TH1F("hVvsTime",";t-#tau_{0} (fm/c);V (fm^{3})",nsteps,0.,tmax);
  TH1F* hRvsTime = new TH1F("hRvsTime",";t-#tau_{0} (fm/c);R_{T} (fm)",nsteps,0.,tmax);
  TH1F* hEtotVsTime = new TH1F("hEtotVsTime",";t-#tau_{0} (fm/c);T^{4}V",nsteps,0.,tmax);
  TH1F* hMeanPhotonsVsTime = new TH1F("hMeanPhotonsVsTime",";t-#tau_{0} (fm/c);<#gamma>/step",nsteps,0.,tmax);
  TH1F* hNPhotonsPerEvent  = new TH1F("hNPhotonsPerEvent",";N_{#gamma}/event;Events",300,0,6000);

  TH2F* hXY = new TH2F("hXY",";x (fm);y (fm)",100,-12,12,100,-12,12);
  TH1F* hR  = new TH1F("hR",";r (fm);Counts",100,0,12);
  TH1F* hZ  = new TH1F("hZ",";z (fm);Counts",100,-12,12);
  TH1F* hEmissionTime = new TH1F("hEmissionTime",";t (fm/c);Counts",nsteps,pmh::tau0,pmh::tau0+tmax);

  TH1F* hSigQinv=new TH1F("hSigQinv","C_{2} vs q_{inv};q_{inv} (GeV);C_{2}",40,0,0.2);
  TH1F* hRefQinv=new TH1F("hRefQinv","",40,0,0.2);
  TH1F* hSigQabs=new TH1F("hSigQabs","C_{2} vs |q|;|q| (GeV);C_{2}",40,0,0.2);
  TH1F* hRefQabs=new TH1F("hRefQabs","",40,0,0.2);
  TH2F* hSig2D=new TH2F("hSig2D","C_{2}(q_{inv},#DeltaE);q_{inv} (GeV);#DeltaE (GeV)",40,0,0.2,40,0,0.5);
  TH2F* hRef2D=new TH2F("hRef2D","",40,0,0.2,40,0,0.5);
  TH3F* hSig3D=new TH3F("hSig3D",";q_{out} (GeV);q_{side} (GeV);q_{long} (GeV)",
                        40,0,0.2, 40,0,0.2, 40,0,0.2);
  TH3F* hRef3D=new TH3F("hRef3D","",40,0,0.2, 40,0,0.2, 40,0,0.2);

  hSigQinv->Sumw2(); hRefQinv->Sumw2();
  hSigQabs->Sumw2(); hRefQabs->Sumw2();
  hSig2D->Sumw2();   hRef2D->Sumw2();
  hSig3D->Sumw2();   hRef3D->Sumw2();

  double gSumP=0, gSumT=0; long gN=0;
  long totalPairs=0;



  //******************************************************
  // Event loop
  //******************************************************

  for(int iEv = 0; iEv < nEvents; ++iEv){

    std::vector<double> vE,vpx,vpy,vpz,vt,vx,vy,vz;
    const bool doQA = (iEv < nQAevents);

    int    nPhotonsEvent = 0;
    double V0 = -1;

    for(int i = 0; i < nsteps; ++i){

      const double dt  = i * timestep;
      const double tau = pmh::tau0 + dt;

      const double RT     = pmh::R0 + pmh::vT * dt;
      const double RL     = RT;
      const double Volume = M_PI*RT*RT*(2.0*RL);
      if(V0 < 0) V0 = Volume;

      const double Temperature = pmh::T0 * pow(V0/Volume, 0.25);   // E ~ T^4 V = const
      if (Temperature < pmh::Tfo){
        if(iEv==0) cout << "freeze-out at step " << i << ", t=" << dt
                        << ", T=" << Temperature << endl;
        break;
      }

      const double meanPhotons = 2.0*pmh::norm*Volume*timestep*pow(Temperature, pmh::n+3.0);
      const int    nPhotons    = gRandom->Poisson(meanPhotons);

      if(iEv==0){
        hTvsTime->SetBinContent(i+1, Temperature);
        hVvsTime->SetBinContent(i+1, Volume);
        hRvsTime->SetBinContent(i+1, RT);
        hEtotVsTime->SetBinContent(i+1, pow(Temperature,4)*Volume);
        hMeanPhotonsVsTime->SetBinContent(i+1, meanPhotons);
      }

      nPhotonsEvent += nPhotons;

    //******************************************************
    // Photon loop
    //******************************************************


      for(int ip = 0; ip < nPhotons; ++ip){

        const double p   = SamplePhotonMomentum(Temperature);   // no TF1, check if it is ok
        const double ct  = gRandom->Uniform(-1,1);
        const double phi = gRandom->Uniform(0,TMath::TwoPi());
        const double st  = sqrt(1-ct*ct);
        const double E=p, px=p*st*cos(phi), py=p*st*sin(phi), pz=p*ct;

        const double rho = RT*sqrt(gRandom->Uniform());
        const double pp  = gRandom->Uniform(0,TMath::TwoPi());
        const double x=rho*cos(pp), y=rho*sin(pp), z=gRandom->Uniform(-RL,RL);

        vE.push_back(E); vpx.push_back(px); vpy.push_back(py); vpz.push_back(pz);
        vt.push_back(tau); vx.push_back(x); vy.push_back(y); vz.push_back(z);

        if(doQA){
          const double pT = sqrt(px*px+py*py);
          hPtVsTime->Fill(pT, dt); hPtVsTemp->Fill(pT, Temperature);
          hP->Fill(p); hPt->Fill(pT); hPhi->Fill(phi);
          if(fabs(p-pz)>1e-9) hEta->Fill(0.5*log((p+pz)/(p-pz)));
          hXY->Fill(x,y); hR->Fill(rho); hZ->Fill(z); hEmissionTime->Fill(tau);
        }
        gSumP += p; gSumT += Temperature; ++gN;
      }
    } // steps

    hNPhotonsPerEvent->Fill(nPhotonsEvent);

    const size_t N = vE.size();
    totalPairs += (long)N*(N-1)/2;
    for(size_t a=0;a<N;++a) for(size_t b=a+1;b<N;++b){

      const double q0=vE[a]-vE[b], qx=vpx[a]-vpx[b], qy=vpy[a]-vpy[b], qz=vpz[a]-vpz[b];
      const double dtt=vt[a]-vt[b], dx=vx[a]-vx[b], dy=vy[a]-vy[b], dz=vz[a]-vz[b];

      const double qdotx = (q0*dtt - qx*dx - qy*dy - qz*dz)/pmh::hbarc;
      const double w = 1.0 + pmh::lambdaInj*cos(qdotx);

      const double q3   = sqrt(qx*qx+qy*qy+qz*qz);
      const double qinv = sqrt(fabs(q3*q3 - q0*q0));
      const double dE   = fabs(q0);
      hSigQinv->Fill(qinv,w); hRefQinv->Fill(qinv,1.0);
      hSigQabs->Fill(q3,  w); hRefQabs->Fill(q3,  1.0);
      hSig2D ->Fill(qinv,dE,w); hRef2D ->Fill(qinv,dE,1.0);

      const double Kx=0.5*(vpx[a]+vpx[b]), Ky=0.5*(vpy[a]+vpy[b]);
      const double Kz=0.5*(vpz[a]+vpz[b]), K0=0.5*(vE[a]+vE[b]);
      const double betaz=Kz/K0, gam=1.0/sqrt(1-betaz*betaz);
      const double qlong=gam*(qz-betaz*q0);
      const double kT=sqrt(Kx*Kx+Ky*Ky);
      if(kT>1e-9){
        const double qout =(qx*Kx+qy*Ky)/kT;       // along pair k_T
        const double qside=(qy*Kx-qx*Ky)/kT;       // perpendicular
        hSig3D->Fill(fabs(qout),fabs(qside),fabs(qlong), w);
        hRef3D->Fill(fabs(qout),fabs(qside),fabs(qlong), 1.0);
      }
    }

    if(iEv%20==0 || iEv==nEvents-1) printProgress((iEv+1.0)/nEvents);
  } // events

  // CF STUff

  // TODO: Something is looking weird for the 3D correlation function. Check if the q_out, q_side, q_long calculation is correct. Maybe the boost to LCMS is not done correctly.


  TH1F* hCqinv=(TH1F*)hSigQinv->Clone("hCqinv"); hCqinv->Divide(hRefQinv);
  TH1F* hCqabs=(TH1F*)hSigQabs->Clone("hCqabs"); hCqabs->Divide(hRefQabs);
  TH2F* hC2D  =(TH2F*)hSig2D ->Clone("hC2D");    hC2D ->Divide(hRef2D);
 
  const int nCut  = 8;              
  const int nFull = 40;             // whole range = the diluted projection
  TH1D* hCout =SlicedC(hSig3D,hRef3D,'x',nCut, "hCout");   // q_side,q_long < 0.04
  TH1D* hCside=SlicedC(hSig3D,hRef3D,'y',nCut, "hCside");
  TH1D* hClong=SlicedC(hSig3D,hRef3D,'z',nCut, "hClong");
  TH1D* hCoutFull =SlicedC(hSig3D,hRef3D,'x',nFull,"hCoutFull");  // diluted, for contrast
  TH1D* hCsideFull=SlicedC(hSig3D,hRef3D,'y',nFull,"hCsideFull");
  TH1D* hClongFull=SlicedC(hSig3D,hRef3D,'z',nFull,"hClongFull");

  TF1* fBEqinv = new TF1("fBEqinv",FitBE,0,0.15,2); fBEqinv->SetParameters(0.5,3.0);
  TF1* fBEqabs = new TF1("fBEqabs",FitBE,0,0.15,2); fBEqabs->SetParameters(0.5,3.0);
  hCqinv->Fit(fBEqinv,"RQ0");
  hCqabs->Fit(fBEqabs,"RQ0");

  TFile* fout = TFile::Open("pmh.root","Recreate");

  TDirectory* dSpectra = fout->mkdir("spectra"); dSpectra->cd();
  hPtVsTime->Write(); hPtVsTemp->Write();
  hP->Write(); hPt->Write(); hEta->Write(); hPhi->Write();

  TDirectory* dFireball = fout->mkdir("fireball"); dFireball->cd();
  hTvsTime->Write(); hVvsTime->Write(); hRvsTime->Write(); hEtotVsTime->Write();
  hMeanPhotonsVsTime->Write(); hNPhotonsPerEvent->Write();

  TDirectory* dGeom = fout->mkdir("source_geometry"); dGeom->cd();
  hXY->Write(); hR->Write(); hZ->Write(); hEmissionTime->Write();

  TDirectory* dRaw = fout->mkdir("hbt_raw"); dRaw->cd();
  hSigQinv->Write(); hRefQinv->Write(); hSigQabs->Write(); hRefQabs->Write();
  hSig2D->Write();   hRef2D->Write();
  hSig3D->Write();   hRef3D->Write();

  TDirectory* dCorr = fout->mkdir("hbt_correlations"); dCorr->cd();
  hCqinv->Write(); hCqabs->Write(); hC2D->Write();
  hCout->Write();      hCside->Write();      hClong->Write();
  hCoutFull->Write();  hCsideFull->Write();  hClongFull->Write();

  fout->Close();
  cout << "wrote pmh.root" << endl;
}