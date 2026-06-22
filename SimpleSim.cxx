
namespace pmh {
  const double hbarc   = 0.1973269804;   // GeV*fm 
 
  // Fireball
  const double T0      = 0.300;          // GeV   Initial temperature
  const double Tc      = 0.155;          // GeV   QGP hadron gas transition
  const double Tfo     = 0.120;          // GeV   Freeze-out
  const double tau0    = 0.4;            // fm/c  
  const double R0      = 3.0;            // fm    initial radius
  const double vT      = 0.5;            // c     expansion rate
 
}


void SimpleSim(){
    for(int i = 0; i < 100; i++){
        cout << "i = " << i << endl;
    }
}