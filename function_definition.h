#include "call_libraries.h"  // call libraries from ROOT and C++
#include "trk_efficiency_correction.h" // track efficiency correction
#include "histogram_definition_new.h" // define histograms
#include "weights.h" // weights applied

/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
int get_Ntrkoff(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr, float* chi2, unsigned char* ndof, unsigned char* nlayer, unsigned char* nhits, int* algo, float* mva){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
		int NDF = (int) ndof[ii];
		int NLayer = (int) nlayer[ii];
		int NHits = (int) nhits[ii];	
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] == false) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		double calomatching = ((pfEcal[ii]+pfHcal[ii])/cosh(eta[ii]))/pt[ii];
		
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){if(pt[ii] <= 0.4) continue;}
		
		if(col_sys=="pp" && yearofdatataking==2017){if(pt[ii] <= 0.5) continue;}
		
		if(col_sys=="XeXe" && col_energy==5440 && yearofdatataking==2017){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.15) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
		}
		
		if(col_sys=="PbPb" && col_energy==5020 && yearofdatataking==2018){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.18) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
			if(algo[ii]==6 && mva[ii]<0.98) continue;
		}
		
		Ntrk_off=Ntrk_off+1;
	
	}
	
	return Ntrk_off;
}

/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
trkeff_file: efficiency file
use_centrality: centrality or multiplicity
mult: multiplicity or centrality bun
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
phi: track phi
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
double get_Ntrkcorr(TFile *trkeff_file, bool use_centrality, int mult, TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, float *phi, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr, float* chi2, unsigned char* ndof, unsigned char* nlayer, unsigned char* nhits, int* algo, float* mva){
	double Ntrk_off = 0.0;
	for(int ii=0; ii<size; ii++){ 
		int NDF = (int) ndof[ii];
		int NLayer = (int) nlayer[ii];
		int NHits = (int) nhits[ii];	
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] == false) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		double calomatching = ((pfEcal[ii]+pfHcal[ii])/cosh(eta[ii]))/pt[ii];
		
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){if(pt[ii] <= 0.4) continue;}
		
		if(col_sys=="pp" && yearofdatataking==2017){if(pt[ii] <= 0.5) continue;}
		
		if(col_sys=="XeXe" && col_energy==5440 && yearofdatataking==2017){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.15) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
		}
		
		if(col_sys=="PbPb" && col_energy==5020 && yearofdatataking==2018){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.18) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
			if(algo[ii]==6 && mva[ii]<0.98) continue;
		}
		
		Ntrk_off = Ntrk_off + getTrkCorrWeight(trkeff_file, use_centrality, col_sys.Data(), yearofdatataking, col_energy, mult, pt[ii], eta[ii], phi[ii]);

	}
	
	return Ntrk_off;
}

/*
Find Ntrk offline for pT >= 1-> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
int get_Ntrkoff_at1(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr, float* chi2, unsigned char* ndof, unsigned char* nlayer, unsigned char* nhits, int* algo, float* mva){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
		int NDF = (int) ndof[ii];
		int NLayer = (int) nlayer[ii];
		int NHits = (int) nhits[ii];	
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] == false) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		double calomatching = ((pfEcal[ii]+pfHcal[ii])/cosh(eta[ii]))/pt[ii];
		
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){if(pt[ii] < 1.0) continue;}
		
		if(col_sys=="pp" && yearofdatataking==2017){if(pt[ii] < 1.0) continue;}
		
		if(col_sys=="XeXe" && col_energy==5440 && yearofdatataking==2017){
			if(pt[ii] < 1.0) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.15) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
		}
		
		if(col_sys=="PbPb" && col_energy==5020 && yearofdatataking==2018){
			if(pt[ii] < 1.0) continue; 
			if((chi2[ii]/NDF)/NLayer >= 0.18) continue;
		 	if(NHits < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
			if(algo[ii]==6 && mva[ii]<0.98) continue;
		}
		
		Ntrk_off=Ntrk_off+1;
	
	}
	
	return Ntrk_off;
}

/*
Find number of tracks for reco -> updated for all systems (and easy to update for future systems)
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
*/
int get_simple_mult_reco(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge){
	float Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){
			if(pt[ii] <= 0.4) continue;
		}else{
			if(pt[ii] <= 0.5) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}
/*
Find number of tracks for gen -> updated for all systems (and easy to update for future systems)
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
*/
int get_simple_mult_gen(TString col_sys, int col_energy, int yearofdatataking, int size, std::vector<float> *eta, std::vector<float> *pt, std::vector<int> *charge){
	float Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){
		if(fabs(eta->at(ii)) > 2.4) continue; 
		if(fabs(charge->at(ii)) == 0)continue;
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){
			if(pt->at(ii) <= 0.4) continue;
		}else{
			if(pt->at(ii) <= 0.5) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}

/*
Calculate jet Aj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float asymmetry(float pt_leading, float pt_subleading){
	float Avariable = (pt_leading - pt_subleading)/(pt_leading + pt_subleading);
	return Avariable;
}

/*
Calculate jet Xj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float xjvar(float pt_leading, float pt_subleading){
	float XJvariable = pt_subleading/pt_leading;
	return XJvariable;
}

/*
Calculate Delta Eta
--> Arguments
eta1: eta of first object
eta2: eta of second object
*/
float deltaeta(float eta1, float eta2){
	float deltaEta = ( eta1 - eta2 );
	return deltaEta;
}

/*
Calculate Delta Phi in the range [-pi,pi]
--> Arguments
phi1: phi of first object
phi2: phi of second object
*/
float deltaphi(float phi1, float phi2){
	if(phi1 < -3.0*TMath::Pi() || phi2 < -3.0*TMath::Pi()) return -999.9;
	float deltaPhi = ( phi1 - phi2 );
	if( deltaPhi >  TMath::Pi() ) deltaPhi  += -2*TMath::Pi(); 
   	if( deltaPhi < -TMath::Pi() ) deltaPhi  +=  2*TMath::Pi();
	return deltaPhi;
}

/*
Calculate Delta Phi in the range [-pi/2 , 3/2 pi]
--> Arguments
phi1: phi of first object
phi2: phi of second object
*/
float deltaphi2PC(float phi1, float phi2){    
	if(phi1 < -3.0*TMath::Pi() || phi2 < -3.0*TMath::Pi()) return -999.9;	
	float deltaPhi = (phi1 - phi2);
	if( deltaPhi >  1.5*TMath::Pi() ) deltaPhi += -2.*TMath::Pi();
	if( deltaPhi < -0.5*TMath::Pi() ) deltaPhi +=  2.*TMath::Pi();
	return deltaPhi;
}

/*
Calculate delta R (distance)
--> Arguments
eta1: eta of first object
phi1: phi of first object
eta2: eta of second object
phi2: phi of second object
*/
float deltaR(float eta1, float phi1, float eta2, float phi2){
	float deltaR = sqrt(deltaeta(eta1,eta2)*deltaeta(eta1,eta2) + deltaphi(phi1,phi2)*deltaphi(phi1,phi2));
	return deltaR;
}

/*
Find the leading and subleading jets, return jet leading and subleading pt, eta and phi
--> Arguments
pt: jet pT
eta: jet Eta
phi: jet Phi
mass: jet Mass
flavor: jet flavour
jetindex: jet index in the loop
leadpt: leading jet pT
leadeta: leading jet Eta
leadphi: leading jet Phi
leadmass: leading jet Mass
leadflavor: leading jet flavour
leadindex: leading jet index in the loop
sublpt: subleading jet pT
subleta: subleading jet Eta
sublphi: subleading jet Phi
sublmass: subleading jet Mass
sublflavor: subleading jet flavour
sublindex: subleading jet index in the loop
thrdpt: third jet pT
thrdeta: third jet Eta
thrdphi: third jet Phi
thrdmass: third jet Mass
thrdflavor: third jet flavour
thrdindex: third jet index in the loop
*/
void find_leading_subleading_third(float pt, float eta, float phi, float mass, float flavor, int jetindex, float &leadpt, float &leadeta, float &leadphi, float &leadmass, float &leadflavor, int &leadindex, float &sublpt, float &subleta, float &sublphi, float &sublmass, float &sublflavor, int &sublindex, float &thrdpt, float &thrdeta, float &thrdphi, float &thrdmass, float &thrdflavor, int &thrdindex, float &fourpt){
    if( pt > leadpt ) {
    	fourpt = thrdpt;
        thrdpt = sublpt; 
        thrdeta = subleta;
        thrdphi = sublphi;
        thrdmass = sublmass;
        thrdflavor = sublflavor;
        thrdindex = sublindex;
        sublpt = leadpt;
        subleta = leadeta;
        sublphi = leadphi;
        sublmass = leadmass;
        sublflavor = leadflavor;
        sublindex = leadindex;
        leadpt = pt;
        leadeta = eta;
        leadphi = phi;
        leadmass = mass;
        leadflavor = flavor;
        leadindex = jetindex;
    } else if( pt > sublpt ) {
    	fourpt = thrdpt;
        thrdpt = sublpt;
        thrdeta = subleta;
        thrdphi = sublphi;
        thrdmass = sublmass;
        thrdflavor = sublflavor;
        thrdindex = sublindex;
        sublpt = pt;
        subleta = eta;
        sublphi = phi;
        sublmass = mass;
        sublflavor = flavor;
        sublindex = jetindex;
    } else if( pt > thrdpt ) {
    	fourpt = thrdpt;
        thrdpt = pt;
        thrdeta = eta;
        thrdphi = phi;
        thrdmass = mass;
        thrdflavor = flavor;
        thrdindex = jetindex;
    } else if ( pt > fourpt ) { fourpt = pt; }
}

/*
Function to compare TVector3 objects based on Pt
vec1: vector 1
vec2: vector 2
*/
bool sortByPt(const TVector3& vec1, const TVector3& vec2) { return vec1.Pt() > vec2.Pt(); } // Sorting in descending order of Pt


// Function to find the closest vector to the reference vector in phi
TVector3 findClosestVectorInPhi(const std::vector<TVector3>& sortedVector) {

    // Define the reference vector (i == 1)
    const TVector3& refVector = sortedVector[1];

    // Initialize variables to store the closest vector and its difference in phi
    TVector3 closestVector;
    double minPhiDiff = std::numeric_limits<double>::max();

    // Loop over the sorted vector of TVector3 objects
    double dphivalue = -999.9;
    int vecindex = -999;
//	cout << "===================================" << endl;   
    for (int i = 0; i < sortedVector.size(); ++i) {
        // Skip the reference vector itself and index 0
        if (i == 0 || i == 1) continue;

        // Calculate the difference in phi
        double phiDiff = std::abs(refVector.DeltaPhi(sortedVector[i]));

		//cout << "i: " << i << "; delta Phi: " << phiDiff << "; pT1st:  " << sortedVector[0].Pt() << "GeV; pT2nd: " << refVector.Pt() << " GeV; pTnext: " << sortedVector[i].Pt() << " GeV;" << endl;		

        // Update the closest vector if phi difference is smaller
        if (phiDiff < minPhiDiff) { minPhiDiff = phiDiff; closestVector = sortedVector[i]; vecindex = i;  dphivalue = phiDiff; }
    }
/*
        cout << "===================================" << endl;
        cout << "===================================" << endl;
        cout << "Closest" << endl;
        cout << "i: " << vecindex << "; delta Phi: " << dphivalue << "; pT1st:  " << sortedVector[0].Pt() << "GeV; pT2nd: " << refVector.Pt() << " GeV; pTnext: " << closestVector.Pt() << " GeV;" << endl;
        cout << "-----------------------------------" << endl;
        cout << "-----------------------------------" << endl;
*/
    return closestVector;
}


/*
Find bin dynamically 
--> Arguments
quant_vec: vector with binning
quant: variable
*/
int find_my_bin(std::vector<double> quant_vec, double quant){
	int bin = -999;
	for(int ii = 0; ii < quant_vec.size()-1; ii++) {if(quant >= quant_vec[ii] && quant < quant_vec[ii+1]){bin = ii;} }
	return bin;
}

/*
Measure the correlation between objects
--> Arguments
jets: vector with jet informations
jet_w: vector with jet weight informations
tracks: vector with track informations
tracks_w: vector with track weight informations
histo_corr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
histjet: multidimentional histogram for jets in correlations {pT, Eta, Phi}
histtrk: multidimentional histogram for tracks in correlations {pT, Eta, Phi}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
do_rotation: true means apply/fill rotation method, otherwise use false
N_rot: number of rotation (only use if do_rotation is true)
histo_rot: histogram using rotation method (only use if do_rotation is true)
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_corr_subeg0: if sube > 0 save in this histogram
JetR: jet radii defined in input_variables.h
histo_injet: in jet multiplicity
flow: true for flow measurement false for jet shapes
*/
void correlation(std::vector<TVector3> jets, std::vector<double> jets_w, std::vector<TVector3> tracks, std::vector<double> tracks_w, THnSparse* histo_corr, THnSparse* histjet, THnSparse* histtrk, float event_weight, int mult, float extra_variable, bool do_rotation, int N_rot, THnSparse* histo_rot, std::vector<int> sube_trk, THnSparse* histo_corr_subeg0, float JetR, THnSparse* histo_injet, bool flow){
	// get correlation histograms
	double multcentbin = (double) mult;
	double extrabin = (double) extra_variable;
	for (int a = 0; a < jets.size(); a++){ // start loop over jets
       	double jet_weight = jets_w[a];
       	int injettrk = 0;
		for (int b = 0; b < tracks.size(); b++){ // start loop over tracks
			double trkpt = tracks[b].Pt();
			double trketa = tracks[b].Eta();
			int subetrk = sube_trk[b];
            // track efficiency correction for reco
            double trk_weight = tracks_w[b];
           	// Find track and multiplicity bins
			double trkbin = (double) trkpt;
			// Fill jet and track quantities
			double x_jet[5]={jets[a].Pt(),jets[a].Eta(),jets[a].Phi(), (double)multcentbin,(double)extrabin}; histjet->Fill(x_jet,jet_weight*event_weight);
			double x_trk[5]={tracks[b].Pt(),tracks[b].Eta(),tracks[b].Phi(), (double)multcentbin,(double)extrabin}; histtrk->Fill(x_trk,trk_weight*event_weight);
			// Fill correlation histograms
			double del_phi = deltaphi2PC(jets[a].Phi(), tracks[b].Phi());
			double del_eta = deltaeta(jets[a].Eta(), tracks[b].Eta());
			double xhisto[5]={del_phi,del_eta,(double)trkbin, (double)multcentbin,(double)extrabin}; 
			if(flow){if(subetrk==0){histo_corr->Fill(xhisto,jet_weight*trk_weight*event_weight);}else{histo_corr_subeg0->Fill(xhisto,jet_weight*trk_weight*event_weight);}
			}else{if(subetrk==0){histo_corr->Fill(xhisto,jet_weight*trk_weight*event_weight*trkpt);}else{histo_corr_subeg0->Fill(xhisto,jet_weight*trk_weight*event_weight*trkpt);}}
			double delta_R = deltaR(jets[a].Eta(), jets[a].Phi(), tracks[b].Eta(), tracks[b].Phi());
			if(delta_R < JetR) injettrk = injettrk + 1;
		}
		double xinjet[3]={(double)injettrk, (double)multcentbin,(double)extrabin}; histo_injet->Fill(xinjet,jet_weight*event_weight);

		// get rotation histograms 
		if(do_rotation){
			for(int c = 0; c < N_rot; c++){
				float alpha = gRandom->Uniform(2.0*TMath::Pi()); // make a random number between 0 and 2pi
				float newphi = jets[a].Phi() + alpha; // add to the jet phi and reorder the axis
				if (newphi > 3.0*TMath::Pi()/2.) newphi -= 2.0*TMath::Pi();
				if (newphi < -TMath::Pi()/2.) newphi += 2.0*TMath::Pi();
				float neweta = -1.0*jets[a].Eta(); // invert the jet eta
				for (int d = 0; d < tracks.size(); d++){ // start loop over tracks using the new jet eta and phi (similar as above)
					double trkpt = tracks[d].Pt();
					double trketa = tracks[d].Eta();
           			// track efficiency correction for reco
           			double trk_weight = tracks_w[d];
            		// Find track and multiplicity bins
					double trkbin = (double) trkpt;
					// Fill correlation histograms     
					double del_phi_rot = deltaphi2PC(newphi, tracks[d].Phi());
					double del_eta_rot = deltaeta(neweta, tracks[d].Eta());
					double x_rot[5]={del_phi_rot,del_eta_rot,(double)trkbin,(double)multcentbin,(double)extrabin}; 
					if(flow){histo_rot->Fill(x_rot,jet_weight*trk_weight*event_weight);}else{histo_rot->Fill(x_rot,jet_weight*trk_weight*event_weight*trkpt);}
				}
			}
		} //end rotation
	} //end jet track loop
}

/*
Function to fill vectors to be used during the mixing (the ones with &)
--> Arguments
similar_events: true for using tracks only if event has one jet within the jet cut or falt for all tracks
nev: histogram with number of events stored to be used in the mixing
jets: vector with jet informations
jet_weight: vector with jet weight informations
tracks: vector with track informations
trk_weight: vector with track weight informations
mult: multiplicity or centrality
vertexz: Z vertex in centimeters
weight: event weight
ev_jet_vector: vector to be used in the mixing with jet information for each event
ev_jet_weight_vector: vector to be used in the mixing with jet weight information for each event
ev_track_vector: vector to be used in the mixing with track information for each event
ev_trk_weight_vector: vector to be used in the mixing with track weight information for each event
multvec: vector to be used in the mixing with event multiplicity or centrality information
vzvec: vector to be used in the mixing with event Z vertex position information
weightvec: vector to be used in the mixing with event weight information
*/
void fillvectors(bool similar_events, TH1* nev, std::vector<TVector3> jets, std::vector<double> jet_weight, std::vector<TVector3> tracks, std::vector<double> trk_weight, int mult, double vertexz, double weight, double extra_variable, std::vector<std::vector<TVector3>> &ev_jet_vector, std::vector<std::vector<double>> &ev_jet_weight_vector, std::vector<std::vector<TVector3>> &ev_track_vector, std::vector<std::vector<double>> &ev_trk_weight_vector, std::vector<int> &multvec, std::vector<double> &vzvec, std::vector<double> &weightvec, std::vector<double> &extravec){
			if(similar_events){
				if(jets.size() > 0 && tracks.size() > 0){
					nev->Fill(0);
					ev_jet_vector.push_back(jets);
					ev_track_vector.push_back(tracks);
					ev_jet_weight_vector.push_back(jet_weight);
					ev_trk_weight_vector.push_back(trk_weight);
					multvec.push_back(mult);
					vzvec.push_back(vertexz);
					weightvec.push_back(weight);
					extravec.push_back(extra_variable);
				}
			}else{
				nev->Fill(0);
				ev_jet_vector.push_back(jets);
				ev_track_vector.push_back(tracks);
				ev_jet_weight_vector.push_back(jet_weight);
				ev_trk_weight_vector.push_back(trk_weight);
				multvec.push_back(mult);
				vzvec.push_back(vertexz);
				weightvec.push_back(weight);
				extravec.push_back(extra_variable);
			}
}


/*
Measure the 2 particle correlation
--> Arguments
tracks: vector with track informations
tracks_w: vector with track weight informations
histo_2pcorr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_2pcorr_subeg0: if sube > 0 save in this histogram
histo_2pcorr_subeg0_cross: cross correlation sub == 0 and sub != 0
*/
void twoparticlecorrelation(std::vector<TVector3> tracks, std::vector<double> tracks_w, THnSparse* histo_2pcorr, float event_weight, int mult, float extra_variable, std::vector<int> sube_trk, THnSparse* histo_2pcorr_subeg0, THnSparse* histo_2pcorr_subeg0_cross){
    // Find track and multiplicity  and extra variable bins
	double multcentbin = (double) mult;
	double extrabin = (double) extra_variable;
	// get correlation histograms
	for (int a = 0; a < tracks.size(); a++){ // start loop over tracks
		double trkpt1 = tracks[a].Pt();
        double trk_weight1 = tracks_w[a];
		int subetrk1 = sube_trk[a];
		int trkbin1 = (int) find_my_bin(trk_pt_bins,trkpt1);
		for (int b = 0; b < tracks.size(); b++){ // start loop over tracks+1
			double trkpt2 = tracks[b].Pt();
			double trk_weight2 = tracks_w[b];
			int subetrk2 = sube_trk[b];
			int trkbin2 = (int) find_my_bin(trk_pt_bins,trkpt2);
			if(trkbin1 != trkbin2) continue; // only same bin to get vn as sqrt of Vn
			double trkbin = (double) trkpt1;
           	// track efficiency correction for reco
            double trk_weight = trk_weight1*trk_weight2;
			// Fill correlation histograms
			double del_phi = deltaphi2PC(tracks[a].Phi(), tracks[b].Phi());
			double del_eta = deltaeta(tracks[a].Eta(), tracks[b].Eta());
			if(del_phi == 0 && del_eta == 0 && trkpt1 == trkpt2) continue; // do not fill histograms if particles are identical
			double x2pc[5]={del_phi,del_eta,(double)trkbin,(double)multcentbin,(double)extrabin}; 
			if(subetrk1==0 && subetrk2==0){histo_2pcorr->Fill(x2pc,trk_weight*event_weight);
			}else if(subetrk1>0 && subetrk2>0){histo_2pcorr_subeg0->Fill(x2pc,trk_weight*event_weight);
			}else{histo_2pcorr_subeg0_cross->Fill(x2pc,trk_weight*event_weight);}
		} // b loop
	} // a loop
}

/*
Get area-based underlying events using rho variable (from Yi)
--> Arguments
etamin: min eta strip
etamax: max eta strip
rho: density
Jet_Eta: jet eta for each jet
Jet_R: jet cone R
*/
double GetUE(std::vector<double> *etamin, vector<double> *etamax, vector<double> *rho, double Jet_Eta, double Jet_R){

   double Result = 0;

   if(etamin == nullptr) return -1;
   if(etamax == nullptr) return -1;
   if(rho == nullptr) return -1;

   int NBin = etamin->size();
   if(NBin == 0) return -1;

   for(int i = 0; i < NBin; i++){
   
      if(etamax->at(i) < Jet_Eta - Jet_R) continue;
      if(etamin->at(i) > Jet_Eta + Jet_R) continue;

      double XMin = (std::max(etamin->at(i), Jet_Eta - Jet_R) - Jet_Eta) / Jet_R;
      double XMax = (std::min(etamax->at(i), Jet_Eta + Jet_R) - Jet_Eta) / Jet_R;

      if(XMin <= -1) XMin = -0.99999;
      if(XMax >= +1) XMax = +0.99999;

      double High = XMax * std::sqrt(1 - XMax * XMax) + std::asin(XMax);
      double Low = XMin * std::sqrt(1 - XMin * XMin) + std::asin(XMin);

      Result = Result + Jet_R * Jet_R * (High - Low) * rho->at(i);
      
   }

   return Result;
}


double TransformToUnfoldingAxis_xjptave(const double xj, const double jetPtAve){
  const int nJetPtAveBins = nPtaveBins;
  const double maxxj = 1.0;
  double transformedxj = xj;
  for(int iJetPtAve = 1; iJetPtAve < nJetPtAveBins+1; iJetPtAve++){
    if(jetPtAve >= PtaveBinsClone[iJetPtAve]){
      transformedxj += maxxj;
    } else {
      return transformedxj;
    }
  }
  // We should never reach this point. If we are here, just return error code -1
  return -1;
}

double TransformToUnfoldingAxis_pt1pt2(const double pt1, const double pt2){
  const int nJetPt2Bins = nPtLSLBins2;
  const double maxpt1 = 2000.0;
  double transformedpt1 = pt1;
  for(int iJetPt = 1; iJetPt < nJetPt2Bins+1; iJetPt++){
    if(pt2 >= PtLSLBins2[iJetPt]){
      transformedpt1 += maxpt1;
    } else {
      return transformedpt1;
    }
  }
  // We should never reach this point. If we are here, just return error code -1
  return -1;
}
