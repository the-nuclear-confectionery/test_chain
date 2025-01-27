// Restored and expanded version of the code with calculations for the original and new histograms
// Author: Willian M. Serenone

#include "TTree.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "smash/pdgcode.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/GenParticle.h"
#include <iostream>
#include <string>

void check_args(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Wrong parameters. Correct call:" << std::endl
                  << argv[0] << " [input-file] [output-file] [(pseudo)rapidity-cut]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void charged_spectra(const std::string& file_in_path, const TString& file_out_path, double rap_cut = 0.5) {
    TFile* file_in = TFile::Open(file_in_path.c_str());
    if (!file_in || file_in->IsZombie()) {
        std::cerr << "Error: Could not open input file " << file_in_path << std::endl;
        exit(EXIT_FAILURE);
    }

    TTree* tree = static_cast<TTree*>(file_in->Get("hepmc3_tree"));
    if (!tree) {
        std::cerr << "Error: Tree 'hepmc3_tree' not found in input file." << std::endl;
        file_in->Close();
        exit(EXIT_FAILURE);
    }

    HepMC3::GenEventData* sample = new HepMC3::GenEventData();
    tree->SetBranchAddress("hepmc3_event", &sample);

    const int nbins_pt = 120;
    const double max_pt = 6.0;
    const int nbins_rap = 54;
    const double max_rap = 5.4;
    const int nphi = 100;

    TFile* file_out = TFile::Open(file_out_path, "RECREATE");
    if (!file_out || file_out->IsZombie()) {
        std::cerr << "Error: Could not create output file " << file_out_path << std::endl;
        file_in->Close();
        exit(EXIT_FAILURE);
    }

    //Qvectors hist
    const int n_qvec_max = 10;
    // Histogram definitions for real and imaginary Q-vectors for charged particles
    std::vector<TH1D> hReQ_charged;
    std::vector<TH1D> hImQ_charged;

    // Histogram definitions for real and imaginary Q-vectors for identified particles
    std::vector<TH1D> hReQ_pions_minus;
    std::vector<TH1D> hImQ_pions_minus;
    std::vector<TH1D> hReQ_pions_plus;
    std::vector<TH1D> hImQ_pions_plus;
    std::vector<TH1D> hReQ_protons_minus;
    std::vector<TH1D> hImQ_protons_minus;
    std::vector<TH1D> hReQ_protons_plus;
    std::vector<TH1D> hImQ_protons_plus;
    std::vector<TH1D> hReQ_lambda_plus;
    std::vector<TH1D> hImQ_lambda_plus;
    std::vector<TH1D> hReQ_lambda_minus;
    std::vector<TH1D> hImQ_lambda_minus;
    std::vector<TH1D> hReQ_sigma_plus;
    std::vector<TH1D> hImQ_sigma_plus;
    std::vector<TH1D> hReQ_sigma_minus;
    std::vector<TH1D> hImQ_sigma_minus;
    std::vector<TH1D> hReQ_omega;
    std::vector<TH1D> hImQ_omega;


    // Histograms for charged particles
    TH1D hpt_charged("hpt_charged", "p_{T} spectra - charged particles", nbins_pt, 0, max_pt);
    TH1D hphi_charged("hphi_charged", "#varphi spectra - charged particles", nphi, -M_PI, M_PI);
    TH1D heta_charged("heta_charged", "pseudorapidity spectra - charged particles", nbins_rap, -max_rap, max_rap);

    // Histograms for identified particles
    TH1D hpt_pions_minus("hpt_pions_minus", "p_{T} spectra - #pi^{-}", nbins_pt, 0, max_pt);
    TH1D hphi_pions_minus("hphi_pions_minus", "#varphi spectra - #pi^{-}", nphi, -M_PI, M_PI);
    TH1D hy_pions_minus("hy_pions_minus", "dN/dy - #pi^{-}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_pions_plus("hpt_pions_plus", "p_{T} spectra - #pi^{+}", nbins_pt, 0, max_pt);
    TH1D hphi_pions_plus("hphi_pions_plus", "#varphi spectra - #pi^{+}", nphi, -M_PI, M_PI);
    TH1D hy_pions_plus("hy_pions_plus", "dN/dy - #pi^{+}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_protons_minus("hpt_protons_minus", "p_{T} spectra - p^{-}", nbins_pt, 0, max_pt);
    TH1D hphi_protons_minus("hphi_protons_minus", "#varphi spectra - p^{-}", nphi, -M_PI, M_PI);
    TH1D hy_protons_minus("hy_protons_minus", "pseudorapidity spectra - p^{-}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_protons_plus("hpt_protons_plus", "p_{T} spectra - p^{+}", nbins_pt, 0, max_pt);
    TH1D hphi_protons_plus("hphi_protons_plus", "#varphi spectra - p^{+}", nphi, -M_PI, M_PI);
    TH1D hy_protons_plus("hy_protons_plus", "pseudorapidity spectra - p^{+}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_lambda_plus("hpt_lambda_plus", "p_{T} spectra - Lambda^{+}", nbins_pt, 0, max_pt);
    TH1D hphi_lambda_plus("hphi_lambda_plus", "#varphi spectra - Lambda^{+}", nphi, -M_PI, M_PI);
    TH1D hy_lambda_plus("hy_lambda_plus", "pseudorapidity spectra - Lambda^{+}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_lambda_minus("hpt_lambda_minus", "p_{T} spectra - Lambda^{-}", nbins_pt, 0, max_pt);
    TH1D hphi_lambda_minus("hphi_lambda_minus", "#varphi spectra - Lambda^{-}", nphi, -M_PI, M_PI);
    TH1D hy_lambda_minus("hy_lambda_minus", "pseudorapidity spectra - Lambda^{-}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_sigma_plus("hpt_sigma_plus", "p_{T} spectra - Sigma^{+}", nbins_pt, 0, max_pt);
    TH1D hphi_sigma_plus("hphi_sigma_plus", "#varphi spectra - Sigma^{+}", nphi, -M_PI, M_PI);
    TH1D hy_sigma_plus("hy_sigma_plus", "pseudorapidity spectra - Sigma^{+}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_sigma_minus("hpt_sigma_minus", "p_{T} spectra - Sigma^{-}", nbins_pt, 0, max_pt);
    TH1D hphi_sigma_minus("hphi_sigma_minus", "#varphi spectra - Sigma^{-}", nphi, -M_PI, M_PI);
    TH1D hy_sigma_minus("hy_sigma_minus", "pseudorapidity spectra - Sigma^{-}", nbins_rap, -max_rap, max_rap);

    TH1D hpt_omega("hpt_omega", "p_{T} spectra - Omega", nbins_pt, 0, max_pt);
    TH1D hphi_omega("hphi_omega", "#varphi spectra - Omega", nphi, -M_PI, M_PI);
    TH1D hy_omega("hy_omega", "pseudorapidity spectra - Omega", nbins_rap, -max_rap, max_rap);
    for (int in = 1; in <= n_qvec_max; ++in) {
        hReQ_charged.push_back(TH1D(TString::Format("ReQ_charged_%d", in), TString::Format("Real part of Q%d - charged particles", in), 100, -5.0, 5.0));
        hImQ_charged.push_back(TH1D(TString::Format("ImQ_charged_%d", in), TString::Format("Imag part of Q%d - charged particles", in), 100, -5.0, 5.0));

        hReQ_pions_minus.push_back(TH1D(TString::Format("ReQ_pions_minus_%d", in), TString::Format("Real part of Q%d - #pi^{-}", in), 100, -5.0, 5.0));
        hImQ_pions_minus.push_back(TH1D(TString::Format("ImQ_pions_minus_%d", in), TString::Format("Imag part of Q%d - #pi^{-}", in), 100, -5.0, 5.0));

        hReQ_pions_plus.push_back(TH1D(TString::Format("ReQ_pions_plus_%d", in), TString::Format("Real part of Q%d - #pi^{+}", in), 100, -5.0, 5.0));
        hImQ_pions_plus.push_back(TH1D(TString::Format("ImQ_pions_plus_%d", in), TString::Format("Imag part of Q%d - #pi^{+}", in), 100, -5.0, 5.0));

        hReQ_protons_minus.push_back(TH1D(TString::Format("ReQ_protons_minus_%d", in), TString::Format("Real part of Q%d - p^{-}", in), 100, -5.0, 5.0));
        hImQ_protons_minus.push_back(TH1D(TString::Format("ImQ_protons_minus_%d", in), TString::Format("Imag part of Q%d - p^{-}", in), 100, -5.0, 5.0));

        hReQ_protons_plus.push_back(TH1D(TString::Format("ReQ_protons_plus_%d", in), TString::Format("Real part of Q%d - p^{+}", in), 100, -5.0, 5.0));
        hImQ_protons_plus.push_back(TH1D(TString::Format("ImQ_protons_plus_%d", in), TString::Format("Imag part of Q%d - p^{+}", in), 100, -5.0, 5.0));

        hReQ_lambda_plus.push_back(TH1D(TString::Format("ReQ_lambda_plus_%d", in), TString::Format("Real part of Q%d - Lambda^{+}", in), 100, -5.0, 5.0));
        hImQ_lambda_plus.push_back(TH1D(TString::Format("ImQ_lambda_plus_%d", in), TString::Format("Imag part of Q%d - Lambda^{+}", in), 100, -5.0, 5.0));

        hReQ_lambda_minus.push_back(TH1D(TString::Format("ReQ_lambda_minus_%d", in), TString::Format("Real part of Q%d - Lambda^{-}", in), 100, -5.0, 5.0));
        hImQ_lambda_minus.push_back(TH1D(TString::Format("ImQ_lambda_minus_%d", in), TString::Format("Imag part of Q%d - Lambda^{-}", in), 100, -5.0, 5.0));

        hReQ_sigma_plus.push_back(TH1D(TString::Format("ReQ_sigma_plus_%d", in), TString::Format("Real part of Q%d - Sigma^{+}", in), 100, -5.0, 5.0));
        hImQ_sigma_plus.push_back(TH1D(TString::Format("ImQ_sigma_plus_%d", in), TString::Format("Imag part of Q%d - Sigma^{+}", in), 100, -5.0, 5.0));

        hReQ_sigma_minus.push_back(TH1D(TString::Format("ReQ_sigma_minus_%d", in), TString::Format("Real part of Q%d - Sigma^{-}", in), 100, -5.0, 5.0));
        hImQ_sigma_minus.push_back(TH1D(TString::Format("ImQ_sigma_minus_%d", in), TString::Format("Imag part of Q%d - Sigma^{-}", in), 100, -5.0, 5.0));

        hReQ_omega.push_back(TH1D(TString::Format("ReQ_omega_%d", in), TString::Format("Real part of Q%d - Omega", in), 100, -5.0, 5.0));
        hImQ_omega.push_back(TH1D(TString::Format("ImQ_omega_%d", in), TString::Format("Imag part of Q%d - Omega", in), 100, -5.0, 5.0));
    }


    int nsamples = tree->GetEntries();
    std::cout << "Number of samples: " << nsamples << std::endl;
    double dn_deta = 0.0;
    double ncharged = 0.0;
    double n_total = 0.0;
    double mean_pt = 0.0;
    for (int isample = 0; isample < nsamples; ++isample) {
        tree->GetEntry(isample);
    
        for (const auto& part : sample->particles) {
            int id = part.pid;
            smash::PdgCode pdgcode(std::to_string(id));
            if (pdgcode == smash::PdgCode::invalid()) continue;
    
            double eta = part.momentum.eta();
            double pt = part.momentum.pt();
            double phi = part.momentum.phi();
    
            if (phi > M_PI) phi -= 2 * M_PI;
            if (phi < -M_PI) phi += 2 * M_PI;
            n_total += 1.0;
            // Charged particle analysis
            bool is_charged = fabs(pdgcode.charge()) > 1.E-4;
            if (is_charged) {
                heta_charged.Fill(eta);
                if (fabs(eta) < rap_cut) {
                    hphi_charged.Fill(phi);
                    hpt_charged.Fill(pt);
                    dn_deta += 1.0;
                    mean_pt += pt;
                }
                ncharged += 1.0;

            }
    
            // Identified particle analysis
            switch (id) {
                case -211:  // Pion-
                    hy_pions_minus.Fill(eta);
                    hpt_pions_minus.Fill(pt);
                    hphi_pions_minus.Fill(phi);
                    break;
                case 211:  // Pion+
                    hy_pions_plus.Fill(eta);
                    hpt_pions_plus.Fill(pt);
                    hphi_pions_plus.Fill(phi);
                    break;
                case 2212:  // Proton+
                    hy_protons_plus.Fill(eta);
                    hpt_protons_plus.Fill(pt);
                    hphi_protons_plus.Fill(phi);
                    break;
                case -2212:  // Proton-
                    hy_protons_minus.Fill(eta);
                    hpt_protons_minus.Fill(pt);
                    hphi_protons_minus.Fill(phi);
                    break;
                case 3122:  // Lambda+
                    hy_lambda_plus.Fill(eta);
                    hpt_lambda_plus.Fill(pt);
                    hphi_lambda_plus.Fill(phi);
                    break;
                case -3122:  // Lambda-
                    hy_lambda_minus.Fill(eta);
                    hpt_lambda_minus.Fill(pt);
                    hphi_lambda_minus.Fill(phi);
                    break;
                case 3222:  // Sigma+
                    hy_sigma_plus.Fill(eta);
                    hpt_sigma_plus.Fill(pt);
                    hphi_sigma_plus.Fill(phi);
                    break;
                case 3112:  // Sigma-
                    hy_sigma_minus.Fill(eta);
                    hpt_sigma_minus.Fill(pt);
                    hphi_sigma_minus.Fill(phi);
                    break;
                case 3334:  // Omega
                    hy_omega.Fill(eta);
                    hpt_omega.Fill(pt);
                    hphi_omega.Fill(phi);
                    break;
                default:
                    break;
            }
            for (int in = 1; in <= n_qvec_max; ++in) {
                if (fabs(pdgcode.charge()) > 1.E-4) {
                    hReQ_charged[in - 1].Fill(eta, cos(in * phi));
                    hImQ_charged[in - 1].Fill(eta, sin(in * phi));
                }
        
                switch (id) {
                    case -211:  // Pion-
                        hReQ_pions_minus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_pions_minus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 211:  // Pion+
                        hReQ_pions_plus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_pions_plus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 2212:  // Proton+
                        hReQ_protons_plus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_protons_plus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case -2212:  // Proton-
                        hReQ_protons_minus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_protons_minus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 3122:  // Lambda+
                        hReQ_lambda_plus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_lambda_plus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case -3122:  // Lambda-
                        hReQ_lambda_minus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_lambda_minus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 3222:  // Sigma+
                        hReQ_sigma_plus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_sigma_plus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 3112:  // Sigma-
                        hReQ_sigma_minus[in - 1].Fill(eta, cos(in * phi));
                        hImQ_sigma_minus[in - 1].Fill(eta, sin(in * phi));
                        break;
                    case 3334:  // Omega
                        hReQ_omega[in - 1].Fill(eta, cos(in * phi));
                        hImQ_omega[in - 1].Fill(eta, sin(in * phi));
                        break;
                    default:
                        break;
                }
            }
        }
    }
    
    file_out->Write();
    file_out->Close();
    file_in->Close();
    std::cout << "Mean pT: " << mean_pt  / nsamples << std::endl;
    std::cout << "dN/deta: " << dn_deta / nsamples << std::endl;
    std::cout << "Ncharged: " << ncharged / nsamples << std::endl;
    std::cout << "N: " << n_total / nsamples << std::endl;
}       


int main(int argc, char** argv)
{
    check_args(argc, argv);
    charged_spectra(argv[1],argv[2],std::stod(argv[3]));
}


























