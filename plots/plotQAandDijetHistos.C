void plotQAandDijetHistograms() {
   //const char* file1 = "histos_MC_procesed_gamma_proton_wtFRGcuts_OCt28_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_TOTAL.root";
   // const char* file2 = "histos_MC_procesed_QCD_wt1FRGcuts_OCt29_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0.root";
  //  const char* file3 = "histos_Data_wZDC_pPb_gamma_proton_wtFRGcuts_OCt28_pPb_Data_8160GeV_Quench_ak4PFJetAnalyzer_nojettrig_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0_TOTAL.root";
    
   //  const char* file1 =
   // "histos_MC_gam_pom_FRGless4_FRGgt2_Hfmin4_Feb18_Total_pPb_MC_8160GeV_Quench_ak4PFJetAnalyzer_nojettrig_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0_ESC_norefsample__jetshape_20250218.root";
   //  const char* file2 = "histos_MC_gam_proton_FRGles4_FRGgt8_Hfmin4_Feb18_TOTAL.root";
    //const char* file3 = "histos_DATA_gam_pom_FRGles4_FRGgt8_Hfmin4_Feb18_TOTAL.root";
    //const char* file4 = "histos_MC_QCD_NOembeded_pthat15_FRGles4_FRGgt8_Hfmin4_Feb18_pPb_MC_8160GeV_Quench_ak4PFJetAnalyzer_nojettrig_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0_ESC_norefsample__jetshape_20250218.root";
   // const char* file4 = "histos_MC_QCD_pthat15_FRGles4_FRGgt8_Hfmin4_Feb18_pPb_MC_8160GeV_Quench_ak4PFJetAnalyzer_nojettrig_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0_ESC_norefsample__jetshape_20250218.root";
    
    
   const char* file1 ="histos_MC_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";
   const char* file2 = "histos_MC_gam_pom_FRGles3_FRGaftjets_jetsAftFRG_March12_TOTAL.root";    //gam_pom MC
    const char* file3 = "histos_DATA_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";
    const char* file4 ="histos_MC_QCD_pthat15To30_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_pPb_MC_8160GeV_Qch_ak4PFJetAnalyzer_nojettrig_Jptmin_20p0_Jptmax_800p0_Jetamin_N2p5_Jetamax_2p5_ESC_tor2_3rdm0_c0p0_rem4thjet0__noref__jtsh_20250228.root";
    
    const char* qa_folder = "QA_histograms";
    const char* qa_hist_name = "hist_FRG_aft_Jetptcuts";
    
    const char* dijet_folder = "dijet_histograms";
    const char* dijet_hist_name = "hist_etaDijet_reco";
    int axes[] = {0, 1, 2, 3, 4, 5, 6, 7};
    
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);
    TFile *f4 = TFile::Open(file4);
    
    if (!f1 || !f2 || !f3) {
        std::cerr << "Error: One or more files could not be opened!" << std::endl;
        return;
    }
    
    TDirectory *qa_dir1 = (TDirectory*)f1->Get(qa_folder);
    TDirectory *qa_dir2 = (TDirectory*)f2->Get(qa_folder);
    TDirectory *qa_dir3 = (TDirectory*)f3->Get(qa_folder);
    TDirectory *qa_dir4 = (TDirectory*)f4->Get(qa_folder);
    
    TH1D *qa_hist1 = dynamic_cast<TH1D*>(qa_dir1->Get(qa_hist_name));
    TH1D *qa_hist2 = dynamic_cast<TH1D*>(qa_dir2->Get(qa_hist_name));
    TH1D *qa_hist3 = dynamic_cast<TH1D*>(qa_dir3->Get(qa_hist_name));
    TH1D *qa_hist4 = dynamic_cast<TH1D*>(qa_dir4->Get(qa_hist_name));
    
    if (!qa_hist1 || !qa_hist2 || !qa_hist3) {
        std::cerr << "Error: QA histograms could not be retrieved!" << std::endl;
        return;
    }
    
    

    // Plot QA histogram
    TCanvas *c1 = new TCanvas("c1", "QA Histogram", 800, 600);
    gPad->SetLogy();
    qa_hist1->Scale(1.2144*0.01);//original 3.9189  //5.28 da un mejor agreement pero el "real es 3.9189" gam-pom calculado con base en el nuevo #events del nuevo MC
    qa_hist2->Scale(0.013796*0.1); //0.1972 // factor from xsec from pythia gam-proton
    qa_hist3->Scale(1.0); 
    qa_hist4->Scale(248.94*0.01); //248.94*0.01 QCD embeded; 915.72 QCD no embeded
    
    // Create a sum histogram
    TH1D *qa_hist_sum = (TH1D*)qa_hist1->Clone("qa_hist_sum");
    qa_hist_sum->Add(qa_hist2);
    qa_hist_sum->Add(qa_hist4);
    
    qa_hist_sum->SetLineColor(kMagenta);
    qa_hist_sum->SetLineStyle(2);
    qa_hist_sum->SetLineWidth(2);
    
    qa_hist1->SetLineColor(kRed);
    qa_hist2->SetLineColor(kBlue);
    qa_hist3->SetMarkerStyle(kFullCircle);
    qa_hist3->SetMarkerColor(kBlack);
    qa_hist3->SetLineColor(kBlack);
    qa_hist3->SetMarkerSize(0.7);
    qa_hist4->SetLineColor(kGreen);

    
    qa_hist1->SetStats(kFALSE);
    qa_hist2->SetStats(kFALSE);
    qa_hist3->SetStats(kFALSE);
    qa_hist4->SetStats(kFALSE);
    qa_hist_sum->SetStats(kFALSE);
    
   // qa_hist1->GetYaxis()->SetRangeUser(0, 100000000); // Set ratio range
    qa_hist1->GetYaxis()->SetTitle("counts");
    qa_hist1->GetYaxis()->SetTitleSize(0.03); 
    //qa_hist1->GetYaxis()->SetTitleOffset(0.35);  // Adjust the Y-axis title offset to align it
    qa_hist1->GetXaxis()->SetTitleOffset(1.2);
    qa_hist1->GetXaxis()->SetTitle("Gap_{photon} ");
    qa_hist1->GetXaxis()->SetTitleSize(0.03);
    
   // qa_hist1->SetMaximum(qa_hist3->GetMaximum());
    qa_hist1->SetMaximum(1e4);
    qa_hist1->SetMinimum(1e1);

        
    qa_hist1->Draw("HIST");
    qa_hist2->Draw("HIST SAME");
    qa_hist4->Draw("HIST SAME");
    qa_hist_sum->Draw("HIST SAME");
    qa_hist3->Draw("P SAME");
    
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.85);
    legend->AddEntry(qa_hist1, "PYTHIA8 photo-ND", "l");
    legend->AddEntry(qa_hist2, "PYTHIA8 photo-SD", "l");
    legend->AddEntry(qa_hist4, "PYTHIA8 ND embedded", "l");
    legend->AddEntry(qa_hist_sum, "MC Total", "l");
    legend->AddEntry(qa_hist3, "Data", "p");
    legend->SetFillStyle(0);  // Set fill style to transparent
    legend->SetBorderSize(0); // Remove the border around the legend
    legend->Draw();
    
    // Add CMS Work In Progress text using TLatex
    TLatex cmsText;
cmsText.SetTextSize(0.025);  // Adjust text size 
cmsText.SetTextFont(61);  // Set font to bold for CMS
cmsText.SetNDC();  // Use normalized coordinates for positioning
cmsText.DrawLatex(0.12, 0.92, "CMS");  // Position CMS on the top left

TLatex workInProgressText;
workInProgressText.SetTextSize(0.025);  // Adjust size for the second line
workInProgressText.SetTextFont(42);  // Set font to normal (non-bold) for Work in Progress
workInProgressText.SetNDC();  // Use normalized coordinates for positioning
workInProgressText.DrawLatex(0.16, 0.92, "Work in Progress");  // Position it next to CMS

// Add pPb 8.16TeV (2016) text using TLatex
TLatex pPbText;
pPbText.SetTextFont(42);  // Normal text font
pPbText.SetTextSize(0.025);  // Adjust text size for pPb label
pPbText.SetNDC();  // Use normalized coordinates for positioning
pPbText.DrawLatex(0.75, 0.92, "pPb 8.16TeV (2016)");  // Position pPb label
    
    c1->SaveAs("QA_Histogram_FRG_sum_March13.pdf");
    delete c1;
    
    
    
    f1->Close();
    f2->Close();
    f3->Close();
    f4->Close();
}
