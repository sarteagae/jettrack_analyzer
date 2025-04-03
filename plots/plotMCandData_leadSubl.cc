void plotNormalizedHistogramsWithSubleading() {
    // Leading jet files and corresponding details (no change)
   // const char* file1 = "histos_MC_procesed_gamma_proton_wtFRGcuts_OCt28_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_TOTAL.root";
   // const char* file2 = "histos_MC_procesed_QCD_wt1FRGcuts_OCt29_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0.root";
   // const char* file3 = "histos_Data_wZDC_pPb_gamma_proton_wtFRGcuts_OCt28_pPb_Data_8160GeV_Quench_ak4PFJetAnalyzer_nojettrig_Jptmin_10p0_Jptmax_8160p0_Jetamin_N1p0_Jetamax_1p0_TOTAL.root";
    
    
    const char* file1 ="histos_MC_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";
   const char* file2 = ".root";
    const char* file3 = "histos_DATA_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";

    
    const char* folder1 = "QA_histograms";
    const char* folder2 = "QA_histograms";
    const char* folder3 = "QA_histograms";
    const char* thn1_name = "hist_reco_leadjet";
    const char* thn2_name = "hist_reco_leadjet";
    const char* thn3_name = "hist_reco_leadjet";
    int axis1 = 0;
    int axis2 = 0;
    int axis3 = 0;

    // Subleading jet files and details (unchanged)
    const char* sub_file1 = "histos_MC_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";
    const char* sub_file2 = ".root";
    const char* sub_file3 = "histos_DATA_gam_proton_FRGles3_FRGaftjets_jetsAftFRGFeb28_Jptmin20_Jetamin2p5_TOTAL.root";
    const char* sub_folder = "QA_histograms";
    const char* sub_thn_name = "hist_reco_subljet";
    int sub_axis = 0;

    // Open ROOT files for leading and subleading jets
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);
    TFile *fsub1 = TFile::Open(sub_file1);
    TFile *fsub2 = TFile::Open(sub_file2);
    TFile *fsub3 = TFile::Open(sub_file3);

    if (!f1 || !f2 || !f3 || !fsub1 || !fsub2 || !fsub3) {
        std::cerr << "Error: One or more files could not be opened!" << std::endl;
        return;
    }

    // Retrieve THnSparse histograms for leading jets
    TDirectory *dir1 = (TDirectory*)f1->Get(folder1);
    TDirectory *dir2 = (TDirectory*)f2->Get(folder2);
    TDirectory *dir3 = (TDirectory*)f3->Get(folder3);
    
    THnSparseD *thn1 = dynamic_cast<THnSparseD*>(dir1->Get(thn1_name));
    THnSparseD *thn2 = dynamic_cast<THnSparseD*>(dir2->Get(thn2_name));
    THnSparseD *thn3 = dynamic_cast<THnSparseD*>(dir3->Get(thn3_name));

    // Retrieve THnSparse histograms for subleading jets
    TDirectory *sub_dir1 = (TDirectory*)fsub1->Get(sub_folder);
    TDirectory *sub_dir2 = (TDirectory*)fsub2->Get(sub_folder);
    TDirectory *sub_dir3 = (TDirectory*)fsub3->Get(sub_folder);
    
    THnSparseD *sub_thn1 = dynamic_cast<THnSparseD*>(sub_dir1->Get(sub_thn_name));
    THnSparseD *sub_thn2 = dynamic_cast<THnSparseD*>(sub_dir2->Get(sub_thn_name));
    THnSparseD *sub_thn3 = dynamic_cast<THnSparseD*>(sub_dir3->Get(sub_thn_name));

    // Check if the histograms were successfully retrieved
    if (!thn1 || !thn2 || !thn3 || !sub_thn1 || !sub_thn2 || !sub_thn3) {
        std::cerr << "Error: One or more histograms could not be retrieved!" << std::endl;
        return;
    }

    // Project the THnSparse histograms to 1D histograms
    TH1D *hist1 = thn1->Projection(axis1);
    TH1D *hist2 = thn2->Projection(axis2);
    TH1D *hist3 = thn3->Projection(axis3);

    TH1D *sub_hist1 = sub_thn1->Projection(sub_axis);
    TH1D *sub_hist2 = sub_thn2->Projection(sub_axis);
    TH1D *sub_hist3 = sub_thn3->Projection(sub_axis);

    // Normalize histograms
    hist1->Scale(1.0); //5.28 //0.028
    hist2->Scale(0.1972);
    hist3->Scale(1.0); // Optional normalization for data histogram
    sub_hist1->Scale(1.0);
    sub_hist2->Scale(0.1972);
    sub_hist3->Scale(1.0); // Optional normalization for data histogram

    // Styling for main histograms
    hist1->SetStats(kFALSE);
    hist2->SetStats(kFALSE);
    hist3->SetStats(kFALSE);
    sub_hist1->SetStats(kFALSE);
    sub_hist2->SetStats(kFALSE);
    sub_hist3->SetStats(kFALSE);
    
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);
    hist3->SetMarkerStyle(kFullCircle);
    hist3->SetMarkerColor(kBlack);
    hist3->SetLineColor(kBlack);
    hist3->SetMarkerSize(0.7);

    sub_hist1->SetLineColor(kRed);
    sub_hist1->SetLineStyle(2);
    sub_hist2->SetLineColor(kBlue);
    sub_hist2->SetLineStyle(2);
    sub_hist3->SetMarkerStyle(kOpenCircle);
    sub_hist3->SetMarkerColor(kBlack);
    sub_hist3->SetLineColor(kBlack);
    sub_hist3->SetMarkerSize(0.7);

    // Create a Canvas with two pads (upper for histograms, lower for ratio)
    TCanvas *c1 = new TCanvas("c1", "Jets Comparison", 800, 600);
    c1->Divide(1, 2); // Divide the canvas into 2 pads: 1 for histograms, 2 for ratio

    // Upper pad for histograms
    c1->cd(1);
    gPad->SetLogy();
    hist1->GetXaxis()->SetRangeUser(0, 180);
    hist2->GetXaxis()->SetRangeUser(0, 180);
    hist3->GetXaxis()->SetRangeUser(0, 180);
    sub_hist1->GetXaxis()->SetRangeUser(0, 180);
    sub_hist2->GetXaxis()->SetRangeUser(0, 180);
    sub_hist3->GetXaxis()->SetRangeUser(0, 180);
    
    hist1->SetTitle(""); // Remove title
    hist2->SetTitle(""); // Remove title
    hist3->SetTitle(""); // Remove title
    sub_hist1->SetTitle(""); // Remove title
    sub_hist2->SetTitle(""); // Remove title
    sub_hist3->SetTitle("");
    
    hist1->GetXaxis()->SetTitle("Pt [GeV]");
    hist1->GetYaxis()->SetTitle("Counts");

    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");
    hist3->Draw("P SAME");
    sub_hist1->Draw("HIST SAME");
    sub_hist2->Draw("HIST SAME");
    sub_hist3->Draw("P SAME");

       
    // Add legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.85);
    legend->AddEntry(hist1, "PYTHIA8 photp-ND ", "l");
    legend->AddEntry(hist2, "PYTHIA ND", "l");
    legend->AddEntry(hist3, "leading jet (Data)", "p");
    legend->AddEntry(sub_hist1, "PYTHIA8 photp-ND", "l");
    legend->AddEntry(sub_hist2, "PYTHIA ND", "l");
    legend->AddEntry(sub_hist3, "Subleading Jet (Data)", "p");
    // Remove the box lines and set transparent background
    legend->SetFillStyle(0);  // Set fill style to transparent
    legend->SetBorderSize(0); // Remove the border around the legend
    legend->Draw();

    // Add CMS Work In Progress text using TLatex
    TLatex cmsText;
cmsText.SetTextSize(0.035);  // Adjust text size 
cmsText.SetTextFont(61);  // Set font to bold for CMS
cmsText.SetNDC();  // Use normalized coordinates for positioning
cmsText.DrawLatex(0.12, 0.92, "CMS");  // Position CMS on the top left

TLatex workInProgressText;
workInProgressText.SetTextSize(0.035);  // Adjust size for the second line
workInProgressText.SetTextFont(42);  // Set font to normal (non-bold) for Work in Progress
workInProgressText.SetNDC();  // Use normalized coordinates for positioning
workInProgressText.DrawLatex(0.16, 0.92, "Work in Progress");  // Position it next to CMS

// Add pPb 8.16TeV (2016) text using TLatex
TLatex pPbText;
pPbText.SetTextFont(42);  // Normal text font
pPbText.SetTextSize(0.035);  // Adjust text size for pPb label
pPbText.SetNDC();  // Use normalized coordinates for positioning
pPbText.DrawLatex(0.75, 0.92, "pPb 8.16TeV (2016)");  // Position pPb label

    // Lower pad for ratio plot
    c1->cd(2);
    gPad->SetPad(0.0, 0.5, 1.0, 0.3); // Set the position and size of the lower pad
    gPad->SetTopMargin(0.04); // Reduce the top margin for ratio plot
    gPad->SetBottomMargin(0.3); // Increase the bottom margin for ratio plot

    // Create ratio histograms
    TH1D *ratio_hist1 = (TH1D*)hist1->Clone("ratio_hist1");
    ratio_hist1->Divide(hist3);
    TH1D *sub_ratio_hist1 = (TH1D*)sub_hist1->Clone("sub_ratio_hist1");
    sub_ratio_hist1->Divide(sub_hist3);

    ratio_hist1->SetStats(kFALSE);
    sub_ratio_hist1->SetStats(kFALSE);

    ratio_hist1->SetMarkerStyle(kFullCircle);
    ratio_hist1->SetMarkerColor(kBlack);
    ratio_hist1->SetLineColor(kBlack);
    ratio_hist1->SetMarkerSize(0.7);
    
    sub_ratio_hist1->SetMarkerStyle(kOpenCircle);
    sub_ratio_hist1->SetMarkerColor(kBlack);
    sub_ratio_hist1->SetLineColor(kBlack);
    sub_ratio_hist1->SetMarkerSize(0.7);
    


    ratio_hist1->Draw("EP");
    sub_ratio_hist1->Draw("EP SAME");

    ratio_hist1->GetYaxis()->SetRangeUser(0, 2.5); // Set ratio range
    ratio_hist1->GetYaxis()->SetTitle("MC/Data");
    ratio_hist1->GetYaxis()->SetTitleSize(0.08); 
    ratio_hist1->GetYaxis()->SetTitleOffset(0.35);  // Adjust the Y-axis title offset to align it
    ratio_hist1->GetXaxis()->SetTitle("Pt [GeV] ");
    ratio_hist1->GetXaxis()->SetTitleSize(0.08);

    // Save plot
    c1->SaveAs("Lead_vs_SubLead_Jets_with_Ratio.pdf");

    // Cleanup
    delete c1;
    delete hist1;
    delete hist2;
    delete hist3;
    delete sub_hist1;
    delete sub_hist2;
    delete sub_hist3;
    delete ratio_hist1;
    delete sub_ratio_hist1;
    delete legend;
    f1->Close();
    f2->Close();
    f3->Close();
    fsub1->Close();
    fsub2->Close();
    fsub3->Close();
}


