#include <Vizard.h>


namespace insur {
    // public
    /**
     * The constructor builds the logical structure within the <i>TGeoManager</i> that is used to display
     * a tracker geometry in ROOT. It assigns different materials and media to the various categories of
     * geometry elements.
     */
    Vizard::Vizard() {
        // internal flag
        geometry_created = false;
        // ROOT geometry manager
        gm = new TGeoManager("display", "Tracker");
        // dummy material definitions for each category
        matvac = new TGeoMaterial("Vacuum", 0, 0, 0);
        matact = new TGeoMaterial("Si", a_silicon, z_silicon, d_silicon);
        matserf = new TGeoMaterial("C ", a_carbon, z_carbon, d_carbon);
        matlazy = new TGeoMaterial("Cu", a_copper, z_copper, d_copper);
        // dummy medium definitions for each category
        medvac = new TGeoMedium("Vacuum", 0, matvac);
        medact = new TGeoMedium("Silicon", 1, matact);
        medserf = new TGeoMedium("Copper", 2, matserf);
        medlazy = new TGeoMedium("Carbon", 3, matlazy);
        // hierarchy definitions to group individual volumes
        barrels = new TGeoVolumeAssembly("Barrels");
        endcaps = new TGeoVolumeAssembly("Endcaps");
        services = new TGeoVolumeAssembly("Services");
        supports = new TGeoVolumeAssembly("Supports");
        active = new TGeoVolumeAssembly("Active Modules");
        inactive = new TGeoVolumeAssembly("Inactive Surfaces");
        // top-level volume definition
        top = gm->MakeBox("WORLD", medvac, outer_radius + top_volume_pad, outer_radius + top_volume_pad, max_length + top_volume_pad);
        // definition of tree hierarchy for visualisation
        active->AddNode(barrels, 0);
        active->AddNode(endcaps, 0);
        inactive->AddNode(services, 0);
        inactive->AddNode(supports, 0);
        top->AddNode(active, 0);
        top->AddNode(inactive, 0);
        // declaration of top volume within ROOT geometry manager
        gm->SetTopVolume(top);
	// Some stylish option
	gStyle->SetOptStat(0);
	const UInt_t numberOfSteps = 5;
	Double_t stops[numberOfSteps] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[numberOfSteps]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[numberOfSteps] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[numberOfSteps]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	Int_t myPalette[temperature_levels];
    
	Int_t colorIndex = TColor::CreateGradientColorTable(numberOfSteps, stops, red, green, blue, temperature_levels);
	for (int i=0;i<temperature_levels;i++) myPalette[i] = colorIndex+i;
	gStyle->SetPalette(temperature_levels, myPalette);
    }
    
    /**
     * The destructor deletes the instance of <i>TGeoManager</i> that was created in the constructor. Since
     * any object that was added to it in the constructor is now owned by it, this is the only step that is necessary
     * to clean up memory: everything else will be deleted from within the <i>TGeoManager</i>.
     */
    Vizard::~Vizard() {
        delete gm;
    }
    
    /**
     * This function turns the abstract representation of the active and inactive surfaces in a tracker geometry into
     * a series of ROOT shapes. Those shapes are added as leaves to the collections of volumes that were previously
     * initialised in the constructor. Once all of them have been added to the volume tree, the geometry manager is
     * closed, making it ready to be displayed or written to file.
     *
     * NOTE: It is highly recommended to use the <i>simplified</i> flag to have layers and discs represented as
     * bounding boxes rather than by individual module. In a typical case, the modules in a tracker number in the
     * thousands; since volume creation and placement in this function is at this point unoptimised, this will cause
     * considerable strain on the resources of whatever is used to visualise the geometry tree later.
     *
     * @param am A reference to the tracker object that contains the collection of active surfaces
     * @param is A reference to the collection of inactive surfaces
     * @param simplified A flag indicating whether to draw bounding boxes around the layers/discs or whether to display each module individually
     */
    void Vizard::buildVisualization(Tracker& am, InactiveSurfaces& is, bool simplified) {
        int c = 0;
        TGeoVolume* vol=NULL;
        TGeoTranslation* trans=NULL;
        TGeoCombiTrans* trafo=NULL;
        Layer* current=NULL;
        std::vector<Module*> templates;
        // barrels
        if (simplified) {
            // layer loop, one tube per layer
            for (unsigned int i = 0; i < am.getBarrelLayers()->size(); i++) {
                current = am.getBarrelLayers()->at(i);
                // short layers
                if ((current->getMinZ() > 0) || (current->getMaxZ() < 0)) {
                    vol = gm->MakeTube("", medact, current->getMinRho(), current->getMaxRho(), (current->getMaxZ() - current->getMinZ()) / 2.0);
                    vol->SetLineColor(kRed);
                    trans = new TGeoTranslation(0, 0, current->getMaxZ() - (current->getMaxZ() - current->getMinZ()) / 2.0);
                    barrels->AddNode(vol, c, trans);
                }
                // regular layers
                else {
                    vol = gm->MakeTube("", medact, current->getMinRho(), current->getMaxRho(), current->getMaxZ());
                    vol->SetLineColor(kRed);
                    barrels->AddNode(vol, c);
                }
                c++;
            }
        }
        else c = detailedModules(am.getBarrelLayers(), vol, trafo, barrels, c);
        // endcaps
        if (simplified) {
            // disc loop, one (very short) tube per disc
            for (unsigned int i = 0; i < am.getEndcapLayers()->size(); i++) {
                current = am.getEndcapLayers()->at(i);
                vol = gm->MakeTube("", medact, current->getMinRho(), current->getMaxRho(),
                        (current->getMaxZ() - current->getMinZ()) / 2.0);
                vol->SetLineColor(kRed);
                trans = new TGeoTranslation(0, 0, current->getMaxZ() - (current->getMaxZ() - current->getMinZ()) / 2.0);
                endcaps->AddNode(vol, c, trans);
                c++;
            }
        }
        else c = detailedModules(am.getEndcapLayers(), vol, trafo, endcaps, c);
        
        // services
        int skip = is.getBarrelServices().size() / 2;
        // barrel services loop using symmetries with respect to z=0
        for (int i = 0; i < skip; i++) {
            vol = gm->MakeTube("", medserf, is.getBarrelServicePart(i).getInnerRadius(),
                    is.getBarrelServicePart(i).getInnerRadius() + is.getBarrelServicePart(i).getRWidth(),
                    is.getBarrelServicePart(i).getZLength() / 2.0);
            vol->SetLineColor(kBlue);
            trans = new TGeoTranslation(0, 0, (is.getBarrelServicePart(i).getZOffset() + is.getBarrelServicePart(i).getZLength() / 2.0));
            services->AddNode(vol, c, trans);
            trans = new TGeoTranslation(0, 0, (is.getBarrelServicePart(i + skip).getZOffset() + is.getBarrelServicePart(i + skip).getZLength() / 2.0));
            services->AddNode(vol, c + skip, trans);
            c++;
        }
        c = c + skip;
        skip = is.getEndcapServices().size() / 2;
        // endcap services loop using symmetries with respect to z=0
        for (int i = 0; i < skip; i++) {
            vol = gm->MakeTube("", medserf, is.getEndcapServicePart(i).getInnerRadius(),
                    is.getEndcapServicePart(i).getInnerRadius() + is.getEndcapServicePart(i).getRWidth(),
                    is.getEndcapServicePart(i).getZLength() / 2.0);
            vol->SetLineColor(kBlue);
            trans = new TGeoTranslation(0, 0, (is.getEndcapServicePart(i).getZOffset() + is.getEndcapServicePart(i).getZLength() / 2.0));
            services->AddNode(vol, c, trans);
            trans = new TGeoTranslation(0, 0, (is.getEndcapServicePart(i + skip).getZOffset() + is.getEndcapServicePart(i + skip).getZLength() / 2.0));
            services->AddNode(vol, c + skip, trans);
            c++;
        }
        c = c + skip;
        
        // supports
        skip = is.getSupports().size();
        // support parts loop, using all entries
        for (int i = 0; i < skip; i++) {
            // process entry if its rightmost point is in z+ - this includes tubes that cross z=0 but not disc supports in z-
            if ((is.getSupportPart(i).getZOffset() + is.getSupportPart(i).getZLength()) > 0) {
                vol = gm->MakeTube("", medlazy, is.getSupportPart(i).getInnerRadius(),
                        is.getSupportPart(i).getInnerRadius() + is.getSupportPart(i).getRWidth(),
                        is.getSupportPart(i).getZLength() / 2.0);
                vol->SetLineColor(kGray);
                trans = new TGeoTranslation(0, 0, (is.getSupportPart(i).getZOffset() + is.getSupportPart(i).getZLength() / 2.0));
                supports->AddNode(vol, c, trans);
                c++;
                // use symmetries with respect to z=0: if volume is completely in z+, it will have a twin in z-
                if (is.getSupportPart(i).getZOffset() > 0) {
                    trans = new TGeoTranslation(0, 0, (0.0 - is.getSupportPart(i).getZOffset() - is.getSupportPart(i).getZLength() / 2.0));
                    supports->AddNode(vol, c, trans);
                    c++;
                }
            }
        }
        // check overlaps, write status to cout
        gm->CloseGeometry();
        geometry_created = true;
    }
    
    /**
     * This function writes the previously created geometry tree (including the geometry manager) to a ROOT
     * file. If it finds that the internal representation using ROOT shapes has not been initialised, it prints an 
     * error message and does nothing.
     * @param rootfilename The name of the output file that will be written to the application's default directory for root files
     */
    void Vizard::display(std::string rootfilename) {
        if (geometry_created) {
            std::string outfilename = default_rootfiledir + "/";
            if (rootfilename.empty()) outfilename = outfilename + default_rootfile;
            else outfilename = outfilename + rootfilename;
            TFile f(outfilename.c_str(), "recreate");
            if (f.IsZombie()) {
                std::cout << root_wrong << std::endl;
            }
            else {
                gm->Write();
                f.Close();
                std::cout << "Geometry written to file '" << outfilename << "'." << std::endl;
            }
            // display top volume after opening file in ROOT with:
            // TFile f("rootfiles/output.root");
            // TGeoManager* gm = (TGeoManager*)f.Get("display");
            // gm->GetMasterVolume()->Draw("ogl");
            // press 'w' for wireframe or 't' for outline view
            // when done: delete gm, f.Close()
        }
        else std::cout << msg_uninitialised << std::endl;
    }
    
    /**
     * This convenience function provides a frame for creation of a geometry tree from a tracker object and a 
     * collection of inactive surfaces and for writing them to a ROOT file in a single step.
     * @param am A reference to the tracker object that contains the collection of active surfaces
     * @param is A reference to the collection of inactive surfaces
     * @param rootfilename The name of the output file that will be written to the application's default directory for ROOT files
     * @param simplified A flag indicating whether to draw bounding boxes around the layers/discs or whether to display each module individually
     */
    void Vizard::display(Tracker& am, InactiveSurfaces& is, std::string rootfilename, bool simplified) {
        buildVisualization(am, is, simplified);
        display(rootfilename);
    }
    
    /**
     * This convenience function writes the feeder/neighbour relations of a collection of inactive surfaces to a
     * default file.
     * @param is A reference to the collection of inactive surfaces
     */
    void Vizard::writeNeighbourGraph(InactiveSurfaces& is) {
        writeNeighbourGraph(is, default_graphfile);
    }
    
    /**
     * This function writes the feeder/neighbour relations in a collection of inactive surfaces to a very simple text
     * file. It essentially lists all edges of the neighbour graph, first those in the barrels, then those in the endcaps,
     * but otherwise in a more or less unordered heap: the more order in the source collection, the more order in
     * the output. If no name is given for the output file, a default filename is used.
     * @param is A reference to the collection of inactive surfaces
     * @param outfile The name of the output file that will be written to the application's default directory for graph files
     */
    void Vizard::writeNeighbourGraph(InactiveSurfaces& is, std::string outfile) {
        std::string filename = default_graphdir + "/";
        if (outfile.empty()) filename = filename + default_graphfile;
        else filename = filename + outfile;
        std::cout << "Preparing to write neighbour graph to " << filename << "..." << std::endl;
        try {
            std::ofstream outstream(filename.c_str());
            if (outstream) {
                // barrel services loop
                outstream << "BARREL SERVICES:" << std::endl << std::endl;
                for (unsigned int i = 0; i < is.getBarrelServices().size(); i++) {
                    outstream << "Barrel element " << i << ": service is ";
                    if (is.getBarrelServicePart(i).isFinal()) outstream << "final and ";
                    else outstream << "not final and ";
                    if (is.getBarrelServicePart(i).isVertical()) outstream << "vertical.";
                    else outstream << "horizontal.";
                    outstream << std::endl << "Feeder type: ";
                    switch (is.getBarrelServicePart(i).getFeederType()) {
                        case InactiveElement::no_in: outstream << "none, ";
                        break;
                        case InactiveElement::tracker: outstream << "tracker, ";
                        break;
                        case InactiveElement::barrel: outstream << "barrel service, ";
                        break;
                        case InactiveElement::endcap: outstream << "endcap service, ";
                        break;
                        default: outstream << "something weird, ";
                    }
                    outstream << "feeder index = " << is.getBarrelServicePart(i).getFeederIndex() << ".";
                    outstream << std::endl << "Neighbour type: ";
                    switch (is.getBarrelServicePart(i).getNeighbourType()) {
                        case InactiveElement::no_in: outstream << "none, ";
                        break;
                        case InactiveElement::tracker: outstream << "tracker, ";
                        break;
                        case InactiveElement::barrel: outstream << "barrel service, ";
                        break;
                        case InactiveElement::endcap: outstream << "endcap service, ";
                        break;
                        default: outstream << "something weird, ";
                    }
                    outstream << "neighbour index = " << is.getBarrelServicePart(i).getNeighbourIndex() << ".";
                    outstream << std::endl << std::endl;
                }
                // endcap services
                outstream << "ENDCAP SERVICES:" << std::endl << std::endl;
                for (unsigned int i = 0; i < is.getEndcapServices().size(); i++) {
                    outstream << "Endcap element " << i << ": service is ";
                    if (is.getEndcapServicePart(i).isFinal()) outstream << "final and ";
                    else outstream << "not final and ";
                    if (is.getEndcapServicePart(i).isVertical()) outstream << "vertical.";
                    else outstream << "horizontal.";
                    outstream << std::endl << "Feeder type: ";
                    switch (is.getEndcapServicePart(i).getFeederType()) {
                        case InactiveElement::no_in: outstream << "none, ";
                        break;
                        case InactiveElement::tracker: outstream << "tracker, ";
                        break;
                        case InactiveElement::barrel: outstream << "barrel service, ";
                        break;
                        case InactiveElement::endcap: outstream << "endcap service, ";
                        break;
                        default: outstream << "something weird, ";
                    }
                    outstream << "feeder index = " << is.getEndcapServicePart(i).getFeederIndex() << ".";
                    outstream << std::endl << "Neighbour type: ";
                    switch (is.getEndcapServicePart(i).getNeighbourType()) {
                        case InactiveElement::no_in: outstream << "none, ";
                        break;
                        case InactiveElement::tracker: outstream << "tracker, ";
                        break;
                        case InactiveElement::barrel: outstream << "barrel service, ";
                        break;
                        case InactiveElement::endcap: outstream << "endcap service, ";
                        break;
                        default: outstream << "something weird, ";
                    }
                    outstream << "neighbour index = " << is.getEndcapServicePart(i).getNeighbourIndex() << ".";
                    outstream << std::endl << std::endl;
                }
                outstream.close();
                std::cout << "Neighbour graph written to " << filename << "." << std::endl;
            }
            else std::cout << graph_wrong << std::endl;
        }
        catch (std::bad_alloc ba) {
            std::cerr << exc_badalloc_graph << graph_nowrite << std::endl;
        }
    }
    
    /**
     * This function is meant to write the feeder/neighbour relations of a given collection of inactive surfaces
     * to a DOT file instead of the quick and dirty text format that is used at the moment.
     *
     * NOTE: This function is currently in DEVELOPMENT HELL. There is no way of knowing if it will ever
     * finished, or when. So for now, the function can be called, but it does NOTHING AT ALL. Don't say you
     * haven't been warned.
     * @param is A reference to the collection of inactive surfaces
     * @param outfile The name of the output file that will be written to the application's default directory for graph files
     */
    void Vizard::dotGraph(InactiveSurfaces& is, std::string outfile) {
        const std::string preamble = "digraph tracker";
        const std::string ori = "rankdir=DU"; // check if this is possible!
        const std::string shape = "node [shape=box]";
        const std::string label = "label=";
        const std::string edge = "->";
    }
    
    /**
     * This function draws some of the histograms that were filled during material budget analysis and
     * embeds the resulting image in an HTML file for easy access. If no name is given for the output file,
     * a default filename is used.
     * @param a A reference to the analysing class that examined the material budget and filled the histograms
     * @param outfilename The name of the output file that will be written to the application's default directory for material budget summaries
     */
  void Vizard::histogramSummary(Analyzer& a, std::string outfilename) {
        THStack rcontainer("rstack", "Radiation Length by Category");
        THStack icontainer("istack", "Interaction Length by Category");
        TH1D *cr = NULL, *ci = NULL, *fr1 = NULL, *fi1 = NULL, *fr2 = NULL, *fi2 = NULL;
        TH1D *acr = NULL, *aci = NULL, *ser = NULL, *sei = NULL, *sur = NULL, *sui = NULL;
        TH2D *ir = NULL, *ii = NULL;
        TVirtualPad* pad;
        // filename gymnastics
        std::string outfile = default_summarypath + "/";
        std::string pngoutfile, pngout, pngpath;
        std::string svgout, svgpath;
        std::string Cout, Cpath;
        if (outfilename.empty()) {
            outfile = outfile + default_summary;
            pngoutfile = default_summary;
        }
        else {
            outfile = outfile + outfilename;
            pngoutfile = outfilename;
        }
        if (pngoutfile.substr(pngoutfile.size() - 4) == "html") {
            pngoutfile = pngoutfile.substr(0, pngoutfile.size() - 4);
            if (pngoutfile.substr(pngoutfile.size() - 1) == ".") pngoutfile = pngoutfile.substr(0, pngoutfile.size() - 1);
        }
        // output initialisation and headers
        std::ostringstream htmlstream;
        std::ofstream outstream(outfile.c_str());
        htmlstream << "<html><title>" << outfilename << "</title><body>";
        htmlstream << "<p><big><b>1D Overview</b></big></p>";
        TCanvas c("matbudgetcanvas", "Material Budgets over Eta", 900, 400);
        c.SetFillColor(color_plot_background);
        c.Divide(2,1);
        pad = c.GetPad(0);
        pad->SetFillColor(color_pad_background);
        pad = c.GetPad(1);
        pad->cd();
        // total tracking volume rlength
        cr = (TH1D*)a.getHistoGlobalR().Clone();
        fr1 = (TH1D*)a.getHistoExtraServicesR().Clone();
        fr2 = (TH1D*)a.getHistoExtraSupportsR().Clone();
        fr1 = (TH1D*)a.getHistoExtraServicesR().Clone();
        fr2 = (TH1D*)a.getHistoExtraSupportsR().Clone();
        cr->Add(fr1);
        cr->Add(fr2);
        cr->SetFillColor(kGray + 2);
        cr->SetNameTitle("rfullvolume", "Radiation Length Over Full Tracker Volume");
        cr->SetXTitle("Eta");
        cr->Draw();
        pad = c.GetPad(2);
        pad->cd();
        // total tracking volume ilength
        ci = (TH1D*)a.getHistoGlobalI().Clone();
        fi1 = (TH1D*)a.getHistoExtraServicesI().Clone();
        fi2 = (TH1D*)a.getHistoExtraSupportsI().Clone();
        fi1 = (TH1D*)a.getHistoExtraServicesI().Clone();
        fi2 = (TH1D*)a.getHistoExtraSupportsI().Clone();
        ci->Add(fi1);
        ci->Add(fi2);
        ci->SetFillColor(kGray + 1);
        ci->SetNameTitle("ifullvolume", "Interaction Length Over Full Tracker Volume");
        ci->SetXTitle("Eta");
        ci->Draw();
        // write total plots to file
        pngout = pngoutfile + ".fullvolume.png";
        pngpath = default_summarypath + "/" + pngout;
        svgout = pngoutfile + ".fullvolume.svg";
        svgpath = default_summarypath + "/" + svgout;
        Cout = pngoutfile + ".fullvolume.C";
        Cpath = default_summarypath + "/" + Cout;
        c.SaveAs(pngpath.c_str());
        c.SaveAs(svgpath.c_str());
        c.SaveAs(Cpath.c_str());
        htmlstream << "<img src=\"" << pngout << "\" /><br><p><small><b>Average radiation length in full volume ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*cr, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in full volume (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*ci, etaMaxAvg) << "</b></small></p>";
        // work area cleanup an re-init
        c.Clear("D");
        if (cr) delete cr;
        if (ci) delete ci;
        pad = c.GetPad(1);
        pad->cd();
        // global plots in tracking volume: radiation length
        cr = (TH1D*)a.getHistoGlobalR().Clone();
        cr->SetFillColor(kGray + 2);
        cr->SetNameTitle("rglobal", "Overall Radiation Length");
        cr->SetXTitle("Eta");
        cr->Draw();
        pad = c.GetPad(2);
        pad->cd();
        // global plots in tracking volume: interaction length
        ci = (TH1D*)a.getHistoGlobalI().Clone();
        ci->SetFillColor(kGray + 1);
        ci->SetNameTitle("iglobal", "Overall Interaction Length");
        ci->SetXTitle("Eta");
        ci->Draw();
        // write global tracking volume plots to file
        pngout = pngoutfile + ".global.png";
        pngpath = default_summarypath + "/" + pngout;
        svgout = pngoutfile + ".global.svg";
        svgpath = default_summarypath + "/" + svgout;
        Cout = pngoutfile + ".global.C";
        Cpath = default_summarypath + "/" + Cout;
        c.SaveAs(pngpath.c_str());
        c.SaveAs(svgpath.c_str());
        c.SaveAs(Cpath.c_str());
        htmlstream << "<br><img src=\"" << pngout << "\" /><br><p><small><b>Average radiation length in tracking volume ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*cr, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in tracking volume (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*ci, etaMaxAvg) << "</b></small></p>";
        // work area cleanup and re-init
        c.Clear("D");
        pad = c.GetPad(1);
        pad->cd();
        // radiation length in tracking volume by active, serving or passive
        sur = (TH1D*)a.getHistoSupportsAllR().Clone();
        sur->SetFillColor(kOrange + 4);
        sur->SetXTitle("Eta");
        rcontainer.Add(sur);
        ser = (TH1D*)a.getHistoServicesAllR().Clone();
        ser->SetFillColor(kBlue);
        ser->SetXTitle("Eta");
        rcontainer.Add(ser);
        acr = (TH1D*)a.getHistoModulesAllR().Clone();
        acr->SetFillColor(kRed);
        acr->SetXTitle("Eta");
        rcontainer.Add(acr);
        rcontainer.Draw();
        // interaction length in tracking volume by active, serving or passive
        pad = c.GetPad(2);
        pad->cd();
        sui = (TH1D*)a.getHistoSupportsAllI().Clone();
        sui->SetFillColor(kOrange + 2);
        sui->SetXTitle("Eta");
        icontainer.Add(sui);
        sei = (TH1D*)a.getHistoServicesAllI().Clone();
        sei->SetFillColor(kAzure - 2);
        sei->SetXTitle("Eta");
        icontainer.Add(sei);
        aci = (TH1D*)a.getHistoModulesAllI().Clone();
        aci->SetFillColor(kRed - 3);
        aci->SetXTitle("Eta");
        icontainer.Add(aci);
        icontainer.Draw();
        // write asl category plots to file
        pngout = pngoutfile + ".asl.png";
        pngpath = default_summarypath + "/" + pngout;
        svgout = pngoutfile + ".asl.svg";
        svgpath = default_summarypath + "/" + svgout;
        Cout = pngoutfile + ".asl.C";
        Cpath = default_summarypath + "/" + Cout;
        c.SaveAs(pngpath.c_str());
        c.SaveAs(svgpath.c_str());
        c.SaveAs(Cpath.c_str());
        htmlstream << "<img src=\"" << pngout << "\" /><br>";
        // average values by active, service and passive
        htmlstream << "<br><p><small><b>Average radiation length in modules ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*acr, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average radiation length in services ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*ser, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average radiation length in supports ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*sur, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<br><p><small><b>Average interaction length in modules (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*aci, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in services (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*sei, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in supports (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*sui, etaMaxAvg) << "</b></small></p>";
        // work area cleanup and re-init
        c.Clear("D");
        pad = c.GetPad(1);
        pad->cd();
        // radiation length in isolines
        ir = (TH2D*)a.getHistoIsoR().Clone();
        ir->SetNameTitle("isor", "Radiation Length Contours");
        ir->SetContour(temperature_levels, NULL);
        ir->SetXTitle("z");
        ir->SetYTitle("r");
        ir->Draw("COLZ");
        pad = c.GetPad(2);
        pad->cd();
        // interaction length in isolines
        ii = (TH2D*)a.getHistoIsoI().Clone();
        ii->SetNameTitle("isoi", "Interaction Length Contours");
        ii->SetContour(temperature_levels, NULL);
        ii->SetXTitle("z");
        ii->SetYTitle("r");
        ii->Draw("COLZ");
        // write isoline plots to file
        pngout = pngoutfile + ".twodee.png";
        pngpath = default_summarypath + "/" + pngout;
        svgout = pngoutfile + ".twodee.svg";
        svgpath = default_summarypath + "/" + svgout;
        Cout = pngoutfile + ".twodee.C";
        Cpath = default_summarypath + "/" + Cout;
        c.SaveAs(pngpath.c_str());
        c.SaveAs(svgpath.c_str());
        c.SaveAs(Cpath.c_str());
        htmlstream << "<br><p><big><b>2D Overview</b></big></p><img src=\"" << pngout << "\" /><br>";
        // general cleanup
        if (cr) delete cr;
        if (ci) delete ci;
        if (fr1) delete fr1;
        if (fi1) delete fi1;
        if (fr2) delete fr2;
        if (fi2) delete fi2;
        if (acr) delete acr;
        if (aci) delete aci;
        if (ser) delete ser;
        if (sei) delete sei;
        if (sur) delete sur;
        if (sui) delete sui;
        // write html frame to file
        htmlstream << "</body></html>";
        outstream << htmlstream.str() << std::endl;
        outstream.close();
        std::cout << "HTML file written to " << outfile << std::endl;
    }

#ifdef USING_ROOTWEB
    /**
     * This function draws some of the histograms that were filled during material budget analysis
     * with the rootweb library
     * @param a A reference to the analysing class that examined the material budget and filled the histograms
     * @param site the RootWSite object for the output
     */
    void Vizard::histogramSummary(Analyzer& a, RootWSite& site) {
      // Initialize the page with the material budget
      RootWPage* myPage;
      RootWContent* myContent;
      RootWTable* myTable;
      RootWImage* myImage;
      TCanvas* myCanvas;
      TVirtualPad* myPad;
      myPage = new RootWPage("Material");
      myPage->setAddress("material.html");
      site.addPage(myPage);

      // 1D Overview
      myContent = new RootWContent("1D Overview");
      myPage->addContent(myContent);

      // Book histograms
      THStack* rcontainer = new THStack("rstack", "Radiation Length by Category");
      THStack* icontainer = new THStack("istack", "Interaction Length by Category");
      TH1D *cr = NULL, *ci = NULL, *fr1 = NULL, *fi1 = NULL, *fr2 = NULL, *fi2 = NULL;
      TH1D *acr = NULL, *aci = NULL, *ser = NULL, *sei = NULL, *sur = NULL, *sui = NULL;
      TH2D *ir = NULL, *ii = NULL;

      // Output initialisation and headers
      myCanvas = new TCanvas("overviewMaterial");
      myCanvas->SetFillColor(color_plot_background);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(color_pad_background);
      myPad = myCanvas->GetPad(1);
      myPad->cd();
      // Total tracking volume rlength
      cr = (TH1D*)a.getHistoGlobalR().Clone();
      fr1 = (TH1D*)a.getHistoExtraServicesR().Clone();
      fr2 = (TH1D*)a.getHistoExtraSupportsR().Clone();
      fr1 = (TH1D*)a.getHistoExtraServicesR().Clone();
      fr2 = (TH1D*)a.getHistoExtraSupportsR().Clone();
      cr->Add(fr1);
      cr->Add(fr2);
      cr->SetFillColor(kGray + 2);
      cr->SetNameTitle("rfullvolume", "Radiation Length Over Full Tracker Volume");
      cr->SetXTitle("Eta");
      cr->Draw();
      myPad = myCanvas->GetPad(2);
      myPad->cd();
      // Total Tracking volume ilength
      ci = (TH1D*)a.getHistoGlobalI().Clone();
      fi1 = (TH1D*)a.getHistoExtraServicesI().Clone();
      fi2 = (TH1D*)a.getHistoExtraSupportsI().Clone();
      fi1 = (TH1D*)a.getHistoExtraServicesI().Clone();
      fi2 = (TH1D*)a.getHistoExtraSupportsI().Clone();
      ci->Add(fi1);
      ci->Add(fi2);
      ci->SetFillColor(kGray + 1);
      ci->SetNameTitle("ifullvolume", "Interaction Length Over Full Tracker Volume");
      ci->SetXTitle("Eta");
      ci->Draw();
      // Put the total plots to the site
      myImage = new RootWImage(myCanvas, 900, 400);
      myImage->setComment("Material in full volume");
      myTable = new RootWTable();
      // TODO: put etaMaxAvg correctly in the string :)
      myTable->setContent(1,1,"Average radiation length in full volume (eta = [0, 2.4])");
      myTable->setContent(2,1,"Average interaction length in full volume (eta = [0, 2.4])");
      myTable->setContent(1,2,averageHistogramValues(*cr, etaMaxAvg), 5);
      myTable->setContent(2,2,averageHistogramValues(*ci, etaMaxAvg), 5);
      myContent->addItem(myTable);
      myContent->addItem(myImage);
      
      // Detailed plots
      myContent = new RootWContent("Tracking volume", false);
      myPage->addContent(myContent);
      // Work area re-init
      myCanvas = new TCanvas("materialInTrackingVolume");
      myCanvas->SetFillColor(color_plot_background);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(color_pad_background);
      myPad = myCanvas->GetPad(1);
      myPad->cd();
      // global plots in tracking volume: radiation length
      cr = (TH1D*)a.getHistoGlobalR().Clone();
      cr->SetFillColor(kGray + 2);
      cr->SetNameTitle("rglobal", "Overall Radiation Length");
      cr->SetXTitle("Eta");
      cr->Draw();
      myPad = myCanvas->GetPad(2);
      myPad->cd();
      // global plots in tracking volume: interaction length
      ci = (TH1D*)a.getHistoGlobalI().Clone();
      ci->SetFillColor(kGray + 1);
      ci->SetNameTitle("iglobal", "Overall Interaction Length");
      ci->SetXTitle("Eta");
      ci->Draw();
      // Write global tracking volume plots to web pag
      myImage = new RootWImage(myCanvas, 900, 400);
      myImage->setComment("Material in tracking volume");
      myTable = new RootWTable();
      myTable->setContent(1,1,"Average radiation length in tracking volume (eta = [0, 2.4])");
      myTable->setContent(2,1,"Average interaction length in tracking volume (eta = [0, 2.4])");
      myTable->setContent(1,2,averageHistogramValues(*cr, etaMaxAvg), 5);
      myTable->setContent(2,2,averageHistogramValues(*ci, etaMaxAvg), 5);
      myContent->addItem(myTable);
      myContent->addItem(myImage);

      // Detailed plots
      myContent = new RootWContent("Detailed", false);
      myPage->addContent(myContent);
      // Work area re-init
      myCanvas = new TCanvas("detailedMaterial");
      myCanvas->SetFillColor(color_plot_background);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(color_pad_background);
      myPad = myCanvas->GetPad(1);
      myPad->cd();
      // radiation length in tracking volume by active, serving or passive
      sur = (TH1D*)a.getHistoSupportsAllR().Clone();
      sur->SetFillColor(kOrange + 4);
      sur->SetXTitle("Eta");
      rcontainer->Add(sur);
      ser = (TH1D*)a.getHistoServicesAllR().Clone();
      ser->SetFillColor(kBlue);
      ser->SetXTitle("Eta");
      rcontainer->Add(ser);
      acr = (TH1D*)a.getHistoModulesAllR().Clone();
      acr->SetFillColor(kRed);
      acr->SetXTitle("Eta");
      rcontainer->Add(acr);
      rcontainer->Draw();
      // interaction length in tracking volume by active, serving or passive
      myPad = myCanvas->GetPad(2);
      myPad->cd();
      sui = (TH1D*)a.getHistoSupportsAllI().Clone();
      sui->SetFillColor(kOrange + 2);
      sui->SetXTitle("Eta");
      icontainer->Add(sui);
      sei = (TH1D*)a.getHistoServicesAllI().Clone();
      sei->SetFillColor(kAzure - 2);
      sei->SetXTitle("Eta");
      icontainer->Add(sei);
      aci = (TH1D*)a.getHistoModulesAllI().Clone();
      aci->SetFillColor(kRed - 3);
      aci->SetXTitle("Eta");
      icontainer->Add(aci);
      icontainer->Draw();

      // Write asl category plots to web page
      myImage = new RootWImage(myCanvas, 900, 400);
      myImage->setComment("Detailed");
      myTable = new RootWTable();
      // Average values by active, service and passive      
      myTable->setContent(0, 0, "Average (eta = [0, 2.4])");
      myTable->setContent(1, 0, "modules");
      myTable->setContent(2, 0, "services");
      myTable->setContent(3, 0, "supports");
      myTable->setContent(0, 1, "Radiation length");
      myTable->setContent(0, 2, "Interaction length");
      myTable->setContent(1, 1, averageHistogramValues(*acr, etaMaxAvg), 5);
      myTable->setContent(2, 1, averageHistogramValues(*ser, etaMaxAvg), 5);
      myTable->setContent(3, 1, averageHistogramValues(*sur, etaMaxAvg), 5);
      myTable->setContent(1, 2, averageHistogramValues(*aci, etaMaxAvg), 5);
      myTable->setContent(2, 2, averageHistogramValues(*sei, etaMaxAvg), 5);
      myTable->setContent(3, 2, averageHistogramValues(*sui, etaMaxAvg), 5);
      myContent->addItem(myTable);
      myContent->addItem(myImage);


      // Work area re-init
      myCanvas = new TCanvas("countourMaterial");
      myCanvas->SetFillColor(color_plot_background);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(color_pad_background);
      myPad = myCanvas->GetPad(1);
      myPad->cd();

      // Countour plots
      myContent = new RootWContent("Contours", true);
      myPage->addContent(myContent);
      // radiation length in isolines
      ir = (TH2D*)a.getHistoIsoR().Clone();
      ir->SetNameTitle("isor", "Radiation Length Contours");
      ir->SetContour(temperature_levels, NULL);
      ir->SetXTitle("z");
      ir->SetYTitle("r");
      ir->Draw("COLZ");
      myPad = myCanvas->GetPad(2);
      myPad->cd();
      // interaction length in isolines
      ii = (TH2D*)a.getHistoIsoI().Clone();
      ii->SetNameTitle("isoi", "Interaction Length Contours");
      ii->SetContour(temperature_levels, NULL);
      ii->SetXTitle("z");
      ii->SetYTitle("r");
      ii->Draw("COLZ");
      // Write isoline plots to web page
      myImage = new RootWImage(myCanvas, 900, 400);
      myImage->setComment("Material 2D distributions");
      myContent->addItem(myImage);
    }

#endif
    
    // private
    /**
     * This function bundles the placement of a collection of individual modules in a ROOT geometry tree for 
     * visualisation. It loops through the provided module vectors, determining their modules' corners and position
     * in space. Using that information, it adds a ROOT shape and a 3D transformation to a volume assembly
     * note from the geometry tree for each module found in the vectors.
     * @param layers A pointer to the list of layers or discs that is to be displayed
     * @param v A pointer to a template volume that is to be adjusted to the module shape
     * @param t A pointer to a transformation object that will describe the template volume's position in space
     * @param a A pointer to the assembly node that the template volume will be added to
     * @param counter The element counter that keeps track of how many volumes have been added to the geometry tree
     * @return The new value of the element counter
     */
    int Vizard::detailedModules(std::vector<Layer*>* layers,
            TGeoVolume* v, TGeoCombiTrans* t, TGeoVolumeAssembly* a, int counter) {
        Layer* current;
        Module* mod;
        if (!layers->empty()) {
            //  init of volume object for modules
            v = gm->MakeArb8("", medact, 0);
            v->SetLineColor(kRed);
            // layer loop
            for (unsigned int i = 0; i < layers->size(); i++) {
                current = layers->at(i);
                // module loop
                for (unsigned int j = 0; j < current->getModuleVector()->size(); j++) {
                    mod = current->getModuleVector()->at(j);
                    // place volume v according to information in module mod
                    t = modulePlacement(mod, v);
                    // add volume v to scene graph using translation t
                    a->AddNode(v, counter, t);
                    counter++;
                }
            }
        }
        else std::cout << "detailedModules(): layers vector is empty." << std::endl;
        return counter;
    }
    
    /**
     * This geometry function computes the transformation matrix that describes the position of a given module
     * in space. It also sets the shape of the provided template volume to correspond to that of the module.
     * @param m A pointer to the module object that needs to be visualised
     * @param v A pointer to the template volume that will represent the module in the visualisation
     * @return A pointer to the finished transformation matrix object
     */
    TGeoCombiTrans* Vizard::modulePlacement(Module* m, TGeoVolume* v) {
        XYZVector ex, ey, ez, b, c, d, p;
        TGeoArb8* arb;
        TGeoRotation* rot;
        TGeoCombiTrans* tr;
        // copy of module placement parameters in Module class
        b = m->getCorner(1) - m->getCorner(0);
        c = m->getCorner(2) - m->getCorner(0);
        d = m->getCorner(3) - m->getCorner(0);
        ex = b / b.R();
        p = (d.Dot(ex) * ex);
        // unit vectors for module coordinate system
        ey = d - p;
        ey = ey / ey.R();
        ez = ex.Cross(ey);
        // set vertices in volume v according to extracted module measurements
        arb = (TGeoArb8*)(v->GetShape());
        for (int i = 0; i < 5; i = i + 4) {
            arb->SetVertex(i, 0, 0);
            arb->SetVertex(i + 1, b.R(), 0);
            arb->SetVertex(i + 2, c.Dot(ex), c.Dot(ey));
            arb->SetVertex(i + 3, d.Dot(ex), d.Dot(ey));
        }
        // set position of module within the tracker volume
        double matrix[9];
        matrix[0] = ex.X();
        matrix[1] = ey.X();
        matrix[2] = ez.X();
        matrix[3] = ex.Y();
        matrix[4] = ey.Y();
        matrix[5] = ez.Y();
        matrix[6] = ex.Z();
        matrix[7] = ey.Z();
        matrix[8] = ez.Z();
        rot = new TGeoRotation();
        rot->SetMatrix(matrix);
        // save position in transformation object
        tr = new TGeoCombiTrans(m->getCorner(0).X(), m->getCorner(0).Y(), m->getCorner(0).Z(), rot);
        return tr;
    }
    
    /**
     * This function computes the average of a range of histogram bins: from the first to the one that includes a
     * cutoff value along the axis.
     * @param histo A reference to the histogram data
     * @param cutoff The cutoff value
     * @return The average value of the bins within range
     */
    double Vizard::averageHistogramValues(TH1D& histo, double cutoff) {
        double avg = 0.0;
        int cobin = 1;
        // find last relevant bin
        while ((cobin < histo.GetNbinsX()) && (histo.GetBinLowEdge(cobin) < cutoff)) cobin++;
        // calculate average
        if (cobin >= histo.GetNbinsX() - 1) avg = histo.GetMean();
        else {
	  for (int i = 1; i <= cobin; i++) avg = avg + histo.GetBinContent(i) / (double)cobin;
        }
        return avg;
    }
  

#ifdef USING_ROOTWEB
  /**
   * This function draws the profile of hits obtained by the analysis of the geometry
   * together with the summaries in tables with the rootweb library. It also actually does a couple of
   * calculations to count modules and such, to put the results in the tables.
   * @param analyzer A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   */
  bool Vizard::geometrySummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
    // Formatting parameters
    int coordPrecision = 0;
    int areaPrecision = 1;
    int occupancyPrecision = 1;
    int pitchPrecision = 0;
    int stripLengthPrecision = 1;
    int millionChannelPrecision = 2;
    int powerPrecision = 1;
    int costPrecision  = 1;
    int powerPerUnitPrecision = 2;
    int costPerUnitPrecision  = 1;
    
    // A bunch of indexes
    std::map<std::string, Module*> typeMap;
    std::map<std::string, int> typeMapCount;
    std::map<std::string, long> typeMapCountChan;
    std::map<std::string, double> typeMapMaxOccupancy;
    std::map<std::string, double> typeMapAveOccupancy;
    std::map<std::string, Module*>::iterator typeMapIt;
    std::map<int, Module*> ringTypeMap;
    std::string aSensorTag;
    LayerVector::iterator layIt;
    ModuleVector::iterator modIt;
    ModuleVector* aLay;
    double totAreaPts=0;
    double totAreaStrips=0;
    int totCountMod=0;
    int totCountSens=0;
    long totChannelStrips=0;
    long totChannelPts=0;

    RootWPage* myPage = new RootWPage("Geometry");
    // TODO: the web site should decide which page to call index.html
    myPage->setAddress("index.html");
    site.addPage(myPage);
    RootWContent* myContent;

    // Grab a list of layers from teh tracker object
    LayerVector& layerSet = tracker.getLayers();
    int nMB = tracker.getNMB();
    ModuleVector& endcapSample = tracker.getEndcapSample();


    //********************************//
    //*                              *//
    //*       Layers and disks       *//
    //*                              *//
    //********************************//
    myContent = new RootWContent("Layers and disks");
    myPage->addContent(myContent);
    RootWTable* layerTable = new RootWTable(); myContent->addItem(layerTable);
    RootWTable* diskTable = new RootWTable(); myContent->addItem(diskTable);
    RootWTable* ringTable = new RootWTable(); myContent->addItem(ringTable);
    
    
    std::vector<std::string> layerNames;
    std::vector<double> layerRho;
    std::vector<std::string> diskNames;
    std::vector<double> diskZ;
    std::vector<std::string> ringNames;
    std::vector<double> ringRho1;
    std::vector<double> ringRho2;
    
    Layer* aLayer;
    BarrelLayer* aBarrelLayer;
    EndcapLayer* anEndcapDisk;
    double aRingRho;
    
    layerTable->setContent(0,0, "Layer");
    layerTable->setContent(1,0, "r");
    diskTable->setContent(0,0, "Disk");
    diskTable->setContent(1,0, "z");
    ringTable->setContent(0,0, "Ring");
    ringTable->setContent(1,0, "r"+subStart+"min"+subEnd);
    ringTable->setContent(2,0, "r"+subStart+"max"+subEnd);
    
    // Build the module type maps
    // with a pointer to a sample module
    // Build the layer summary BTW
    int nBarrelLayers=0;
    int nDisks=0;
    for (layIt=layerSet.begin(); layIt!=layerSet.end(); layIt++) {
      aLayer = (*layIt);
      if ( (aBarrelLayer=dynamic_cast<BarrelLayer*>(aLayer)) ) {
	if (aBarrelLayer->getMaxZ(+1)>0) {
	  ++nBarrelLayers;
	  //std::cerr << "Layer number " << nBarrelLayers << std::endl;
	  layerTable->setContent(0,nBarrelLayers, aBarrelLayer->getName());
	  layerTable->setContent(1,nBarrelLayers, aBarrelLayer->getAverageRadius(), coordPrecision);
	}
      }
      if ( (anEndcapDisk=dynamic_cast<EndcapLayer*>(aLayer)) ) {
	if (anEndcapDisk->getAverageZ()>0) {
	  ++nDisks;
	  diskTable->setContent(0,nDisks, anEndcapDisk->getName());
	  diskTable->setContent(1,nDisks, anEndcapDisk->getAverageZ(), coordPrecision);
	}
      }
      aLay = (*layIt)->getModuleVector();
      for (modIt=aLay->begin(); modIt!=aLay->end(); modIt++) {
	aSensorTag=(*modIt)->getSensorTag();
	typeMapCount[aSensorTag]++;
	typeMapCountChan[aSensorTag]+=(*modIt)->getNChannels();
	if (((*modIt)->getOccupancyPerEvent()*nMB)>typeMapMaxOccupancy[aSensorTag]) {
	  typeMapMaxOccupancy[aSensorTag]=(*modIt)->getOccupancyPerEvent()*nMB;
	}
	typeMapAveOccupancy[aSensorTag]+=(*modIt)->getOccupancyPerEvent()*nMB;
	totCountMod++;
	totCountSens+=(*modIt)->getNFaces();
	if ((*modIt)->getReadoutType()==Module::Strip) {
	  totChannelStrips+=(*modIt)->getNChannels();
	  totAreaStrips+=(*modIt)->getArea()*(*modIt)->getNFaces();
	}
	if ((*modIt)->getReadoutType()==Module::Pt) {
	  totChannelPts+=(*modIt)->getNChannels();
	  totAreaPts+=(*modIt)->getArea()*(*modIt)->getNFaces();
	}
	if (typeMap.find(aSensorTag)==typeMap.end()){
	  // We have a new sensor geometry
	  typeMap[aSensorTag]=(*modIt);
	}
      }
    }
    
    EndcapModule* anEC;
    int aRing;
    // Look into the endcap sample in order to indentify and measure rings
    for (ModuleVector::iterator moduleIt=endcapSample.begin(); moduleIt!=endcapSample.end(); moduleIt++) {
      if ( (anEC=dynamic_cast<EndcapModule*>(*moduleIt)) ) {
	aRing=anEC->getRing();
	if (ringTypeMap.find(aRing)==ringTypeMap.end()){
	  // We have a new sensor geometry
	  ringTypeMap[aRing]=(*moduleIt);
	}
      } else {
	std::cout << "ERROR: found a non-Endcap module in the map of ring types" << std::endl;
      }
    }
    
    std::ostringstream myName;
    for (std::map<int, Module*>::iterator typeIt = ringTypeMap.begin();
	 typeIt!=ringTypeMap.end(); typeIt++) {
      if ( (anEC=dynamic_cast<EndcapModule*>((*typeIt).second)) ) {
	aRing=(*typeIt).first;
	ringTable->setContent(0,aRing,aRing);
	ringTable->setContent(1,aRing, aRingRho = anEC->getDist(), coordPrecision);
	ringTable->setContent(2,aRing, aRingRho = anEC->getDist()+anEC->getHeight(), coordPrecision);
      } else {
	std::cout << "ERROR: found a non-Endcap module in the map of ring types (twice...)" << std::endl;
      }
    }
    

    // A bit of variables
    std::vector<std::string> names;
    std::vector<std::string> tags;
    std::vector<std::string> types;
    std::vector<std::string> areastrips;
    std::vector<std::string> areapts;
    std::vector<std::string> occupancies;
    std::vector<std::string> pitchpairs;
    std::vector<std::string> striplengths;
    std::vector<std::string> segments;
    std::vector<std::string> nstrips;
    std::vector<std::string> numbermods;
    std::vector<std::string> numbersens;
    std::vector<std::string> channelstrips;
    std::vector<std::string> channelpts;
    std::vector<std::string> powers;
    std::vector<std::string> costs;
    
    double totalPower=0;
    double totalCost=0;
    
    std::ostringstream aName;
    std::ostringstream aTag;
    std::ostringstream aType;
    std::ostringstream anArea;
    std::ostringstream anOccupancy;
    std::ostringstream aPitchPair;
    std::ostringstream aStripLength;
    std::ostringstream aSegment;
    std::ostringstream anNstrips;
    std::ostringstream aNumberMod;
    std::ostringstream aNumberSens;
    std::ostringstream aChannel;
    std::ostringstream aPower;
    std::ostringstream aCost;
    int barrelCount=0;
    int endcapCount=0;
    Module* aModule;
    


    //********************************//
    //*                              *//
    //*       Modules                *//
    //*                              *//
    //********************************//
    myContent = new RootWContent("Modules", false);
    myPage->addContent(myContent);
    RootWTable* moduleTable = new RootWTable(); myContent->addItem(moduleTable);

    static const int tagRow = 1;
    static const int typeRow = 2;
    static const int areastripRow = 3;
    static const int areaptRow = 4;
    static const int occupancyRow = 5;
    static const int pitchpairsRow = 6;
    static const int striplengthRow = 7;
    static const int segmentsRow = 8;
    static const int nstripsRow = 9;
    static const int numbermodsRow = 10;
    static const int numbersensRow = 11;
    static const int channelstripRow = 12;
    static const int channelptRow = 13;
    static const int powerRow = 14;
    static const int costRow = 15;

    // Row names
    moduleTable->setContent(tagRow,0,"Tag");
    moduleTable->setContent(typeRow,0,"Type");
    moduleTable->setContent(areastripRow,0,"Area (mm"+superStart+"2"+superEnd+")");
    moduleTable->setContent(areaptRow,0,"Area (mm"+superStart+"2"+superEnd+")");
    moduleTable->setContent(occupancyRow,0,"Occup (max/av)");
    moduleTable->setContent(pitchpairsRow,0,"Pitch (min/max)");
    moduleTable->setContent(striplengthRow,0,"Strip length");
    moduleTable->setContent(segmentsRow,0,"Segments x Chips");
    moduleTable->setContent(nstripsRow,0,"Chan/Sensor");
    moduleTable->setContent(numbermodsRow,0,"N. mod");
    moduleTable->setContent(numbersensRow,0,"N. sens");
    moduleTable->setContent(channelstripRow,0,"Channels (M)");
    moduleTable->setContent(channelptRow,0,"Channels (M)");
    moduleTable->setContent(powerRow,0,"Power (kW)");
    moduleTable->setContent(costRow,0,"Cost (MCHF)");
    
    int loPitch;
    int hiPitch;
    
    // Summary cycle: prepares the rows cell by cell
    int iType=0;
    for (typeMapIt=typeMap.begin(); typeMapIt!=typeMap.end(); typeMapIt++) {
      ++iType;
      // Name
      aName.str("");
      aModule=(*typeMapIt).second;
      if (dynamic_cast<BarrelModule*>(aModule)) {
	aName << std::dec << "B" << subStart << ++barrelCount << subEnd;
      }
      if (dynamic_cast<EndcapModule*>(aModule)) {
	aName << std::dec << "E" << subStart << ++endcapCount << subEnd;
      }
      // Tag
      aTag.str("");
      aTag << smallStart << aModule->getTag() << smallEnd;
      // Type
      aType.str("");
      aType << (*typeMapIt).second->getType();
      // Area
      anArea.str("");
      anArea << std::dec << std::fixed << std::setprecision(areaPrecision) << (*typeMapIt).second->getArea();
      if ((*typeMapIt).second->getArea()<0) { anArea << "XXX"; }
      // Occupancy
      anOccupancy.str("");
      anOccupancy << std::dec << std::fixed << std::setprecision(occupancyPrecision) <<  typeMapMaxOccupancy[(*typeMapIt).first]*100<< "/" <<typeMapAveOccupancy[(*typeMapIt).first]*100/typeMapCount[(*typeMapIt).first] ; // Percentage
      // Pitches
      aPitchPair.str("");
      loPitch=int((*typeMapIt).second->getLowPitch()*1e3);
      hiPitch=int((*typeMapIt).second->getHighPitch()*1e3);
      if (loPitch==hiPitch) {
	aPitchPair << std::dec << std::fixed << std::setprecision(pitchPrecision) << loPitch;
      } else {
	aPitchPair << std::dec << std::fixed << std::setprecision(pitchPrecision)<< loPitch
		   << "/" << std::fixed << std::setprecision(pitchPrecision) << hiPitch;
      }
      // Strip Lengths
      aStripLength.str("");
      aStripLength << std::fixed << std::setprecision(stripLengthPrecision)
		   << (*typeMapIt).second->getHeight()/(*typeMapIt).second->getNSegments();
      // Segments
      aSegment.str("");
      aSegment << std::dec << (*typeMapIt).second->getNSegments()
	       << "x" << int( (*typeMapIt).second->getNStripAcross() / 128. );
      // Nstrips
      anNstrips.str("");
      anNstrips << std::dec << (*typeMapIt).second->getNChannelsPerFace();
      // Number Mod
      aNumberMod.str("");
      aNumberMod << std::dec << typeMapCount[(*typeMapIt).first];
      // Number Sensor
      aNumberSens.str("");
      aNumberSens << std::dec << typeMapCount[(*typeMapIt).first]*((*typeMapIt).second->getNFaces());
      // Channels
      aChannel.str("");
      aChannel << std::fixed << std::setprecision(millionChannelPrecision)
	       << typeMapCountChan[(*typeMapIt).first] / 1e6 ;
      // Power and cost
      aPower.str("");
      aCost.str("");
      aPower << std::fixed << std::setprecision(powerPrecision) <<
        typeMapCountChan[(*typeMapIt).first] *           // number of channels in type
        1e-3 *                                           // conversion from W to kW
        tracker.getPower((*typeMapIt).second->getReadoutType()); // power consumption in W/channel
      totalPower += typeMapCountChan[(*typeMapIt).first] * 1e-3 * tracker.getPower((*typeMapIt).second->getReadoutType());
      aCost  << std::fixed << std::setprecision(costPrecision) <<
        (*typeMapIt).second->getArea() * 1e-2 *          // area in cm^2
        (*typeMapIt).second->getNFaces() *               // number of faces
        tracker.getCost((*typeMapIt).second->getReadoutType()) * // price in CHF*cm^-2
        1e-6 *                                           // conversion CHF-> MCHF
        typeMapCount[(*typeMapIt).first];                // Number of modules
      totalCost +=(*typeMapIt).second->getArea() * 1e-2 * (*typeMapIt).second->getNFaces() * tracker.getCost((*typeMapIt).second->getReadoutType()) * 1e-6 * typeMapCount[(*typeMapIt).first];
      

      moduleTable->setContent(0,iType,aName.str());
      moduleTable->setContent(tagRow,iType,aTag.str());
      moduleTable->setContent(typeRow,iType,aType.str());
      moduleTable->setContent(occupancyRow,iType,anOccupancy.str());
      moduleTable->setContent(pitchpairsRow,iType,aPitchPair.str());
      moduleTable->setContent(striplengthRow,iType,aStripLength.str());
      moduleTable->setContent(segmentsRow,iType,aSegment.str());
      moduleTable->setContent(nstripsRow,iType,anNstrips.str());
      moduleTable->setContent(numbermodsRow,iType,aNumberMod.str());
      moduleTable->setContent(numbersensRow,iType,aNumberSens.str());
      moduleTable->setContent(powerRow,iType,aPower.str());
      moduleTable->setContent(costRow,iType,aCost.str());
      
      if ((*typeMapIt).second->getReadoutType()==Module::Strip) {
	moduleTable->setContent(channelstripRow,iType,aChannel.str());
	moduleTable->setContent(areastripRow,iType,anArea.str());
	moduleTable->setContent(channelptRow,iType,"--");
	moduleTable->setContent(areaptRow,iType,"--");
      } else {
	moduleTable->setContent(channelstripRow,iType,"--");
	moduleTable->setContent(areastripRow,iType,"--");
	moduleTable->setContent(channelptRow,iType,aChannel.str());
	moduleTable->setContent(areaptRow,iType,anArea.str());
      }
    }

    // Score totals
    ++iType;
    moduleTable->setContent(0,iType,"Total");
    moduleTable->setContent(tagRow,iType,"");
    moduleTable->setContent(typeRow,iType,"");
    anArea.str("");
    anArea << emphStart << std::fixed << std::setprecision(areaPrecision) << totAreaStrips/1e6
	   << "(m" << superStart << "2" << superEnd << ")" << emphEnd;
    moduleTable->setContent(areastripRow,iType,anArea.str());
    anArea.str("");
    anArea << emphStart << std::fixed << std::setprecision(areaPrecision) << totAreaPts/1e6
	   << "(m" << superStart << "2" << superEnd << ")" << emphEnd;
    moduleTable->setContent(areaptRow,iType,anArea.str());
    moduleTable->setContent(occupancyRow,iType,"");
    moduleTable->setContent(pitchpairsRow,iType,"");
    moduleTable->setContent(striplengthRow,iType,"");
    moduleTable->setContent(segmentsRow,iType,"");
    moduleTable->setContent(nstripsRow,iType,"");
    aNumberMod.str("");
    aNumberMod << emphStart << totCountMod << emphEnd;
    aNumberSens.str("");
    aNumberSens << emphStart << totCountSens << emphEnd;
    moduleTable->setContent(numbermodsRow,iType,aNumberMod.str());
    moduleTable->setContent(numbersensRow,iType,aNumberSens.str());
    aChannel.str("");
    aChannel << emphStart << std::fixed
	     << std::setprecision(millionChannelPrecision)
	     << totChannelStrips / 1e6 << emphEnd;
    moduleTable->setContent(channelstripRow,iType,aChannel.str());
    aChannel.str("");
    aChannel << emphStart << std::fixed
	     << std::setprecision(millionChannelPrecision)
	     << totChannelPts / 1e6 << emphEnd;
    moduleTable->setContent(channelptRow,iType,aChannel.str());
    aPower.str("");
    aCost.str("");
    aPower   << std::fixed << std::setprecision(powerPrecision) << totalPower;
    aCost    << std::fixed << std::setprecision(costPrecision) << totalCost;
    moduleTable->setContent(powerRow,iType,aPower.str());
    moduleTable->setContent(costRow,iType,aCost.str());

    //********************************//
    //*                              *//
    //*       Plots                  *//
    //*                              *//
    //********************************//
    RootWImage* myImage;
    TCanvas *summaryCanvas = NULL;
    TCanvas *YZCanvas = NULL;
    TCanvas *myCanvas = NULL;
    createSummaryCanvas(tracker.getMaxL(), tracker.getMaxR(), analyzer, summaryCanvas, YZCanvas);
    //TVirtualPad* myPad;
    myContent = new RootWContent("Plots");
    myPage->addContent(myContent);

    if (summaryCanvas) {
      myImage = new RootWImage(summaryCanvas, 800, 800);
      myImage->setComment("Tracker summary: modules position in XY (endcap and barrel), YZ and number of hits vs. eta");
      myContent->addItem(myImage);
    }

    if (YZCanvas) {
      myImage = new RootWImage(YZCanvas, 800, 800);
      myImage->setComment("YZ Section of the tracker barrel");
      myContent->addItem(myImage);
    }

    /*
    myCanvas = new TCanvas("XYViewBarrel", "XYViewBarrel", 800, 800);
    myCanvas->cd();
    myPad = summaryCanvas->GetPad(padXY);
    if (myPad) {
      myPad->DrawClonePad();
      myImage = new RootWImage(myCanvas, 800, 800);
      myImage->setComment("XY Section of the tracker barrel");
      myContent->addItem(myImage);
    }

    myCanvas = new TCanvas("XYViewEndcap", "XYViewEndcap", 800, 800);
    myCanvas->cd();
    myPad = summaryCanvas->GetPad(padEC);
    if (myPad) {
      myPad->DrawClonePad();
      myImage = new RootWImage(myCanvas, 800, 800);
      myImage->setComment("XY View of the tracker endcap");
      myContent->addItem(myImage);
    }
    */

    myCanvas = new TCanvas("EtaProfile", "Eta profile", 800, 800);
    myCanvas->cd();
    analyzer.getEtaProfileCanvas().DrawClonePad();
    myImage = new RootWImage(myCanvas, 800, 800);
    myImage->setComment("Hit coverage in eta");
    myContent->addItem(myImage);

    TCanvas* hitMapCanvas = new TCanvas("hitmapcanvas", "Hit Map", 800, 800);
    int prevStat = gStyle->GetOptStat();
    gStyle->SetOptStat(0);
    hitMapCanvas->cd();
    gStyle->SetPalette(1);
    hitMapCanvas->SetFillColor(color_plot_background);
    hitMapCanvas->SetBorderMode(0);
    hitMapCanvas->SetBorderSize(0);
    analyzer.getMapPhiEta().Draw("colz");
    hitMapCanvas->Modified();
    gStyle->SetOptStat(prevStat);
    myImage = new RootWImage(hitMapCanvas, 800, 800);
    myImage->setComment("Hit coverage in eta, phi");
    myContent->addItem(myImage);

    // TODO: bandwidth
    // (also todo: handle this properly)
    /*
    if (bandWidthCanvas_) {
      pngFileName = fileName+"_bandwidth.png";
      svgFileName = fileName+"_bandwidth.svg";
      bandWidthCanvas_->DrawClonePad();
      bandWidthCanvas_->SaveAs(pngFileName.c_str());
      bandWidthCanvas_->SaveAs(svgFileName.c_str());
    }
    myfile << "<h5>Bandwidth useage estimate:</h5>" << std::endl;
    myfile << "<p>(Pt modules: ignored)</p>" << std::endl
	   << "<p>Sparsified (binary) bits/event: 23 bits/chip + 9 bit/hit</p>" << std::endl
	   << "<p>Unsparsified (binary) bits/event: 16 bits/chip + 1 bit/channel</p>" << std::endl
	   << "<p>100 kHz trigger, " << nMB_ << " minimum bias events assumed</p>" << std::endl;
    myfile << "<img src=\"summaryPlots_bandwidth.png\" />" << std::endl;
    myfile << "</body></html>" << std::endl; */

    // TODO: Conditions and files
    /* 
       myfile << clearStart << emphStart << "Geometry configuration file:     " << emphEnd
       << "<a href=\".//" << configFile << "\">" << configFile << "</a>" << clearEnd << "<br/>" << std::endl;
       myfile << clearStart << emphStart << "Module types configuration file: " << emphEnd
       << "<a href=\".//" << dressFile << "\">" << dressFile << "</a>" << clearEnd << "<br/>" << std::endl;
       myfile << clearStart << emphStart << "Minimum bias per bunch crossing: " << emphEnd
       << nMB_ << clearEnd << std::endl;
       if (barrelModuleCoordinatesFile!="")
       myfile << "<br/>" << clearStart << emphStart << "Barrel modules coordinate file:  " << emphEnd
       << "<a href=\".//" << barrelModuleCoordinatesFile << "\">" << barrelModuleCoordinatesFile << "</a>" << clearEnd << std::endl;
       if (endcapModuleCoordinatesFile!="")
       myfile << "<br/>" << clearStart << emphStart << "End-cap modules coordinate file:  " << emphEnd
       << "<a href=\".//" << endcapModuleCoordinatesFile << "\">" << endcapModuleCoordinatesFile << "</a>" << clearEnd << std::endl;
       
       // TODO: make an object that handles this properly:
       myfile << "<h5>Cost estimates:</h5>" << std::endl;
       myfile << "<p>Pt modules: "
       << std::fixed << std::setprecision(costPerUnitPrecision)
       << getCost(Module::Pt) << " CFH/cm2 - Strip modules: "
       << std::fixed << std::setprecision(costPerUnitPrecision)
       << getCost(Module::Strip) << " CHF/cm2</p>" << std::endl;
       myfile << "<h5>Power estimates:</h5>" << std::endl;
       myfile << "<p>Pt modules: "
       << std::fixed << std::setprecision(powerPerUnitPrecision)
       << getPower(Module::Pt)*1e3 << " mW/chan - Strip modules: "
       << std::fixed << std::setprecision(powerPerUnitPrecision)
       << getPower(Module::Strip)*1e3 << " mW/chan</p>" << std::endl;
  */


    //summaryCanvas->SaveAs(epsFileName.c_str());
    //summaryCanvas->SaveAs(gifFileName.c_str());

    // TODO: make this meaningful!
    return true;
    
  }

#endif


  // private
  // Draws tickmarks on 3d canvases
  // @param myView the TView where to draw ticks
  // @param maxL maximum tracker length in z
  // @param maxR maximum tracker radius in rho
  // @param noAxis number of the axis: Enumerate sections by axis
  //        index normal to draw plane (if x=1, y=2, z=3)
  // @param spacing grid tick spacing
  // @param option the options to pass to the Draw() method
  void Vizard::drawTicks(TView* myView, double maxL, double maxR, int noAxis/*=1*/, double spacing /*= 100.*/, Option_t* option /*= "same"*/) {
    TPolyLine3D* aLine;
    Color_t gridColor_hard = color_hard_grid;
    int gridStyle_solid = 1;
    std::string theOption(option);
    
    int i;
    int j;
    int k;
    
    double topMax = (maxL > maxR) ? maxL : maxR;
    topMax = ceil(topMax/spacing)*spacing;
    
    double aValue[3];
    double minValue[3];
    double maxValue[3];
    
    i=(noAxis)%3;
    j=(noAxis+1)%3;
    k=(noAxis+2)%3;
    
    maxL *= 1.1;
    maxR *= 1.1;
    
    if (noAxis==1) {
        minValue[0]=0;
        maxValue[0]=+maxR;
        minValue[1]=0;
        maxValue[1]=+maxR;
        minValue[2]=0;
        maxValue[2]=+maxL;
    } else {
        minValue[0]=-maxR;
        maxValue[0]=+maxR;
        minValue[1]=-maxR;
        maxValue[1]=+maxR;
        minValue[2]=0;
        maxValue[2]=+maxL;
    }
    
    aValue[k]=-topMax;
    
    if (noAxis==1) {
      double etaStep=.2;
      double etaMax = 2.1;
      // Add the eta ticks
      double theta;
      double tickLength = 2 * spacing;
      double tickDistance = spacing;
      double startR = maxR + tickDistance;
      double startL = maxL + tickDistance;
      double endR = maxR + tickDistance + tickLength;
      double endL = maxL + tickDistance + tickLength;
      XYZVector startTick;
      XYZVector endTick;
      Double_t pw[3];
      Double_t pn[3];
      TText* aLabel;
      char labelChar[10];
      for (double eta=0; eta<etaMax; eta+=etaStep) {
	aLine = new TPolyLine3D(2);
	theta = 2 * atan(exp(-eta));
	startTick = XYZVector(0, sin(theta), cos(theta));
	startTick *= startR/startTick.Rho();
	endTick = startTick / startTick.Rho() * endR;
	if (startTick.Z()>startL) {
	  startTick *= startL/startTick.Z();
	  endTick *=  endL/endTick.Z();
	}
	pw[0]=0.;
	pw[1]=endTick.Y();
	pw[2]=endTick.Z();
	myView->WCtoNDC(pw, pn);
	sprintf(labelChar, "%.01f", eta);
	aLabel = new TText(pn[0], pn[1], labelChar);
	aLabel->SetTextSize(aLabel->GetTextSize()*.6);
	aLabel->SetTextAlign(21);
	aLabel->Draw(theOption.c_str());
	theOption="same";
	endTick = (endTick+startTick)/2.;
	aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
	aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
	aLine->SetLineStyle(gridStyle_solid);
	aLine->SetLineColor(gridColor_hard);
	aLine->Draw("same");
      }
      aLine = new TPolyLine3D(2);
      theta = 2 * atan(exp(-2.5));
      startTick = XYZVector(0, sin(theta), cos(theta));
      startTick *= startR/startTick.Rho();
      endTick = startTick / startTick.Rho() * endR;
      if (startTick.Z()>startL) {
	startTick *= startL/startTick.Z();
	endTick *=  endL/endTick.Z();
      }
      pw[0]=0.;
      pw[1]=endTick.Y();
      pw[2]=endTick.Z();
      myView->WCtoNDC(pw, pn);
      sprintf(labelChar, "%.01f", 2.5);
      aLabel = new TText(pn[0], pn[1], labelChar);
      aLabel->SetTextSize(aLabel->GetTextSize()*.8);
      aLabel->SetTextAlign(21);
      aLabel->Draw(theOption.c_str());
      theOption="same";
      endTick = (endTick+startTick)/2.;
      aLine->SetPoint(0, 0., 0., 0.);
      aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
      aLine->SetLineStyle(gridStyle_solid);
      aLine->SetLineColor(gridColor_hard);
      aLine->Draw("same");
      
      
      
      for (double z=0; z<=maxL ; z+=(4*spacing)) {
	aLine = new TPolyLine3D(2);
	startTick = XYZVector(0, 0, z);
	endTick = XYZVector(0, -(tickLength/2), z);
	aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
	aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
	pw[0]=0.;
	pw[1]=-tickLength;
	pw[2]=endTick.Z();
	myView->WCtoNDC(pw, pn);
	sprintf(labelChar, "%.0f", z);
	aLabel = new TText(pn[0], pn[1], labelChar);
	aLabel->SetTextSize(aLabel->GetTextSize()*.6);
	aLabel->SetTextAlign(23);
	aLabel->Draw(theOption.c_str());
	theOption="same";
	aLine->SetLineStyle(gridStyle_solid);
	aLine->SetLineColor(gridColor_hard);
	aLine->Draw("same");
      }
      
      for (double y=0; y<=maxR ; y+=(2*spacing)) {
	aLine = new TPolyLine3D(2);
	startTick = XYZVector(0, y, 0);
	endTick = XYZVector(0, y, -(tickLength/2));
	aLine->SetPoint(0, 0., startTick.Y(), startTick.Z());
	aLine->SetPoint(1, 0., endTick.Y(), endTick.Z());
	pw[0]=0.;
	pw[1]=endTick.Y();
	pw[2]=-tickLength;
	myView->WCtoNDC(pw, pn);
	sprintf(labelChar, "%.0f", y);
	aLabel = new TText(pn[0], pn[1], labelChar);
	aLabel->SetTextSize(aLabel->GetTextSize()*.6);
	aLabel->SetTextAlign(32);
	aLabel->Draw(theOption.c_str());
	theOption="same";
	aLine->SetLineStyle(gridStyle_solid);
	aLine->SetLineColor(gridColor_hard);
	aLine->Draw("same");
      }
    }
  }

  // private
  // Draws a grid on the current canvas
  // @param maxL maximum tracker length in z
  // @param maxR maximum tracker radius in rho
  // @param noAxis number of the axis: Enumerate sections by axis
  //        index normal to draw plane (if x=1, y=2, z=3)
  // @param spacing grid tick spacing
  // @param option the options to pass to the Draw() method
  void Vizard::drawGrid(double maxL, double maxR, int noAxis/*=1*/, double spacing /*= 100.*/, Option_t* option /*= "same"*/) {
    TPolyLine3D* aLine;
    Color_t gridColor = color_grid;
    Color_t gridColor_hard = color_hard_grid;
    Color_t thisLineColor;
    
    std::string theOption(option);
    
    int i;
    int j;
    int k;
    
    double topMax = (maxL > maxR) ? maxL : maxR;
    topMax = ceil(topMax/spacing)*spacing;
    
    double aValue[3];
    double minValue[3];
    double maxValue[3];
    double runValue;
    int thisLineStyle;
    
    i=(noAxis)%3;
    j=(noAxis+1)%3;
    k=(noAxis+2)%3;
    
    maxL *= 1.1;
    maxR *= 1.1;
    
    if (noAxis==1) {
        minValue[0]=0;
        maxValue[0]=+maxR;
        minValue[1]=0;
        maxValue[1]=+maxR;
        minValue[2]=0;
        maxValue[2]=+maxL;
    } else {
        minValue[0]=-maxR;
        maxValue[0]=+maxR;
        minValue[1]=-maxR;
        maxValue[1]=+maxR;
        minValue[2]=0;
        maxValue[2]=+maxL;
    }
    
    aValue[k]=-topMax;
    for(runValue = -topMax; runValue<=topMax; runValue+=spacing) {
      
      // Special line for axis
      if (runValue==0) {
	thisLineStyle=1;
	thisLineColor=gridColor_hard;
      } else {
	thisLineStyle=2;
	thisLineColor=gridColor;
      }
      
      // Parallel to j
      if ((runValue<=maxValue[i])&&(runValue>=minValue[i])) {
	aValue[i] = runValue;
	aLine = new TPolyLine3D(2);
	aValue[j] = minValue[j];
	aLine->SetPoint(0, aValue[0], aValue[1], aValue[2]);
	aValue[j] = maxValue[j];
	aLine->SetPoint(1, aValue[0], aValue[1], aValue[2]);
	aLine->SetLineStyle(thisLineStyle);
	aLine->SetLineColor(thisLineColor);
	aLine->Draw(theOption.c_str());
	theOption="same";
      };
      
      // Parallel to i
      if ((runValue<=maxValue[j])&&(runValue>=minValue[j])) {
	aValue[j] = runValue;
	aLine = new TPolyLine3D(2);
	aValue[i] = minValue[i];
	aLine->SetPoint(0, aValue[0], aValue[1], aValue[2]);
	aValue[i] = maxValue[i];
	aLine->SetPoint(1, aValue[0], aValue[1], aValue[2]);
	aLine->SetLineStyle(thisLineStyle);
	aLine->SetLineColor(thisLineColor);
	aLine->Draw(theOption.c_str());
	theOption="same";
      };
      
    } 
  }
  
  // private
  // Creates a new 4-pad canvas with XY and YZ views with all the useful details, like the axis ticks
  // and the eta reference. The fourth pad contains a miniature of the eta profile coverage
  // if you need any of these you can get them with GetPad()
  // @param maxZ maximum tracker's Z coordinate to be shown
  // @param maxRho maximum tracker's Rho coordinate to be shown
  // @param analyzer A reference to the analysing class that examined the material budget and filled the histograms 
  // @return a pointer to the new TCanvas
  void Vizard::createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer, TCanvas* summaryCanvas, TCanvas* YZCanvas) {
    Int_t irep;
    TVirtualPad* myPad;

    YZCanvas = new TCanvas("YZCanvas", "YZView Canvas", 800, 800 );
    summaryCanvas = new TCanvas("summaryCanvas", "Summary Canvas", 800, 800);
    summaryCanvas->SetFillColor(color_pad_background);
    summaryCanvas->Divide(2, 2);

    for (int i=1; i<=4; i++) { myPad=summaryCanvas->GetPad(i); myPad->SetFillColor(color_plot_background);  }

    // First pad
    // YZView
    myPad = summaryCanvas->GetPad(padYZ);
    myPad->SetFillColor(color_plot_background);
    myPad->cd();
    if (analyzer.getGeomLiteYZ()) {
      drawGrid(maxZ, maxRho, ViewSectionYZ);
      analyzer.getGeomLiteYZ()->DrawClonePad();
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
      myPad->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
      drawTicks(myPad->GetView(), maxZ, maxRho, ViewSectionYZ);      

      YZCanvas->cd();
      myPad = YZCanvas->GetPad(0);
      myPad->SetFillColor(color_plot_background);
      drawGrid(maxZ, maxRho, ViewSectionYZ);
      analyzer.getGeomLiteYZ()->DrawClonePad();
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
      myPad->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
      drawTicks(myPad->GetView(), maxZ, maxRho, ViewSectionYZ);
    }

    // Second pad
    // XYView (barrel)
    myPad = summaryCanvas->GetPad(padXY);
    myPad->cd();
    myPad->SetFillColor(color_plot_background);
    if (analyzer.getGeomLiteXY()) {
      drawGrid(maxZ, maxRho, ViewSectionXY);
      analyzer.getGeomLiteXY()->DrawClonePad();
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
      myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
    }
    
    // Third pad
    // Plots
    myPad = summaryCanvas->GetPad(padProfile);
    myPad->cd();
    myPad->SetFillColor(color_plot_background);
    analyzer.getEtaProfileCanvas().DrawClonePad();
    
    // Fourth pad
    // XYView (EndCap)
    myPad = summaryCanvas->GetPad(padEC);
    myPad->cd();
    myPad->SetFillColor(color_plot_background);
    if (analyzer.getGeomLiteEC()) {
      drawGrid(maxZ, maxRho, ViewSectionXY);
      analyzer.getGeomLiteEC()->DrawClonePad();
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
      myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
    }
    
    for (int i=1; i<=4; i++) { myPad=summaryCanvas->GetPad(i); myPad->SetBorderMode(0); }
    summaryCanvas->Modified();   
    //return summaryCanvas;
  }
  
  

}

