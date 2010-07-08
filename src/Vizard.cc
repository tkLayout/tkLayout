#include <Vizard.h>
#include <TStyle.h>
#include <TColor.h>

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
        c.SetFillColor(kWhite);
        c.Divide(2,1);
        pad = c.GetPad(0);
        pad->SetFillColor(kGray);
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

    /**
     * This function draws some of the histograms that were filled during material budget analysis and
     * embeds the resulting image in an HTML file for easy access. If no name is given for the output file,
     * a default filename is used.
     * @param a A reference to the analysing class that examined the material budget and filled the histograms
     * @param outfilename The name of the output file that will be written to the application's default directory for material budget summaries
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
      myCanvas->SetFillColor(kWhite);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(kGray);
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
      myCanvas->SetFillColor(kWhite);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(kGray);
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
      myCanvas->SetFillColor(kWhite);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(kGray);
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
      myCanvas->SetFillColor(kWhite);
      myCanvas->Divide(2,1);
      myPad = myCanvas->GetPad(0);
      myPad->SetFillColor(kGray);
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
}
