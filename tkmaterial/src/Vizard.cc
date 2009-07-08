#include <Vizard.h>
namespace insur {
    // public
    /**
     * The constructor builds the logical structure within the <i>TGeoManager</i> that is used to display
     * a tracker geometry in ROOT. It assigns different materials and media to the various categories of
     * geometry elements.
     */
    Vizard::Vizard() {
        geometry_created = false;
        gm = new TGeoManager("display", "Tracker");
        matvac = new TGeoMaterial("Vacuum", 0, 0, 0);
        matact = new TGeoMaterial("Si", a_silicon, z_silicon, d_silicon);
        matserf = new TGeoMaterial("C ", a_carbon, z_carbon, d_carbon);
        matlazy = new TGeoMaterial("Cu", a_copper, z_copper, d_copper);
        medvac = new TGeoMedium("Vacuum", 0, matvac);
        medact = new TGeoMedium("Silicon", 1, matact);
        medserf = new TGeoMedium("Copper", 2, matserf);
        medlazy = new TGeoMedium("Carbon", 3, matlazy);
        barrels = new TGeoVolumeAssembly("Barrels");
        endcaps = new TGeoVolumeAssembly("Endcaps");
        services = new TGeoVolumeAssembly("Services");
        supports = new TGeoVolumeAssembly("Supports");
        active = new TGeoVolumeAssembly("Active Modules");
        inactive = new TGeoVolumeAssembly("Inactive Surfaces");
        top = gm->MakeBox("WORLD", medvac, outer_radius + top_volume_pad, outer_radius + top_volume_pad, max_length + top_volume_pad);
        active->AddNode(barrels, 0);
        active->AddNode(endcaps, 0);
        inactive->AddNode(services, 0);
        inactive->AddNode(supports, 0);
        top->AddNode(active, 0);
        top->AddNode(inactive, 0);
        gm->SetTopVolume(top);
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
        TGeoVolume* vol;
        TGeoTranslation* trans;
        TGeoCombiTrans* trafo;
        Layer* current;
        std::vector<Module*> templates;
        // barrels
        if (simplified) {
            for (unsigned int i = 0; i < am.getBarrelLayers()->size(); i++) {
                current = am.getBarrelLayers()->at(i);
                if ((current->getMinZ() > 0) || (current->getMaxZ() < 0)) {
                    vol = gm->MakeTube("", medact, current->getMinRho(), current->getMaxRho(), (current->getMaxZ() - current->getMinZ()) / 2.0);
                    vol->SetLineColor(kRed);
                    trans = new TGeoTranslation(0, 0, current->getMaxZ() - (current->getMaxZ() - current->getMinZ()) / 2.0);
                    barrels->AddNode(vol, c, trans);
                }
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
        
        // supports => inner and outer support tube will not be displayed
        skip = is.getSupports().size() - 2;
        for (int i = 2; i < skip + 2; i++) {
            if ((is.getSupportPart(i).getZOffset() + is.getSupportPart(i).getZLength()) > 0) {
                vol = gm->MakeTube("", medlazy, is.getSupportPart(i).getInnerRadius(),
                        is.getSupportPart(i).getInnerRadius() + is.getSupportPart(i).getRWidth(),
                        is.getSupportPart(i).getZLength() / 2.0);
                vol->SetLineColor(kGray);
                trans = new TGeoTranslation(0, 0, (is.getSupportPart(i).getZOffset() + is.getSupportPart(i).getZLength() / 2.0));
                supports->AddNode(vol, c, trans);
                c++;
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
        std::string outfile = default_summarypath + "/";
        std::string pngoutfile, pngout, pngpath;
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
        pngout = pngoutfile + ".fullvolume.png";
        pngpath = default_summarypath + "/" + pngout;
        c.SaveAs(pngpath.c_str());
        htmlstream << "<img src=\"" << pngout << "\" /><br><p><small><b>Average radiation length in full volume ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*cr, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in full volume (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*ci, etaMaxAvg) << "</b></small></p>";
        c.Clear("D");
        if (cr) delete cr;
        if (ci) delete ci;
        pad = c.GetPad(1);
        pad->cd();
        cr = (TH1D*)a.getHistoGlobalR().Clone();
        cr->SetFillColor(kGray + 2);
        cr->SetNameTitle("rglobal", "Overall Radiation Length");
        cr->SetXTitle("Eta");
        cr->Draw();
        pad = c.GetPad(2);
        pad->cd();
        ci = (TH1D*)a.getHistoGlobalI().Clone();
        ci->SetFillColor(kGray + 1);
        ci->SetNameTitle("iglobal", "Overall Interaction Length");
        ci->SetXTitle("Eta");
        ci->Draw();
        pngout = pngoutfile + ".global.png";
        pngpath = default_summarypath + "/" + pngout;
        c.SaveAs(pngpath.c_str());
        htmlstream << "<br><img src=\"" << pngout << "\" /><br><p><small><b>Average radiation length in tracking volume ";
        htmlstream << "(eta = [0, 2.4]): " << averageHistogramValues(*cr, etaMaxAvg) << "</b></small></p>";
        htmlstream << "<p><small><b>Average interaction length in tracking volume (eta = [0, 2.4]): ";
        htmlstream << averageHistogramValues(*ci, etaMaxAvg) << "</b></small></p>";
        c.Clear("D");
        pad = c.GetPad(1);
        pad->cd();
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
        icontainer.Add(aci);
        aci->SetXTitle("Eta");
        icontainer.Draw();
        pngout = pngoutfile + ".asl.png";
        pngpath = default_summarypath + "/" + pngout;
        c.SaveAs(pngpath.c_str());
        htmlstream << "<img src=\"" << pngout << "\" /><br>";
        //TODO: add information about average rlength and ilength per category a, s or l
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
        c.Clear("D");
        pad = c.GetPad(1);
        pad->cd();
        ir = (TH2D*)a.getHistoIsoR().Clone();
        ir->SetNameTitle("isor", "Radiation Length Contours");
        ir->SetContour(temperature_levels);
        ir->SetXTitle("z");
        ir->SetYTitle("r");
        ir->Draw("COLZ");
        pad = c.GetPad(2);
        pad->cd();
        ii = (TH2D*)a.getHistoIsoI().Clone();
        ii->SetNameTitle("isoi", "Interaction Length Contours");
        ii->SetContour(temperature_levels);
        ii->SetXTitle("z");
        ii->SetYTitle("r");
        ii->Draw("COLZ");
        pngout = pngoutfile + ".twodee.png";
        pngpath = default_summarypath + "/" + pngout;
        c.SaveAs(pngpath.c_str());
        htmlstream << "<br><p><big><b>2D Overview</b></big></p><img src=\"" << pngout << "\" /><br>";
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
        htmlstream << "</body></html>";
        outstream << htmlstream.str() << std::endl;
        outstream.close();
        std::cout << "HTML file written to " << outfile << std::endl;
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
            std::cout << "detailedModules(): layers vector is not empty." << std::endl;
            v = gm->MakeArb8("", medact, 0);
            v->SetLineColor(kRed);
            for (unsigned int i = 0; i < layers->size(); i++) {
                std::cout << "detailedModules(): layer " << i << std::endl;
                current = layers->at(i);
                for (unsigned int j = 0; j < current->getModuleVector()->size(); j++) {
                    mod = current->getModuleVector()->at(j);
                    t = modulePlacement(mod, v);
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
        b = m->getCorner(1) - m->getCorner(0);
        c = m->getCorner(2) - m->getCorner(0);
        d = m->getCorner(3) - m->getCorner(0);
        ex = b / b.R();
        p = (d.Dot(ex) * ex);
        ey = d - p;
        ey = ey / ey.R();
        ez = ex.Cross(ey);
        arb = (TGeoArb8*)(v->GetShape());
        for (int i = 0; i < 5; i = i + 4) {
            arb->SetVertex(i, 0, 0);
            arb->SetVertex(i + 1, b.R(), 0);
            arb->SetVertex(i + 2, c.Dot(ex), c.Dot(ey));
            arb->SetVertex(i + 3, d.Dot(ex), d.Dot(ey));
        }
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
        while ((cobin < histo.GetNbinsX()) && (histo.GetBinLowEdge(cobin) < cutoff)) cobin++;
        if (cobin >= histo.GetNbinsX() - 1) avg = histo.GetMean();
        else {
            for (int i = 1; i <= cobin; i++) avg = avg + histo.GetBinContent(i) / (double)cobin;
        }
        return avg;
    }
}
