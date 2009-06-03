#include <Vizard.h>
namespace insur {
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
    
    Vizard::~Vizard() {
        delete gm;
    }
    
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
    
    void Vizard::display(Tracker& am, InactiveSurfaces& is, std::string rootfilename, bool simplified) {
        buildVisualization(am, is, simplified);
        display(rootfilename);
    }
    
    void Vizard::writeNeighbourGraph(InactiveSurfaces& is) {
        writeNeighbourGraph(is, default_graphfile);
    }
    
    void Vizard::writeNeighbourGraph(InactiveSurfaces& is, std::string outfile) {
        if (geometry_created) {
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
        else std::cout << msg_uninitialised << std::endl;
    }
    
    void Vizard::dotGraph(InactiveSurfaces& is, std::string outfile) {
        const std::string preamble = "digraph tracker";
        const std::string ori = "rankdir=DU"; // check if this is possible!
        const std::string shape = "node [shape=box]";
        const std::string label = "label=";
        const std::string edge = "->";
    }
    
    void Vizard::histogramSummary(Analyzer& a, std::string outfilename) {
        THStack rcontainer("rstack", "Radiation Length by Category");
        THStack icontainer("istack", "Interaction Length by Category");
        TH1D *cr = NULL, *ci = NULL, *acr = NULL, *aci = NULL, *ser = NULL, *sei = NULL, *sur = NULL, *sui = NULL;
        TVirtualPad* pad;
        std::string outfile = default_summarypath + "/";
        std::string pngoutfile, pngout;
        if (outfilename.empty()) {
            outfile = outfile + default_summary;
            pngoutfile = default_summary + ".png";
        }
        else {
            outfile = outfile + outfilename;
            pngoutfile = outfilename + ".png";
        }
        pngout = default_summarypath + "/" + pngoutfile;
        std::string htmlcontents;
        std::ofstream outstream(outfile.c_str());
        TCanvas c("matbudgetcanvas", "Material Budgets over Eta", 800, 800);
        c.SetFillColor(kWhite);
        c.Divide(2, 2);
        pad = c.GetPad(0);
        pad->SetFillColor(kGray);
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
        pad = c.GetPad(3);
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
        pad = c.GetPad(4);
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
        c.SaveAs(pngout.c_str());
        if (cr) delete cr;
        if (ci) delete ci;
        if (acr) delete acr;
        if (aci) delete aci;
        if (ser) delete ser;
        if (sei) delete sei;
        if (sur) delete sur;
        if (sui) delete sui;
        htmlcontents = "<html><title>" + outfilename + "</title><body>";
        htmlcontents = htmlcontents + "<img src=\"" + pngoutfile + "\" />" + "</body></html>";
        outstream << htmlcontents << std::endl;
        outstream.close();
        std::cout << "HTML file written to " << outfile << ", image written to " << pngout << std::endl;
    }
    
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
}
