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
        Module* mod;
        Layer* current;
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
        else {
            if (!am.getBarrelLayers()->empty() && !am.getBarrelLayers()->at(0)->getModuleVector()->empty()) {
                mod = am.getBarrelLayers()->at(0)->getModuleVector()->at(0);
                vol = gm->MakeArb8("", medact, 0);
                vol->SetLineColor(kRed);
                for (unsigned int i = 0; i < am.getBarrelLayers()->size(); i++) {
                    current = am.getBarrelLayers()->at(i);
                    for (unsigned int j = 0; j < current->getModuleVector()->size(); j++) {
                        mod = current->getModuleVector()->at(j);
                        trafo = modulePlacement(mod, vol);
                        barrels->AddNode(vol, c, trafo);
                        c++;
                    }
                }
            }
        }
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
        else {
            if (!am.getEndcapLayers()->empty() && !am.getEndcapLayers()->at(0)->getModuleVector()->empty()) {
                // TODO: for each endcap:
                // for each ring:
                // define module shape
                // create volume("", shape, medact)
                // for all modules in this ring, in every disc:
                // define transform trafo as a combination of trans and rot
                // endcaps->AddNode(vol, some_index, trafo);
                // c++;
            }
        }
        
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
                std::cout << "Something went wrong creating output file.";
                std::cout << " Existing geometry was not written to file." << std::endl;
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
        else std::cout << "Vizard::buildVisualization(am, is) needs to be called first to build the visual geometry objects." << std::endl;
    }
    
    void Vizard::display(Tracker& am, InactiveSurfaces& is, std::string rootfilename, bool simplified) {
        buildVisualization(am, is, simplified);
        display(rootfilename);
    }
    
    void Vizard::writeNeighbourGraph(InactiveSurfaces& is) {
        writeNeighbourGraph(is, default_graphfile);
    }
    
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
            else std::cout << "File stream reported error state: neighbour graph not written to file." << std::endl;
        }
        catch (std::bad_alloc ba) {
            std::cerr << "Error: caught bad_alloc exception in Vizard::writeNeighbourGraph(). ";
            std::cerr << "Neighbour graph was not written to file." << std::endl;
        }
    }
    
    void Vizard::dotGraph(InactiveSurfaces& is, std::string outfile) {
        const std::string preamble = "digraph tracker";
        const std::string ori = "rankdir=DU"; // check if this is possible!
        const std::string shape = "node [shape=box]";
        const std::string label = "label=";
        const std::string edge = "->";
    }
    
    void Vizard::histogramSummary(Analyzer& a, std::string outfilename) {
        TH1D container;
        std::string outfile = default_summarypath + "/";
        if (outfilename.empty()) outfile = outfile + default_summary;
        else outfile = outfile + outfilename;
        std::string pngout = outfile + ".png";
        std::string htmlcontents;
        std::ofstream outstream(outfile.c_str());
        TCanvas c("matbudgetcanvas", "Material Budgets over Eta", 800, 800);
        c.Divide(2, 2);
        c.cd(0);
        a.getHistoGlobalR().Draw();
        c.cd(1);
        a.getHistoGlobalI().Draw();
        c.cd(2);
        // TODO: add correctly
        container.Add(&(a.getHistoSupportsAllR()));
        container.Add(&(a.getHistoServicesAllR()));
        container.Add(&(a.getHistoModulesAllR()));
        container.DrawCopy();
        container.Reset();
        container.Add(&(a.getHistoSupportsAllR()));
        container.Add(&(a.getHistoServicesAllR()));
        container.DrawCopy();
        container.Reset();
        container.Add(&(a.getHistoSupportsAllR()));
        container.DrawCopy();
        c.cd(3);
        // TODO: add correctly
        container.Reset();
        container.Add(&(a.getHistoSupportsAllI()));
        container.Add(&(a.getHistoServicesAllI()));
        container.Add(&(a.getHistoModulesAllI()));
        container.DrawCopy();
        container.Reset();
        container.Add(&(a.getHistoSupportsAllI()));
        container.Add(&(a.getHistoServicesAllI()));
        container.DrawCopy();
        container.Reset();
        container.Add(&(a.getHistoSupportsAllI()));
        container.DrawCopy();
        c.SaveAs(pngout.c_str());
        htmlcontents = "<html><title>" + outfilename + "</title><body>";
        htmlcontents = htmlcontents + "<img src=\"" + pngout + "\" />" + "</body></html>";
        outstream << htmlcontents << std::endl;
        outstream.close();
        std::cout << "HTML file written to " << outfile << ", image written to " << pngout << std::endl;
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
