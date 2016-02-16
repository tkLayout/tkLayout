/**
 * @file Vizard.cc
 * @brief This class takes care of visualisation for both geometry and analysis results
 */

#include <Vizard.h>
#include <TPolyLine.h>
#include <Units.h>

#include "mainConfigHandler.h"

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
    matact = new TGeoMaterial("Si", mat_a_silicon, mat_z_silicon, mat_d_silicon);
    matserf = new TGeoMaterial("C ", mat_a_carbon, mat_z_carbon, mat_d_carbon);
    matlazy = new TGeoMaterial("Cu", mat_a_copper, mat_z_copper, mat_d_copper);
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
    top = gm->MakeBox("WORLD", medvac, geom_max_radius + geom_top_volume_pad, geom_max_radius + geom_top_volume_pad, geom_max_length + geom_top_volume_pad);
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
    Int_t myPalette[vis_temperature_levels];
    gStyle->SetNumberContours(vis_temperature_levels);

    Int_t colorIndex = TColor::CreateGradientColorTable(numberOfSteps, stops, red, green, blue, vis_temperature_levels);
    for (int i=0;i<vis_temperature_levels;i++) myPalette[i] = colorIndex+i;
    gStyle->SetPalette(vis_temperature_levels, myPalette);
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
/*    int c = 0;
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
    geometry_created = true;*/
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


  void Vizard::writeNeighbourGraph(InactiveSurfaces& is, std::string outfile) {
    std::string filename = default_graphdir + "/";
    if (outfile.empty()) filename = filename + default_graphfile;
    else filename = filename + outfile;
    std::cout << "Preparing to write neighbour graph to " << filename << "..." << std::endl;
    std::ofstream outstream(filename.c_str());
    writeNeighbourGraph(is, outstream);
    outstream.close();    
    std::cout << "Neighbour graph written to " << filename << "." << std::endl;
  }
  /**
   * This function writes the feeder/neighbour relations in a collection of inactive surfaces to a very simple text
   * file. It essentially lists all edges of the neighbour graph, first those in the barrels, then those in the endcaps,
   * but otherwise in a more or less unordered heap: the more order in the source collection, the more order in
   * the output. If no name is given for the output file, a default filename is used.
   * @param is A reference to the collection of inactive surfaces
   * @param outfile The name of the output file that will be written to the application's default directory for graph files
   */
  void Vizard::writeNeighbourGraph(InactiveSurfaces& is, std::ostream& outstream) {
    try {
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



  void Vizard::histogramSummary(Analyzer& a, MaterialBudget& materialBudget, bool debugServices, RootWSite& site) {
    histogramSummary(a, materialBudget, debugServices, site, "outer");
  }

  /**
   * This function writes the tables of the weight summaries for the modules
   * with the rootweb library
   * @param a A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   * @param name a qualifier that goes in parenthesis in the title (outer or strip, for example)
   */ 
 
  // TODO: if weightGrid is actually unused, then remove it
  void Vizard::weigthSummart(Analyzer& a, WeightDistributionGrid& weightGrid, RootWSite& site, std::string name) {
    RootWContent* myContent;

    // Initialize the page with the material budget
    std::string pageTitle="Weights";
    if (name!="") pageTitle+=" (" +name+")";
    std::string pageAddress="weights"+name+".html";
    RootWPage& myPage = site.addPage(pageTitle);
    myPage.setAddress(pageAddress);

    // weight plot
    myContent = new RootWContent("Overview plot", true);
    myPage.addContent(myContent);

    std::map<std::string, SummaryTable>* summaryTables;

    // Write the summary for barrel first and endcap second
    for (int i=0; i<2; ++i) {
      if (i==0) summaryTables = &a.getBarrelWeightSummary();
      else summaryTables = &a.getEndcapWeightSummary();

      std::map<std::string, SummaryTable>::iterator it;
      for (it=summaryTables->begin(); it!=summaryTables->end(); ++it) {
        // Create one content per layer
        RootWContent& myContent = myPage.addContent(it->first, false);
        RootWTable& myTable = myContent.addTable();
        myTable.setContent(it->second.getContent());
      }
    }

    // Write the category summary for barrel first and endcap second
    for (int i=0; i<2; ++i) {
      if (i==0) summaryTables = &a.getBarrelWeightComponentSummary();
      else summaryTables = &a.getEndcapWeightComponentSummary();

      std::map<std::string, SummaryTable>::iterator it;
      for (it=summaryTables->begin(); it!=summaryTables->end(); ++it) {
        // Create one content per layer
        RootWContent& myContent = myPage.addContent(it->first + " - components", false);
        RootWTable& myTable = myContent.addTable();
        myTable.setContent(it->second.getContent());
      }
    }
  }



  /**
   * This function draws some of the histograms that were filled during material budget analysis
   * with the rootweb library
   * @param a A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   * @param name a qualifier that goes in parenthesis in the title (outer or strip, for example)
   */
  void Vizard::histogramSummary(Analyzer& a, MaterialBudget& materialBudget, bool debugServices, RootWSite& site, std::string name) {
    // Initialize the page with the material budget
    RootWPage* myPage;
    RootWContent* myContent;
    RootWTable* myTable;
    RootWImage* myImage;
    TCanvas* myCanvas;
    TVirtualPad* myPad;
    std::string pageTitle="Material";
    if (name!="") pageTitle+=" (" +name+")";
    myPage = new RootWPage(pageTitle);
    std::string pageAddress="material"+name+".html";
    myPage->setAddress(pageAddress);
    site.addPage(myPage);

    std::string name_overviewMaterial = std::string("overviewMaterial") + name ;
    std::string name_materialInTrackingVolume = std::string("materialInTrackingVolume") + name ;
    std::string name_detailedMaterial = std::string("detailedMaterial") + name ;
    std::string name_countourMaterial = std::string("countourMaterial") + name ;
    std::string name_mapMaterialRadiation = std::string("mapMaterialRadiation") + name ;
    std::string name_mapMaterialInteraction = std::string("mapMaterialInteraction") + name ;
    std::string name_hadronsHitsNumber = std::string("hadronsHitsNumber") + name ;
    std::string name_hadronsTracksFraction = std::string("hadronsTracksFraction") + name ;
    std::string name_hadTrackRanger = std::string("hadTrackRanger") + name ;


    // 1D Overview
    myContent = new RootWContent("1D Overview");
    myPage->addContent(myContent);

    ostringstream label;
    std::map<int, std::vector<double> > averages;

    // Book histograms
    THStack* rcontainer = new THStack("rstack", "Radiation Length by Category");
    THStack* icontainer = new THStack("istack", "Interaction Length by Category");
    TH1D *cr = NULL, *ci = NULL, *fr1 = NULL, *fi1 = NULL, *fr2 = NULL, *fi2 = NULL;
    TH1D *acr = NULL, *aci = NULL, *ser = NULL, *sei = NULL, *sur = NULL, *sui = NULL;
#ifdef MATERIAL_SHADOW
    TH2D *ir = NULL, *ii = NULL;
#endif
    TH2D *mapRad = NULL, *mapInt = NULL;
    TProfile *ciProf, *crProf;

    // Output initialisation and headers
    myCanvas = new TCanvas(name_overviewMaterial.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
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
    crProf = newProfile(cr);
    crProf->SetNameTitle("rfullvolume", "Radiation Length Over Full Tracker Volume");
    crProf->SetXTitle("#eta");
    crProf->Rebin(materialRebin);
    crProf->Draw("hist");
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
    ciProf = newProfile(ci);
    ciProf->SetNameTitle("ifullvolume", "Interaction Length Over Full Tracker Volume");
    ciProf->SetXTitle("#eta");
    ciProf->Rebin(materialRebin);
    ciProf->Draw("hist");
    // Put the total plots to the site
    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in full volume");
    myImage->setName("matFull");
    myTable = new RootWTable();
    std::ostringstream aStringStream;
    aStringStream.str("");
    aStringStream << "Average radiation length in full volume (eta = [0, ";
    aStringStream << std::dec << std::fixed
                << std::setprecision(1) << a.getEtaMaxMaterial();
    aStringStream << "])";
    myTable->setContent(1, 1, aStringStream.str().c_str());
    aStringStream.str("");
    aStringStream << "Average interaction length in full volume (eta = [0, ";
    aStringStream << std::dec << std::fixed
      << std::setprecision(1) << a.getEtaMaxMaterial();
    aStringStream << "])";
    myTable->setContent(2, 1, aStringStream.str().c_str());
    myTable->setContent(1, 2, averageHistogramValues(*cr, a.getEtaMaxMaterial()), 5);
    myTable->setContent(2, 2, averageHistogramValues(*ci, a.getEtaMaxMaterial()), 5);
    myContent->addItem(myTable);
    myContent->addItem(myImage);

    // Material summary table
    RootWTable* materialSummaryTable = new RootWTable();
      {
        double averageValue;
        materialSummaryTable->setContent(0,0,"Material");
        materialSummaryTable->setContent(1,0,"Rad. length");
        materialSummaryTable->setContent(2,0,"Int. length");
        materialSummaryTable->setContent(3,0,"Photon conversion");
        for (unsigned int j=1; j< geom_name_eta_regions.size(); ++j) {
          // Column: the cut name
          materialSummaryTable->setContent(0,j, geom_name_eta_regions[j]);

          // First row: the radiation length
          averageValue = averageHistogramValues(*cr, geom_range_eta_regions[j-1], geom_range_eta_regions[j]);
          materialSummaryTable->setContent(1,j, averageValue ,4);
          addSummaryLabelElement("radiation length ("+geom_name_eta_regions[j]+") for "+name);
          addSummaryElement(averageValue);

          // Third row: the photon conversion probability
          averageValue *= -7./9.;
          averageValue = 1 - exp(averageValue);
          materialSummaryTable->setContent(3,j, averageValue ,4);
          addSummaryLabelElement("photon conversion probability ("+geom_name_eta_regions[j]+") for "+name);
          addSummaryElement(averageValue);

          // Second row: the interaction length
          averageValue = averageHistogramValues(*ci, geom_range_eta_regions[j-1], geom_range_eta_regions[j]);
          materialSummaryTable->setContent(2,j, averageValue ,4);
          addSummaryLabelElement("interaction length ("+geom_name_eta_regions[j]+") for "+name);
          addSummaryElement(averageValue);
        }
      }

    // Detailed plots
    myContent = new RootWContent("Tracking volume", false);
    myPage->addContent(myContent);
    // Work area re-init
    myCanvas = new TCanvas(name_materialInTrackingVolume.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // global plots in tracking volume: radiation length
    cr = (TH1D*)a.getHistoGlobalR().Clone();
    cr->SetFillColor(kGray + 2);
    crProf = newProfile(cr);
    crProf->SetNameTitle("rglobal", "Overall Radiation Length");
    crProf->SetXTitle("#eta");
    crProf->Rebin(materialRebin);
    crProf->Draw("hist");
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    // global plots in tracking volume: interaction length
    ci = (TH1D*)a.getHistoGlobalI().Clone();
    ci->SetFillColor(kGray + 1);
    ciProf = newProfile(ci);
    ciProf->SetNameTitle("iglobal", "Overall Interaction Length");
    ciProf->SetXTitle("#eta");
    ciProf->Rebin(materialRebin);
    ciProf->Draw("hist");
    // Write global tracking volume plots to web pag
    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in tracking volume");
    myImage->setName("matTrack");
    myTable = new RootWTable();
    char titleString[256];
    sprintf(titleString, "Average radiation length in tracking volume (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(1, 1, titleString);
    sprintf(titleString, "Average interaction length in tracking volume (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(2, 1, titleString);
    myTable->setContent(1, 2, averageHistogramValues(*cr, a.getEtaMaxMaterial()), 5);
    myTable->setContent(2, 2, averageHistogramValues(*ci, a.getEtaMaxMaterial()), 5);
    myContent->addItem(myTable);
    myContent->addItem(myImage);

    // Detailed plots
    myContent = new RootWContent("Detailed", false);
    myPage->addContent(myContent);
    // Work area re-init
    myCanvas = new TCanvas(name_detailedMaterial.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // radiation length in tracking volume by active, serving or passive
    sur = (TH1D*)a.getHistoSupportsAllR().Clone();
    sur->SetFillColor(kOrange + 4);
    sur->SetXTitle("#eta");
    rcontainer->Add(sur);
    ser = (TH1D*)a.getHistoServicesAllR().Clone();
    ser->SetFillColor(kBlue);
    ser->SetXTitle("#eta");
    rcontainer->Add(ser);
    acr = (TH1D*)a.getHistoModulesAllR().Clone();
    acr->SetFillColor(kRed);
    acr->SetXTitle("#eta");
    rcontainer->Add(acr);
    rcontainer->Draw();
    // interaction length in tracking volume by active, serving or passive
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    sui = (TH1D*)a.getHistoSupportsAllI().Clone();
    sui->SetFillColor(kOrange + 2);
    sui->SetXTitle("#eta");
    icontainer->Add(sui);
    sei = (TH1D*)a.getHistoServicesAllI().Clone();
    sei->SetFillColor(kAzure - 2);
    sei->SetXTitle("#eta");
    icontainer->Add(sei);
    aci = (TH1D*)a.getHistoModulesAllI().Clone();
    aci->SetFillColor(kRed - 3);
    aci->SetXTitle("#eta");
    icontainer->Add(aci);
    icontainer->Draw();

    // Write asl category plots to web page
    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Detailed");
    myImage->setName("matTrackDet");
    myTable = new RootWTable();
    // Average values by active, service and passive
    sprintf(titleString, "Average (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(1, 0, "modules");
    myTable->setContent(2, 0, "services");
    myTable->setContent(3, 0, "supports");
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");
    myTable->setContent(1, 1, averageHistogramValues(*acr, a.getEtaMaxMaterial()), 5);
    myTable->setContent(2, 1, averageHistogramValues(*ser, a.getEtaMaxMaterial()), 5);
    myTable->setContent(3, 1, averageHistogramValues(*sur, a.getEtaMaxMaterial()), 5);
    myTable->setContent(1, 2, averageHistogramValues(*aci, a.getEtaMaxMaterial()), 5);
    myTable->setContent(2, 2, averageHistogramValues(*sei, a.getEtaMaxMaterial()), 5);
    myTable->setContent(3, 2, averageHistogramValues(*sui, a.getEtaMaxMaterial()), 5);
    myContent->addItem(myTable);
    myContent->addItem(myImage);





    myContent = new RootWContent("Module components detail", false);
    myPage->addContent(myContent);

    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");


    THStack* rCompStack = new THStack("rcompstack", "Radiation Length by Component");
    THStack* iCompStack = new THStack("icompstack", "Interaction Length by Component");

    TLegend* compLegend = new TLegend(0.1,0.6,0.35,0.9);

    myCanvas = new TCanvas(("moduleComponentsRI"+name).c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    std::map<std::string, TH1D*>& rActiveComps = a.getHistoActiveComponentsR();
    int compIndex = 1;
    for (std::map<std::string, TH1D*>::iterator it = rActiveComps.begin(); it != rActiveComps.end(); ++it) {
      it->second->SetLineColor(Palette::color(compIndex));
      it->second->SetFillColor(Palette::color(compIndex));
      it->second->SetXTitle("#eta");
      it->second->SetTitle(it->first.c_str());
      compLegend->AddEntry(it->second, it->first.c_str());
      rCompStack->Add(it->second);
      myTable->setContent(compIndex, 0, it->first);
      myTable->setContent(compIndex++, 1, averageHistogramValues(*it->second, a.getEtaMaxMaterial()), 5);
    }
    rCompStack->Draw();
    compLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    std::map<std::string, TH1D*>& iActiveComps = a.getHistoActiveComponentsI();
    compIndex = 1;
    for (std::map<std::string, TH1D*>::iterator it = iActiveComps.begin(); it != iActiveComps.end(); ++it) {
      it->second->SetLineColor(Palette::color(compIndex));
      it->second->SetFillColor(Palette::color(compIndex));
      it->second->SetXTitle("#eta");
      it->second->SetTitle(it->first.c_str());
      iCompStack->Add(it->second);
      myTable->setContent(compIndex++, 2, averageHistogramValues(*it->second, a.getEtaMaxMaterial()), 5);
    }
    iCompStack->Draw();
    compLegend->Draw();

    myContent->addItem(myTable);

    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by component type in modules");
    myImage->setName("riDistrComp");
    myContent->addItem(myImage);



    // Work area re-init
    myCanvas = new TCanvas(name_countourMaterial.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // Countour plots
    myContent = new RootWContent("Material distribution", true);
    myPage->addContent(myContent);

#ifdef MATERIAL_SHADOW        
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
    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material 2D distributions");
    myImage->setName("matShadow");
    myContent->addItem(myImage);
#endif // MATERIAL_SHADOW

    // Radiation length plot
    myCanvas = new TCanvas(name_mapMaterialRadiation.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->cd();
    mapRad = (TH2D*)a.getHistoMapRadiation().Clone();
    mapRad->SetContour(vis_temperature_levels, NULL);
    //myCanvas->SetLogz();
    mapRad->Draw("COLZ");
    myImage = new RootWImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation length material map");
    myImage->setName("matMapR");
    myContent->addItem(myImage);

    // Interaction length plot
    myCanvas = new TCanvas(name_mapMaterialInteraction.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->cd();
    mapInt = (TH2D*)a.getHistoMapInteraction().Clone();
    mapInt->SetContour(vis_temperature_levels, NULL);
    mapInt->Draw("COLZ");
    myImage = new RootWImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Interaction length material map");
    myImage->setName("matMapI");
    myContent->addItem(myImage);

    // Nuclear interactions
    myContent = new RootWContent("Nuclear interactions", true);
    myPage->addContent(myContent);

    // Number of hits
    myCanvas = new TCanvas(name_hadronsHitsNumber.c_str());
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    TGraph* hadronTotalHitsGraph = new TGraph(a.getHadronTotalHitsGraph());
    TGraph* hadronAverageHitsGraph = new TGraph(a.getHadronAverageHitsGraph());
    hadronTotalHitsGraph->SetMarkerStyle(8);
    hadronTotalHitsGraph->SetMarkerColor(kBlack);
    hadronTotalHitsGraph->SetMinimum(0);
    hadronTotalHitsGraph->Draw("alp");
    hadronAverageHitsGraph->SetMarkerStyle(8);
    hadronAverageHitsGraph->SetMarkerColor(kRed);
    hadronAverageHitsGraph->Draw("same lp");

    // Number of hits
    std::vector<TGraph> hadronGoodTracksFraction=a.getHadronGoodTracksFraction();
    std::vector<double> hadronNeededHitsFraction=a.getHadronNeededHitsFraction();
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    TLegend* myLegend = new TLegend(0.75, 0.16, .95, .40);
    // Old-style palette by Stefano, with custom-generated colors
    // Palette::prepare(hadronGoodTracksFraction.size()); // there was a 120 degree phase here
    // Replaced by the libreOffice-like palette
    TH1D* ranger = new TH1D(name_hadTrackRanger.c_str(),"", 100, 0, a.getEtaMaxMaterial());
    ranger->SetMaximum(1.);
    TAxis* myAxis;
    myAxis = ranger->GetXaxis();
    myAxis->SetTitle("#eta");
    myAxis = ranger->GetYaxis();
    myAxis->SetTitle("Tracks fraction");
    ranger->Draw();
    ostringstream tempSS;
    std::map<int, std::string> fractionTitles;
    for (unsigned int i=0;
         i<hadronGoodTracksFraction.size();
         ++i) {
      TGraph& myGraph = hadronGoodTracksFraction.at(i);
      //std::cerr << "Good Hadrons fractions at (" << i <<") has " << myGraph.GetN() << " points" << std::endl;
      //double xx, yy;
      //myGraph.GetPoint(myGraph.GetN()-1, xx, yy);
      //std::cerr << "Last point (x,y) = ("<< xx <<", " << yy <<")" << std::endl;
      averages[i] = Analyzer::average(myGraph, geom_range_eta_regions);
      closeGraph(myGraph);
      myGraph.SetFillColor(Palette::color(i+1));
      myGraph.Draw("same F");
      tempSS.str("");
      if (hadronNeededHitsFraction.at(i)!=Analyzer::ZeroHitsRequired) {
        if (hadronNeededHitsFraction.at(i)==Analyzer::OneHitRequired)
          tempSS << "1 hit required";
        else
          tempSS << int(hadronNeededHitsFraction.at(i)*100)
            << "%% hits required";
        fractionTitles[i]=tempSS.str();
        myLegend->AddEntry(&myGraph, fractionTitles[i].c_str(), "F");
      }
    }
    ranger->Draw("sameaxis");
    myLegend->Draw();
    myImage = new RootWImage(myCanvas, 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hits occupancy & track efficiency for hadrons");
    myImage->setName("hadHitsTracks");
    myContent->addItem(myImage);

    //if (name=="outer") {    
    // Summary table
    RootWContent& summaryContent = myPage->addContent("Summary", false);
    RootWTable& cutsTable = summaryContent.addTable();
    RootWTable& summaryTable = summaryContent.addTable();
    summaryTable.setContent(0,0,"Clean pions");

    // Table explaining the cuts
    cutsTable.setContent(0,0,"Region");
    cutsTable.setContent(1,0,"etaMin");
    cutsTable.setContent(2,0,"etaMax");

    myTable = &cutsTable;
    for (unsigned int iBorder=0; iBorder<geom_name_eta_regions.size()-1; ++iBorder) {
      myTable->setContent(0,iBorder+1,geom_name_eta_regions[iBorder+1]);
      label.str(""); label << std::fixed << std::setprecision(1) << geom_range_eta_regions[iBorder];
      myTable->setContent(1,iBorder+1,label.str());
      label.str(""); label << std::fixed << std::setprecision(1) << geom_range_eta_regions[iBorder+1];
      myTable->setContent(2,iBorder+1,label.str());
    }

    int delta=1;
    for (unsigned int i=0;
         i<hadronGoodTracksFraction.size();
         ++i) {
      if (fractionTitles[i]!="") {
        summaryTable.setContent(i+delta,0,fractionTitles[i]);
        int j=1;
        double myValue;
        for ( std::vector<double>::iterator it = averages[i].begin();
              it != averages[i].end(); ++it ) {
          //std::cout << "average: " << (*it) << std::endl;
          myValue = 100 * (1 - (*it));
          summaryTable.setContent(i+delta,j, myValue,4);
          summaryTable.setContent(0,j, geom_name_eta_regions[j]);
          addSummaryLabelElement(fractionTitles[i]+"("+geom_name_eta_regions[j]+") for "+name);
          addSummaryElement(myValue);
          j++;
        }
      } else delta--;
    }
    summaryContent.addItem(materialSummaryTable);
    //} else {
    //  delete materialSummaryTable;
    //}

    if (debugServices) {
        drawInactiveSurfacesSummary(materialBudget, *myPage);
    }

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
/*    Layer* current;
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
    else std::cout << "detailedModules(): layers vector is empty." << std::endl;*/
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
/*    b = m->getCorner(1) - m->getCorner(0);
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
    tr = new TGeoCombiTrans(m->getCorner(0).X(), m->getCorner(0).Y(), m->getCorner(0).Z(), rot); */
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
   // if (cobin >= histo.GetNbinsX() - 1) avg = histo.GetMean();
   // else {
      for (int i = 1; i <= cobin; i++) avg = avg + histo.GetBinContent(i) / (double)cobin;
   // }
    return avg;
  }

  /**
   * This function computes the average of a range of histogram bins: from the first to the one that includes a
   * cutoff value along the axis.
   * @param histo A reference to the histogram data
   * @param cutoffStart The lower cutoff value
   * @param cutoffEnd The upper cutoff value
   * @return The average value of the bins within range
   */
  double Vizard::averageHistogramValues(TH1D& histo, double cutoffStart, double cutoffEnd) {
    double avg = 0.0;
    int coBinStart = 1;
    int coBinEnd = 1;
    if (cutoffStart >= cutoffEnd) return 0;
    // find first relevant bin
    while ((coBinStart < histo.GetNbinsX()) && (histo.GetBinLowEdge(coBinStart) < cutoffStart)) coBinStart++;
    coBinEnd=coBinStart;
    // find last relevant bin
    while ((coBinEnd < histo.GetNbinsX()) && (histo.GetBinLowEdge(coBinEnd) < cutoffEnd)) coBinEnd++;
    double coBinN=coBinEnd-coBinStart+1; // TODO: IMPORTANT check this
    // calculate average
    if (coBinStart> histo.GetNbinsX() - 1) coBinStart= histo.GetNbinsX() - 1;
    if (coBinEnd > histo.GetNbinsX() - 1) coBinEnd= histo.GetNbinsX() - 1;
    for (int i = coBinStart; i <= coBinEnd; i++) avg = avg + histo.GetBinContent(i) / coBinN;
    return avg;
  }


  std::string Vizard::getSummaryString() {
    return summaryCsv_;
  }

  std::string Vizard::getSummaryLabelString() {
    return summaryCsvLabels_;
  }

  void Vizard::setSummaryString(std::string myString) {
    summaryCsv_ = myString;
  }

  void Vizard::setSummaryLabelString(std::string myString) {
    summaryCsvLabels_ = myString;
  }

  void Vizard::addSummaryElement(std::string element, bool first /*= false*/ ) {
    if (!first) summaryCsv_+=csv_separator;
    summaryCsv_+=element;
  }

  void Vizard::addSummaryLabelElement(std::string element, bool first /*= false*/ ) {
    if (!first) summaryCsvLabels_+=csv_separator;
    summaryCsvLabels_+=element;
  }

  void Vizard::addSummaryElement(double element, bool first /*= false*/ ) {
    std::ostringstream myElement;
    std::string myString;
    myElement.str("");
    myElement << element;
    addSummaryElement(myElement.str(), first);
  }

  void Vizard::addOccupancyElement(std::string element) {
    // 2 is a magic adjustment factoroccupancyCsv_ += element;
    occupancyCsv_ += csv_separator;
  }

  void Vizard::addOccupancyElement(double element) {
    std::ostringstream myElement;
    std::string myString;
    myElement.str("");
    myElement << element;
    addOccupancyElement(myElement.str());
  }

  void Vizard::addOccupancyEOL() {
    occupancyCsv_ += "\n";
  }

  /**
   * This function draws the profile of hits obtained by the analysis of the geometry
   * together with the summaries in tables with the rootweb library. It also actually does a couple of
   * calculations to count modules and such, to put the results in the tables.
   * @param analyzer A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   */
  bool Vizard::geometrySummary(Analyzer& analyzer, Tracker& tracker, SimParms& simparms, InactiveSurfaces* inactive, RootWSite& site, bool& debugResolution, std::string name) {
    trackers_.push_back(&tracker);

    std::map<std::string, double>& tagMapWeight = analyzer.getTagWeigth();
    

    std::string pageTitle = "Geometry";
    if (name!="") pageTitle+=" (" +name+")";
    
    RootWPage* myPage = new RootWPage(pageTitle);
    // TODO: the web site should decide which page to call index.html

    std::string pageAddress="index"+name+".html";
    myPage->setAddress(pageAddress);

    site.addPage(myPage, 100);
    RootWContent* myContent;


    // Inactive surfaces
    double inactiveSurfacesTotalMass;
    if (inactive) {
      std::vector<InactiveElement>& inactiveBarrelServices = inactive->getBarrelServices();
      std::vector<InactiveElement>& inactiveEndcapServices = inactive->getEndcapServices();
      std::vector<InactiveElement>& inactiveSupports = inactive->getSupports();
      std::vector<InactiveElement> allInactives;
      allInactives.reserve( inactiveBarrelServices.size() + inactiveEndcapServices.size() + inactiveSupports.size() );
      allInactives.insert(allInactives.end(), inactiveBarrelServices.begin(), inactiveBarrelServices.end() );
      allInactives.insert(allInactives.end(), inactiveEndcapServices.begin(), inactiveEndcapServices.end() );
      allInactives.insert(allInactives.end(), inactiveSupports.begin(), inactiveSupports.end() );
      inactiveSurfacesTotalMass = 0;
      for (const auto elem : allInactives ) {
        if (elem.getTotalMass()>0) inactiveSurfacesTotalMass += elem.getTotalMass();
      }
    }

    //********************************//
    //*                              *//
    //*       Layers and disks       *//
    //*                              *//
    //********************************//
    myContent = new RootWContent("Layers and disks");
    myPage->addContent(myContent);


    // Build the module type maps
    // with a pointer to a sample module
    // Build the layer summary BTW

    class LayerDiskSummaryVisitor : public ConstGeometryVisitor {
    public:
      RootWTable* layerTable = new RootWTable();
      RootWTable* diskTable = new RootWTable();
      RootWTable* ringTable = new RootWTable();
      std::map<std::string, std::set<std::string> > tagMapPositions;
      std::map<std::string, int> tagMapCount;
      std::map<std::string, long> tagMapCountChan;
      std::map<std::string, double> tagMapMaxStripOccupancy;
      std::map<std::string, double> tagMapAveStripOccupancy;
      std::map<std::string, double> tagMapMaxHitOccupancy;
      std::map<std::string, double> tagMapAveHitOccupancy;
      std::map<std::string, double> tagMapAveRphiResolution;
      std::map<std::string, double> tagMapAveRphiResolutionRmse;
      std::map<std::string, double> tagMapSumXResolution;
      std::map<std::string, double> tagMapSumSquaresXResolution;
      std::map<std::string, double> tagMapCountXResolution;
      std::map<std::string, double> tagMapIsParametrizedXResolution;
      std::map<std::string, double> tagMapAveYResolution;
      std::map<std::string, double> tagMapAveYResolutionRmse;
      std::map<std::string, double> tagMapSumYResolution;
      std::map<std::string, double> tagMapSumSquaresYResolution;
      std::map<std::string, double> tagMapCountYResolution;
      std::map<std::string, double> tagMapIsParametrizedYResolution;
      std::map<std::string, double> tagMapAveRphiResolutionTrigger;
      std::map<std::string, double> tagMapAveYResolutionTrigger;
      std::map<std::string, double> tagMapSensorPowerAvg;
      std::map<std::string, double> tagMapSensorPowerMax;
      std::map<std::string, const DetectorModule*> tagMap;
      std::map<int, const EndcapModule*> ringTypeMap;

      int nBarrelLayers=0;
      int nDisks=0;
      int totalBarrelModules = 0;
      int totalEndcapModules = 0;

      double totArea = 0;
      int totCountMod = 0;
      int totCountSens = 0;
      long totChannel = 0;
      double totalSensorPower = 0;

      double nMB;

      void preVisit() {
        layerTable->setContent(0, 0, "Layer");
        layerTable->setContent(1, 0, "r");
        layerTable->setContent(2, 0, "z_max");
        layerTable->setContent(3, 0, "# mod");
        layerTable->setContent(4, 0, "# rods");
        diskTable->setContent(0, 0, "Disk");
        diskTable->setContent(1, 0, "z");
        diskTable->setContent(2, 0, "# mod");
        ringTable->setContent(0, 0, "Ring");
        ringTable->setContent(1, 0, "r"+subStart+"min"+subEnd);
        ringTable->setContent(2, 0, "r"+subStart+"low"+subEnd);
        ringTable->setContent(3, 0, "r"+subStart+"high"+subEnd);
        ringTable->setContent(4, 0, "r"+subStart+"max"+subEnd);
      }

      void visit(const SimParms& s) override { nMB = s.numMinBiasEvents(); }

      void visit(const Layer& l) override {
        if (l.maxZ() < 0.) return;
        ++nBarrelLayers;
        int nModules = l.totalModules();
        totalBarrelModules += nModules;
        layerTable->setContent(0, nBarrelLayers, l.myid());
        layerTable->setContent(1, nBarrelLayers, l.placeRadius(), coordPrecision);
        layerTable->setContent(2, nBarrelLayers, l.maxZ(), coordPrecision);
        layerTable->setContent(3, nBarrelLayers, nModules);
        layerTable->setContent(4, nBarrelLayers, l.numRods());
      }

      void visit(const Disk& d) override {
        if (d.averageZ() < 0.) return;
        ++nDisks;
        int nModules = d.totalModules();
        totalEndcapModules += nModules;
        diskTable->setContent(0, nDisks, d.myid());
        diskTable->setContent(1, nDisks, d.averageZ(), coordPrecision);
        diskTable->setContent(2, nDisks, nModules);
      }

      void visit(const Module& m) override {
        TagMaker tmak(m);

        std::string aSensorTag = tmak.sensorGeoTag;
        tagMapPositions[aSensorTag].insert(tmak.posTag);
        tagMapCount[aSensorTag]++;
        tagMapCountChan[aSensorTag] += m.totalChannels();
        tagMapMaxStripOccupancy[aSensorTag] = MAX(m.stripOccupancyPerEvent()*nMB, tagMapMaxStripOccupancy[aSensorTag]);
        tagMapMaxHitOccupancy[aSensorTag] = MAX(m.hitOccupancyPerEvent()*nMB, tagMapMaxHitOccupancy[aSensorTag]);
        tagMapAveStripOccupancy[aSensorTag] += m.stripOccupancyPerEvent()*nMB;
        tagMapAveHitOccupancy[aSensorTag] += m.hitOccupancyPerEvent()*nMB;
	// modules' spatial resolution along the local X axis is not parametrized
	if (!m.hasAnyResolutionLocalXParam()) {
	  tagMapIsParametrizedXResolution[aSensorTag] = false;
	  tagMapAveRphiResolution[aSensorTag] += m.nominalResolutionLocalX(); 
	  tagMapAveRphiResolutionRmse[aSensorTag] += 0.; 
	}
	// modules' spatial resolution along the local X axis is parametrized
	else {
	  tagMapIsParametrizedXResolution[aSensorTag] = true;
	  if (boost::accumulators::count(m.rollingParametrizedResolutionLocalX) > 0) { // calculation only on hit modules
	    tagMapSumXResolution[aSensorTag] += sum(m.rollingParametrizedResolutionLocalX);
	    tagMapSumSquaresXResolution[aSensorTag] += moment<2>(m.rollingParametrizedResolutionLocalX) * boost::accumulators::count(m.rollingParametrizedResolutionLocalX);	  
	    tagMapCountXResolution[aSensorTag] += double(boost::accumulators::count(m.rollingParametrizedResolutionLocalX));
	  }
	}
	// modules' spatial resolution along the local Y axis is not parametrized
        if (!m.hasAnyResolutionLocalYParam()) { 
	  tagMapIsParametrizedYResolution[aSensorTag] = false;
	  tagMapAveYResolution[aSensorTag] += m.nominalResolutionLocalY(); 
	  tagMapAveYResolutionRmse[aSensorTag] += 0.; 
	}
	// modules' spatial resolution along the local Y axis is parametrized
	else {
	  tagMapIsParametrizedYResolution[aSensorTag] = true;
	  if (boost::accumulators::count(m.rollingParametrizedResolutionLocalY) > 0) { // calculation only on hit modules
	    tagMapSumYResolution[aSensorTag] += sum(m.rollingParametrizedResolutionLocalY);
	    tagMapSumSquaresYResolution[aSensorTag] += moment<2>(m.rollingParametrizedResolutionLocalY) * boost::accumulators::count(m.rollingParametrizedResolutionLocalY);	  													     
	    tagMapCountYResolution[aSensorTag] += double(boost::accumulators::count(m.rollingParametrizedResolutionLocalY));
	  }
	}
        //tagMapAveRphiResolutionTrigger[aSensorTag] += m.resolutionRPhiTrigger();
        //tagMapAveYResolutionTrigger[aSensorTag] += m.resolutionYTrigger();
        tagMapSensorPowerAvg[aSensorTag] += m.irradiationPower();
        if (tagMapSensorPowerMax[aSensorTag] < m.irradiationPower()) tagMapSensorPowerMax[aSensorTag] = m.irradiationPower();
        totCountMod++;
        totCountSens += m.numSensors();
        totChannel += m.totalChannels();
        totArea += m.area()*m.numSensors();
        totalSensorPower += m.irradiationPower();
        if (tagMap.find(aSensorTag)==tagMap.end()){
          // We have a new sensor geometry
          tagMap[aSensorTag] = &m;
        }
      }

      void visit(const EndcapModule& m) override {
        if (m.disk() != 1 && m.side() != 1) return;
        if (ringTypeMap.find(m.ring())==ringTypeMap.end()){
          // We have a new sensor geometry
          ringTypeMap[m.ring()] = &m;
        }

      }

      void postVisit() {
        layerTable->setContent(0, nBarrelLayers+1, "Total");
        layerTable->setContent(3, nBarrelLayers+1, totalBarrelModules);
        diskTable->setContent(0, nDisks+1, "Total");
        diskTable->setContent(2, nDisks+1, totalEndcapModules*2);

        std::ostringstream myName;
        for (auto typeIt = ringTypeMap.begin();
             typeIt!=ringTypeMap.end(); typeIt++) {
          auto* anEC = (*typeIt).second;
          int aRing=(*typeIt).first;
          ringTable->setContent(0, aRing, aRing);
          ringTable->setContent(1, aRing, anEC->minR(), coordPrecision);
          ringTable->setContent(2, aRing, sqrt(pow(anEC->minR(),2)+pow(anEC->minWidth()/2.,2)), coordPrecision); // Ugly, this should be accessible as a method
          ringTable->setContent(3, aRing, anEC->minR()+anEC->length(), coordPrecision);
          ringTable->setContent(4, aRing, anEC->maxR(), coordPrecision);
        }
      }
    };

    LayerDiskSummaryVisitor v;
    v.preVisit();
    simparms.accept(v);
    tracker.accept(v);
    v.postVisit();

    myContent->addItem(v.layerTable);
    myContent->addItem(v.diskTable);
    myContent->addItem(v.ringTable);





    double totalPower=0; 
    double totalCost=0;
    double moduleTotalWeight=0;

    std::ostringstream aName;
    std::ostringstream aTag;
    std::ostringstream aType;
    std::ostringstream aThickness;
    std::ostringstream aModuleArea;
    std::ostringstream aTotalArea;
    std::ostringstream aStripOccupancy;
    std::ostringstream aHitOccupancy;
    std::ostringstream anRphiResolution;
    std::ostringstream anRphiResolutionRmse;
    std::ostringstream aYResolution;
    std::ostringstream aYResolutionRmse;
    std::ostringstream anRphiResolutionTrigger;
    std::ostringstream aYResolutionTrigger;
    std::ostringstream aPitchPair;
    std::ostringstream aStripLength;
    std::ostringstream aSegment;
    std::ostringstream anNstrips;
    std::ostringstream aNumberMod;
    std::ostringstream aNumberSens;
    std::ostringstream aChannel;
    std::ostringstream aPower;
    std::ostringstream aPowerPerModule;
    std::ostringstream aSensorPower;
    std::ostringstream aSensorPowerPerModuleAvg;
    std::ostringstream aSensorPowerPerModuleMax;
    std::ostringstream aCost;
    std::ostringstream aWeight;
    int barrelCount=0;
    int endcapCount=0;

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
    static const int thicknessRow = 3;
    static const int moduleAreaRow = 4;
    static const int totalAreaRow = 5;
    static const int numbermodsRow = 6;
    static const int numbersensRow = 7;
    static const int channelRow = 8;
    static const int nstripsRow = 9;
    static const int segmentsRow = 10;
    static const int striplengthRow = 11;
    static const int pitchpairsRow = 12;
    static const int rphiResolutionRow = 13;
    static const int rphiResolutionRmseRow = 14;
    static const int yResolutionRow = 15;
    static const int yResolutionRmseRow = 16;
    static const int rphiResolutionTriggerRow = 17;
    static const int yResolutionTriggerRow = 18;
    static const int stripOccupancyRow = 19;
    static const int hitOccupancyRow = 20;
    static const int powerPerModuleRow = 21;
    static const int sensorPowerPerModuleAvgRow = 22;
    static const int sensorPowerPerModuleMaxRow = 23;
    static const int powerRow = 24;
    static const int sensorPowerRow = 25;
    static const int costRow = 26;
    static const int moduleWeightRow = 27;
    static const int inactiveWeightRow = 28;
    static const int totalWeightRow = 29;

    // Row names
    moduleTable->setContent(tagRow, 0, "Tag");
    moduleTable->setContent(typeRow, 0, "Type");
    moduleTable->setContent(thicknessRow, 0, "Sensor spacing");
    moduleTable->setContent(moduleAreaRow, 0, "Sensor area (mm"+superStart+"2"+superEnd+")");
    moduleTable->setContent(totalAreaRow, 0, "Total area (m"+superStart+"2"+superEnd+")");
    moduleTable->setContent(stripOccupancyRow, 0, "Strip Occ (max/av)");
    moduleTable->setContent(hitOccupancyRow, 0, "Hit Occ (max/av)");
    moduleTable->setContent(rphiResolutionRow, 0, "Local X axis resolution : " + muLetter + " ("+muLetter+"m)");
    moduleTable->setContent(rphiResolutionRmseRow, 0, "Local X axis resolution : " + sigmaLetter + " ("+muLetter+"m)");
    moduleTable->setContent(yResolutionRow, 0, "Local Y axis resolution : " + muLetter + " ("+muLetter+"m)");
    moduleTable->setContent(yResolutionRmseRow, 0, "Local Y axis resolution : " + sigmaLetter + " ("+muLetter+"m)");
    moduleTable->setContent(rphiResolutionTriggerRow, 0, "R/Phi resolution [pt] ("+muLetter+"m)");
    moduleTable->setContent(yResolutionTriggerRow, 0, "Y resolution [pt] ("+muLetter+"m)");
    moduleTable->setContent(pitchpairsRow, 0, "Pitch (min/max) ("+muLetter+"m)");
    moduleTable->setContent(striplengthRow, 0, "Strip length (mm)");
    moduleTable->setContent(segmentsRow, 0, "Segments x Chips");
    moduleTable->setContent(nstripsRow, 0, "Chan/Sensor");
    moduleTable->setContent(numbermodsRow, 0, "N. mod");
    moduleTable->setContent(numbersensRow, 0, "N. sens");
    moduleTable->setContent(channelRow, 0, "Channels (M)");
    moduleTable->setContent(powerRow, 0, "FE Power (kW)");
    moduleTable->setContent(sensorPowerRow, 0, "Sensor power (kW)");
    moduleTable->setContent(powerPerModuleRow, 0, "FE Power/mod (mW)");
    moduleTable->setContent(sensorPowerPerModuleAvgRow, 0, "Agerage sensor power/mod (mW)");
    moduleTable->setContent(sensorPowerPerModuleMaxRow, 0, "Max sensor power/mod (mW)");
    moduleTable->setContent(costRow, 0, "Cost (MCHF)");
    moduleTable->setContent(moduleWeightRow, 0, "Weight (av, g)");
    moduleTable->setContent(inactiveWeightRow, 0, "Service Weight");
    moduleTable->setContent(totalWeightRow, 0, "Total Weight");

    int loPitch;
    int hiPitch;


    setOccupancyString("");

    addOccupancyElement(tracker.myid());
    addOccupancyElement("");
    addOccupancyElement("");
    addOccupancyEOL();
    addOccupancyElement("radius");
    addOccupancyElement("occupancy [%]");
    addOccupancyElement("rphi resolution [um]");

    // Summary cycle: prepares the rows cell by cell
    int iType=0;

    for (auto tagMapIt=v.tagMap.begin(); tagMapIt!=v.tagMap.end(); tagMapIt++) {
      ++iType;
      // Name
      aName.str("");
      auto aModule=(*tagMapIt).second;
      if (dynamic_cast<const BarrelModule*>(aModule)) {
        aName << std::dec << "B" << subStart << ++barrelCount << subEnd;
      }
      if (dynamic_cast<const EndcapModule*>(aModule)) {
        aName << std::dec << "E" << subStart << ++endcapCount << subEnd;
      }
      // Tag
      aTag.str("");
      //aTag << smallStart << aModule->getTag() << smallEnd;
      aTag << smallStart;
      for (std::set<std::string>::iterator strIt = v.tagMapPositions[(*tagMapIt).first].begin();
           strIt!=v.tagMapPositions[(*tagMapIt).first].end(); ++strIt) 
        aTag << (*strIt) << "<br/> ";
      aTag << smallEnd;
      // Type
      aType.str("");
      aType << (*tagMapIt).second->moduleType();
      // Thickness
      aThickness.str("");
      aThickness << (*tagMapIt).second->dsDistance();
      // Area
      aModuleArea.str("");
      aModuleArea << std::dec << std::fixed << std::setprecision(areaPrecision) << (*tagMapIt).second->area();
      aTotalArea.str("");
      aTotalArea << std::dec << std::fixed << std::setprecision(areaPrecision) << (*tagMapIt).second->area() *
        (*tagMapIt).second->numSensors() * v.tagMapCount[(*tagMapIt).first] * 1e-6;
      // if ((*tagMapIt).second->getArea()<0) { anArea << "XXX"; } // TODO: what's this?
      // Occupancy
      aStripOccupancy.str("");
      aHitOccupancy.str("");
      aStripOccupancy << std::dec << std::fixed << std::setprecision(occupancyPrecision) <<  v.tagMapMaxStripOccupancy[(*tagMapIt).first]*100<< "/" <<v.tagMapAveStripOccupancy[(*tagMapIt).first]*100/v.tagMapCount[(*tagMapIt).first] ; // Percentage
      aHitOccupancy << std::dec << std::fixed << std::setprecision(occupancyPrecision) <<  v.tagMapMaxHitOccupancy[(*tagMapIt).first]*100<< "/" <<v.tagMapAveHitOccupancy[(*tagMapIt).first]*100/v.tagMapCount[(*tagMapIt).first] ; // Percentage

      addOccupancyEOL();
      addOccupancyElement((aModule->minR() + aModule->maxR())/2);
      addOccupancyElement(v.tagMapAveStripOccupancy[(*tagMapIt).first]*100/v.tagMapCount[(*tagMapIt).first]);
      addOccupancyElement(v.tagMapAveHitOccupancy[(*tagMapIt).first]*100/v.tagMapCount[(*tagMapIt).first]);

      // Formulae used for mean and RMSE in the parametric case :
      // mean = S / N
      // rmse = sqrt( (N*Q - S^2) / N^2 )
      // with the following notation :
      // N : count
      // S : sum
      // Q : sum of squares

      // RphiResolution
      anRphiResolution.str("");
      // modules' spatial resolution along the local X axis is not parametrized
      if (!v.tagMapIsParametrizedXResolution[(*tagMapIt).first]) {
	anRphiResolution << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << v.tagMapAveRphiResolution[(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um
      }
      // modules' spatial resolution along the local X axis is parametrized
      else {
	if (v.tagMapCountXResolution[(*tagMapIt).first] != 0) {
	  anRphiResolution << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << (v.tagMapSumXResolution[(*tagMapIt).first]) / (v.tagMapCountXResolution[(*tagMapIt).first]) / Units::um; // mm -> um
	}
	else {
	  anRphiResolution << "n/a"; // resolution is parametrized but no run with -r or -R
	}
      }

      // RphiResolution Rmse
      anRphiResolutionRmse.str("");
      // modules' spatial resolution along the local X axis is not parametrized
      if (!v.tagMapIsParametrizedXResolution[(*tagMapIt).first]) {
	anRphiResolutionRmse << std::dec << std::fixed << std::setprecision(rphiResolutionRmsePrecision) << v.tagMapAveRphiResolutionRmse[(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um
      }
      // modules' spatial resolution along the local X axis is parametrized
      else {
	if (v.tagMapCountXResolution[(*tagMapIt).first] != 0) {
	  anRphiResolutionRmse << std::dec << std::fixed << std::setprecision(rphiResolutionRmsePrecision) << sqrt(fabs((v.tagMapCountXResolution[(*tagMapIt).first] * v.tagMapSumSquaresXResolution[(*tagMapIt).first] - pow(v.tagMapSumXResolution[(*tagMapIt).first], 2)) / pow(v.tagMapCountXResolution[(*tagMapIt).first], 2))) / Units::um; // mm -> um
	}
	else {
	  anRphiResolutionRmse << "n/a"; // resolution is parametrized but no run with -r or -R
	}
      }

      // YResolution
      aYResolution.str("");
      // modules' spatial resolution along the local Y axis is not parametrized
      if (!v.tagMapIsParametrizedYResolution[(*tagMapIt).first]) {
	aYResolution << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << v.tagMapAveYResolution[(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um
      }
      // modules' spatial resolution along the local Y axis is parametrized
      else {
	if (v.tagMapCountYResolution[(*tagMapIt).first] != 0) {
	  aYResolution << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << (v.tagMapSumYResolution[(*tagMapIt).first]) / (v.tagMapCountYResolution[(*tagMapIt).first]) / Units::um; // mm -> um
	}
	else {
	  aYResolution << "n/a"; // resolution is parametrized but no run with -r or -R
	}
      }

      // YResolution Rmse
      aYResolutionRmse.str("");
      // modules' spatial resolution along the local Y axis is not parametrized
      if (!v.tagMapIsParametrizedYResolution[(*tagMapIt).first]) {
	aYResolutionRmse << std::dec << std::fixed << std::setprecision(rphiResolutionRmsePrecision) << v.tagMapAveYResolutionRmse[(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um
      }
      // modules' spatial resolution along the local Y axis is parametrized
      else {
	if (v.tagMapCountYResolution[(*tagMapIt).first] != 0) {
	  aYResolutionRmse << std::dec << std::fixed << std::setprecision(rphiResolutionRmsePrecision) << sqrt(fabs((v.tagMapCountYResolution[(*tagMapIt).first] * v.tagMapSumSquaresYResolution[(*tagMapIt).first] - pow(v.tagMapSumYResolution[(*tagMapIt).first], 2)) / pow(v.tagMapCountYResolution[(*tagMapIt).first], 2))) / Units::um; // mm -> um
	}
	else {
	  aYResolutionRmse << "n/a"; // resolution is parametrized but no run with -r or -R
	}
      }

      // RphiResolution (trigger)
      anRphiResolutionTrigger.str("");
      if ( v.tagMapAveRphiResolutionTrigger[(*tagMapIt).first] != v.tagMapAveRphiResolution[(*tagMapIt).first] )
        anRphiResolutionTrigger << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << v.tagMapAveRphiResolutionTrigger[(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um
      // YResolution (trigger)
      aYResolutionTrigger.str("");
      if ( v.tagMapAveYResolutionTrigger[(*tagMapIt).first] != v.tagMapAveYResolution[(*tagMapIt).first] )
        aYResolutionTrigger << std::dec << std::fixed << std::setprecision(rphiResolutionPrecision) << v.tagMapAveYResolutionTrigger [(*tagMapIt).first] / v.tagMapCount[(*tagMapIt).first] / Units::um; // mm -> um

      // Pitches
      aPitchPair.str("");
      loPitch=int((*tagMapIt).second->outerSensor().minPitch() / Units::um); // mm -> um
      hiPitch=int((*tagMapIt).second->outerSensor().maxPitch() / Units::um); // mm -> um
      addOccupancyElement((loPitch+hiPitch)/2);

      if (loPitch==hiPitch) {
        aPitchPair << std::dec << std::fixed << std::setprecision(pitchPrecision) << loPitch;
      } else {
        aPitchPair << std::dec << std::fixed << std::setprecision(pitchPrecision)<< loPitch
          << "/" << std::fixed << std::setprecision(pitchPrecision) << hiPitch;
      }

      // Strip Lengths and segmentation
      aStripLength.str("");
      aSegment.str("");
      // One number only if all the same
      if ((*tagMapIt).second->minSegments() == (*tagMapIt).second->maxSegments()) {
        // Strip length
        aStripLength << std::fixed << std::setprecision(stripLengthPrecision)
          << (*tagMapIt).second->length()/(*tagMapIt).second->minSegments();  // CUIDADO!!!! what happens with single sided modules????
        // Segments
        aSegment << std::dec << (*tagMapIt).second->minSegments()
          << "x" << (*tagMapIt).second->outerSensor().numROCX();
      } else { // They are different
        for (int iFace=0; iFace<(*tagMapIt).second->numSensors(); ++iFace) {
          // Strip length
          aStripLength << std::fixed << std::setprecision(stripLengthPrecision)
            << (*tagMapIt).second->length()/(*tagMapIt).second->sensors().at(iFace).numSegmentsEstimate();
          // Segments
          aSegment << std::dec << (*tagMapIt).second->sensors().at(iFace).numSegmentsEstimate()
            << "x" << (*tagMapIt).second->sensors().at(iFace).numROCX();
          if (iFace<(*tagMapIt).second->numSensors() - 1) {
            aStripLength << ", ";
            aSegment << ", ";
          }
        }
      }

      // Nstrips
      anNstrips.str("");
      if ( (*tagMapIt).second->minChannels() == (*tagMapIt).second->maxChannels()) {
        anNstrips << std::dec << (*tagMapIt).second->minChannels();
      } else {
        for (int iFace=0; iFace<(*tagMapIt).second->numSensors(); ++iFace) {
          anNstrips << std::dec << (*tagMapIt).second->sensors().at(iFace).numChannels();
          if (iFace<(*tagMapIt).second->numSensors()-1) anNstrips << ", ";
        }
      }


      // Number Mod
      aNumberMod.str("");
      aNumberMod << std::dec << v.tagMapCount[(*tagMapIt).first];
      // Number Sensor
      aNumberSens.str("");
      aNumberSens << std::dec << v.tagMapCount[(*tagMapIt).first]*((*tagMapIt).second->numSensors());
      // Channels
      aChannel.str("");
      aChannel << std::fixed << std::setprecision(millionChannelPrecision)
        << v.tagMapCountChan[(*tagMapIt).first] / 1e6 ;

      // Power (per module and total)
      aPower.str("");
      double powerPerModule =  tagMapIt->second->totalPower(); // power [mW] of a module with this # strips // CUIDADO needs to take into account numChannels
      aPower << std::fixed << std::setprecision(totalPowerPrecision) << powerPerModule * v.tagMapCount[tagMapIt->first] * 1e-6; // converted from mW to kW 
      // number of modules of this type

      aPowerPerModule.str("");
      aPowerPerModule << std::fixed << std::setprecision(modulePowerPrecision) << powerPerModule;
      totalPower += powerPerModule * v.tagMapCount[tagMapIt->first];
      // Power in sensors (per module and total)
      aSensorPower.str("");
      aSensorPowerPerModuleAvg.str("");
      aSensorPowerPerModuleMax.str("");
      double totalSensorPowerTag = v.tagMapSensorPowerAvg[tagMapIt->first];
      if (totalSensorPowerTag > 1e-6) { // non-zero checking with double
        aSensorPower << std::fixed << std::setprecision(totalPowerPrecision)
          << totalSensorPowerTag * 1e-3; // converted from W to kW
        aSensorPowerPerModuleAvg << std::fixed << std::setprecision(modulePowerPrecision)
          << totalSensorPowerTag / v.tagMapCount[tagMapIt->first] * 1e3; // converted from W to mW
        aSensorPowerPerModuleMax << std::fixed << std::setprecision(modulePowerPrecision)
          << v.tagMapSensorPowerMax[tagMapIt->first] * 1e3; // converted from W to mW
      } else { // if sensor power is 0, it means we haven't run the power analysis (-p)
        aSensorPower << "n/a";
        aSensorPowerPerModuleAvg << "n/a";
        aSensorPowerPerModuleMax << "n/a";
      }

      // Cost
      aCost.str("");
      aCost  << std::fixed << std::setprecision(costPrecision) <<
        (*tagMapIt).second->area() * 1e-2 *          // area in cm^2
        (*tagMapIt).second->numSensors() *               // number of faces
        simparms.calcCost((*tagMapIt).second->readoutType()) * // price in CHF*cm^-2
        1e-6 *                                           // conversion CHF-> MCHF
        v.tagMapCount[(*tagMapIt).first];                // Number of modules
      totalCost +=(*tagMapIt).second->area() * 1e-2 * (*tagMapIt).second->numSensors() * simparms.calcCost((*tagMapIt).second->readoutType()) * 1e-6 * v.tagMapCount[(*tagMapIt).first];

      // Weight
      aWeight.str("");
      TagMaker tmak(*aModule);
      if (tagMapWeight[tmak.sensorGeoTag] > 1e-6) { // non-zero check for double
        aWeight << std::fixed << std::setprecision(weightPrecision) <<
          tagMapWeight[tmak.sensorGeoTag] / v.tagMapCount[(*tagMapIt).first];
      } else {
        aWeight << "n/a";
      }
      moduleTotalWeight += tagMapWeight[tmak.sensorGeoTag];

      moduleTable->setContent(0, iType, aName.str());
      moduleTable->setContent(tagRow, iType, aTag.str());
      moduleTable->setContent(typeRow, iType, aType.str());
      moduleTable->setContent(stripOccupancyRow, iType, aStripOccupancy.str());
      moduleTable->setContent(hitOccupancyRow, iType, aHitOccupancy.str());
      moduleTable->setContent(rphiResolutionRow, iType, anRphiResolution.str());
      moduleTable->setContent(rphiResolutionRmseRow, iType, anRphiResolutionRmse.str());
      moduleTable->setContent(yResolutionRow, iType, aYResolution.str());
      moduleTable->setContent(yResolutionRmseRow, iType, aYResolutionRmse.str());
      moduleTable->setContent(rphiResolutionTriggerRow, iType, anRphiResolutionTrigger.str());
      moduleTable->setContent(yResolutionTriggerRow, iType, aYResolutionTrigger.str());
      moduleTable->setContent(pitchpairsRow, iType, aPitchPair.str());
      moduleTable->setContent(striplengthRow, iType, aStripLength.str());
      moduleTable->setContent(segmentsRow, iType, aSegment.str());
      moduleTable->setContent(nstripsRow, iType, anNstrips.str());
      moduleTable->setContent(numbermodsRow, iType, aNumberMod.str());
      moduleTable->setContent(numbersensRow, iType, aNumberSens.str());
      moduleTable->setContent(powerRow, iType, aPower.str());
      moduleTable->setContent(powerPerModuleRow, iType, aPowerPerModule.str());
      moduleTable->setContent(sensorPowerRow, iType, aSensorPower.str());
      moduleTable->setContent(sensorPowerPerModuleAvgRow, iType, aSensorPowerPerModuleAvg.str());
      moduleTable->setContent(sensorPowerPerModuleMaxRow, iType, aSensorPowerPerModuleMax.str());
      moduleTable->setContent(costRow, iType, aCost.str());
      moduleTable->setContent(moduleWeightRow, iType, aWeight.str());

      moduleTable->setContent(thicknessRow, iType, aThickness.str());
      moduleTable->setContent(moduleAreaRow, iType, aModuleArea.str());
      moduleTable->setContent(totalAreaRow, iType, aTotalArea.str());
      moduleTable->setContent(channelRow, iType, aChannel.str());
      // moduleTable->setContent(areaRow, iType, anArea.str());

    }

    // Summary in short
    setSummaryString(tracker.myid());
    setSummaryLabelString("Name");
    addSummaryElement(v.totArea/1e6);
    addSummaryElement(v.totCountMod);
    addSummaryElement(v.totCountSens);
    addSummaryElement(v.totChannel / 1e6);
    addSummaryElement(totalPower);
    addSummaryElement(totalCost);
    addSummaryElement(moduleTotalWeight/1.e3);
    addSummaryElement(inactiveSurfacesTotalMass/1.e3);
    addSummaryElement((moduleTotalWeight+inactiveSurfacesTotalMass)/1.e3);

    addSummaryLabelElement("Area (total) mq");
    addSummaryLabelElement("Modules");
    addSummaryLabelElement("Sensors");
    addSummaryLabelElement("Standard channels (M)");
    addSummaryLabelElement("pt channels (M)");
    addSummaryLabelElement("Power (kW)");
    addSummaryLabelElement("Cost (MCHF)");
    addSummaryLabelElement("Modules weight (kg)");
    addSummaryLabelElement("Services weight (kg)");
    addSummaryLabelElement("Total weight (kg)");

    // Score totals
    ++iType;
    moduleTable->setContent(0, iType, "Total");
    moduleTable->setContent(tagRow, iType, "");
    moduleTable->setContent(typeRow, iType, "");
    aTotalArea.str("");
    aTotalArea << emphStart << std::fixed << std::setprecision(areaPrecision) << v.totArea/1e6 << emphEnd;
    //moduleTable->setContent(areaRow, iType, anArea.str());
    //anArea.str("");
    //anArea << emphStart << std::fixed << std::setprecision(areaPrecision) << totArea/1e6
    //<< "(m" << superStart << "2" << superEnd << ")" << emphEnd;
    // moduleTable->setContent(totalAreaRow, iType, aTotalArea.str());
    moduleTable->setContent(totalAreaRow, iType, aTotalArea.str());
    moduleTable->setContent(stripOccupancyRow, iType, "");
    moduleTable->setContent(hitOccupancyRow, iType, "");
    moduleTable->setContent(rphiResolutionRow, iType, "");
    moduleTable->setContent(yResolutionRow, iType, "");
    moduleTable->setContent(pitchpairsRow, iType, "");
    moduleTable->setContent(striplengthRow, iType, "");
    moduleTable->setContent(segmentsRow, iType, "");
    moduleTable->setContent(nstripsRow, iType, "");
    aNumberMod.str("");
    aNumberMod << emphStart << v.totCountMod << emphEnd;
    aNumberSens.str("");
    aNumberSens << emphStart << v.totCountSens << emphEnd;
    moduleTable->setContent(numbermodsRow, iType, aNumberMod.str());
    moduleTable->setContent(numbersensRow, iType, aNumberSens.str());
    aChannel.str("");
    aChannel << emphStart << std::fixed
      << std::setprecision(millionChannelPrecision)
      << v.totChannel / 1e6 << emphEnd;
    moduleTable->setContent(channelRow, iType, aChannel.str());
    // aChannel.str("");
    // aChannel << emphStart << std::fixed
    // 	 << std::setprecision(millionChannelPrecision)
    // 	 << totChannelPts / 1e6 << emphEnd;
    // moduleTable->setContent(channelptRow, iType, aChannel.str());
    aPower.str("");
    aPowerPerModule.str("");
    aSensorPower.str("");
    aSensorPowerPerModuleAvg.str("");
    aSensorPowerPerModuleMax.str("");
    aCost.str("");
    aPower << std::fixed << std::setprecision(totalPowerPrecision) << totalPower * 1e-6;
    if (v.totalSensorPower > 1e-6) { // non-zero check for double
      aSensorPower << std::fixed << std::setprecision(totalPowerPrecision) << v.totalSensorPower * 1e-3;
    } else {
      aSensorPower << "n/a";
    }
    aCost    << std::fixed << std::setprecision(costPrecision) << totalCost;
    moduleTable->setContent(powerRow, iType, aPower.str());
    moduleTable->setContent(powerPerModuleRow, iType, aPowerPerModule.str());
    moduleTable->setContent(sensorPowerRow, iType, aSensorPower.str());
    moduleTable->setContent(sensorPowerPerModuleAvgRow, iType, aSensorPowerPerModuleAvg.str());
    moduleTable->setContent(sensorPowerPerModuleMaxRow, iType, aSensorPowerPerModuleMax.str());
    moduleTable->setContent(costRow, iType, aCost.str());
    aWeight.str("");
    if (moduleTotalWeight > 1e-6) { // non-zero check for double
      aWeight << std::fixed << std::setprecision(weightPrecision) << moduleTotalWeight/1.e3 << " (kg)";
    } else {
      aWeight << "n/a";
    }
    moduleTable->setContent(moduleWeightRow, iType, aWeight.str());
    aWeight.str("");
    if (inactiveSurfacesTotalMass > 1e-6) {
      aWeight << std::fixed << std::setprecision(weightPrecision) << inactiveSurfacesTotalMass/1.e3 << " (kg)";
    } else {
      aWeight << "n/a";
    }
    moduleTable->setContent(inactiveWeightRow, iType, aWeight.str());
    aWeight.str("");
    if (inactiveSurfacesTotalMass+moduleTotalWeight > 1e-6) {
      aWeight << std::fixed << std::setprecision(weightPrecision) << (moduleTotalWeight+inactiveSurfacesTotalMass)/1.e3 << " (kg)";
    } else {
      aWeight << "n/a";
    }
    moduleTable->setContent(totalWeightRow, iType, aWeight.str());

    //********************************//
    //*                              *//
    //*       Plots                  *//
    //*                              *//
    //********************************//
    RootWImage* myImage;
    TCanvas *summaryCanvas = NULL;
    TCanvas *RZCanvas = NULL;
    TCanvas *XYCanvas = NULL;
    TCanvas *XYCanvasEC = NULL;
    TCanvas *myCanvas = NULL;
    //createSummaryCanvas(tracker.getMaxL(), tracker.getMaxR(), analyzer, summaryCanvas, YZCanvas, XYCanvas, XYCanvasEC);
    createSummaryCanvasNicer(tracker, RZCanvas, XYCanvas, XYCanvasEC);
    if (name=="pixel") {
      logINFO("PIXEL HACK for beam pipe");
      TPolyLine* beampipe  = new TPolyLine();
      beampipe->SetPoint(0, 0, 45/2.);
      beampipe->SetPoint(1, 2915/2., 45/2.);
      beampipe->SetPoint(2, 3804/2., 56.6/2.);
      beampipe->SetPoint(3, 3804/2.+1164, 91/2.);
      XYCanvasEC->cd();
      drawCircle(22.5, true, 18); // "grey18"
      XYCanvas->cd();
      drawCircle(22.5, true, 18); // "grey18"
      RZCanvas->cd();
      beampipe->Draw("same");

      TPolyLine* etafour  = new TPolyLine();
      etafour->SetPoint(0, 0, 0);
      etafour->SetPoint(1, 2700, 98.9376398798);
      etafour->Draw("same");
    }
    // createColorPlotCanvas(tracker, 1, RZCanvas);


    //TVirtualPad* myPad;
    myContent = new RootWContent("Plots");
    myPage->addContent(myContent);

    if (summaryCanvas) {
      myImage = new RootWImage(summaryCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Tracker summary: modules position in XY (endcap and barrel), YZ and number of hits vs. eta");
      myContent->addItem(myImage);
    }

    if (RZCanvas) {
      myImage = new RootWImage(RZCanvas, RZCanvas->GetWindowWidth(), RZCanvas->GetWindowHeight() );
      myImage->setComment("RZ position of the modules");
      myContent->addItem(myImage);
    }
    if (XYCanvas) {
      myImage = new RootWImage(XYCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("XY Section of the tracker barrel");
      myContent->addItem(myImage);
    }
    if (XYCanvasEC) {
      myImage = new RootWImage(XYCanvasEC, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("XY Projection of the tracker endcap(s)");
      myContent->addItem(myImage);
    }

    /*
     * myCanvas = new TCanvas("XYViewBarrel", "XYViewBarrel", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
     * myCanvas->cd();
     * myPad = summaryCanvas->GetPad(padXY);
     * if (myPad) {
     * myPad->DrawClonePad();
     * myImage = new RootWImage(myCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
     * myImage->setComment("XY Section of the tracker barrel");
     * myContent->addItem(myImage);
     * }
     *
     * myCanvas = new TCanvas("XYViewEndcap", "XYViewEndcap", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
     * myCanvas->cd();
     * myPad = summaryCanvas->GetPad(padEC);
     * if (myPad) {
     * myPad->DrawClonePad();
     * myImage = new RootWImage(myCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
     * myImage->setComment("XY View of the tracker endcap");
     * myContent->addItem(myImage);
     * }
     */

    // Eta profile big plot
    myCanvas = new TCanvas("EtaProfileHits", "Eta profile (Hit Modules)", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    drawEtaProfiles(*myCanvas, analyzer);
    myImage = new RootWImage(myCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit modules across eta");
    myContent->addItem(myImage);

    myCanvas = new TCanvas("EtaProfileSensors", "Eta profile (Hits)", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    drawEtaProfilesSensors(*myCanvas, analyzer);
    myImage = new RootWImage(myCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit coverage across eta");
    myContent->addItem(myImage);

    myCanvas = new TCanvas("EtaProfileStubs", "Eta profile (Stubs)", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    drawEtaProfilesStubs(*myCanvas, analyzer);
    myImage = new RootWImage(myCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Stub coverage across eta");
    myContent->addItem(myImage);

    if (name != "pixel") totalEtaProfileSensors_ = &analyzer.getTotalEtaProfileSensors();
    else totalEtaProfileSensorsPixel_ = &analyzer.getTotalEtaProfileSensors();

    TCanvas* hitMapCanvas = new TCanvas("hitmapcanvas", "Hit Map", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    hitMapCanvas->cd();
    //gStyle->SetPalette(1);
    hitMapCanvas->SetFillColor(color_plot_background);
    hitMapCanvas->SetBorderMode(0);
    hitMapCanvas->SetBorderSize(0);
    analyzer.getMapPhiEta().Draw("colz");
    analyzer.getMapPhiEta().SetStats(0);
    hitMapCanvas->Modified();
    myImage = new RootWImage(hitMapCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit coverage in eta, phi");
    myContent->addItem(myImage);

    drawEtaCoverage(*myPage, analyzer);
    drawEtaCoverageStubs(*myPage, analyzer);

    // Add detailed geometry info here
    RootWContent* filesContent = new RootWContent("Geometry files", false);
    myPage->addContent(filesContent);
    RootWTextFile* myTextFile;
    // Barrel coordinates
    myTextFile = new RootWTextFile(Form("barrelCoordinates%s.csv", name.c_str()), "Barrel modules coordinate file");
    myTextFile->addText(createBarrelModulesCsv(tracker));
    filesContent->addItem(myTextFile);
    // Endcap coordinates
    myTextFile = new RootWTextFile(Form("endcapCoordinates%s.csv", name.c_str()), "Endcap modules coordinate file");
    myTextFile->addText(createEndcapModulesCsv(tracker));
    filesContent->addItem(myTextFile);
    // All coordinates
    myTextFile = new RootWTextFile(Form("allCoordinates%s.csv", name.c_str()), "Complete coordinate file");
    myTextFile->addText(createAllModulesCsv(tracker));
    filesContent->addItem(myTextFile);

    // If debug requested, add modules' parametrized spatial resolution maps and distributions to website resolution tabs
    if (debugResolution) {

      std::string tag;
      if (name == "" ) tag = "tracker";
      else tag = name;

      RootWContent& parametrizedResolutionContent  = myPage->addContent("Modules parametrized spatial resolution");

      // Modules' parametrized spatial resolution maps
      std::map<std::string, TH2D>& parametrizedResolutionLocalXBarrelMap = analyzer.getParametrizedResolutionLocalXBarrelMap();
      std::map<std::string, TH2D>& parametrizedResolutionLocalYBarrelMap = analyzer.getParametrizedResolutionLocalYBarrelMap();
      std::map<std::string, TH2D>& parametrizedResolutionLocalXEndcapsMap = analyzer.getParametrizedResolutionLocalXEndcapsMap();
      std::map<std::string, TH2D>& parametrizedResolutionLocalYEndcapsMap = analyzer.getParametrizedResolutionLocalYEndcapsMap();

      // Modules' parametrized spatial resolution distributions
      std::map<std::string, TH1D>& parametrizedResolutionLocalXBarrelDistribution = analyzer.getParametrizedResolutionLocalXBarrelDistribution();
      std::map<std::string, TH1D>& parametrizedResolutionLocalYBarrelDistribution = analyzer.getParametrizedResolutionLocalYBarrelDistribution();
      std::map<std::string, TH1D>& parametrizedResolutionLocalXEndcapsDistribution = analyzer.getParametrizedResolutionLocalXEndcapsDistribution();
      std::map<std::string, TH1D>& parametrizedResolutionLocalYEndcapsDistribution = analyzer.getParametrizedResolutionLocalYEndcapsDistribution();

      if (parametrizedResolutionLocalXBarrelMap[tag].GetEntries() == 0 && parametrizedResolutionLocalYBarrelMap[tag].GetEntries() == 0 && parametrizedResolutionLocalXEndcapsMap[tag].GetEntries() == 0 && parametrizedResolutionLocalYEndcapsMap[tag].GetEntries() == 0) {
	parametrizedResolutionContent.addText(Form("Spatial resolution is not parametrized for any module. To get spatial resolution values, please have a look at modules table."));
      }

      // If maps not empty, add modules' parametrized spatial resolution maps and corresponding distributions
      else {
	RootWTable* myTable = new RootWTable();
	parametrizedResolutionContent.addItem(myTable);
	myTable->setContent(1, 1, "These plots are for all modules types, over full barrel or endcaps volumes.");
	myTable->setContent(2, 1, "To get statistics per module type, please have a look at modules table.");
	myTable->setContent(3, 1, " ");
	myTable->setContent(4, 1, " ");

	gStyle->SetTitleX(0.54);
	gStyle->SetTitleW(1);
	gStyle->SetOptStat("emr");
	if (parametrizedResolutionLocalXBarrelMap[tag].GetEntries() != 0) {
	  TCanvas resoXBarCanvas;
	  resoXBarCanvas.SetFillColor(color_plot_background);
	  resoXBarCanvas.Divide(2,1);
	  TVirtualPad* myPad;
	  myPad = resoXBarCanvas.GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoXBarCanvas.GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalXBarrelMap[tag].Draw();
	  myPad = resoXBarCanvas.GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalXBarrelDistribution[tag].SetStats(1);
	  parametrizedResolutionLocalXBarrelDistribution[tag].DrawNormalized();
	  RootWImage& resoXBarImage = parametrizedResolutionContent.addImage(resoXBarCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoXBarImage.setComment(Form("Resolution on local X coordinate for %s barrel modules", tag.c_str()));
	  resoXBarImage.setName(Form("Resolution on local X coordinate for %s barrel modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalYBarrelMap[tag].GetEntries() != 0) {
	  TCanvas resoYBarCanvas;
	  resoYBarCanvas.SetFillColor(color_plot_background);
	  resoYBarCanvas.Divide(2,1);
	  TVirtualPad* myPad;
	  myPad = resoYBarCanvas.GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoYBarCanvas.GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalYBarrelMap[tag].Draw();
	  myPad = resoYBarCanvas.GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalYBarrelDistribution[tag].SetStats(1);
	  parametrizedResolutionLocalYBarrelDistribution[tag].DrawNormalized();
	  RootWImage& resoYBarImage = parametrizedResolutionContent.addImage(resoYBarCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoYBarImage.setComment(Form("Resolution on local Y coordinate for %s barrel modules", tag.c_str()));
	  resoYBarImage.setName(Form("Resolution on local Y coordinate for %s barrel modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalXEndcapsMap[tag].GetEntries() != 0) {
	  TCanvas resoXEndCanvas;
	  resoXEndCanvas.SetFillColor(color_plot_background);
	  resoXEndCanvas.Divide(2,1);
	  TVirtualPad* myPad;
	  myPad = resoXEndCanvas.GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoXEndCanvas.GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalXEndcapsMap[tag].Draw();
	  myPad = resoXEndCanvas.GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalXEndcapsDistribution[tag].SetStats(1);
	  parametrizedResolutionLocalXEndcapsDistribution[tag].DrawNormalized();
	  RootWImage& resoXEndImage = parametrizedResolutionContent.addImage(resoXEndCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoXEndImage.setComment(Form("Resolution on local X coordinate for %s endcaps modules", tag.c_str()));
	  resoXEndImage.setName(Form("Resolution on local X coordinate for %s endcaps modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalYEndcapsMap[tag].GetEntries() != 0) {
	  TCanvas resoYEndCanvas;
	  resoYEndCanvas.SetFillColor(color_plot_background);
	  resoYEndCanvas.Divide(2,1);
	  TVirtualPad* myPad;
	  myPad = resoYEndCanvas.GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoYEndCanvas.GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalYEndcapsMap[tag].Draw();
	  myPad = resoYEndCanvas.GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalYEndcapsDistribution[tag].SetStats(1);
	  parametrizedResolutionLocalYEndcapsDistribution[tag].DrawNormalized();
	  RootWImage& resoYEndImage = parametrizedResolutionContent.addImage(resoYEndCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoYEndImage.setComment(Form("Resolution on local Y coordinate for %s endcaps modules", tag.c_str()));
	  resoYEndImage.setName(Form("Resolution on local Y coordinate for %s endcaps modules", tag.c_str()));
	}
      }
    } // debugResolution

    return true;
  }

  // Draws all the profile plots present in the analyzer into the given TCanvas
  // @param myPad the target TPad
  // @param analyzer the plot data container
  bool Vizard::drawEtaProfiles(TVirtualPad& myPad, Analyzer& analyzer) {
    myPad.cd();
    myPad.SetFillColor(color_plot_background);
    TProfile& totalEtaProfile = analyzer.getTotalEtaProfile();
    std::vector<TProfile>& etaProfiles = analyzer.getTypeEtaProfiles();
    return drawEtaProfilesAny(totalEtaProfile, etaProfiles);
  }

  bool Vizard::drawEtaProfilesSensors(TVirtualPad& myPad, Analyzer& analyzer, bool total/*=true*/) {
    myPad.cd();
    myPad.SetFillColor(color_plot_background);
    TProfile& totalEtaProfileSensors = analyzer.getTotalEtaProfileSensors();
    std::vector<TProfile>& etaProfilesSensors = analyzer.getTypeEtaProfilesSensors();
    return drawEtaProfilesAny(totalEtaProfileSensors, etaProfilesSensors, total);
  }

  bool Vizard::drawEtaProfilesStubs(TVirtualPad& myPad, Analyzer& analyzer) {
    myPad.cd();
    myPad.SetFillColor(color_plot_background);
    TProfile& totalEtaProfileStubs = analyzer.getTotalEtaProfileStubs();
    std::vector<TProfile>& etaProfilesStubs = analyzer.getTypeEtaProfilesStubs();
    return drawEtaProfilesAny(totalEtaProfileStubs, etaProfilesStubs);
  }

  bool Vizard::drawEtaProfilesAny(TProfile& totalEtaProfile, std::vector<TProfile>& etaProfiles, bool total/*=true*/) {
    std::vector<TProfile>::iterator etaProfileIterator;
    //totalEtaProfile.SetMaximum(15); // TODO: make this configurable
    totalEtaProfile.SetMinimum(0); // TODO: make this configurable

    if (total) totalEtaProfile.Draw();
    for (etaProfileIterator=etaProfiles.begin();
         etaProfileIterator!=etaProfiles.end();
         ++etaProfileIterator) {
      etaProfileIterator->SetMinimum(0);
      (*etaProfileIterator).Draw("same");
    }
    return true; // TODO: make this meaningful
  }
  // Draws all the profile plots present in the analyzer into the given TCanvas
  // @param myCanvas the target TCanvas
  // @param analyzer the plot data container
  bool Vizard::drawEtaProfiles(TCanvas& myCanvas, Analyzer& analyzer) {
    TVirtualPad* myVirtualPad = myCanvas.GetPad(0);
    if (!myVirtualPad) return false;
    return drawEtaProfiles(*myVirtualPad, analyzer);
  }

  bool Vizard::drawEtaProfilesSensors(TCanvas& myCanvas, Analyzer& analyzer, bool total/*=true*/) {
    TVirtualPad* myVirtualPad = myCanvas.GetPad(0);
    if (!myVirtualPad) return false;
    return drawEtaProfilesSensors(*myVirtualPad, analyzer, total);
  }

  bool Vizard::drawEtaProfilesStubs(TCanvas& myCanvas, Analyzer& analyzer) {
    TVirtualPad* myVirtualPad = myCanvas.GetPad(0);
    if (!myVirtualPad) return false;
    return drawEtaProfilesStubs(*myVirtualPad, analyzer);
  }


  bool Vizard::drawEtaCoverage(RootWPage& myPage, Analyzer& analyzer) {
    return drawEtaCoverageAny(myPage, analyzer.getLayerEtaCoverageProfiles(), "Hits");
  }

  bool Vizard::drawEtaCoverageStubs(RootWPage& myPage, Analyzer& analyzer) {
    return drawEtaCoverageAny(myPage, analyzer.getLayerEtaCoverageProfilesStubs(), "Stubs");
  }

  bool Vizard::drawEtaCoverageAny(RootWPage& myPage, std::map<std::string, TProfile>& layerEtaCoverage, const std::string& type) {
    if (layerEtaCoverage.size()==0) return false;

    TCanvas* myCanvas;
    RootWContent* myContent = new RootWContent("Layer coverage (" + type + ")", false);
    myPage.addContent(myContent);

    int layerCount = 0;
    for (std::map<std::string, TProfile>::iterator it = layerEtaCoverage.begin(); it!= layerEtaCoverage.end(); ++it) {
      TProfile& aProfile = it->second;
      layerCount++;
      myCanvas = new TCanvas(Form("LayerCoverage%s%s", it->first.c_str(), type.c_str()), ("Layer eta coverage (" + type + ")").c_str(), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myCanvas->cd();


      TPad* upperPad = new TPad(Form("%s_upper", myCanvas->GetName()), "upper", 0, 0.4, 1, 1);
      TPad* lowerPad = new TPad(Form("%s_lower", myCanvas->GetName()), "upper", 0, 0, 1, 0.4);
      myCanvas->cd();
      upperPad->Draw();
      lowerPad->Draw();
      aProfile.SetMinimum(0);
      aProfile.SetMaximum(1.05);
      aProfile.SetMarkerColor(Palette::color(1));
      aProfile.SetLineColor(Palette::color(1));
      aProfile.SetMarkerStyle(1);

      TProfile* zoomedProfile = (TProfile*) aProfile.Clone();
      zoomedProfile->SetMinimum(0.9);
      zoomedProfile->SetMaximum(1.01);
      zoomedProfile->SetTitle("");
      upperPad->cd();
      aProfile.Draw();
      lowerPad->cd();
      zoomedProfile->Draw();

      RootWImage* myImage = new RootWImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Layer coverage in eta for " + type + " (multiple occurrences in the same layer are counted once here)");
      myContent->addItem(myImage);
    }
    return true;
  }

  TCanvas* Vizard::drawFullLayout() {

    TCanvas* result = nullptr;
    std::string aClass;
    PlotDrawer<YZ, Type> yzDrawer;
    //double maxR = 1100;
    //double maxL = 2800;
    //double scaleFactor = maxR / 600;

    for (unsigned int i=0; i< trackers_.size(); ++i) {
      Tracker& tracker = *(trackers_[i]);
      yzDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    }

    int rzCanvasX = vis_max_canvas_sizeX; //int(maxL/scaleFactor);
    int rzCanvasY = vis_min_canvas_sizeY; //int(maxR/scaleFactor);
    result = new TCanvas("FullRZCanvas", "RZView Canvas (full layout)", rzCanvasX, rzCanvasY );
    result->cd();
    yzDrawer.drawFrame<SummaryFrameStyle>(*result);
    yzDrawer.drawModules<ContourStyle>(*result);

    return result;
  }
  

  bool Vizard::additionalInfoSite(const std::string& settingsfile,
                                  Analyzer& analyzer, Analyzer& pixelAnalyzer, Tracker& tracker, SimParms& simparms, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Info");
    myPage->setAddress("info.html");
    site.addPage(myPage);
    RootWContent *simulationContent, *summaryContent, *fullLayoutContent;
    RootWBinaryFile* myBinaryFile;
    std::string trackerName = tracker.myid();

    int materialTracksUsed = analyzer.getMaterialTracksUsed();
    int geometryTracksUsed = analyzer.getGeometryTracksUsed();

    //********************************//
    //*                              *//
    //*  Simulation and files        *//
    //*                              *//
    //********************************//
    
    // Detector full layout
    TCanvas* aLayout = drawFullLayout();
    if (aLayout) {
      fullLayoutContent = new RootWContent("Full layout", true);
      RootWImage* anImage = new RootWImage(aLayout, aLayout->GetWindowWidth(), aLayout->GetWindowHeight() );
      anImage->setComment("RZ position of the modules (full layout)");
      fullLayoutContent->addItem(anImage);
    }

    // Define web-page sections
    if (aLayout) myPage->addContent(fullLayoutContent);
    simulationContent = new RootWContent("Simulation parameters");
    myPage->addContent(simulationContent);
    summaryContent = new RootWContent("Summary");
    myPage->addContent(summaryContent);
     
    THStack* totalEtaStack = new THStack();
    if (totalEtaProfileSensors_) totalEtaStack->Add(totalEtaProfileSensors_->ProjectionX());
    if (totalEtaProfileSensorsPixel_) totalEtaStack->Add(totalEtaProfileSensorsPixel_->ProjectionX());
    TCanvas* totalEtaProfileFull = new TCanvas("TotalEtaProfileFull", "Full eta profile (Hits)", vis_std_canvas_sizeX, vis_std_canvas_sizeY);
    totalEtaProfileFull->cd();
    ((TH1I*)totalEtaStack->GetStack()->Last())->SetMarkerStyle(8);
    ((TH1I*)totalEtaStack->GetStack()->Last())->SetMarkerSize(1);
    ((TH1I*)totalEtaStack->GetStack()->Last())->SetMinimum(0.);
    totalEtaStack->GetStack()->Last()->Draw();
    // add profile for types here...##### 
    drawEtaProfilesSensors(*totalEtaProfileFull, analyzer, false);
    drawEtaProfilesSensors(*totalEtaProfileFull, pixelAnalyzer, false);
    totalEtaStack->GetStack()->Last()->Draw("same");
    RootWImage* myImage = new RootWImage(totalEtaProfileFull, vis_std_canvas_sizeX, vis_std_canvas_sizeY);
    myImage->setComment("Full hit coverage across eta");
    fullLayoutContent->addItem(myImage);


    RootWInfo* cmdLineInfo;
    cmdLineInfo = new RootWInfo("Command line arguments");
    cmdLineInfo->setValue(commandLine_);
    simulationContent->addItem(cmdLineInfo);

    std::string destinationFilename;

    // Let's add the included files from the main configuration
    mainConfigHandler& mainConfig = mainConfigHandler::instance();
    std::vector<std::string> origSet;
    std::vector<std::string> destSet;
    mainConfig.createEncodedFileList(origSet, destSet);
    RootWBinaryFileList* myBinaryFileList = new RootWBinaryFileList(destSet.begin(), destSet.end(), "Geometry configuration file(s)", origSet.begin(), origSet.end());
    simulationContent->addItem(myBinaryFileList);

    RootWInfo* myInfo;
    myInfo = new RootWInfo("Minimum bias per bunch crossing");
    myInfo->setValue(simparms.numMinBiasEvents(), minimumBiasPrecision);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo("Number of tracks used for material");
    myInfo->setValue(materialTracksUsed);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo("Number of tracks used for geometry");
    myInfo->setValue(geometryTracksUsed);
    simulationContent->addItem(myInfo);

    // TODO: make an object that handles this properly:
    RootWTable* myTable = new RootWTable();
    myTable->setContent(1, 0, "CHF/cm"+superStart+"2"+superEnd);
    //myTable->setContent(2, 0, "mW/channel");
    myTable->setContent(0, 1, "Pt modules");
    myTable->setContent(0, 2, "Strip modules");
    myTable->setContent(1, 1, simparms.calcCost(READOUT_PT), costPerUnitPrecision);
    myTable->setContent(1, 2, simparms.calcCost(READOUT_STRIP), costPerUnitPrecision);
    simulationContent->addItem(myTable);

    RootWTable& typesTable = simulationContent->addTable();
    typesTable.setContent(1,0,"mW / channel [chip]");
    typesTable.setContent(2,0,"mW / channel [opto]");
    typesTable.setContent(3,0,"mW / channel [total]");
    typesTable.setContent(4,0,"mW / module [chip]");
    typesTable.setContent(5,0,"mW / module [opto]");
    typesTable.setContent(6,0,"mW / module [total]");
    int iType=1;
    struct ModuleTypeVisitor : public ConstGeometryVisitor {
      std::map<std::string, const Module*> typeMap;
      void visit(const Module& m) { if (!typeMap.count(m.moduleType())) typeMap[m.moduleType()] = &m; }
    };
    ModuleTypeVisitor v;
    tracker.accept(v);
    for (auto it = v.typeMap.begin(); it != v.typeMap.end(); ++it) {
      typesTable.setContent(0,iType, it->first);
      typesTable.setContent(1,iType,it->second->powerStripChip(),2);
      typesTable.setContent(2,iType,it->second->powerStripOptical(),2);
      typesTable.setContent(3,iType,it->second->totalPowerStrip(), 2);
      typesTable.setContent(4,iType,it->second->powerModuleChip(),2);
      typesTable.setContent(5,iType,it->second->powerModuleOptical(),2);
      typesTable.setContent(6,iType,it->second->totalPowerModule(), 2);
      iType++;
    }

    //********************************//
    //*                              *//
    //*  Summary files               *//
    //*                              *//
    //********************************//

    RootWTextFile* myTextFile;

    // Summary of layout and performance
    myTextFile = new RootWTextFile("summary.csv", "Summary variables csv file");
    myTextFile->addText(getSummaryLabelString()+"\n");
    myTextFile->addText(getSummaryString());
    summaryContent->addItem(myTextFile);

    // Occupancy vs. radius
    myTextFile = new RootWTextFile("occupancy.csv", "Occupancy vs. radius");
    myTextFile->addText(occupancyCsv_);
    summaryContent->addItem(myTextFile);

    // Bill of materials
    myTextFile = new RootWTextFile("materials.csv", "Bill of materials");
    myTextFile->addText(analyzer.getBillOfMaterials());
    summaryContent->addItem(myTextFile);

    createTriggerSectorMapCsv(analyzer.getTriggerSectorMap());
    myTextFile = new RootWTextFile("trigger_sector_map.csv", "Trigger Towers to Modules connections");
    myTextFile->addText(triggerSectorMapCsv_);
    summaryContent->addItem(myTextFile);

    createModuleConnectionsCsv(analyzer.getModuleConnectionMap());
    myTextFile = new RootWTextFile("module_connections.csv", "Modules to Trigger Towers connections");
    myTextFile->addText(moduleConnectionsCsv_);
    summaryContent->addItem(myTextFile);

    RootWGraphViz* myGv = new RootWGraphViz("include_graph.gv", "Include structure");
    myGv->addText(mainConfigHandler::instance().createGraphVizFile());
    summaryContent->addItem(myGv);

    return true; // TODO: make this meaningful
  }


  bool Vizard::bandwidthSummary(Analyzer& analyzer, Tracker& tracker, SimParms& simparms, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Bandwidth");
    myPage->setAddress("bandwidth.html");
    site.addPage(myPage);
    RootWContent* myContent;

    //********************************//
    //*                              *//
    //*       Bandwidth              *//
    //*                              *//
    //********************************//
    // (also todo: handle this properly: with a not-hardcoded model)
    myContent = new RootWContent("Distributions and models");
    myPage->addContent(myContent);
    TCanvas* bandWidthCanvas = new TCanvas("ModuleBandwidthC", "Modules needed bandwidthC", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    TCanvas* moduleHitCanvas = new TCanvas("ModuleHitC", "Module hit countC", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    bandWidthCanvas->SetLogy(1);
    moduleHitCanvas->SetLogy(1);

    bandWidthCanvas->cd();
    TH1D& bandwidthDistribution = analyzer.getBandwidthDistribution();
    TH1D& bandwidthDistributionSparsified = analyzer.getBandwidthDistributionSparsified();
    bandwidthDistribution.Draw();
    bandwidthDistributionSparsified.Draw("same");
    TLegend* myLegend = new TLegend(0.75, 0.5, 1, .75);
    myLegend->AddEntry(&bandwidthDistribution, "Unsparsified", "l");
    myLegend->AddEntry(&bandwidthDistributionSparsified, "Sparsified", "l");
    myLegend->Draw();
    RootWImage* myImage = new RootWImage(bandWidthCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Module bandwidth distribution in the sparsified and unsparsified model");
    myContent->addItem(myImage);

    moduleHitCanvas->cd();
    TH1D& chanHitDistribution = analyzer.getChanHitDistribution();
    chanHitDistribution.Draw();
    myImage = new RootWImage(moduleHitCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Distribution of number of hits per bunch crossing (each sensor is counted separately)");
    myContent->addItem(myImage);

    RootWText* myDescription = new RootWText();
    myContent->addItem(myDescription);
    myDescription->addText( "Bandwidth useage estimate:<br/>");
    myDescription->addText( "(Pt modules: ignored)<br/>");
    myDescription->addText( "Sparsified (binary) bits/event: 23 bits/chip + 9 bit/hit<br/>");
    myDescription->addText( "Unsparsified (binary) bits/event: 16 bits/chip + 1 bit/channel<br/>");
    ostringstream aStringStream; aStringStream.str("100 kHz trigger, "); aStringStream << simparms.numMinBiasEvents();
    aStringStream <<" minimum bias events assumed</br>";
    myDescription->addText( aStringStream.str() );


    std::map<std::string, SummaryTable>& particleSummaries = analyzer.getTriggerFrequencyInterestingSummaries();
    std::map<std::string, SummaryTable>& trueSummaries = analyzer.getTriggerFrequencyTrueSummaries();
    std::map<std::string, SummaryTable>& fakeSummaries = analyzer.getTriggerFrequencyFakeSummaries();
    std::map<std::string, SummaryTable>& misfilteredSummaries = analyzer.getTriggerFrequencyMisfilteredSummaries();
    std::map<std::string, SummaryTable>& combinatorialSummaries = analyzer.getTriggerFrequencyCombinatorialSummaries();
    std::map<std::string, SummaryTable>& rateSummaries = analyzer.getTriggerRateSummaries();
    std::map<std::string, SummaryTable>& efficiencySummaries = analyzer.getTriggerEfficiencySummaries();
    std::map<std::string, SummaryTable>& puritySummaries = analyzer.getTriggerPuritySummaries();
    std::map<std::string, SummaryTable>& dataBandwidthSummaries = analyzer.getTriggerDataBandwidthSummaries();
    std::map<std::string, SummaryTable>& stripOccupancySummaries = analyzer.getStripOccupancySummaries();
    std::map<std::string, SummaryTable>& hitOccupancySummaries = analyzer.getHitOccupancySummaries();

    for (std::map<std::string, SummaryTable>::iterator it = particleSummaries.begin(); it != particleSummaries.end(); ++it) {
      myPage->addContent(std::string("Interesting particle rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = trueSummaries.begin(); it != trueSummaries.end(); ++it) {
      myPage->addContent(std::string("True stub rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = misfilteredSummaries.begin(); it != misfilteredSummaries.end(); ++it) {
      myPage->addContent(std::string("Low-pT over-threshold stub rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = combinatorialSummaries.begin(); it != combinatorialSummaries.end(); ++it) {
      myPage->addContent(std::string("Combinatorial stub rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = fakeSummaries.begin(); it != fakeSummaries.end(); ++it) {
      myPage->addContent(std::string("Fake (over-threshold + combinatorial) stub rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = rateSummaries.begin(); it != rateSummaries.end(); ++it) {
      myPage->addContent(std::string("Total (true + fake) stub rate (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = efficiencySummaries.begin(); it != efficiencySummaries.end(); ++it) {
      myPage->addContent(std::string("Stub efficiency (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = puritySummaries.begin(); it != puritySummaries.end(); ++it) {
      myPage->addContent(std::string("Stub purity (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = dataBandwidthSummaries.begin(); it != dataBandwidthSummaries.end(); ++it) {
      myPage->addContent(std::string("Trigger data bandwidth Gbps (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = stripOccupancySummaries.begin(); it != stripOccupancySummaries.end(); ++it) {
      myPage->addContent(std::string("Strip occupancy (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
    for (std::map<std::string, SummaryTable>::iterator it = hitOccupancySummaries.begin(); it != hitOccupancySummaries.end(); ++it) {
      myPage->addContent(std::string("Hit occupancy (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }
   

    myContent = &myPage->addContent("Trigger bandwidth and frequency maps", true);

    TCanvas triggerDataBandwidthCanvas;
    TCanvas triggerFrequencyPerEventCanvas;

    PlotDrawer<YZ, Type, Max> yzbwDrawer(0, 0); // we take the MAX because the Analyzer only sweeps across the first quadrant (up to PI/2),
    PlotDrawer<YZ, Type, Max> yztfDrawer(0, 0); // so there's plenty modules in Phi which don't have their property set, but Max disregards all the 0's

    yzbwDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);
    yztfDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);

    yzbwDrawer.drawFrame<HistogramFrameStyle>(triggerDataBandwidthCanvas);
    yztfDrawer.drawFrame<HistogramFrameStyle>(triggerFrequencyPerEventCanvas);

    yzbwDrawer.drawModules<ContourStyle>(triggerDataBandwidthCanvas);
    yztfDrawer.drawModules<ContourStyle>(triggerFrequencyPerEventCanvas);

    RootWImage& triggerDataBandwidthImage = myContent->addImage(triggerDataBandwidthCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    triggerDataBandwidthImage.setComment("Map of the bandwidth for trigger data in Gbps");
    triggerDataBandwidthImage.setName("triggerDataBandwidthMap");

    RootWImage& triggerFrequencyPerEventImage = myContent->addImage(triggerFrequencyPerEventCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    triggerFrequencyPerEventImage.setComment("Map of the trigger frequency per event (a.k.a. stubs per event)");
    triggerFrequencyPerEventImage.setName("triggerFrequencyPerEventMap");




    StubRateHistos& totalHistos = analyzer.getTotalStubRateHistos();
    StubRateHistos& trueHistos = analyzer.getTrueStubRateHistos();

    std::string currCntName = "";
    for (StubRateHistos::iterator totalIt = totalHistos.begin(), trueIt = trueHistos.begin(); totalIt != totalHistos.end(); ++totalIt, ++trueIt) {
      std::pair<std::string, int> layer = totalIt->first;
      TH1D* totalHisto = totalIt->second;
      TH1D* trueHisto = trueIt->second;
      if (layer.first != currCntName) {
        myContent = &myPage->addContent(std::string("Stub rate plots (") + layer.first + ")", false);
        currCntName = layer.first;
      }
      TCanvas graphCanvas;
      graphCanvas.cd();
      totalHisto->SetLineColor(1);
      totalHisto->SetMinimum(trueHisto->GetMinimum()*.9 < 0.1 ? 0 : trueHisto->GetMinimum()*.9);
      totalHisto->Draw();
      trueHisto->SetLineColor(2);
      trueHisto->Draw("SAME");
      RootWImage& graphImage = myContent->addImage(graphCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      graphImage.setComment("Stub rate for layer " + any2str(layer.second) + " (MHz/cm^2)");
      graphImage.setName("stubRate" + layer.first + any2str(layer.second));

    }


    return true;
  }


  bool Vizard::triggerProcessorsSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Trigger CPUs");
    myPage->setAddress("trigger_cpus.html");
    site.addPage(myPage);

    SummaryTable& processorSummary = analyzer.getProcessorConnectionSummary(); 
    SummaryTable& processorCommonSummary = analyzer.getProcessorCommonConnectionSummary();
    //std::map<std::string, SummaryTable>& moduleSummaries = analyzer.getModuleConnectionSummaries();

    myPage->addContent("Processor inbound connections").addTable().setContent(processorSummary.getContent());
    RootWContent& sharedConnContent = myPage->addContent("Processor shared inbound connections", false);
    TCanvas sharedConnCanvas;
    sharedConnCanvas.cd();
    TH2I& sharedConnMap = analyzer.getProcessorCommonConnectionMap();
    sharedConnMap.GetXaxis()->LabelsOption("v");
    sharedConnMap.GetXaxis()->SetLabelSize(0.03);
    sharedConnMap.GetYaxis()->SetLabelSize(0.03);
    sharedConnMap.Draw("colz");
    RootWImage& sharedConnImage = sharedConnContent.addImage(sharedConnCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    sharedConnImage.setComment("Map of the shared processor connections (on the diagonal unshared connections are reported)");
    sharedConnImage.setName("sharedConnMap");

    sharedConnContent.addTable().setContent(processorCommonSummary.getContent());
    sharedConnContent.addText("Columns and rows both report trigger towers, in the format 't Eta# , Phi#'. Each table cell contains the number of connections the TT on the column shares with the TT on the corresponding row. On the diagonal the number of unshared (i.e. belonging to a single TT) connections for each TT is reported.");

    SummaryTable& processorBandwidthSummary = analyzer.getProcessorInboundBandwidthSummary(); 
    SummaryTable& processorStubSummary = analyzer.getProcessorInboundStubPerEventSummary(); 
    myPage->addContent("Processor inbound bandwidth Gbps").addTable().setContent(processorBandwidthSummary.getContent());
    myPage->addContent("Processor inbound stubs per event").addTable().setContent(processorStubSummary.getContent());

    // Some helper string objects
    ostringstream tempSS;
    std::string tempString;    

    RootWContent& myContent = myPage->addContent("Module outbound connection maps", true);

    TCanvas moduleConnectionEtaCanvas;
    TCanvas moduleConnectionPhiCanvas;
    TCanvas moduleConnectionEndcapPhiCanvas;

    struct EtaConnections {
      const ModuleConnectionMap& mm_;
      EtaConnections(const ModuleConnectionMap& mm) : mm_(mm) {}
      double operator()(const Module& m) { return mm_.at(&m).etaCpuConnections(); }
    };
    struct PhiConnections {
      const ModuleConnectionMap& mm_;
      PhiConnections(const ModuleConnectionMap& mm) : mm_(mm) {}
      double operator()(const Module& m) { return mm_.at(&m).phiCpuConnections(); }
    };
    PlotDrawer<YZFull, EtaConnections, Max> yzDrawer(0, 0, EtaConnections(analyzer.getModuleConnectionMap())); //(2*getDrawAreaZ(tracker), getDrawAreaR(tracker));
    PlotDrawer<XY, PhiConnections, Max> xyDrawer(0, 0, PhiConnections(analyzer.getModuleConnectionMap())); //(getDrawAreaX(tracker), getDrawAreaY(tracker));
    PlotDrawer<XY, PhiConnections, Max> xyecDrawer(0, 0, PhiConnections(analyzer.getModuleConnectionMap())); //(getDrawAreaX(tracker), getDrawAreaY(tracker));

    yzDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    //xyDrawer.addModules<CheckSection<Layer::XYSection> >(tracker.getLayers());
    xyDrawer.addModules<CheckType<BARREL>>(tracker.modules().begin(), tracker.modules().end());
    xyecDrawer.addModules<CheckType<ENDCAP>>(tracker.modules().begin(), tracker.modules().end());

    yzDrawer.drawFrame<HistogramFrameStyle>(moduleConnectionEtaCanvas);
    xyDrawer.drawFrame<HistogramFrameStyle>(moduleConnectionPhiCanvas);
    std::pair<Circle, Circle> petal = analyzer.getSampleTriggerPetal();
    TArc a1(petal.first.x0, petal.first.y0, petal.first.r, (XYPoint(petal.first.x0, petal.first.y0)).Phi()*180./M_PI + 180.);
    TArc a2(petal.second.x0, petal.second.y0, petal.second.r, 0., (XYPoint(petal.second.x0, petal.second.y0)).Phi()*180./M_PI + 180.);
    a1.SetFillStyle(0);
    a2.SetFillStyle(0);
    moduleConnectionPhiCanvas.cd();
    a1.Draw("only");
    a2.Draw("only");

    xyecDrawer.drawFrame<HistogramFrameStyle>(moduleConnectionEndcapPhiCanvas);
    moduleConnectionEndcapPhiCanvas.cd();
    a1.Draw("only");
    a2.Draw("only");

    yzDrawer.drawModules<ContourStyle>(moduleConnectionEtaCanvas);
    xyDrawer.drawModules<ContourStyle>(moduleConnectionPhiCanvas);
    xyecDrawer.drawModules<ContourStyle>(moduleConnectionEndcapPhiCanvas);



    /*
       moduleConnectionEtaCanvas.SetFillColor(color_plot_background);
       moduleConnectionPhiCanvas.SetFillColor(color_plot_background);
    //moduleConnectionEndcapPhiCanvas.SetFillColor(color_plot_background);

    moduleConnectionEtaCanvas.cd();
    moduleConnectionEtaMap.Draw("colz");
    */  
    RootWImage& moduleConnectionEtaImage = myContent.addImage(moduleConnectionEtaCanvas, vis_max_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionEtaImage.setComment("Map of the number of connections to trigger processors per module (eta section)");
    moduleConnectionEtaImage.setName("moduleConnectionEtaMap");
    /*
       moduleConnectionPhiCanvas.cd();
       moduleConnectionPhiMap.Draw("colz");
       */  
    RootWImage& moduleConnectionPhiImage = myContent.addImage(moduleConnectionPhiCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionPhiImage.setComment("Map of the number of connections to trigger processors per barrel module (phi section)");
    moduleConnectionPhiImage.setName("moduleConnectionPhiMap");

    // moduleConnectionEndcapPhiCanvas.cd();
    // moduleConnectionEndcapPhiMap.Draw("colz");

    RootWImage& moduleConnectionEndcapPhiImage = myContent.addImage(moduleConnectionEndcapPhiCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionEndcapPhiImage.setComment("Map of the number of connections to trigger processors per endcap module (phi section)");
    moduleConnectionEndcapPhiImage.setName("moduleConnectionEndcapPhiMap");

    //    myContent = myPage->addContent("Module Connections distribution", true);

    TCanvas moduleConnectionsCanvas("ModuleConnectionsC", "Modules connectionsC", vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionsCanvas.cd();
    TH1I& moduleConnectionsDistribution = analyzer.getModuleConnectionsDistribution();
    moduleConnectionsDistribution.SetFillColor(Palette::color(2));
    moduleConnectionsDistribution.Draw();
    RootWImage& myImage = myContent.addImage(moduleConnectionsCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage.setComment("Module connections distribution");

    return true;
  }

  bool Vizard::irradiatedPowerSummary(Analyzer& a, Tracker& tracker, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Power");
    myPage->setAddress("power.html");
    site.addPage(myPage);

    std::map<std::string, SummaryTable>& powerSummaries = a.getIrradiatedPowerConsumptionSummaries();
    for (std::map<std::string, SummaryTable>::iterator it = powerSummaries.begin(); it != powerSummaries.end(); ++it) {
      myPage->addContent(std::string("Power in irradiated sensors (") + it->first + ")", false).addTable().setContent(it->second.getContent());
    }

    // Some helper string objects
    ostringstream tempSS;
    std::string tempString;    

    // mapBag myMapBag = a.getMapBag();
    //TH2D& irradiatedPowerMap = myMapBag.getMaps(mapBag::irradiatedPowerConsumptionMap)[mapBag::dummyMomentum];
    // TH2D& totalPowerMap = myMapBag.getMaps(mapBag::totalPowerConsumptionMap)[mapBag::dummyMomentum];
    //
    //
    
    struct IrradiationPower {
      double operator()(const Module& m) { return m.irradiationPower(); }
    };

    struct TotalIrradiationPower {
      double operator()(const Module& m) { return m.irradiationPower() + m.totalPower()/1000; }
    };


    PlotDrawer<YZ, IrradiationPower, Average> yzPowerDrawer(0, 0);
    PlotDrawer<YZ, TotalIrradiationPower, Average> yzTotalPowerDrawer(0, 0);

    yzPowerDrawer.addModules<CheckType<BARREL | ENDCAP>>(tracker.modules().begin(), tracker.modules().end());
    yzTotalPowerDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);

    RootWContent& myContent = myPage->addContent("Power maps", true);

    TCanvas irradiatedPowerCanvas;
    TCanvas totalPowerCanvas;

    yzPowerDrawer.drawFrame<HistogramFrameStyle>(irradiatedPowerCanvas);
    yzPowerDrawer.drawModules<ContourStyle>(irradiatedPowerCanvas);


    yzTotalPowerDrawer.drawFrame<HistogramFrameStyle>(totalPowerCanvas);
    yzTotalPowerDrawer.drawModules<ContourStyle>(totalPowerCanvas);

    RootWImage& irradiatedPowerImage = myContent.addImage(irradiatedPowerCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    irradiatedPowerImage.setComment("Map of power dissipation in irradiated modules (W)");
    irradiatedPowerImage.setName("irradiatedPowerMap");
    RootWImage& totalPowerImage = myContent.addImage(totalPowerCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    totalPowerImage.setComment("Map of power dissipation in irradiated modules (W)");
    totalPowerImage.setName("totalPowerMap");


    return true;
  }

  bool Vizard::errorSummary(Analyzer& a, RootWSite& site, std::string additionalTag, bool isTrigger) {

    //********************************//
    //*                              *//
    //*    Resolution estimate       *//
    //*                              *//
    //********************************//

    // Here you should check if the TGraph
    // list is empty => maybe not?
    if (!(a.getRhoGraphs(false, isTrigger).empty() && a.getDGraphs(false, isTrigger).empty() && a.getPhiGraphs(false, isTrigger).empty())) {
      // Create a page for the errors
      std::string pageTitle = "Resolution";
      std::string additionalSummaryTag;
      double verticalScale=1;
      if (additionalTag!="") {
        pageTitle += " ("+additionalTag+")";
        additionalSummaryTag = "_"+additionalTag+"_";
        verticalScale = 10;
      } else {
        additionalSummaryTag = "";
      }
      std::string pageAddress = "errors" + additionalTag + ".html";
      RootWPage& myPage = site.addPage(pageTitle);
      myPage.setAddress(pageAddress);

      // Create the contents
      RootWContent& resolutionContent = myPage.addContent("Track resolution");
      RootWContent& idealResolutionContent = myPage.addContent("Track resolution (without material)");

      std::string scenarioStr;
      for (int scenario=0; scenario<2; ++scenario) {
        bool idealMaterial;
        RootWContent* myContent;
        if (scenario==0) {
          idealMaterial=false;
          myContent = &resolutionContent;

        } else {
          idealMaterial=true;
          myContent = &idealResolutionContent;
          scenarioStr = "noMS";
        }

        TCanvas linearMomentumCanvas;
        TCanvas momentumCanvas;
        TCanvas distanceCanvas;
        TCanvas angleCanvas;
        TCanvas ctgThetaCanvas;
        TCanvas z0Canvas;
        TCanvas pCanvas;

        int myColor=0;
        int nRebin = 2;
        int markerStyle = 21;
        double markerSize = 1.;
        double lineWidth = 2.;

        linearMomentumCanvas.SetGrid(1,1);
        momentumCanvas.SetGrid(1,1);
        distanceCanvas.SetGrid(1,1);
        angleCanvas.SetGrid(1,1);
        ctgThetaCanvas.SetGrid(1,1);
        z0Canvas.SetGrid(1,1);
        pCanvas.SetGrid(1,1);
        std::string plotOption = "";
        std::map<int, TGraph>::iterator g_iter, g_guard;
        // momentum canvas loop
        g_guard = a.getRhoGraphs(idealMaterial, isTrigger).end();
        gStyle->SetGridStyle(style_grid);
        gStyle->SetGridColor(color_hard_grid);
        for (g_iter = a.getRhoGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& momentumGraph = g_iter->second;
          TProfile& momentumProfile = newProfile(momentumGraph, 0, a.getEtaMaxTracker(), nRebin);

          if (idealMaterial) {
            momentumProfile.SetMinimum(vis_min_dPtOverPt);//1E-5*100);
            momentumProfile.SetMaximum(vis_max_dPtOverPt);//.11*100*verticalScale);
          } else {
            momentumProfile.SetMinimum(vis_min_dPtOverPt);//4E-3*100);
            momentumProfile.SetMaximum(vis_max_dPtOverPt);//.5*100*verticalScale);
          }
          linearMomentumCanvas.SetLogy(0);        
          momentumCanvas.SetLogy(1);
          momentumProfile.SetLineColor(momentumColor(myColor));
          momentumProfile.SetMarkerColor(momentumColor(myColor));
          momentumProfile.SetLineWidth(lineWidth);
          myColor++;
          momentumProfile.SetMarkerStyle(markerStyle);
          momentumProfile.SetMarkerSize(markerSize);
          momentumCanvas.SetFillColor(color_plot_background);
          linearMomentumCanvas.SetFillColor(color_plot_background);
          if (momentumGraph.GetN()>0) {
            momentumCanvas.cd();
            momentumProfile.Draw(plotOption.c_str());
            linearMomentumCanvas.cd();
            momentumProfile.Draw(plotOption.c_str());
            plotOption = "same";
          }
        }
        plotOption = "";
        myColor=0;
        // distance canvas loop
        g_guard = a.getDGraphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getDGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& distanceGraph = g_iter->second;
          TProfile& distanceProfile = newProfile(distanceGraph, 0, a.getEtaMaxTracker(), nRebin);
          if (idealMaterial) {
            distanceProfile.SetMinimum(vis_min_dD0);//4*1e-4);
            distanceProfile.SetMaximum(vis_max_dD0);//4E2*1e-4*verticalScale);
          } else {
            distanceProfile.SetMinimum(vis_min_dD0);//4*1e-4);
            distanceProfile.SetMaximum(vis_max_dD0);//4E2*1e-4*verticalScale);
          }
          distanceCanvas.SetLogy();
          distanceProfile.SetLineColor(momentumColor(myColor));
          distanceProfile.SetMarkerColor(momentumColor(myColor));
          distanceProfile.SetLineWidth(lineWidth);
          myColor++;
          distanceProfile.SetMarkerStyle(markerStyle);
          distanceProfile.SetMarkerSize(markerSize);
          distanceCanvas.SetFillColor(color_plot_background);
          if (distanceGraph.GetN()>0) {
            distanceCanvas.cd();
            distanceProfile.Draw(plotOption.c_str());
            plotOption = "same";
          }
        }
        plotOption = "";
        myColor=0;
        // angle canvas loop
        g_guard = a.getPhiGraphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getPhiGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& angleGraph = g_iter->second;
          TProfile& angleProfile = newProfile(angleGraph, 0, a.getEtaMaxTracker(), nRebin);
          if (idealMaterial) {
            angleProfile.SetMinimum(vis_min_dPhi);//1E-5);
            angleProfile.SetMaximum(vis_max_dPhi);//0.01*verticalScale);
          } else {
            angleProfile.SetMinimum(vis_min_dPhi);//1E-5);
            angleProfile.SetMaximum(vis_max_dPhi);//0.01*verticalScale);
          }
          angleCanvas.SetLogy();
          angleProfile.SetLineColor(momentumColor(myColor));
          angleProfile.SetMarkerColor(momentumColor(myColor));
          angleProfile.SetLineWidth(lineWidth);
          myColor++;
          angleProfile.SetMarkerStyle(markerStyle);
          angleProfile.SetMarkerSize(markerSize);
          angleCanvas.SetFillColor(color_plot_background);
          if (angleGraph.GetN() > 0) {
            angleCanvas.cd();
            angleProfile.Draw(plotOption.c_str());
            plotOption = "same";
          }
        }
        plotOption = "";
        myColor=0;
        // ctgTheta canvas loop
        g_guard = a.getCtgThetaGraphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getCtgThetaGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& ctgThetaGraph = g_iter->second;
          TProfile& ctgThetaProfile = newProfile(ctgThetaGraph, 0, a.getEtaMaxTracker(), nRebin);
          ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);//1E-5);
          ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);//0.1*verticalScale);
          ctgThetaCanvas.SetLogy();
          ctgThetaProfile.SetLineColor(momentumColor(myColor));
          ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          ctgThetaProfile.SetMarkerStyle(markerStyle);
          ctgThetaProfile.SetMarkerSize(markerSize);
          ctgThetaCanvas.SetFillColor(color_plot_background);
          if (ctgThetaGraph.GetN() > 0) {
            ctgThetaCanvas.cd();
            ctgThetaProfile.Draw(plotOption.c_str());
            plotOption = "same";
          }
        }
        plotOption = "";
        myColor=0;
        // z0 canvas loop
        g_guard = a.getZ0Graphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getZ0Graphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& z0Graph = g_iter->second;
          TProfile& z0Profile = newProfile(z0Graph, 0, a.getEtaMaxTracker(), nRebin);
          z0Profile.SetMinimum(vis_min_dZ0);//1E-5);
          z0Profile.SetMaximum(vis_max_dZ0);//1*verticalScale);
          z0Canvas.SetLogy();
          z0Profile.SetLineColor(momentumColor(myColor));
          z0Profile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          z0Profile.SetMarkerStyle(markerStyle);
          z0Profile.SetMarkerSize(markerSize);
          z0Canvas.SetFillColor(color_plot_background);
          if (z0Graph.GetN() > 0) {
            z0Canvas.cd();
            z0Profile.Draw(plotOption.c_str());
            plotOption = "p same";
          }
        }
        plotOption = "";
        myColor=0;
        // p canvas loop
        g_guard = a.getPGraphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getPGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& pGraph = g_iter->second;
          TProfile& pProfile = newProfile(pGraph, 0, a.getEtaMaxTracker(), nRebin);
          if (idealMaterial) {
            pProfile.SetMinimum(vis_min_dPtOverPt);//1E-5*100);
            pProfile.SetMaximum(vis_max_dPtOverPt);//.11*100*verticalScale);
          } else {
            pProfile.SetMinimum(vis_min_dPtOverPt);//4E-3*100);
            pProfile.SetMaximum(vis_max_dPtOverPt);//.11*100*verticalScale);
          }
          pCanvas.SetLogy();
          pProfile.SetLineColor(momentumColor(myColor));
          pProfile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          pProfile.SetMarkerStyle(markerStyle);
          pProfile.SetMarkerSize(markerSize);
          pCanvas.SetFillColor(color_plot_background);
          if (pGraph.GetN() > 0) {
            pCanvas.cd();
            pProfile.Draw(plotOption.c_str());
            plotOption = "p same";
          }
        }
        RootWImage& linearMomentumImage = myContent->addImage(linearMomentumCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        linearMomentumImage.setComment("Transverse momentum resolution vs. eta (linear scale)");
        linearMomentumImage.setName(Form("linptres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& momentumImage = myContent->addImage(momentumCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        momentumImage.setComment("Transverse momentum resolution vs. eta");
        momentumImage.setName(Form("ptres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& distanceImage = myContent->addImage(distanceCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        distanceImage.setComment("Distance of closest approach resolution vs. eta");
        distanceImage.setName(Form("dxyres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& angleImage = myContent->addImage(angleCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        angleImage.setComment("Angle resolution vs. eta");
        angleImage.setName(Form("phires_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& ctgThetaImage = myContent->addImage(ctgThetaCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        ctgThetaImage.setComment("CtgTheta resolution vs. eta");
        ctgThetaImage.setName(Form("cotThetares_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& z0Image = myContent->addImage(z0Canvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        z0Image.setComment("z0 resolution vs. eta");
        z0Image.setName(Form("dzres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& pImage = myContent->addImage(pCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        pImage.setComment("Momentum resolution vs. eta");
        pImage.setName(Form("pres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
      }

      // Check that the ideal and real have the same pts
      // Otherwise the table cannot be prepared

      RootWContent& summaryContent = myPage.addContent("Summary", false);
      RootWTable& cutsTable = summaryContent.addTable();
      std::vector<std::string> plotNames;
      std::map<std::string, RootWTable*> tableMap;
      std::map<std::string, RootWTable*>::iterator tableMapIt;
      plotNames.push_back("pt");
      plotNames.push_back("d");
      plotNames.push_back("phi");
      plotNames.push_back("ctg(theta)");
      plotNames.push_back("z0");
      plotNames.push_back("p");
      for (std::vector<std::string>::iterator it=plotNames.begin();
           it!=plotNames.end(); ++it) {
        tableMap[(*it)] = &(summaryContent.addTable());
        tableMap[(*it)]->setContent(0,0,(*it));
      }

      // Prepare the cuts for the averages
      ostringstream label;
      std::string name;      
      RootWTable* myTable;

      // Table explaining the cuts
      cutsTable.setContent(0,0,"Region");
      cutsTable.setContent(1,0,"etaMin");
      cutsTable.setContent(2,0,"etaMax");
      myTable = &cutsTable;
      for (unsigned int iBorder=0; iBorder<geom_name_eta_regions.size()-1; ++iBorder) {
        myTable->setContent(0,iBorder+1,geom_name_eta_regions[iBorder+1]);
        label.str(""); label << std::fixed << std::setprecision(1) << geom_range_eta_regions[iBorder];
        myTable->setContent(1,iBorder+1,label.str());
        label.str(""); label << std::fixed << std::setprecision(1) << geom_range_eta_regions[iBorder+1];
        myTable->setContent(2,iBorder+1,label.str());
      }

      std::map<graphIndex, TGraph*> myPlotMap;
      graphIndex myIndex;

      fillPlotMap(plotNames[0], myPlotMap, &a, &Analyzer::getRhoGraphs, isTrigger);
      fillPlotMap(plotNames[1], myPlotMap, &a, &Analyzer::getDGraphs, isTrigger);
      fillPlotMap(plotNames[2], myPlotMap, &a, &Analyzer::getPhiGraphs, isTrigger);
      fillPlotMap(plotNames[3], myPlotMap, &a, &Analyzer::getCtgThetaGraphs, isTrigger);
      fillPlotMap(plotNames[4], myPlotMap, &a, &Analyzer::getZ0Graphs, isTrigger);
      fillPlotMap(plotNames[5], myPlotMap, &a, &Analyzer::getPGraphs, isTrigger);

      // Cycle over the different measurements
      for (std::vector<std::string>::iterator plotNameIt = plotNames.begin();
           plotNameIt!=plotNames.end(); ++plotNameIt) {

        //std::cerr << "tableMap[\""<< *plotNameIt <<"\"] = " << tableMap[*plotNameIt] << std::endl; // debug
        myTable = tableMap[*plotNameIt];
        if (!myTable) continue;

        // Count the realistic plots' momenta
        std::vector<double> momentum;
        std::vector<double>::iterator momentumIt;

        for (std::map<graphIndex, TGraph*>::iterator myPlotMapIt = myPlotMap.begin();
             myPlotMapIt!=myPlotMap.end(); ++myPlotMapIt) {
          myIndex =  (*myPlotMapIt).first;
          //std::cerr << "Check3: myIndex.name = " << myIndex.name << std::endl; // debug
          if (myIndex.name==(*plotNameIt)) {
            //std::cerr << "found momentum " << myIndex.p <<std::endl; // debug
            momentumIt = std::find(momentum.begin(), momentum.end(), myIndex.p);
            if (momentumIt == momentum.end()) momentum.push_back(myIndex.p);
          }
        }

        std::sort(momentum.begin(), momentum.end());
        //std::cerr << "momentum.size() = " << momentum.size() <<std::endl; // debug

        // Fill the table with the values
        // First the heading of momentum
        int baseColumn;
        std::vector<double> averagesReal;
        std::vector<double> averagesIdeal;
        TGraph* myGraph;
        int myColor = kBlack;
        myIndex.name=(*plotNameIt);
        std::ostringstream myLabel;
        for (unsigned int i=0; i<momentum.size(); ++i) {
          baseColumn = (geom_name_eta_regions.size()-1)*i + 1;
          myTable->setContent(0, baseColumn, momentum[i],0);
          myIndex.p=momentum[i];
          myIndex.ideal = false;
          myGraph = myPlotMap[myIndex];
          myTable->setContent(2, 0, "Real");
          myTable->setContent(3, 0, "Ideal");
          myTable->setContent(4, 0, "Loss");
          if (myGraph) {
            averagesReal=Analyzer::average(*myGraph, geom_range_eta_regions);
            myColor = myGraph->GetMarkerColor();
            myTable->setColor(0, baseColumn, myColor);
          }
          myIndex.ideal = true;
          myGraph = myPlotMap[myIndex];
          if (myGraph) averagesIdeal=Analyzer::average(*myGraph, geom_range_eta_regions);
          for (unsigned int j=0; j<(geom_name_eta_regions.size()-1); ++j) {
            myTable->setContent(1, baseColumn+j, geom_name_eta_regions[j+1]);
            myTable->setColor(1, baseColumn+j, myColor);
            if (averagesReal.size() > j) {
              myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionHigh);
              myTable->setColor(2, baseColumn+j, myColor);
            }
            if (averagesIdeal.size() > j) {
              myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionHigh);
              myTable->setColor(3, baseColumn+j, myColor);
            }
            if ((averagesReal.size() > j)&&(averagesIdeal.size() > j)) {
              myTable->setContent(4, baseColumn+j,averagesReal[j]/averagesIdeal[j],1);
              myTable->setColor(4, baseColumn+j, myColor);
            }
            myLabel.str("");
            myLabel << myIndex.name
              << std::dec << std::fixed << std::setprecision(0) 
              << myIndex.p << "(" << geom_name_eta_regions[j+1] << ")";
            addSummaryLabelElement(myLabel.str()+additionalSummaryTag+"_Real");
            addSummaryLabelElement(myLabel.str()+additionalSummaryTag+"_Ideal");
            addSummaryElement(averagesReal[j]);
            addSummaryElement(averagesIdeal[j]);
          }
        }
      }
      return true;
    }
    return false;
  }

  bool Vizard::taggedErrorSummary(Analyzer& analyzer, RootWSite& site) {

    //********************************//
    //*                              *//
    //*    Resolution estimate       *//
    //*                              *//
    //********************************//

    GraphBag& gb = analyzer.getGraphBag();

    for (auto tag : gb.getTagSet()) {
    
      std::string pageTitle = "Resolution";
      std::string additionalSummaryTag;
      double verticalScale=1;

      pageTitle              += " ("+tag+")";
      additionalSummaryTag    = "_"+tag+"_";
      verticalScale           = 10;
      std::string pageAddress = "errors" + tag + ".html";

      RootWPage* myPage = new RootWPage(pageTitle);
      myPage->setAddress(pageAddress);
      site.addPage(myPage);

      // Create the contents
      RootWContent& resolutionContent_Pt      = myPage->addContent("Track resolution for const Pt across "  +etaLetter+" (material)");
      RootWContent& idealResolutionContent_Pt = myPage->addContent("Track resolution for const Pt across "  +etaLetter+" (no material)", false);
      RootWContent& resolutionContent_P       = myPage->addContent("Track resolution for const P across "   +etaLetter+" (material)"   , false);
      RootWContent& idealResolutionContent_P  = myPage->addContent("Track resolution for const P across "   +etaLetter+" (no material)", false);
  
      // Create a page for the errors - scenarios with/without multiple scattering (active+pasive or just active material), extra scenario includes dipole magnet
      std::string scenarioStr="";
      for (int scenario=0; scenario<2; ++scenario) {
        int idealMaterial;
        RootWContent* myContent;
  
        // Draw case I with const Pt across eta
        if (!gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::RhoGraph_Pt     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::DGraph_Pt       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::Z0Graph_Pt      , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PhiGraph_Pt     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::CtgthetaGraph_Pt, tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PGraph_Pt       , tag).empty()) {
  
          // Set link to myContent
          if (scenario==0) {
            idealMaterial=GraphBag::RealGraph;
            myContent = &resolutionContent_Pt;
            scenarioStr = "MS_Pt";
          } else {
            idealMaterial=GraphBag::IdealGraph;
            myContent = &idealResolutionContent_Pt;
            scenarioStr = "noMS_Pt";
          }
  
          TCanvas linMomCanvas_Pt;
          TCanvas logMomCanvas_Pt;
          TCanvas d0Canvas_Pt;
          TCanvas phiCanvas_Pt;
          TCanvas ctgThetaCanvas_Pt;
          TCanvas z0Canvas_Pt;
          TCanvas pCanvas_Pt;
  
          // Default attributes
          int myColor            = 0;
          int nBins              = insur::vis_n_bins;
          int markerStyle        = 21;
          double markerSize      = 1.;
          double lineWidth       = 2.;
          std::string plotOption = "";
  
          linMomCanvas_Pt.SetGrid(1,1);
          logMomCanvas_Pt.SetGrid(1,1);
          d0Canvas_Pt.SetGrid(1,1);
          phiCanvas_Pt.SetGrid(1,1);
          ctgThetaCanvas_Pt.SetGrid(1,1);
          z0Canvas_Pt.SetGrid(1,1);
          pCanvas_Pt.SetGrid(1,1);
  
          gStyle->SetGridStyle(style_grid);
          gStyle->SetGridColor(color_hard_grid);
  
          // Draw pt
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::RhoGraph_Pt | idealMaterial, tag)) {
  
            const TGraph& momentumGraph = mapel.second;
            TProfile& momentumProfile   = newProfile(momentumGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              momentumProfile.SetMinimum(insur::vis_min_dPtOverPt); //1E-5*100);
              momentumProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            } else {
              momentumProfile.SetMinimum(insur::vis_min_dPtOverPt); //4E-3*100);
              momentumProfile.SetMaximum(insur::vis_max_dPtOverPt); //.5*100*verticalScale);
            }
            linMomCanvas_Pt.SetLogy(0);
            logMomCanvas_Pt.SetLogy(1);
            linMomCanvas_Pt.SetFillColor(color_plot_background);
            logMomCanvas_Pt.SetFillColor(color_plot_background);
  
            momentumProfile.SetLineColor(momentumColor(myColor));
            momentumProfile.SetMarkerColor(momentumColor(myColor));
            momentumProfile.SetLineWidth(lineWidth);
            myColor++;
            momentumProfile.SetMarkerStyle(markerStyle);
            momentumProfile.SetMarkerSize(markerSize);
  
            if (momentumGraph.GetN()>0) {
              linMomCanvas_Pt.cd();
              momentumProfile.Draw(plotOption.c_str());
              logMomCanvas_Pt.cd();
              momentumProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw p
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::PGraph_Pt | idealMaterial, tag)) {
  
            const TGraph& pGraph = mapel.second;
            TProfile& pProfile   = newProfile(pGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              pProfile.SetMinimum(insur::vis_min_dPtOverPt); //1E-5*100);
              pProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            } else {
              pProfile.SetMinimum(insur::vis_min_dPtOverPt); //4E-3*100);
              pProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            }
            pCanvas_Pt.SetLogy();
            pCanvas_Pt.SetFillColor(color_plot_background);
  
            pProfile.SetLineColor(momentumColor(myColor));
            pProfile.SetMarkerColor(momentumColor(myColor));
            pProfile.SetLineWidth(lineWidth);
            myColor++;
            pProfile.SetMarkerStyle(markerStyle);
            pProfile.SetMarkerSize(markerSize);
  
            if (pGraph.GetN() > 0) {
              pCanvas_Pt.cd();
              pProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }
          // Draw d0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::DGraph_Pt | idealMaterial, tag)) {
  
            const TGraph& distanceGraph = mapel.second;
            TProfile& distanceProfile   = newProfile(distanceGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              distanceProfile.SetMinimum(vis_min_dD0);
              distanceProfile.SetMaximum(vis_max_dD0);//*verticalScale);
            } else {
              distanceProfile.SetMinimum(vis_min_dD0);
              distanceProfile.SetMaximum(vis_max_dD0);//*verticalScale);
            }
            d0Canvas_Pt.SetLogy();
            d0Canvas_Pt.SetFillColor(color_plot_background);
  
            distanceProfile.SetLineColor(momentumColor(myColor));
            distanceProfile.SetMarkerColor(momentumColor(myColor));
            distanceProfile.SetLineWidth(lineWidth);
            myColor++;
            distanceProfile.SetMarkerStyle(markerStyle);
            distanceProfile.SetMarkerSize(markerSize);
  
            if (distanceGraph.GetN()>0) {
              d0Canvas_Pt.cd();
              distanceProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw phi angle
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::PhiGraph_Pt | idealMaterial, tag)) {
  
            const TGraph& angleGraph = mapel.second;
            TProfile& angleProfile   = newProfile(angleGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
            angleProfile.SetMinimum(vis_min_dPhi);
            angleProfile.SetMaximum(vis_max_dPhi);//*verticalScale);
            } else {
              angleProfile.SetMinimum(vis_min_dPhi);
              angleProfile.SetMaximum(vis_max_dPhi);//*verticalScale);
            }
            phiCanvas_Pt.SetLogy();
            phiCanvas_Pt.SetFillColor(color_plot_background);
  
            angleProfile.SetLineColor(momentumColor(myColor));
            angleProfile.SetMarkerColor(momentumColor(myColor));
            angleProfile.SetLineWidth(lineWidth);
            myColor++;
            angleProfile.SetMarkerStyle(markerStyle);
            angleProfile.SetMarkerSize(markerSize);
  
            if (angleGraph.GetN() > 0) {
              phiCanvas_Pt.cd();
              angleProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw ctgTheta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::CtgthetaGraph_Pt | idealMaterial, tag)) {
  
            const TGraph& ctgThetaGraph = mapel.second;
            TProfile& ctgThetaProfile   = newProfile(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);
            ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);//*verticalScale);
            ctgThetaCanvas_Pt.SetLogy();
            ctgThetaCanvas_Pt.SetFillColor(color_plot_background);
  
            ctgThetaProfile.SetLineColor(momentumColor(myColor));
            ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
            ctgThetaProfile.SetLineWidth(lineWidth);
            myColor++;
            ctgThetaProfile.SetMarkerStyle(markerStyle);
            ctgThetaProfile.SetMarkerSize(markerSize);
  
            if (ctgThetaGraph.GetN() > 0) {
              ctgThetaCanvas_Pt.cd();
              ctgThetaProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw z0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::Z0Graph_Pt | idealMaterial, tag)) {
  
            const TGraph& z0Graph = mapel.second;
            TProfile& z0Profile   = newProfile(z0Graph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            z0Profile.SetMinimum(vis_min_dZ0);
            z0Profile.SetMaximum(vis_max_dZ0);//*verticalScale);
            z0Canvas_Pt.SetLogy();
            z0Canvas_Pt.SetFillColor(color_plot_background);
  
            z0Profile.SetLineColor(momentumColor(myColor));
            z0Profile.SetMarkerColor(momentumColor(myColor));
            z0Profile.SetLineWidth(lineWidth);
            myColor++;
            z0Profile.SetMarkerStyle(markerStyle);
            z0Profile.SetMarkerSize(markerSize);
            z0Canvas_Pt.SetFillColor(color_plot_background);
  
            if (z0Graph.GetN() > 0) {
              z0Canvas_Pt.cd();
              z0Profile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }
          plotOption = "";
          myColor=0;
  
          RootWImage& linMomImage_Pt = myContent->addImage(linMomCanvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          linMomImage_Pt.setComment("Transverse momentum resolution vs. "+etaLetter+" (linear scale) - const Pt across "+etaLetter);
          linMomImage_Pt.setName(Form("linptres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& logMomImage_Pt = myContent->addImage(logMomCanvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          logMomImage_Pt.setComment("Transverse momentum resolution vs. "+etaLetter+" (log scale) - const Pt across "+etaLetter);
          logMomImage_Pt.setName(Form("ptres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& pImage_Pt = myContent->addImage(pCanvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          pImage_Pt.setComment("Momentum resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          pImage_Pt.setName(Form("pres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& d0Image_Pt = myContent->addImage(d0Canvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          d0Image_Pt.setComment("d0 resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          d0Image_Pt.setName(Form("dxyres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& z0Image_Pt = myContent->addImage(z0Canvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          z0Image_Pt.setComment("z0 resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          z0Image_Pt.setName(Form("dzres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& phiImage_Pt = myContent->addImage(phiCanvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          phiImage_Pt.setComment("Angle resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          phiImage_Pt.setName(Form("phires_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& ctgThetaImage_Pt = myContent->addImage(ctgThetaCanvas_Pt, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          ctgThetaImage_Pt.setComment("Ctg("+thetaLetter+") resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          ctgThetaImage_Pt.setName(Form("cotThetares_%s_%s", tag.c_str(), scenarioStr.c_str()));
        }

        // Draw case II with const P across eta
        if (!gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::RhoGraph_P     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::DGraph_P       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::Z0Graph_P      , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PhiGraph_P     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::CtgthetaGraph_P, tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PGraph_P       , tag).empty()) {
  
          // Set link to myContent
          if (scenario==0) {
            idealMaterial=GraphBag::RealGraph;
            myContent = &resolutionContent_P;
            scenarioStr = "MS_P";
          } else {
            idealMaterial=GraphBag::IdealGraph;
            myContent = &idealResolutionContent_P;
            scenarioStr = "noMS_P";
          }
  
          TCanvas linMomCanvas_P;
          TCanvas logMomCanvas_P;
          TCanvas d0Canvas_P;
          TCanvas phiCanvas_P;
          TCanvas ctgThetaCanvas_P;
          TCanvas z0Canvas_P;
          TCanvas pCanvas_P;
  
          // Default attributes
          int myColor            = 0;
          int nBins              = insur::vis_n_bins;
          int markerStyle        = 21;
          double markerSize      = 1.;
          double lineWidth       = 2.;
          std::string plotOption = "";
  
          linMomCanvas_P.SetGrid(1,1);
          logMomCanvas_P.SetGrid(1,1);
          d0Canvas_P.SetGrid(1,1);
          phiCanvas_P.SetGrid(1,1);
          ctgThetaCanvas_P.SetGrid(1,1);
          z0Canvas_P.SetGrid(1,1);
          pCanvas_P.SetGrid(1,1);
  
          gStyle->SetGridStyle(style_grid);
          gStyle->SetGridColor(color_hard_grid);
  
          // Draw pt
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::RhoGraph_P | idealMaterial, tag)) {
  
            const TGraph& momentumGraph = mapel.second;
            TProfile& momentumProfile   = newProfile(momentumGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              momentumProfile.SetMinimum(insur::vis_min_dPtOverPt); //1E-5*100);
              momentumProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            } else {
              momentumProfile.SetMinimum(insur::vis_min_dPtOverPt); //4E-3*100);
              momentumProfile.SetMaximum(insur::vis_max_dPtOverPt); //.5*100*verticalScale);
            }
            linMomCanvas_P.SetLogy(0);
            logMomCanvas_P.SetLogy(1);
            linMomCanvas_P.SetFillColor(color_plot_background);
            logMomCanvas_P.SetFillColor(color_plot_background);
  
            momentumProfile.SetLineColor(momentumColor(myColor));
            momentumProfile.SetMarkerColor(momentumColor(myColor));
            momentumProfile.SetLineWidth(lineWidth);
            myColor++;
            momentumProfile.SetMarkerStyle(markerStyle);
            momentumProfile.SetMarkerSize(markerSize);
  
            if (momentumGraph.GetN()>0) {
              linMomCanvas_P.cd();
              momentumProfile.Draw(plotOption.c_str());
              logMomCanvas_P.cd();
              momentumProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw p
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::PGraph_P | idealMaterial, tag)) {
  
            const TGraph& pGraph = mapel.second;
            TProfile& pProfile   = newProfile(pGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              pProfile.SetMinimum(insur::vis_min_dPtOverPt); //1E-5*100);
              pProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            } else {
              pProfile.SetMinimum(insur::vis_min_dPtOverPt); //4E-3*100);
              pProfile.SetMaximum(insur::vis_max_dPtOverPt); //.11*100*verticalScale);
            }
            pCanvas_P.SetLogy();
            pCanvas_P.SetFillColor(color_plot_background);
  
            pProfile.SetLineColor(momentumColor(myColor));
            pProfile.SetMarkerColor(momentumColor(myColor));
            pProfile.SetLineWidth(lineWidth);
            myColor++;
            pProfile.SetMarkerStyle(markerStyle);
            pProfile.SetMarkerSize(markerSize);
  
            if (pGraph.GetN() > 0) {
              pCanvas_P.cd();
              pProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }
          // Draw d0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::DGraph_P | idealMaterial, tag)) {
  
            const TGraph& distanceGraph = mapel.second;
            TProfile& distanceProfile   = newProfile(distanceGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              distanceProfile.SetMinimum(vis_min_dD0);
              distanceProfile.SetMaximum(vis_max_dD0);//*verticalScale);
            } else {
                distanceProfile.SetMinimum(vis_min_dD0);
                distanceProfile.SetMaximum(vis_max_dD0);//*verticalScale);
            }
            d0Canvas_P.SetLogy();
            d0Canvas_P.SetFillColor(color_plot_background);
  
            distanceProfile.SetLineColor(momentumColor(myColor));
            distanceProfile.SetMarkerColor(momentumColor(myColor));
            distanceProfile.SetLineWidth(lineWidth);
            myColor++;
            distanceProfile.SetMarkerStyle(markerStyle);
            distanceProfile.SetMarkerSize(markerSize);
  
            if (distanceGraph.GetN()>0) {
              d0Canvas_P.cd();
              distanceProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw phi angle
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::PhiGraph_P | idealMaterial, tag)) {
  
            const TGraph& angleGraph = mapel.second;
            TProfile& angleProfile   = newProfile(angleGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            if (idealMaterial == GraphBag::IdealGraph) {
              angleProfile.SetMinimum(vis_min_dPhi);
              angleProfile.SetMaximum(vis_max_dPhi);//*verticalScale);
            } else {
              angleProfile.SetMinimum(vis_min_dPhi);
              angleProfile.SetMaximum(vis_max_dPhi);//*verticalScale);
            }
            phiCanvas_P.SetLogy();
            phiCanvas_P.SetFillColor(color_plot_background);
  
            angleProfile.SetLineColor(momentumColor(myColor));
            angleProfile.SetMarkerColor(momentumColor(myColor));
            angleProfile.SetLineWidth(lineWidth);
            myColor++;
            angleProfile.SetMarkerStyle(markerStyle);
            angleProfile.SetMarkerSize(markerSize);
  
            if (angleGraph.GetN() > 0) {
              phiCanvas_P.cd();
              angleProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw ctgTheta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::CtgthetaGraph_P | idealMaterial, tag)) {
  
            const TGraph& ctgThetaGraph = mapel.second;
            TProfile& ctgThetaProfile   = newProfile(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);
            ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);//*verticalScale);
            ctgThetaCanvas_P.SetLogy();
            ctgThetaCanvas_P.SetFillColor(color_plot_background);
  
            ctgThetaProfile.SetLineColor(momentumColor(myColor));
            ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
            ctgThetaProfile.SetLineWidth(lineWidth);
            myColor++;
            ctgThetaProfile.SetMarkerStyle(markerStyle);
            ctgThetaProfile.SetMarkerSize(markerSize);
  
            if (ctgThetaGraph.GetN() > 0) {
              ctgThetaCanvas_P.cd();
              ctgThetaProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw z0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::Z0Graph_P | idealMaterial, tag)) {
  
            const TGraph& z0Graph = mapel.second;
            TProfile& z0Profile   = newProfile(z0Graph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
  
            z0Profile.SetMinimum(vis_min_dZ0);
            z0Profile.SetMaximum(vis_max_dZ0);//*verticalScale);
            z0Canvas_P.SetLogy();
            z0Canvas_P.SetFillColor(color_plot_background);
  
            z0Profile.SetLineColor(momentumColor(myColor));
            z0Profile.SetMarkerColor(momentumColor(myColor));
            z0Profile.SetLineWidth(lineWidth);
            myColor++;
            z0Profile.SetMarkerStyle(markerStyle);
            z0Profile.SetMarkerSize(markerSize);
            z0Canvas_P.SetFillColor(color_plot_background);
  
            if (z0Graph.GetN() > 0) {
              z0Canvas_P.cd();
              z0Profile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }
          plotOption = "";
          myColor=0;
            
          RootWImage& linMomImage_P = myContent->addImage(linMomCanvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          linMomImage_P.setComment("Transverse momentum resolution vs. "+etaLetter+" (linear scale) - const P across "+etaLetter);
          linMomImage_P.setName(Form("linptres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& logMomImage_P = myContent->addImage(logMomCanvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          logMomImage_P.setComment("Transverse momentum resolution vs. "+etaLetter+" (log scale) - const P across "+etaLetter);
          logMomImage_P.setName(Form("ptres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& pImage_P = myContent->addImage(pCanvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          pImage_P.setComment("Momentum resolution vs. "+etaLetter+" - const P across "+etaLetter);
          pImage_P.setName(Form("pres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& d0Image_P = myContent->addImage(d0Canvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          d0Image_P.setComment("d0 resolution vs. "+etaLetter+" - const P across "+etaLetter);
          d0Image_P.setName(Form("dxyres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& z0Image_P = myContent->addImage(z0Canvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          z0Image_P.setComment("z0 resolution vs. "+etaLetter+" - const P across "+etaLetter);
          z0Image_P.setName(Form("dzres_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& phiImage_P = myContent->addImage(phiCanvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          phiImage_P.setComment("Angle resolution vs. "+etaLetter+" - const P across "+etaLetter);
          phiImage_P.setName(Form("phires_%s_%s", tag.c_str(), scenarioStr.c_str()));
  
          RootWImage& ctgThetaImage_P = myContent->addImage(ctgThetaCanvas_P, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          ctgThetaImage_P.setComment("Ctg("+thetaLetter+") resolution vs. "+etaLetter+" - const P across "+etaLetter);
          ctgThetaImage_P.setName(Form("cotThetares_%s_%s", tag.c_str(), scenarioStr.c_str()));
        }
      } // Scenarios
 
      // Set Summary content - case I
      //  check that the ideal and real have the same pts
      //  otherwise the table cannot be prepared
      RootWContent& summaryContent_Pt = myPage->addContent("Summary - const Pt across "+etaLetter);
      RootWTable&   cutsSummaryTable  = summaryContent_Pt.addTable();
      RootWTable&   momSummaryTable   = summaryContent_Pt.addTable();
  
      std::map<std::string, RootWTable*>           tableMap_Pt;
      std::map<std::string, RootWTable*>::iterator tableMapIt;
  
      std::vector<std::string> plotNames;
      plotNames.push_back(deltaLetter+"pt/pt [%]:           ");
      plotNames.push_back(deltaLetter+"p/p [%]:             ");
      plotNames.push_back(deltaLetter+"d0 ["+muLetter+"m]:  ");
      plotNames.push_back(deltaLetter+"z0 ["+muLetter+"m]:  ");
      plotNames.push_back(deltaLetter+phiLetter+":          ");
      plotNames.push_back(deltaLetter+"ctg("+thetaLetter+"):");
  
      for (std::vector<std::string>::iterator it=plotNames.begin();
           it!=plotNames.end(); ++it) {
        tableMap_Pt[(*it)] = &(summaryContent_Pt.addTable());
        tableMap_Pt[(*it)]->setContent(0,0,(*it));
      }
        
      // Prepare the cuts for the averages
      ostringstream label;
      std::string name;
      RootWTable* myTable;

      // Table explaining the cuts
      cutsSummaryTable.setContent(0,0,"Region: ");
      cutsSummaryTable.setContent(1,0,"Min "+etaLetter+":");
      cutsSummaryTable.setContent(2,0,"Max "+etaLetter+":");
  
      myTable = &cutsSummaryTable;
      for (unsigned int iBorder=0; iBorder<geom_name_eta_regions.size()-1; ++iBorder) {
        myTable->setContent(0,iBorder+1,geom_name_eta_regions[iBorder+1]);
        label.str(""); label << geom_range_eta_regions[iBorder];
        myTable->setContent(1,iBorder+1,label.str());
        label.str(""); label << geom_range_eta_regions[iBorder+1];
        myTable->setContent(2,iBorder+1,label.str());     
      }

      // Table explaining momenta
      std::vector<double> momentum;
      std::map<graphIndex, TGraph*> myPlotMap_Pt;
          
      fillTaggedPlotMap(gb, plotNames[0], GraphBag::RhoGraph_Pt     , tag, myPlotMap_Pt);
      fillTaggedPlotMap(gb, plotNames[1], GraphBag::PGraph_Pt       , tag, myPlotMap_Pt);
      fillTaggedPlotMap(gb, plotNames[2], GraphBag::DGraph_Pt       , tag, myPlotMap_Pt);
      fillTaggedPlotMap(gb, plotNames[3], GraphBag::Z0Graph_Pt      , tag, myPlotMap_Pt);
      fillTaggedPlotMap(gb, plotNames[4], GraphBag::PhiGraph_Pt     , tag, myPlotMap_Pt);
      fillTaggedPlotMap(gb, plotNames[5], GraphBag::CtgthetaGraph_Pt, tag, myPlotMap_Pt);

                
      std::vector<std::string>::iterator plotNameIt = plotNames.begin();
      std::vector<double>::iterator      momentumIt;
      graphIndex myIndex;
  
      for (std::map<graphIndex, TGraph*>::iterator myPlotMapIt = myPlotMap_Pt.begin();
           myPlotMapIt!=myPlotMap_Pt.end(); ++myPlotMapIt) {
        myIndex =  (*myPlotMapIt).first;
        if (myIndex.name==(*plotNameIt)) {
        momentumIt = std::find(momentum.begin(), momentum.end(), myIndex.p);
        if (momentumIt == momentum.end()) momentum.push_back(myIndex.p);
        }
      }
  
      std::sort(momentum.begin(), momentum.end());
      momSummaryTable.setContent(0,0,"Particle momenta in GeV:  " );
      myTable = &momSummaryTable;
      for (unsigned int iMom=0; iMom<momentum.size(); ++iMom) {
        label.str("");
  
        std::string color = "Unknow";
        if (iMom<insur::color_names.size()) color = color_names[iMom];
  
        if (iMom!=momentum.size()-1) label << momentum[iMom]/Units::GeV << " (" << color << "),";
        else                         label << momentum[iMom]/Units::GeV << " (" << color << ").";
        myTable->setContent(0,iMom+1,label.str());
      }
  
      // Cycle over the different measurements
      for (std::vector<std::string>::iterator plotNameIt = plotNames.begin();
           plotNameIt!=plotNames.end(); ++plotNameIt) {
  
        myTable = tableMap_Pt[*plotNameIt];
        if (!myTable) continue;
  
        // Fill the table with the values, first the heading of momentum
        int baseColumn;
        std::vector<double> averagesReal;
        std::vector<double> averagesIdeal;
        TGraph* myGraph;
        int myColor = kBlack;
        myIndex.name=(*plotNameIt);
        std::ostringstream myLabel;
  
        // Put units & better formatting
        for (unsigned int i=0; i<momentum.size(); ++i) {
          baseColumn = (geom_name_eta_regions.size()-1)*i + 1;
          myTable->setContent(0, baseColumn, momentum[i]/Units::GeV,0);
          myIndex.p=momentum[i];
          myIndex.ideal = false;
          myGraph = myPlotMap_Pt[myIndex];
          myTable->setContent(2, 0, "Real:      ");
          myTable->setContent(3, 0, "Ideal:     ");
          myTable->setContent(4, 0, "Real/Ideal:");
          if (myGraph) {
            averagesReal=Analyzer::average(*myGraph, geom_range_eta_regions);
            myColor = myGraph->GetMarkerColor();
            myTable->setColor(0, baseColumn, myColor);
          }
          myIndex.ideal = true;
          myGraph = myPlotMap_Pt[myIndex];
          if (myGraph) averagesIdeal=Analyzer::average(*myGraph, geom_range_eta_regions);
  
         // Fill resolution for different eta regions
          for (unsigned int j=0; j<(geom_name_eta_regions.size()-1); ++j) {
            myTable->setContent(1, baseColumn+j, geom_name_eta_regions[j+1]);
            myTable->setColor(1, baseColumn+j, myColor);
            if (averagesReal.size() > j) {
  
              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3) myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionStd);
              else                                              myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionHigh);
              myTable->setColor(2, baseColumn+j, myColor);
            }
            if (averagesIdeal.size() > j) {
  
              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3) myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionStd);
              else                                              myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionHigh);
              myTable->setColor(3, baseColumn+j, myColor);
            }
            if ((averagesReal.size() > j)&&(averagesIdeal.size() > j)) {
              myTable->setContent(4, baseColumn+j,averagesReal[j]/averagesIdeal[j],2);
              myTable->setColor(4, baseColumn+j, myColor);
  
            }
          }
        }
      }

      //
      // Set Summary content - case II
      RootWContent& summaryContent_P = myPage->addContent("Summary - const P across "+etaLetter, false);
  
      std::map<std::string, RootWTable*> tableMap_P;
  
      // Table explaining momenta
      std::map<graphIndex, TGraph*> myPlotMap_P;
  
      fillTaggedPlotMap(gb, plotNames[0], GraphBag::RhoGraph_P     , tag, myPlotMap_P);
      fillTaggedPlotMap(gb, plotNames[1], GraphBag::PGraph_P       , tag, myPlotMap_P);
      fillTaggedPlotMap(gb, plotNames[2], GraphBag::DGraph_P       , tag, myPlotMap_P);
      fillTaggedPlotMap(gb, plotNames[3], GraphBag::Z0Graph_P      , tag, myPlotMap_P);
      fillTaggedPlotMap(gb, plotNames[4], GraphBag::PhiGraph_P     , tag, myPlotMap_P);
      fillTaggedPlotMap(gb, plotNames[5], GraphBag::CtgthetaGraph_P, tag, myPlotMap_P);
  
  
      for (std::vector<std::string>::iterator it=plotNames.begin();
           it!=plotNames.end(); ++it) {
        tableMap_P[(*it)] = &(summaryContent_P.addTable());
        tableMap_P[(*it)]->setContent(0,0,(*it));
      }
  
      // Cycle over the different measurements
      for (std::vector<std::string>::iterator plotNameIt = plotNames.begin();
           plotNameIt!=plotNames.end(); ++plotNameIt) {
  
        myTable = tableMap_P[*plotNameIt];
        if (!myTable) continue;
  
        // Fill the table with the values, first the heading of momentum
        int baseColumn;
        std::vector<double> averagesReal;
        std::vector<double> averagesIdeal;
        TGraph* myGraph;
        int myColor = kBlack;
        myIndex.name=(*plotNameIt);
        std::ostringstream myLabel;
  
        // Put units & better formatting
        for (unsigned int i=0; i<momentum.size(); ++i) {
          baseColumn = (geom_name_eta_regions.size()-1)*i + 1;
          myTable->setContent(0, baseColumn, momentum[i]/Units::GeV,0);
          myIndex.p=momentum[i];
          myIndex.ideal = false;
          myGraph = myPlotMap_P[myIndex];
          myTable->setContent(2, 0, "Real:      ");
          myTable->setContent(3, 0, "Ideal:     ");
          myTable->setContent(4, 0, "Real/Ideal:");
          if (myGraph) {
            averagesReal=Analyzer::average(*myGraph, geom_range_eta_regions);
            myColor = myGraph->GetMarkerColor();
            myTable->setColor(0, baseColumn, myColor);
          }
          myIndex.ideal = true;
          myGraph = myPlotMap_P[myIndex];
          if (myGraph) averagesIdeal=Analyzer::average(*myGraph, geom_range_eta_regions);
   
          for (unsigned int j=0; j<(geom_name_eta_regions.size()-1); ++j) {
            myTable->setContent(1, baseColumn+j, geom_name_eta_regions[j+1]);
            myTable->setColor(1, baseColumn+j, myColor);
            if (averagesReal.size() > j) {
  
              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3) myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionStd);
              else                                              myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionHigh);
              myTable->setColor(2, baseColumn+j, myColor);
            }
            if (averagesIdeal.size() > j) {
  
              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3) myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionStd);
              else                                              myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionHigh);
              myTable->setColor(3, baseColumn+j, myColor);
            }
            if ((averagesReal.size() > j)&&(averagesIdeal.size() > j)) {
              myTable->setContent(4, baseColumn+j,averagesReal[j]/averagesIdeal[j],2);
              myTable->setColor(4, baseColumn+j, myColor);
  
            }
          }
        }
      }
    } // For tags
    return true;
  }

  bool Vizard::triggerSummary(Analyzer& a, Tracker& tracker, RootWSite& site, bool extended) {
    //********************************//
    //*                              *//
    //*   Page with the trigger      *//
    //*   summary                    *//
    //*                              *//
    //********************************//
    bool somethingFound = false;

    // Create a page for the errors
    std::string pageTitle = "Trigger";
    std::string pageAddress = "triggerPerf.html";
    RootWPage& myPage = site.addPage(pageTitle);
    myPage.setAddress(pageAddress);  

    // Some helper string objects
    ostringstream tempSS;
    std::string tempString;    

    //********************************//
    //*                              *//
    //*   Eta plot for the trigger   *//
    //*   (Again, with TProfile)     *//
    //*                              *//
    //********************************//

    profileBag aProfileBag = a.getProfileBag();
    std::map<double, TProfile>& triggerProfiles = aProfileBag.getProfiles(profileBag::TriggerProfile | profileBag::TriggeredProfile);
    std::map<double, TProfile>& triggerFractionProfiles = aProfileBag.getProfiles(profileBag::TriggerProfile | profileBag::TriggeredFractionProfile);
    std::map<double, TProfile>& triggerPurityProfiles = aProfileBag.getProfiles(profileBag::TriggerProfile | profileBag::TriggerPurityProfile);

    // Check if profiles exist at all
    if (!triggerProfiles.empty()) {
      somethingFound = true;

      // std::cerr << "found " << triggerProfiles.size() <<" profiles for trigger" << std::endl; // debug

      std::string plotOption;
      int myColor;
      double miny, maxy;

      // Create the contents
      RootWContent& myContent = myPage.addContent("Overall trigger");
      TCanvas pointsCanvas;
      pointsCanvas.SetGrid(1,1);
      plotOption = "E1"; // or "E6"

      // Strings according to the content
      tempString="";
      tempSS.str(""); tempSS << "Number of triggered and triggerable points vs. eta for pT = ";


      // momentum canvas loop
      myColor=0;
      // Style things
      gStyle->SetGridStyle(style_grid);
      gStyle->SetGridColor(color_hard_grid);

      // Loop over the plots and draw on the canvas
      for (std::map<double, TProfile>::iterator plot_iter = triggerProfiles.begin();
           plot_iter != triggerProfiles.end();
           ++plot_iter) {
        const double& myPt = plot_iter->first;

        TProfile& npointsProfile = plot_iter->second;

        miny = npointsProfile.GetBinContent(npointsProfile.GetMinimumBin());
        maxy = npointsProfile.GetBinContent(npointsProfile.GetMaximumBin());
        if ((miny==0)&&(maxy==0)) continue;

        if (myPt!=0) {
          tempSS << tempString.c_str() << myPt; tempString = ", ";
        }

        npointsProfile.SetMinimum(1E-2);
        //npointsProfile.GetXaxis()->SetLimits(0, 2.4);
        npointsProfile.SetLineColor(Palette::color(myColor));
        npointsProfile.SetMarkerColor(Palette::color(myColor));
        npointsProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        npointsProfile.SetMarkerStyle(8);
        pointsCanvas.SetFillColor(color_plot_background);

        pointsCanvas.cd();
        // std::cerr << "About to draw plot " << myPt << std::endl; // debug
        npointsProfile.Draw(plotOption.c_str());
        //plotOption = "E6 same";
        plotOption = "E1 same";
        //plotOption = "same";
      }

      RootWImage& npointsImage = myContent.addImage(pointsCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      npointsImage.setComment(tempSS.str().c_str());
      npointsImage.setName("ntrigpoints");

      // std::cerr << "now to log scale..." << std::endl; // debug

      pointsCanvas.SetLogy();
      RootWImage& npointsLogImage = myContent.addImage(pointsCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS << " (log scale)";
      npointsLogImage.setComment(tempSS.str().c_str());
      npointsLogImage.setName("ntrigpointsLog");

      // std::cerr << "done..." << std::endl; // debug

      TCanvas fractionCanvas;
      fractionCanvas.SetGrid(1,1);
      plotOption = "E1";

      // Strings according to the content
      tempString="";
      tempSS.str(""); tempSS << "Average trigger efficiency vs. eta for pT = ";

      // momentum canvas loop
      myColor=1;
      // Style things
      gStyle->SetGridStyle(style_grid);
      gStyle->SetGridColor(color_hard_grid);

      // Loop over the plots and draw on the canvas
      //miny=1000;
      //maxy=0;
      for (std::map<double, TProfile>::iterator plot_iter = triggerFractionProfiles.begin();
           plot_iter != triggerFractionProfiles.end();
           ++plot_iter) {
        const double& myPt = plot_iter->first;
        TProfile& fractionProfile = plot_iter->second;

        miny = fractionProfile.GetBinContent(fractionProfile.GetMinimumBin());
        maxy = fractionProfile.GetBinContent(fractionProfile.GetMaximumBin());
        //std::cerr << "miny = " << miny << std::endl; // debug
        //std::cerr << "maxy = " << maxy << std::endl; // debug

        if (myPt!=0) {
          tempSS << tempString.c_str() << myPt; tempString = ", ";
        }

        fractionProfile.SetMinimum(1E-2);
        fractionProfile.SetMaximum(100);
        fractionProfile.SetLineColor(Palette::color(myColor));
        fractionProfile.SetMarkerColor(Palette::color(myColor));
        fractionProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        fractionProfile.SetMarkerStyle(8);
        fractionCanvas.SetFillColor(color_plot_background);

        fractionCanvas.cd();
        // std::cerr << "About to draw fraction plot " << myPt << std::endl; // debug
        fractionProfile.Draw(plotOption.c_str());
        plotOption = "E1 same";
        //aValue = fractionProfile.GetBinContent(fractionProfile.GetMaximumBin());
        //if (aValue>maxy) maxy=aValue;
        //aValue = fractionProfile.GetBinContent(fractionProfile.GetMinimumBin());
        //if (aValue<miny) miny=aValue;
        //std::cerr << "Fraction plots between " << miny << " and " << maxy << std::endl;
      }


      RootWImage& fractionImage = myContent.addImage(fractionCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      fractionImage.setComment(tempSS.str().c_str());
      fractionImage.setName("fractiontrigpoints");
      fractionCanvas.SetLogy();
      RootWImage& fractionLogImage = myContent.addImage(fractionCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS << " (log scale)";
      fractionLogImage.setComment(tempSS.str().c_str());
      fractionLogImage.setName("fractiontrigpointsLog");


      TCanvas purityCanvas;
      purityCanvas.SetGrid(1,1);
      plotOption = "E1";

      // Strings according to the content
      tempString="";
      tempSS.str(""); tempSS << "Average stub purity vs. eta for pT > ";

      // momentum canvas loop
      myColor=1;
      // Style things
      gStyle->SetGridStyle(style_grid);
      gStyle->SetGridColor(color_hard_grid);

      // Loop over the plots and draw on the canvas
      //miny=1000;
      //maxy=0;
      for (std::map<double, TProfile>::iterator plot_iter = triggerPurityProfiles.begin();
           plot_iter != triggerPurityProfiles.end();
           ++plot_iter) {
        const double& myPt = plot_iter->first;
        TProfile& purityProfile = plot_iter->second;

        miny = purityProfile.GetBinContent(purityProfile.GetMinimumBin());
        maxy = purityProfile.GetBinContent(purityProfile.GetMaximumBin());
        //std::cerr << "miny = " << miny << std::endl; // debug
        //std::cerr << "maxy = " << maxy << std::endl; // debug

        if (myPt!=0) {
          tempSS << tempString.c_str() << myPt; tempString = ", ";
        }

        purityProfile.SetMinimum(1E-2);
        purityProfile.SetMaximum(100);
        purityProfile.SetLineColor(Palette::color(myColor));
        purityProfile.SetMarkerColor(Palette::color(myColor));
        purityProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        purityProfile.SetMarkerStyle(8);
        purityCanvas.SetFillColor(color_plot_background);

        purityCanvas.cd();
        // std::cerr << "About to draw purity plot " << myPt << std::endl; // debug
        purityProfile.Draw(plotOption.c_str());
        plotOption = "E1 same";
        //aValue = purityProfile.GetBinContent(purityProfile.GetMaximumBin());
        //if (aValue>maxy) maxy=aValue;
        //aValue = purityProfile.GetBinContent(purityProfile.GetMinimumBin());
        //if (aValue<miny) miny=aValue;
        //std::cerr << "Purity plots between " << miny << " and " << maxy << std::endl;
      }


      RootWImage& purityImage = myContent.addImage(purityCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      purityImage.setComment(tempSS.str().c_str());
      purityImage.setName("puritytrigpoints");
      purityCanvas.SetLogy();
      RootWImage& purityLogImage = myContent.addImage(purityCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS << " (log scale)";
      purityLogImage.setComment(tempSS.str().c_str());
      purityLogImage.setName("puritytrigpointsLog");


    } else { // There are no profiles to plot here...
      std::cerr << "ERROR: no trigger performance profile plot to show here" << std::endl;
    }

    //********************************//
    //*                              *//
    //*   Trigger efficiency maps    *//
    //*                              *//
    //********************************//
    mapBag myMapBag = a.getMapBag();
    std::map<double, TH2D>& efficiencyMaps = myMapBag.getMaps(mapBag::efficiencyMap);
    double maxPt = -1;
    // Check if the maps exist at all
    if (!efficiencyMaps.empty()) {
      // std::cerr << "Found " << efficiencyMaps.size() << " efficiency maps"<< std::endl; // debug
      somethingFound = true;
      // Create the content holder for these maps
      RootWContent& myContent = myPage.addContent("Efficiency maps", false);
      for (std::map<double, TH2D>::iterator it = efficiencyMaps.begin();
           it != efficiencyMaps.end(); ++it) {

        // One canvas per map
        TCanvas myCanvas;
        double myPt = it->first;
        if (myPt>maxPt) maxPt=myPt;
        TH2D& myMap = it->second;
        myCanvas.SetFillColor(color_plot_background);
        myCanvas.cd();

        // Actually plot the map
        myMap.SetMinimum(0);
        if (myPt<1.5) { // TODO: make this 1.5 a global constant (also in Analyzer)
          myMap.SetMaximum(0.025);
        } else {
          myMap.SetMaximum(1.0);
        }
        myMap.Draw("colz");

        // Create the image object
        RootWImage& myImage = myContent.addImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        tempSS.str(""); tempSS << "Trigger efficiency map for pT = " << myPt << " GeV/c"; tempString = tempSS.str();
        myImage.setComment(tempString.c_str());
        tempSS.str(""); tempSS << "TriggerEfficiency_" << myPt; tempString = tempSS.str();
        myImage.setName(tempString.c_str());
      }
    } else {
      std::cerr << "ERROR: no trigger efficiency map to show here" << std::endl;
    }

    //********************************//
    //*                              *//
    //*   Stub Efficiency Coverage   *//
    //*                              *//
    //********************************//

/*    std::map<std::string, std::map<std::string, TH1I*>>& profiles = a.getStubEfficiencyCoverageProfiles();
    if (!profiles.empty()) {
      RootWContent& myContent = myPage.addContent("Stub efficiency coverage", false);
      for (const auto& lmel : profiles) {
        TCanvas* myCanvas = new TCanvas(Form("StubEfficiencyCoverageCanvas%s", lmel.first.c_str()), "Stub efficiency eta coverage", vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        //myCanvas.SetFillColor(color_plot_background);
        myCanvas->cd();
        std::vector<std::string> momenta;
        int myColor = 1;
        std::string drawOpts = "";
        for (const auto& mmel : lmel.second) {
          momenta.push_back(mmel.first);
          mmel.second->SetLineColor(Palette::color(myColor));
          mmel.second->SetMarkerColor(Palette::color(myColor));
          //mmel.second->SetFillColor(Palette::color(myColor));
          mmel.second->SetMarkerStyle(1);
          mmel.second->Draw(drawOpts.c_str());
          myCanvas->cd();
          myColor++;
          drawOpts = "same";
          break;
        }
        RootWImage* myImage = new RootWImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        myImage->setComment("Stub efficiency coverage in eta for pT = " + join(momenta, ","));
        myContent.addItem(myImage);
      }
    }*/

    //********************************//
    //*                              *//
    //*   Trigger threshold maps     *//
    //*                              *//
    //********************************//
    std::map<double, TH2D>& thresholdMaps = myMapBag.getMaps(mapBag::thresholdMap);
    // Check if the maps exist at all
    if (!thresholdMaps.empty()) {
      somethingFound = true;
      // std::cerr << "Found " << thresholdMaps.size() << " threshold maps"<< std::endl; // debug

      // Create the content holder for these maps
      RootWContent& myContent = myPage.addContent("Threshold maps", false);
      for (std::map<double, TH2D>::iterator it = thresholdMaps.begin();
           it != thresholdMaps.end(); ++it) {

        // One canvas per map
        TCanvas myCanvas;
        double myEfficiency = it->first;
        TH2D& myMap = it->second;
        myCanvas.SetFillColor(color_plot_background);
        myCanvas.cd();

        // Actually plot the map
        myMap.SetMinimum(0);
        if (maxPt>0) myMap.SetMaximum(maxPt);
        else myMap.SetMaximum(10);
        myMap.Draw("colz");

        // Create the image object
        RootWImage& myImage = myContent.addImage(myCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        tempSS.str(""); tempSS << "Trigger threshold map for eff = " << myEfficiency * 100 << " %";
        tempString = tempSS.str();
        myImage.setComment(tempString.c_str());
        tempSS.str(""); tempSS << "TriggerThreshold_" << myEfficiency; tempString = tempSS.str();
        myImage.setName(tempString.c_str());
      }
    } else {
      std::cerr << "ERROR: no threshold map to show here" << std::endl;
    }


    //********************************//
    //*                              *//
    //*   Configuration maps         *//
    //*                              *//
    //********************************//
    TH2D& suggestedSpacingMap = myMapBag.getMaps(mapBag::suggestedSpacingMap)[mapBag::dummyMomentum];
    TH2D& suggestedSpacingMapAW = myMapBag.getMaps(mapBag::suggestedSpacingMapAW)[mapBag::dummyMomentum];
    TH2D& nominalCutMap = myMapBag.getMaps(mapBag::nominalCutMap)[mapBag::dummyMomentum];


    // Create the content holder for these maps
    RootWContent& myContent = myPage.addContent("Module configuration maps", false);

    // One canvas per map
    TCanvas thickCanvas;
    TCanvas windowCanvas;
    TCanvas suggestedSpacingCanvas;
    TCanvas suggestedSpacingAWCanvas;
    TCanvas nominalCutCanvas;
    thickCanvas.SetFillColor(color_plot_background);
    windowCanvas.SetFillColor(color_plot_background);
    suggestedSpacingCanvas.SetFillColor(color_plot_background);
    suggestedSpacingAWCanvas.SetFillColor(color_plot_background);
    nominalCutCanvas.SetFillColor(color_plot_background);

    struct Spacing { double operator()(const Module& m) { return m.dsDistance(); } };
    PlotDrawer<YZ, Spacing> thicknessDrawer;
    thicknessDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    thicknessDrawer.drawFrame<HistogramFrameStyle>(thickCanvas);
    thicknessDrawer.drawModules<ContourStyle>(thickCanvas);

    struct TriggerWindow { double operator()(const Module& m) { return m.triggerWindow(); } };
    PlotDrawer<YZ, TriggerWindow> windowDrawer;
    windowDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    windowDrawer.drawFrame<HistogramFrameStyle>(windowCanvas);
    windowDrawer.drawModules<ContourStyle>(windowCanvas);


    // Actually plot the maps
    //thickCanvas.cd();
    //thicknessMap.Draw("colz");
    //windowCanvas.cd();
    //windowMap.Draw("colz");
    if (extended) {
      suggestedSpacingCanvas.cd();
      suggestedSpacingMap.Draw("colz");
      suggestedSpacingAWCanvas.cd();
      suggestedSpacingMapAW.Draw("colz");
      nominalCutCanvas.cd();
      //struct PtCut { double operator()(const Module& m) { return PtErrorAdapter(m).getPtCut(); } };
      //PlotDrawer<YZ, PtCut> cutDrawer;
      //cutDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
      //cutDrawer.drawFrame<HistogramFrameStyle>(nominalCutCanvas);
      //cutDrawer.drawModules<ContourStyle>(nominalCutCanvas);
      nominalCutCanvas.SetLogz();
      nominalCutMap.Draw("colz");
    }

    // Create the image objects
    RootWImage& thicknessImage = myContent.addImage(thickCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    thicknessImage.setComment("Map of sensor distances");
    thicknessImage.setName("ThicknessMap");
    RootWImage& windowImage = myContent.addImage(windowCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    windowImage.setComment("Map of selection windows");
    windowImage.setName("WindowMap");
    if (extended) {
      RootWImage& suggestedSpacingImage = myContent.addImage(suggestedSpacingCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      suggestedSpacingImage.setComment("Map of selection suggestedSpacings [default window]");
      suggestedSpacingImage.setName("SuggestedSpacingMap");
      RootWImage& suggestedSpacingAWImage = myContent.addImage(suggestedSpacingAWCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      suggestedSpacingAWImage.setComment("Map of selection suggestedSpacings [selected windows]");
      suggestedSpacingAWImage.setName("SuggestedSpacingMapAW");
      RootWImage& nominalCutImage = myContent.addImage(nominalCutCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      nominalCutImage.setComment("Map of nominal pT cut");
      nominalCutImage.setName("NominalCutMap");
    }


    if (!extended) return somethingFound;

    //********************************//
    //*                              *//
    //* Sensor spacing tuning plots  *//
    //*                              *//
    //********************************//

    std::vector<std::string> profileNames = aProfileBag.getProfileNames(profileBag::TriggerProfileName);

    if (profileNames.size()!=0) {
      somethingFound = true;

      // std::cerr << "Found " << profileNames.size() << " spacing tuning profiles" << std::endl; // debug

      // Create the content
      RootWContent& spacingSummaryContent = myPage.addContent("Sensor spacing tuning summary", true);

      //********************************//
      //*                              *//
      //* Spacing tuning summary       *//
      //*                              *//
      //********************************//

      int nBins = profileNames.size();

      TH1D myFrame("myFrame", "", nBins, 0, nBins);
      myFrame.SetYTitle("Optimal distance range [mm]");
      myFrame.SetMinimum(0);
      myFrame.SetMaximum(6);
      TAxis* xAxis = myFrame.GetXaxis();

      TGraphErrors rangeGraphBad;
      TGraphErrors rangeGraph;
      int rangeGraphPoints;

      std::map<int, TGraphErrors>& spacingTuningGraphs = a.getSpacingTuningGraphs();
      TGraphErrors& availableSpacings = spacingTuningGraphs[-1];


      for (unsigned int i=0; i<profileNames.size(); ++i) {
        double min = a.getTriggerRangeLowLimit(profileNames[i]);
        double max = a.getTriggerRangeHighLimit(profileNames[i]);
        tempString = profileNames[i];
        tempString.substr(profileBag::TriggerProfileName.size(), tempString.size()-profileBag::TriggerProfileName.size());
        xAxis->SetBinLabel(i+1, tempString.c_str());
        if (min<max) {
          rangeGraphPoints=rangeGraph.GetN();
          rangeGraph.SetPoint(rangeGraphPoints, i+0.5, (min+max)/2.);
          rangeGraph.SetPointError(rangeGraphPoints, 0.25, (max-min)/2.);
        } else {
          rangeGraphPoints=rangeGraphBad.GetN();
          rangeGraphBad.SetPoint(rangeGraphPoints, i+0.5, (min+max)/2.);
          rangeGraphBad.SetPointError(rangeGraphPoints, 0.25, (min-max)/2.);
        }
      }


      TCanvas rangeCanvas;
      rangeCanvas.SetFillColor(color_plot_background);
      rangeCanvas.SetGrid(0,1);
      myFrame.Draw();
      rangeGraph.SetFillColor(Palette::color(1));
      rangeGraph.Draw("same 2");
      rangeGraphBad.SetFillColor(Palette::color(2));
      rangeGraphBad.Draw("same 2");
      availableSpacings.SetMarkerStyle(0);
      availableSpacings.Draw("p same");

      RootWImage& RangeImage = spacingSummaryContent.addImage(rangeCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS.str(""); tempSS << "Sensor distance range tuning";
      RangeImage.setComment(tempSS.str());
      tempSS.str(""); tempSS << "TriggerRangeTuning";
      RangeImage.setName(tempSS.str());

      TCanvas tuningCanvas;
      tuningCanvas.SetFillColor(color_plot_background);
      tuningCanvas.SetGrid(0,1);
      //std::map<int, TGraphErrors>& spacingTuningGraphs = a.getSpacingTuningGraphs();
      std::map<int, TGraphErrors>& spacingTuningGraphsBad = a.getSpacingTuningGraphsBad();
      TH1D& spacingTuningFrame = a.getSpacingTuningFrame();
      spacingTuningFrame.Draw();

      for (std::map<int, TGraphErrors>::iterator it = spacingTuningGraphs.begin();
           it!=spacingTuningGraphs.end(); ++it) {
        const int& spacingTuningCounter= it->first;
        TGraphErrors& tuningGraph = it->second;
        if (spacingTuningCounter>0) {
          tuningGraph.SetFillColor(Palette::color(spacingTuningCounter+1));
          tuningGraph.SetFillStyle(1001);
          tuningGraph.Draw("same 2");
        } else {
          tuningGraph.SetMarkerStyle(0);
          tuningGraph.Draw("p same");
        }
      }
      for (std::map<int, TGraphErrors>::iterator it = spacingTuningGraphsBad.begin();
           it!=spacingTuningGraphsBad.end(); ++it) {
        const int& spacingTuningCounter= it->first;
        TGraphErrors& tuningGraph = it->second;
        tuningGraph.SetFillColor(Palette::color(spacingTuningCounter+1));
        tuningGraph.SetFillStyle(3007);
        tuningGraph.Draw("same 2");
      }

      RootWImage& tuningImage = spacingSummaryContent.addImage(tuningCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS.str(""); tempSS << "Sensor distance range tuning for different windows";
      tuningImage.setComment(tempSS.str());
      tempSS.str(""); tempSS << "TriggerRangeTuningWindows";
      tuningImage.setName(tempSS.str());

      TH1D& spacingDistribution = a.getHistoOptimalSpacing(false);
      TCanvas spacingCanvas;
      spacingCanvas.SetFillColor(color_plot_background);
      spacingCanvas.cd();
      spacingDistribution.SetFillColor(Palette::color(1));
      spacingDistribution.Draw();
      RootWImage& spacingImage = spacingSummaryContent.addImage(spacingCanvas, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      spacingImage.setComment("Distribution of minimal spacing for low pT rejection @ standard window");
      spacingImage.setName("SpacingDistribution");
      TH1D& spacingDistributionAW = a.getHistoOptimalSpacing(true);
      TCanvas spacingCanvasAW;
      spacingCanvasAW.SetFillColor(color_plot_background);
      spacingCanvasAW.cd();
      spacingDistributionAW.SetFillColor(Palette::color(1));
      spacingDistributionAW.Draw();
      RootWImage& spacingImageAW = spacingSummaryContent.addImage(spacingCanvasAW, vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      spacingImageAW.setComment("Distribution of minimal spacing for low pT rejection @ selected window");
      spacingImageAW.setName("SpacingDistributionAW");



      //********************************//
      //*                              *//
      //* Spacing tuning details       *//
      //*                              *//
      //********************************//

      // Create the content
      RootWContent& spacingDetailedContent = myPage.addContent("Sensor spacing tuning (detailed)", false);

      for (std::vector<std::string>::const_iterator itName=profileNames.begin(); itName!=profileNames.end(); ++itName) {
        std::map<double, TProfile>& tuningProfiles = aProfileBag.getNamedProfiles(*itName);

        int myColor = 1;
        TCanvas tuningCanvas;
        tuningCanvas.SetFillColor(color_plot_background);
        tuningCanvas.cd();

        std::string plotOption = "E1";
        for (std::map<double, TProfile>::iterator itProfile = tuningProfiles.begin() ; itProfile!= tuningProfiles.end(); ++itProfile) {
          TProfile& tuningProfile = itProfile->second;
          tuningProfile.SetMaximum(100);
          tuningProfile.SetMinimum(0);
          //tuningProfile.SetMaximum(10);
          //tuningProfile.SetMinimum(-10);
          tuningProfile.SetLineColor(Palette::color(myColor));
          tuningProfile.SetFillColor(Palette::color(myColor));
          tuningProfile.SetMarkerColor(Palette::color(myColor));
          myColor++;
          //tuningProfile.SetMarkerStyle(8);
          tuningProfile.Draw(plotOption.c_str());
          plotOption = "same E1";
        }

        RootWImage& tuningImage = spacingDetailedContent.addImage(tuningCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        tempString = (*itName);
        tempString = tempString.substr(profileBag::TriggerProfileName.size(), tempString.size()-profileBag::TriggerProfileName.size());
        tempSS.str(""); tempSS << "Sensor distance tuning for " << tempString.c_str();
        tuningImage.setComment(tempSS.str());
        tempSS.str(""); tempSS << "TriggerTuning" << tempString.c_str();
        tuningImage.setName(tempSS.str());
      }



    }

    //********************************//
    //*                              *//
    //* Turn-on curves plots         *//
    //*                              *//
    //********************************//

    std::vector<std::string> turnonNames = aProfileBag.getProfileNames(profileBag::TurnOnCurveName);

    if (turnonNames.size()!=0) {
      somethingFound = true;

      //********************************//
      //*                              *//
      //* Per-module type turn-on      *//
      //*                              *//
      //********************************//

      // Create the content
      RootWContent& turnOnDetailedContent = myPage.addContent("Modules turnon curves (detailed)", false);

      for (std::vector<std::string>::const_iterator itName=turnonNames.begin(); itName!=turnonNames.end(); ++itName) {
        std::map<double, TProfile>& turnonProfiles = aProfileBag.getNamedProfiles(*itName);

        int myColor = 1;
        TCanvas turnonCanvas;
        turnonCanvas.SetFillColor(color_plot_background);
        turnonCanvas.cd();

        std::string plotOption = "E1";
        tempSS.str("");
        std::string windowsStringSeparator = "";
        for (std::map<double, TProfile>::iterator itProfile = turnonProfiles.begin() ; itProfile!= turnonProfiles.end(); ++itProfile) {
          TProfile& turnonProfile = itProfile->second;
          turnonProfile.SetMaximum(100);
          turnonProfile.SetMinimum(0);
          turnonProfile.SetLineColor(Palette::color(myColor));
          turnonProfile.SetFillColor(Palette::color(myColor));
          turnonProfile.SetMarkerColor(Palette::color(myColor));
          myColor++;
          turnonProfile.Draw(plotOption.c_str());
          plotOption = "same E1";
          tempSS << windowsStringSeparator << std::setprecision(0) << itProfile->first;
          windowsStringSeparator = ", ";
        }
        std::string windowList = tempSS.str();

        tempString = (*itName);
        tempString = tempString.substr(profileBag::TurnOnCurveName.size(), tempString.size()-profileBag::TurnOnCurveName.size());
        RootWImage& turnonImage = turnOnDetailedContent.addImage(turnonCanvas, vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        tempSS.str(""); tempSS << "Sensor turnon curve for " << tempString.c_str() << " with windows of " << windowList;
        turnonImage.setComment(tempSS.str());
        tempSS.str(""); tempSS << "TriggerTurnon" << tempString.c_str();
        turnonImage.setName(tempSS.str());}
    }

    return somethingFound;
  }



  bool Vizard::neighbourGraphSummary(InactiveSurfaces& is, RootWSite& site) {
    std::stringstream ss;
    writeNeighbourGraph(is, ss);

    RootWPage& myPage = site.addPage("Neighbours");
    RootWContent& newContent = myPage.addContent("Neighbour graph", true);
    newContent.addText("<pre>"+ss.str()+"</pre>");

    return true; 
  }


  
  void Vizard::fillTaggedPlotMap(GraphBag& gb,
                                 const string& plotName,
                                 int graphType,
                                 const string& tag,
                                 std::map<graphIndex, TGraph*>& myPlotMap) {
    graphIndex myIndex;
    double p;
    TGraph* myGraph;
    std::vector<std::string> plotNames;

    myIndex.name=plotName;
    //std::cerr << "myIndex.name=" << myIndex.name << std::endl; // debug
    for (int i=0; i<2; ++i) {
      if (i==0) myIndex.ideal=false;
      else myIndex.ideal=true;
      std::map<int, TGraph>& ptGraphsIdeal = gb.getTaggedGraphs((myIndex.ideal ? GraphBag::IdealGraph : GraphBag::RealGraph) | graphType, tag);
      std::map<int, TGraph>::iterator graphsIterator;
      for (graphsIterator=ptGraphsIdeal.begin();
           graphsIterator!=ptGraphsIdeal.end();
           ++graphsIterator) {
        myGraph = &(*graphsIterator).second;
        p = (*graphsIterator).first;
        myIndex.p = p;
        myPlotMap[myIndex] = myGraph;
        //std::cerr << "myIndex.name=" << myIndex.name << std::endl; // debug
      }
    }

  }



  // TODO: describe this here, if it ever worked
  void Vizard::fillPlotMap(std::string& plotName, 
                           std::map<graphIndex, TGraph*>& myPlotMap,
                           Analyzer *a,
                           std::map<int, TGraph>& (Analyzer::*retriveFunction)(bool, bool),
                           bool isTrigger) {
    graphIndex myIndex;
    double p;
    TGraph* myGraph;
    std::vector<std::string> plotNames;

    myIndex.name=plotName;
    //std::cerr << "myIndex.name=" << myIndex.name << std::endl; // debug
    for (int i=0; i<2; ++i) {
      if (i==0) myIndex.ideal=false;
      else myIndex.ideal=true;
      std::map<int, TGraph>& ptGraphsIdeal = (a->*retriveFunction)(myIndex.ideal, isTrigger);
      std::map<int, TGraph>::iterator graphsIterator;
      for (graphsIterator=ptGraphsIdeal.begin();
           graphsIterator!=ptGraphsIdeal.end();
           ++graphsIterator) {
        myGraph = &(*graphsIterator).second;
        p = (*graphsIterator).first;
        myIndex.p = p;
        myPlotMap[myIndex] = myGraph;
        //std::cerr << "myIndex.name=" << myIndex.name << std::endl; // debug
      }
    }

  }


  // public
  // creates a page with all the logs taken from the messagelogger objects
  // @param site a reference to the site we want to work onto
  // @param loggerVector a vector of references to some messageLogger objects
  // @return true if any log was written
  bool Vizard::makeLogPage(RootWSite& site) {
    bool anythingFound=false;
    RootWPage& myPage = site.addPage("Log page");
    if (!MessageLogger::hasEmptyLog(MessageLogger::ERROR))
      myPage.setAlert(1);
    else if (!MessageLogger::hasEmptyLog(MessageLogger::WARNING))
      myPage.setAlert(0.5);
    for (int iLevel=0; iLevel < MessageLogger::NumberOfLevels; ++iLevel) {
      if (!MessageLogger::hasEmptyLog(iLevel)) {
        bool defaultOpen=false;
        if (iLevel<=MessageLogger::WARNING) defaultOpen=true;
        anythingFound=true;
        RootWContent& newContent = myPage.addContent(MessageLogger::getLevelName(iLevel), defaultOpen);
        newContent.addText("<pre>"+MessageLogger::getLatestLog(iLevel)+"</pre>");
        //MessageLogger::getLatestLog(iLevel);
      }
    }
    return anythingFound;
  }



  // private
  // Draws tickmarks on 3d canvases
  // @param myView the TView where to draw ticks
  // @param maxL maximum tracker length in z
  // @param maxR maximum tracker radius in rho
  // @param noAxis number of the axis: Enumerate sections by axis
  //        index normal to draw plane (if x=1, y=2, z=3)
  // @param spacing grid tick spacing
  // @param option the options to pass to the Draw() method
  void Vizard::drawTicks(Analyzer& analyzer, TView* myView, double maxL, double maxR, int noAxis/*=1*/, double spacing /*= 100.*/, Option_t* option /*= "same"*/) {
    TPolyLine3D* aLine;
    Color_t gridColor_hard = color_hard_grid;
    int gridStyle_solid = 1;
    std::string theOption(option);


    double topMax = (maxL > maxR) ? maxL : maxR;
    topMax = ceil(topMax/spacing)*spacing;

    maxL *= 1.1;
    maxR *= 1.1;

    if (noAxis==1) {
      double etaStep=.2;
      double etaMax = analyzer.getEtaMaxGeometry() - etaStep/2.;
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
      double eta;
      for (eta=0; eta<etaMax+etaStep; eta+=etaStep) {
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
      theta = 2 * atan(exp(-analyzer.getEtaMaxGeometry()));
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
      sprintf(labelChar, "%.01f", analyzer.getEtaMaxGeometry());
      aLabel = new TText(pn[0], pn[1], labelChar);
      aLabel->SetTextSize(aLabel->GetTextSize()*.8);
      aLabel->SetTextAlign(21);
      aLabel->Draw("same");
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
  // Creates 4 new canvas with XY and YZ views with all the useful details, like the axis ticks
  // and the eta reference.
  // @param maxZ maximum tracker's Z coordinate to be shown
  // @param maxRho maximum tracker's Rho coordinate to be shown
  // @param analyzer A reference to the analysing class that examined the material budget and filled the histograms
  // @return a pointer to the new TCanvas
  void Vizard::createSummaryCanvas(double maxZ, double maxRho, Analyzer& analyzer,
                                   TCanvas *&YZCanvas, TCanvas *&XYCanvas,
                                   TCanvas *&XYCanvasEC) {
    Int_t irep;
    TVirtualPad* myPad;

    YZCanvas = new TCanvas("YZCanvas", "YZView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY );
    XYCanvas = new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY );
    XYCanvasEC = new TCanvas("XYCanvasEC", "XYView Canvas (Endcap)", vis_min_canvas_sizeX, vis_min_canvas_sizeY );

    // YZView
    if (analyzer.getGeomLiteYZ()) {
      YZCanvas->cd();
      myPad = YZCanvas->GetPad(0);
      drawGrid(maxZ, maxRho, ViewSectionYZ);
      analyzer.getGeomLiteYZ()->DrawClonePad();
      myPad->SetBorderMode(0);
      myPad->SetFillColor(color_plot_background);
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(0, 0, 0, maxZ, maxZ, maxZ);
      myPad->GetView()->SetView(0 /*long*/, 270/*lat*/, 270/*psi*/, irep);
      drawTicks(analyzer, myPad->GetView(), maxZ, maxRho, ViewSectionYZ);
    }

    // XYView (barrel)
    if (analyzer.getGeomLiteXY()) {
      XYCanvas->cd();
      myPad = XYCanvas->GetPad(0);
      drawGrid(maxZ, maxRho, ViewSectionXY);
      analyzer.getGeomLiteXY()->DrawClonePad();
      myPad->SetFillColor(color_plot_background);
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
      myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
    }

    // XYView (EndCap)
    if (analyzer.getGeomLiteEC()) {
      XYCanvasEC->cd();
      myPad = XYCanvasEC->GetPad(0);
      drawGrid(maxZ, maxRho, ViewSectionXY);
      analyzer.getGeomLiteEC()->DrawClonePad();
      myPad->SetFillColor(color_plot_background);
      myPad->GetView()->SetParallel();
      myPad->GetView()->SetRange(-maxRho, -maxRho, -maxRho, maxRho, maxRho, maxRho);
      myPad->GetView()->SetView(0 /*long*/, 0/*lat*/, 270/*psi*/, irep);
    }

    //return summaryCanvas;
  }    

  void Vizard::createSummaryCanvasNicer(Tracker& tracker,
                                        TCanvas *&RZCanvas, TCanvas *&XYCanvas,
                                        TCanvas *&XYCanvasEC) {

    double scaleFactor = tracker.maxR()/600;

    int rzCanvasX = insur::vis_max_canvas_sizeX;//int(tracker.maxZ()/scaleFactor);
    int rzCanvasY = insur::vis_min_canvas_sizeX;//int(tracker.maxR()/scaleFactor);

    RZCanvas = new TCanvas("RZCanvas", "RZView Canvas", rzCanvasX, rzCanvasY );
    RZCanvas->cd();

    PlotDrawer<YZ, Type> yzDrawer;
    yzDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);
    yzDrawer.drawFrame<SummaryFrameStyle>(*RZCanvas);
    yzDrawer.drawModules<ContourStyle>(*RZCanvas);


    XYCanvas = new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY );
    XYCanvas->cd();
    PlotDrawer<XY, Type> xyBarrelDrawer;
    xyBarrelDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL);
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas);
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas);

    XYCanvasEC = new TCanvas("XYCanvasEC", "XYView Canvas (Endcap)", vis_min_canvas_sizeX, vis_min_canvas_sizeY );
    XYCanvasEC->cd();
    PlotDrawer<XY, Type> xyEndcapDrawer; 
    xyEndcapDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), ENDCAP);
    xyEndcapDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasEC);
    xyEndcapDrawer.drawModules<ContourStyle>(*XYCanvasEC);

  }



  /*
   * Returns always the same color for a given momentum index
   * @param iMomentum index of the momentum
   * @return the color index in ROOT
   */
  int Vizard::momentumColor(int iMomentum) {
    if (iMomentum==0) return kBlack;
    if (iMomentum==1) return kBlue;
    if (iMomentum==2) return kRed;   
    if (iMomentum==3) return kGreen;   
    return iMomentum+1;
  }

  /*
   * Modifies a TGraph, so that it looks like a
   * histogram (can be filled)
   * @param myGraph a reference to the TGraph to be modified
   */  
  void Vizard::closeGraph(TGraph& myGraph) {
    double x, y, x0, y0;
    myGraph.GetPoint(myGraph.GetN()-1, x, y);
    myGraph.GetPoint(0, x0, y0);
    myGraph.SetPoint(myGraph.GetN(), x,0);
    myGraph.SetPoint(myGraph.GetN(), x0,0);
  }

  // Helper function to convert a histogram into a TProfile
  TProfile* Vizard::newProfile(TH1D* sourceHistogram) {
    TProfile* resultProfile;
    resultProfile = new TProfile(Form("%s_profile",sourceHistogram->GetName()),
                                 sourceHistogram->GetTitle(),
                                 sourceHistogram->GetNbinsX(),
                                 sourceHistogram->GetXaxis()->GetXmin(),
                                 sourceHistogram->GetXaxis()->GetXmax());
    for (int i=1; i<=sourceHistogram->GetNbinsX(); ++i) {
      resultProfile->Fill(sourceHistogram->GetBinCenter(i), sourceHistogram->GetBinContent(i));
    } 
    resultProfile->SetLineColor(sourceHistogram->GetLineColor());
    resultProfile->SetLineWidth(sourceHistogram->GetLineWidth());
    resultProfile->SetLineStyle(sourceHistogram->GetLineStyle());
    resultProfile->SetFillColor(sourceHistogram->GetFillColor());
    resultProfile->SetFillStyle(sourceHistogram->GetFillStyle());
    return resultProfile;
  }

  TProfile& Vizard::newProfile(const TGraph& sourceGraph, double xlow, double xup, int rebin /* = 1 */, int nBins) {
    TProfile* resultProfile;
    int nPoints = sourceGraph.GetN();
    // Rebin by factor 1 or user defined factor
    if (nBins==0) nPoints /= rebin;
    // Or set new number of bins
    else if (nBins <= nPoints) nPoints = nBins;
    resultProfile = new TProfile(Form("%s_profile", sourceGraph.GetName()), sourceGraph.GetTitle(), nPoints, xlow, xup);
    double x, y;

    for (int i=0; i<sourceGraph.GetN(); ++i) {
      sourceGraph.GetPoint(i, x, y);
      resultProfile->Fill(x, y);
    }

    return (*resultProfile);
  }
  
  void Vizard::createTriggerSectorMapCsv(const TriggerSectorMap& tsm) {
    triggerSectorMapCsv_.clear();
    triggerSectorMapCsv_ = "eta_idx, phi_idx, module_list" + csv_eol; 
    for (TriggerSectorMap::const_iterator it = tsm.begin(); it != tsm.end(); ++it) {
      triggerSectorMapCsv_ += any2str(it->first.first) + csv_separator + any2str(it->first.second);
      for (std::set<int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        triggerSectorMapCsv_ += csv_separator + any2str(*it2);
      }
      triggerSectorMapCsv_ += csv_eol;
    }
  }

  void Vizard::createModuleConnectionsCsv(const ModuleConnectionMap& moduleConnections) {
    std::stringstream ss;
    ss << "cnt, z, rho, phi, sebid, tt_list" << csv_eol;
    for (const auto& mapel : moduleConnections) {
      auto pos = mapel.first->posRef();
      ss << pos.cnt << csv_separator << pos.z << csv_separator << pos.rho << csv_separator << pos.phi << csv_separator << mapel.second.sebCoords;
      for (const auto& conn : mapel.second.connectedProcessors) {
        ss << csv_separator << 't' << conn.first << '_' << conn.second;
      }
      ss << csv_eol;
    }
    moduleConnectionsCsv_ = ss.str();
  }

  std::string Vizard::createAllModulesCsv(const Tracker& t) {
    class TrackerVisitor : public ConstGeometryVisitor {
      std::stringstream output_;
      string sectionName_;
      int layerId_;
    public:
      void preVisit() {
        output_ << "Section/C:Layer/I:Ring/I:r_mm/D:z_mm/D:phi_rad/D:sensorSpacing_mm/D" <<  std::endl;
        output_ << "Section/C, Layer/I, Ring/I, r_mm/D, z_mm/D, phi_rad/D, sensorSpacing_mm/D" << std::endl;
      }
      void visit(const Barrel& b) { sectionName_ = b.myid(); }
      void visit(const Endcap& e) { sectionName_ = e.myid(); }
      void visit(const Layer& l)  { layerId_ = l.myid(); }
      void visit(const Disk& d)  { layerId_ = d.myid(); }
      void visit(const Module& m) {
        output_ << sectionName_ << ", "
          << layerId_ << ", "
          << m.moduleRing() << ", "
          << std::fixed << std::setprecision(6)
          << m.center().Rho() << ", "
          << m.center().Z() << ", "
          << m.center().Phi() << ", "
          << m.dsDistance()
          << std::endl;
      }

      std::string output() const { return output_.str(); }
    };

    TrackerVisitor v;
    v.preVisit();
    t.accept(v);
    return v.output();
  }

  std::string Vizard::createBarrelModulesCsv(const Tracker& t) {
    class BarrelVisitor : public ConstGeometryVisitor {
      std::stringstream output_;
      string barName_;
      int layId_;
      int numRods_;
    public:
      void preVisit() {
        output_ << "Barrel-Layer name, r(mm), z(mm), ss(mm), num mods" << std::endl;
      }
      void visit(const Barrel& b) {
        barName_ = b.myid();
      }
      void visit(const Layer& l) {
        layId_ = l.myid();
        numRods_ = l.numRods();
      }
      void visit(const BarrelModule& m) {
        if (m.posRef().phi > 2) return;
        output_ << barName_ << "-L" << layId_ << ", " << std::fixed << std::setprecision(3) << m.center().Rho() << ", " << m.center().Z() << ", " << m.dsDistance() << ", " << numRods_/2. << std::endl;
      }

      std::string output() const { return output_.str(); }
    };

    BarrelVisitor v;
    v.preVisit();
    t.accept(v);
    return v.output();
  }
  
  std::string Vizard::createEndcapModulesCsv(const Tracker& t) {
    class EndcapVisitor : public ConstGeometryVisitor {
      double minZ_;
    public:
      std::stringstream output;
      void preVisit() {
        output << "Ring, r(mm), phi(deg), z(mm), base_inner(mm), base_outer(mm), height(mm)" <<std::endl;
      }
      void visit(const Endcap& e) {
        minZ_ = e.minZ();
      }
      void visit(const EndcapModule& m) {
        if (m.disk() != 1 || m.minZ() < 0.) return;

        // Print the data in fixed-precision
        // Limit the precision to one micron for lengths and 1/1000 degree for angles
        output << std::fixed;
        output << m.ring() << ", " 
               << std::fixed << std::setprecision(3) << m.center().Rho() << ", "
               << std::fixed << std::setprecision(3) << m.center().Phi()/M_PI*180. << ", "
               << std::fixed << std::setprecision(3) << m.center().Z() - minZ_ << ", "
               << std::fixed << std::setprecision(3) << m.minWidth() << ", "
               << std::fixed << std::setprecision(3) << m.maxWidth() << ", "
               << std::fixed << std::setprecision(3) << m.length() << std::endl;
      }
    };
    EndcapVisitor v;
    v.preVisit(); 
    t.accept(v);
    return v.output.str();
  }

  void Vizard::drawCircle(double radius, bool full, int color/*=kBlack*/) {
    TEllipse* myEllipse = new TEllipse(0,0,radius);
    if (full) {
      myEllipse->SetFillColor(color);
      myEllipse->SetFillStyle(1001);
    } else {
      myEllipse->SetFillStyle(0);
    }
    myEllipse->Draw();
  }

  
  void Vizard::drawInactiveSurfacesSummary(MaterialBudget& materialBudget, RootWPage& myPage) {
    Tracker& myTracker = materialBudget.getTracker();
    std::string myTrackerName = myTracker.myid();

    RootWContent& myContent = myPage.addContent("Service details");

    auto& barrelServices = materialBudget.getInactiveSurfaces().getBarrelServices();
    auto& endcapServices = materialBudget.getInactiveSurfaces().getEndcapServices();
    auto& supports = materialBudget.getInactiveSurfaces().getSupports();

    // We put all services inside the same container
    std::vector<InactiveElement> allServices;
    allServices.reserve( barrelServices.size() + endcapServices.size() + supports.size() ); // preallocate memory
    allServices.insert( allServices.end(), barrelServices.begin(), barrelServices.end() );
    allServices.insert( allServices.end(), endcapServices.begin(), endcapServices.end() );
    allServices.insert( allServices.end(), supports.begin(), supports.end() );

    // Counting services with an ad-hoc index
    int serviceId = 0;
    double z1, z2, r1, r2, length, il, rl;
    double mass;
    std::stringstream myStringStream;

    // Graphic representation of the services in the rz plane
    double maxR = myTracker.maxR()*1.2;
    double maxZ = myTracker.maxZ()*1.2;
    TCanvas* servicesCanvas = new TCanvas("servicesCanvas", "servicesCanvas"); // TODO Factory for canvases?!
    servicesCanvas->cd();
    TH2D* aServicesFrame = new TH2D("aServicesFrame", ";z [mm];r [mm]", 200, -maxZ, maxZ, 100, 0, maxR);
    maxZ=0; maxR=0;
    aServicesFrame->Draw();
    TBox* myBox;
    TText* myText;

    myStringStream << "serviceID/I,elementID/I,z1/D,z2/D,r1/D,r2/D,Element/C,mass/D,mass_per_length/D,rl/D,il/D,local/I" << std::endl;

    for (auto& iter : allServices) {
      z1 = iter.getZOffset();
      z2 = iter.getZOffset()+iter.getZLength();
      r1 = iter.getInnerRadius();
      r2 = iter.getInnerRadius()+iter.getRWidth();
      length = iter.getLength();
      rl = iter.getRadiationLength();
      il = iter.getInteractionLength();

      // Update the maxZ and maxR with respect to the inactive surfaces
      if (fabs(z1)>maxZ) maxZ=fabs(z1);
      if (fabs(z2)>maxZ) maxZ=fabs(z2);
      if (fabs(r1)>maxR) maxR=fabs(r1);
      if (fabs(r2)>maxR) maxR=fabs(r2);

      bool isEmpty = true;

      const std::map<std::string, double>& localMasses = iter.getLocalMasses();

      int elementId=0;
      for (auto& massIt : localMasses) {
	mass = massIt.second;
	if (mass!=0) isEmpty=false;
	myStringStream << serviceId << ","
                       << elementId++ << ","
		       << z1 << ","
		       << z2 << ","
		       << r1 << ","
		       << r2 << ","
		       << massIt.first << ","
		       << mass << ","
		       << mass/length << ","
                       << rl << ","
                       << il << "," 
                       << "1" << std::endl;
      }

      myBox = new TBox(z1, r1, z2, r2);
      myBox->SetLineColor(kBlack);
      myBox->SetFillStyle(3003);
      if (isEmpty) myBox->SetFillColor(kRed);
      else myBox->SetFillColor(kGray);
      myBox->Draw("l");
	
      myText = new TText((z1+z2)/2, (r1+r2)/2, Form("%d", serviceId));
      myText->SetTextAlign(22);
      myText->SetTextSize(2e-2);
      if (isEmpty) myText->SetTextColor(kRed);
      else myText->SetTextColor(kBlack);
      myText->Draw();
    
      serviceId++;
    }
    
    aServicesFrame->GetXaxis()->SetRangeUser(-maxZ, maxZ);
    aServicesFrame->GetYaxis()->SetRangeUser(0, maxR);

    RootWImage& servicesImage = myContent.addImage(servicesCanvas, vis_max_canvas_sizeX, vis_min_canvas_sizeY);
    servicesImage.setComment("Display of the rz positions of the service volumes. Ignoring services with no material.");
    servicesImage.setName("InactiveSurfacesPosition");

    RootWTextFile* myTextFile = new RootWTextFile(Form("inactiveSurfacesMaterials_%s.csv", myTrackerName.c_str()), "file containing all the materials");
    myTextFile->addText(myStringStream.str());
    myContent.addItem(myTextFile);

  }
  //create an extra tab for xml file linking
    bool Vizard::createXmlSite(RootWSite& site,std::string xmldir,std::string layoutdir) {
    RootWPage* myPage = new RootWPage("XML");
    myPage->setAddress("xml.html");
    site.addPage(myPage);

    std::vector<std::string> pixelxmlfilenames,trackerxmlfilenames;
    boost::filesystem::path xmlDirectory(xmldir);
    boost::filesystem::directory_iterator end_iter;
    if ( boost::filesystem::exists(xmlDirectory) && boost::filesystem::is_directory(xmlDirectory)) {
      for( boost::filesystem::directory_iterator dir_iter(xmlDirectory) ; dir_iter != end_iter ; ++dir_iter) {
         if ( boost::filesystem::is_regular_file( dir_iter->path() ) ) {
            // assign current file name to current_file and echo it out to the console.
            std::string current_file =dir_iter->path().filename().string();
            std::cout << current_file << "\tpath=" << dir_iter->path().string() << std::endl;
             if( current_file.find(".xml") != std::string::npos ) {
              boost::filesystem::copy_file( dir_iter->path(),
                                            layoutdir + current_file,
                                            boost::filesystem::copy_option::overwrite_if_exists);
              if( current_file.find("pixel") != std::string::npos )
                pixelxmlfilenames.push_back(current_file);
              else 
                trackerxmlfilenames.push_back(current_file);
             }
        }
      }
    }
    else std::cerr << "XML directory does not exist" << std::endl;

    RootWBinaryFileList* pxFileList = new RootWBinaryFileList(pixelxmlfilenames.begin(), pixelxmlfilenames.end(), 
                                          "xml for Inner Pixel",pixelxmlfilenames.begin(), pixelxmlfilenames.end());
 
    RootWBinaryFileList* tkFileList = new RootWBinaryFileList(trackerxmlfilenames.begin(), trackerxmlfilenames.end(), 
                                          "xml for Outer Tracker",trackerxmlfilenames.begin(), trackerxmlfilenames.end());
    RootWContent* content = new RootWContent("xml files");
    content->addItem(pxFileList);
    content->addItem(tkFileList);
    myPage->addContent(content);
    
  }

}
