/**
 * @file Vizard.cc
 * @brief This class takes care of visualisation for both geometry and analysis results
 */

#include <Vizard.hh>
#include <Units.hh>

#include <TPolyLine.h>
#include <TPaveText.h>
#include <ReportModuleCount.hh>

#include <boost/filesystem.hpp>

#include "MainConfigHandler.hh"


//#include "MaterialTab.hh"

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

    // Each title centered and specified width
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleW(1.);
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
    TGeoVolume* vol=nullptr;
    TGeoTranslation* trans=nullptr;
    TGeoCombiTrans* trafo=nullptr;
    Layer* current=nullptr;
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



  void Vizard::histogramSummary(Analyzer& a, MaterialBudget& materialBudget, RootWSite& site) {
    histogramSummary(a, materialBudget, site, "outer");
  }

  /**
   * This function writes the tables of the weight summaries for the modules
   * with the rootweb library
   * @param a A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   * @param name a qualifier that goes in parenthesis in the title (outer or strip, for example)
   */

  // TODO: if weightGrid is actually unused, then remove it
  void Vizard::weigthSummary(Analyzer& a, MaterialBudget& materialBudget, WeightDistributionGrid& weightGrid, RootWSite& site, std::string name) {

    //RootWContent* myContent;

    // Initialize the page with the material budget
    std::string pageTitle="Weights";
    if (name!="") pageTitle+=" (" +name+")";
    std::string pageAddress="weights"+name+".html";
    RootWPage& myPage = site.addPage(pageTitle);
    myPage.setAddress(pageAddress);


    const WeightsPerSubdetector& weightsPerSubdetector = computeDetailedWeights(materialBudget, myPage);

    // weight plot
    //myContent = new RootWContent("Overview plot", true);
    //myPage.addContent(myContent);

    /*
    std::map<std::string, SummaryTable>* summaryTables;

    // Write the summary for barrel first and endcap second
    for (int i=0; i<2; ++i) {
      if (i==0) summaryTables = &a.getBarrelWeightSummary();
      else summaryTables = &a.getEndcapWeightSummary();

      std::map<std::string, SummaryTable>::iterator it;
      for (it=summaryTables->begin(); it!=summaryTables->end(); ++it) {
        // Create one content per layer
        RootWContent& myContent = myPage.addContent(it->first + " - modules only (chemical elements)", false);
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
        RootWContent& myContent = myPage.addContent(it->first + " - modules only (components)", false);
        RootWTable& myTable = myContent.addTable();
        myTable.setContent(it->second.getContent());
      }
    }
    */


    // TOTAL WEIGHT
    //const WeightsPerSubdetector& weightsPerSubdetectorAndComponent = a.getWeightBySubdetector();

    double totalWeight = 0.;

    for (const auto& subdetectorIt : weightsPerSubdetector) {
      const std::string subdetectorName = subdetectorIt.first;
      RootWContent& myContent = myPage.addContent(subdetectorName, true);
      RootWTable& myTable = myContent.addTable();

      const WeightsPerMechanicalCategory& weightsPerMechanicalCategory = subdetectorIt.second;

      double totalWeightInSubdetector = 0.;
      int rowCounter = 0;
      const bool boldCell = true;
      const int weightPrecision = 1;

      for (const auto& mechanicalCategoryIt : weightsPerMechanicalCategory) {

	const std::string mechanicalCategory = mechanicalCategoryIt.first;

	const WeightsPerComponent& weightsPerComponent = mechanicalCategoryIt.second;

	double totalWeightInMechanicalCategory = 0.;

	myTable.setContent(rowCounter, 0, mechanicalCategory, boldCell);
	myTable.setContent(rowCounter, 1, "mass (kg)", boldCell);
	rowCounter++;
 
	for (const auto& componentIt : weightsPerComponent) {
	  const std::string componentName = componentIt.first;
	  const double mass = componentIt.second;
	  std::ostringstream massStream;
	  const double massInKg = mass / 1000.;  // kg
	  totalWeightInMechanicalCategory += massInKg;

	  myTable.setContent(rowCounter, 0, componentName);
	  myTable.setContent(rowCounter, 1, massInKg, weightPrecision);
	  rowCounter++;
	} // component

	const int mechanicalCategoryTotalColor = kBlue;

	myTable.setContent(rowCounter, 0, "TOTAL " + mechanicalCategory, !boldCell, mechanicalCategoryTotalColor);
	myTable.setContent(rowCounter, 1, totalWeightInMechanicalCategory, weightPrecision, !boldCell, mechanicalCategoryTotalColor);
	rowCounter++;
	myTable.setContent(rowCounter, 0, " ");
	myTable.setContent(rowCounter, 1, " ");
	rowCounter++;
	myTable.setContent(rowCounter, 0, " ");
	myTable.setContent(rowCounter, 1, " ");
	rowCounter++;
	totalWeightInSubdetector += totalWeightInMechanicalCategory;	
      } // mechanical category
      
      myTable.setContent(rowCounter, 0, "TOTAL " + subdetectorName, boldCell);
      myTable.setContent(rowCounter, 1, totalWeightInSubdetector, weightPrecision, boldCell);
      totalWeight += totalWeightInSubdetector;
    } // subdetector

    const std::string grandTotal = "GRAND TOTAL (kg): " + any2str(totalWeight);
    myPage.addContent(grandTotal, false);

    //RootWContent& totalContent = myPage.addContent(grandTotal, false);
    //RootWTable& totalTable = totalContent.addTable();
    //totalTable.setContent(0, 0, totalWeight);
  }



/**
 * This function draws some of the histograms that were filled during material budget analysis
   * with the rootweb library
   * @param a A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   * @param name a qualifier that goes in parenthesis in the title (outer or strip, for example)
   */
  void Vizard::histogramSummary(Analyzer& a, MaterialBudget& materialBudget, RootWSite& site, std::string name) {
    materialBudgets_.push_back(&materialBudget);
    
    // Initialize the page with the material budget
    RootWPage* myPage;
    RootWContent* myContent;
    RootWContent* myContentDetails;
    RootWTable* myTable;
    RootWImage* myImage;
    TVirtualPad* myPad;
    std::string pageTitle="Material";
    if (name!="") pageTitle+=" (" +name+")";
    myPage = new RootWPage(pageTitle);
    std::string pageAddress="material"+name+".html";
    myPage->setAddress(pageAddress);
    site.addPage(myPage);

    std::string name_overviewMaterial = std::string("overviewMaterial") + name ;
    std::string name_categoriesMaterial = std::string("categoriesMaterial") + name ;
    std::string name_overviewMaterialTrackingVolume = std::string("overviewMaterialTrackingVolume") + name ;    
    std::string name_countourMaterial = std::string("countourMaterial") + name ;
    std::string name_mapMaterialRadiation = std::string("mapMaterialRadiation") + name ;
    std::string name_mapMaterialInteraction = std::string("mapMaterialInteraction") + name ;
    std::string name_hadronsHitsNumber = std::string("hadronsHitsNumber") + name ;
    std::string name_hadronsTracksFraction = std::string("hadronsTracksFraction") + name ;
    std::string name_hadTrackRanger = std::string("hadTrackRanger") + name ;


    // Overview
    myContent = new RootWContent("Overview (Full volume)");
    myPage->addContent(myContent);

    ostringstream label;
    std::map<int, std::vector<double> > averages;

    // Book histograms
    TH1D *cr = nullptr, *ci = nullptr;
    TH1D *acr = nullptr, *aci = nullptr, *ser = nullptr, *sei = nullptr, *sur = nullptr, *sui = nullptr;
#ifdef MATERIAL_SHADOW
    TH2D *ir = nullptr, *ii = nullptr;
#endif
    TH2D *mapRad = nullptr, *mapInt = nullptr;
    TProfile *crProf, *ciProf;

    // 1D OVERVIEW (FULL VOLUME)
    std::unique_ptr<TCanvas> myCanvas(new TCanvas(name_overviewMaterial.c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // Full volume rlength
    cr = (TH1D*)a.getHistoGlobalR().Clone();
    crProf = newProfile(cr, 0., a.getEtaMaxMaterial(), materialNBins);
    //crProf->Rebin(10);
    crProf->SetTitle("Radiation Length Over Full Tracker Volume; #eta; x/X_{0}");
    crProf->SetFillColor(kGray + 2);
    crProf->SetLineColor(kBlue);
    crProf->Draw("hist");
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    // Full volume ilength
    ci = (TH1D*)a.getHistoGlobalI().Clone();
    ciProf = newProfile(ci, 0., a.getEtaMaxMaterial(), materialNBins);
    ciProf->SetTitle("Interaction Length Over Full Tracker Volume; #eta; #lambda/#lambda_{0}");
    ciProf->SetFillColor(kGray + 2);
    ciProf->SetLineColor(kBlue);
    ciProf->Draw("hist");

    // Put the full volume materials plots to the site
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in full volume");
    myImage->setName("matOverviewFull");
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

    // Bill of materials ##########3
    RootWTextFile* myTextFile = new RootWTextFile(Form("materials_%s.csv", name.c_str()), "Bill of materials");
    myTextFile->addText(a.getBillOfMaterials());
    myContent->addItem(myTextFile);

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

          // Third row: the photon conversion probability
          averageValue *= -7./9.;
          averageValue = 1 - exp(averageValue);
          materialSummaryTable->setContent(3,j, averageValue ,4);

          // Second row: the interaction length
          averageValue = averageHistogramValues(*ci, geom_range_eta_regions[j-1], geom_range_eta_regions[j]);
          materialSummaryTable->setContent(2,j, averageValue ,4);
        }
      }



    // CATEGORIES DETAILS (FULL VOLUME)
    myContent = new RootWContent("Categories details (Full volume)", false);
    myPage->addContent(myContent);
    // Work area re-init
    myCanvas.reset(new TCanvas(name_categoriesMaterial.c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // radiation length in tracking volume by active, serving or passive
    THStack* rcontainer = new THStack("rstack", "Radiation Length by Category");
    TProfile* surProf = newProfile((TH1D*)a.getHistoSupportsAllR().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    sur = surProf->ProjectionX();
    sur->SetLineColor(kOrange + 4);
    sur->SetFillColor(kOrange + 4);
    rcontainer->Add(sur);
    TProfile* serProf = newProfile((TH1D*)a.getHistoServicesAllR().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    ser = serProf->ProjectionX();
    ser->SetLineColor(kBlue);
    ser->SetFillColor(kBlue);
    rcontainer->Add(ser);
    TProfile* acrProf = newProfile((TH1D*)a.getHistoModulesAllR().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    acr = acrProf->ProjectionX();
    acr->SetLineColor(kRed);
    acr->SetFillColor(kRed);
    rcontainer->Add(acr);
    rcontainer->Draw("hist");
    //rcontainer->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();

    // interaction length in tracking volume by active, serving or passive
    THStack* icontainer = new THStack("istack", "Interaction Length by Category");
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    TProfile* suiProf = newProfile((TH1D*)a.getHistoSupportsAllI().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    sui = suiProf->ProjectionX();
    sui->SetLineColor(kOrange + 2);
    sui->SetFillColor(kOrange + 2);
    icontainer->Add(sui);
    TProfile* seiProf = newProfile((TH1D*)a.getHistoServicesAllI().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    sei = seiProf->ProjectionX();
    sei->SetLineColor(kAzure - 2);
    sei->SetFillColor(kAzure - 2);
    icontainer->Add(sei);
    TProfile* aciProf = newProfile((TH1D*)a.getHistoModulesAllI().Clone(), 0., a.getEtaMaxMaterial(), materialNBins);
    aci = aciProf->ProjectionX();
    aci->SetLineColor(kRed - 3);
    aci->SetFillColor(kRed - 3);
    icontainer->Add(aci);
    icontainer->Draw("hist");
    //icontainer->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();

    // Write asl category plots to web page
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Detailed");
    myImage->setName("matCategoriesFull");
    myTable = new RootWTable();
    // Average values by active, service and passive
    char titleString[256];
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



    // COMPONENTS DETAILS (FULL VOLUME)
    myContent = new RootWContent("Components details (Full volume)", false);
    myPage->addContent(myContent);

    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");


    THStack* rCompStack = new THStack("rcompstack", "Radiation Length by Component");
    THStack* iCompStack = new THStack("icompstack", "Interaction Length by Component");

    TLegend* compLegend = new TLegend(0.1,0.6,0.35,0.9);

    myCanvas.reset(new TCanvas(("componentsRI"+name).c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    std::map<std::string, TH1D*>& rActiveComps = a.getHistoActiveComponentsR();
    int compIndex = 1;
    TProfile* prof;
    TH1D* histo;
    for (std::map<std::string, TH1D*>::iterator it = rActiveComps.begin(); it != rActiveComps.end(); ++it) {
      prof = newProfile((TH1D*)it->second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();   
      histo->SetLineColor(Palette::color(compIndex));
      histo->SetFillColor(Palette::color(compIndex));
      histo->SetTitle(it->first.c_str());
      compLegend->AddEntry(histo, it->first.c_str());
      rCompStack->Add(histo);
      myTable->setContent(compIndex, 0, it->first);
      myTable->setContent(compIndex++, 1, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    rCompStack->Draw("hist");
    //rCompStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    compLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    std::map<std::string, TH1D*>& iActiveComps = a.getHistoActiveComponentsI();
    compIndex = 1;
    for (std::map<std::string, TH1D*>::iterator it = iActiveComps.begin(); it != iActiveComps.end(); ++it) {
      prof = newProfile((TH1D*)it->second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();   
      histo->SetLineColor(Palette::color(compIndex));
      histo->SetFillColor(Palette::color(compIndex));
      histo->SetTitle(it->first.c_str());
      iCompStack->Add(histo);
      myTable->setContent(compIndex++, 2, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    iCompStack->Draw("hist");
    //iCompStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    compLegend->Draw();

    myContent->addItem(myTable);

    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in full volume as a function of &eta; (excluding beam pipe and including services, supports and module material (split by component)");
    myImage->setName("matComponentsFull");
    myContent->addItem(myImage);



    // SERVICES DETAILS (FULL VOLUME)
    myContent = new RootWContent("Services details (Full volume)", false);
    myPage->addContent(myContent);

    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");

    THStack* rServicesCompStack = new THStack("rservicescompstack", "Radiation Length by Component");
    THStack* iServicesCompStack = new THStack("iservicescompstack", "Interaction Length by Component");

    TLegend* servicesCompLegend = new TLegend(0.1,0.6,0.35,0.9);

    myCanvas.reset(new TCanvas(("ServicesComponentsRI"+name).c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    std::map<std::string, TH1D*>& rServicesComps = a.getHistoServicesDetailsR();
    int servicesCompIndex = 1;
    for (std::map<std::string, TH1D*>::iterator it = rServicesComps.begin(); it != rServicesComps.end(); ++it) {
      prof = newProfile((TH1D*)it->second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();     
      histo->SetLineColor(Palette::color(servicesCompIndex));
      histo->SetFillColor(Palette::color(servicesCompIndex));
      histo->SetTitle(it->first.c_str());
      servicesCompLegend->AddEntry(histo, it->first.c_str());
      rServicesCompStack->Add(histo);
      myTable->setContent(servicesCompIndex, 0, it->first);
      myTable->setContent(servicesCompIndex++, 1, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    rServicesCompStack->Draw("hist");  
    //rServicesCompStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    servicesCompLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    std::map<std::string, TH1D*>& iServicesComps = a.getHistoServicesDetailsI();
    servicesCompIndex = 1;
    for (std::map<std::string, TH1D*>::iterator it = iServicesComps.begin(); it != iServicesComps.end(); ++it) {
      prof = newProfile((TH1D*)it->second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(servicesCompIndex));
      histo->SetFillColor(Palette::color(servicesCompIndex));
      histo->SetTitle(it->first.c_str());
      iServicesCompStack->Add(histo);
      myTable->setContent(servicesCompIndex++, 2, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    iServicesCompStack->Draw("hist");
    //rServicesCompStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    servicesCompLegend->Draw();

    myContent->addItem(myTable);
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by component type in services");
    myImage->setName("matServicesComponentsFull");
    myContent->addItem(myImage);



    // COMPONENTS DETAILS (TRACKING VOLUME)
    myContentDetails = new RootWContent("Components details (Tracking volume)", false);

    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");


    THStack* rCompTrackingVolumeStack = new THStack("rcomptrackingvolumestack", "Radiation Length by Component in tracking volume");
    THStack* iCompTrackingVolumeStack = new THStack("icomptrackingvolumestack", "Interaction Length by Component in tracking volume");

    TLegend* compLegendTrackingVolume = new TLegend(0.1,0.6,0.35,0.9);

    myCanvas.reset(new TCanvas(("componentsTrackingVolumeRI"+name).c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    std::map<std::string, TH1D*> rCompsTrackingVolume;
    if (name == "outer") {
      rCompsTrackingVolume = a.getHistoOuterTrackingVolumeR();
      rCompsPixelTrackingVolume_ = a.getHistoPixelTrackingVolumeR();
    }
    else rCompsTrackingVolume = rCompsPixelTrackingVolume_;
    int compIndexTrackingVolume = 1;

    for (const auto& it : rCompsTrackingVolume) {
      prof = newProfile((TH1D*)it.second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(compIndexTrackingVolume));
      histo->SetFillColor(Palette::color(compIndexTrackingVolume));
      histo->SetTitle(it.first.c_str());
      compLegendTrackingVolume->AddEntry(histo, it.first.c_str());
      rCompTrackingVolumeStack->Add(histo);
      myTable->setContent(compIndexTrackingVolume, 0, it.first);
      myTable->setContent(compIndexTrackingVolume++, 1, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    rCompTrackingVolumeStack->Draw("hist");
    compLegendTrackingVolume->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    std::map<std::string, TH1D*> iCompsTrackingVolume;
    if (name == "outer") {
      iCompsTrackingVolume = a.getHistoOuterTrackingVolumeI();
      iCompsPixelTrackingVolume_ = a.getHistoPixelTrackingVolumeI();
    }
    else iCompsTrackingVolume = iCompsPixelTrackingVolume_;
    compIndexTrackingVolume = 1;

    for (const auto& it : iCompsTrackingVolume) {
      prof = newProfile((TH1D*)it.second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(compIndexTrackingVolume));
      histo->SetFillColor(Palette::color(compIndexTrackingVolume));
      histo->SetTitle(it.first.c_str());
      iCompTrackingVolumeStack->Add(histo);
      myTable->setContent(compIndexTrackingVolume++, 2, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    iCompTrackingVolumeStack->Draw("hist");
    compLegendTrackingVolume->Draw();

    myContentDetails->addItem(myTable);

    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by component type in total tracking volume");
    myImage->setName("matComponentsTrackingVolume");
    myContentDetails->addItem(myImage);




    // 1D OVERVIEW (TRACKING VOLUME)
    TH1D *rTrackingVolume = nullptr, *iTrackingVolume = nullptr;
    TProfile *rTrackingVolumeProf = nullptr, *iTrackingVolumeProf = nullptr;
    myContent = new RootWContent("1D Overview (Tracking volume)", false);
    myPage->addContent(myContent);
    // Work area re-init
    myCanvas.reset(new TCanvas(name_overviewMaterialTrackingVolume.c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // global plots in tracking volume: radiation length      
    if (rCompTrackingVolumeStack->GetHists()) {
      rTrackingVolume = (TH1D*)rCompTrackingVolumeStack->GetStack()->Last()->Clone();
      rTrackingVolumeProf = newProfile(rTrackingVolume, 0., a.getEtaMaxMaterial(), materialNBins);
      rTrackingVolumeProf->SetFillColor(kGray + 2);
      rTrackingVolumeProf->SetLineColor(kBlue);
      rTrackingVolumeProf->SetTitle("Radiation Length within Tracking Volume; #eta; x/X_{0}");
      rTrackingVolumeProf->Draw("hist");
    }
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    // global plots in tracking volume: interaction length
    if (iCompTrackingVolumeStack->GetHists()) {
      iTrackingVolume = (TH1D*)iCompTrackingVolumeStack->GetStack()->Last()->Clone();
      iTrackingVolumeProf = newProfile(iTrackingVolume, 0., a.getEtaMaxMaterial(), materialNBins);
      iTrackingVolumeProf->SetFillColor(kGray + 2);
      iTrackingVolumeProf->SetLineColor(kBlue);
      iTrackingVolumeProf->SetTitle("Interaction Length within Tracking Volume; #eta; #lambda/#lambda_{0}");
      iTrackingVolumeProf->Draw("hist");
    }
    // Write global tracking volume plots to web pag
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in tracking volume");
    myImage->setName("matOverviewTrackingVolume");
    myTable = new RootWTable();
    sprintf(titleString, "Average radiation length in tracking volume (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(1, 1, titleString);
    sprintf(titleString, "Average interaction length in tracking volume (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(2, 1, titleString);
    if (rTrackingVolume) myTable->setContent(1, 2, averageHistogramValues(*rTrackingVolume, a.getEtaMaxMaterial()), 5);
    if (iTrackingVolume) myTable->setContent(2, 2, averageHistogramValues(*iTrackingVolume, a.getEtaMaxMaterial()), 5);
    myContent->addItem(myTable);
    myContent->addItem(myImage);

    myPage->addContent(myContentDetails);



    // SERVICES DETAILS (TRACKING VOLUME)
    myContent = new RootWContent("Services details (Tracking volume)", false);
    myPage->addContent(myContent);

    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f])", a.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");

    THStack* rServicesDetailsTrackingVolumeStack = new THStack("rservicescomptrackingvolumestack", "Radiation Length by Component");
    THStack* iServicesDetailsTrackingVolumeStack = new THStack("iservicescomptrackingvolumestack", "Interaction Length by Component");

    TLegend* servicesTrackingVolumeLegend = new TLegend(0.1,0.6,0.35,0.9);

    myCanvas.reset(new TCanvas(("ServicesDetailsTrackingVolumeRI"+name).c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    std::map<std::string, TH1D*> rServicesDetailsTrackingVolume;
    if (name == "outer") {
      rServicesDetailsTrackingVolume = a.getHistoServicesDetailsOuterTrackingVolumeR();
      rServicesDetailsPixelTrackingVolume_ = a.getHistoServicesDetailsPixelTrackingVolumeR();
    }
    else rServicesDetailsTrackingVolume = rServicesDetailsPixelTrackingVolume_;
    int servicesTrackingVolumeIndex = 1;
    for (const auto& it : rServicesDetailsTrackingVolume) {
      prof = newProfile((TH1D*)it.second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();     
      histo->SetLineColor(Palette::color(servicesTrackingVolumeIndex));
      histo->SetFillColor(Palette::color(servicesTrackingVolumeIndex));
      histo->SetTitle(it.first.c_str());
      servicesTrackingVolumeLegend->AddEntry(histo, it.first.c_str());
      rServicesDetailsTrackingVolumeStack->Add(histo);
      myTable->setContent(servicesTrackingVolumeIndex, 0, it.first);
      myTable->setContent(servicesTrackingVolumeIndex++, 1, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    rServicesDetailsTrackingVolumeStack->Draw("hist");  
    //rServicesDetailsTrackingVolumeStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    servicesTrackingVolumeLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    std::map<std::string, TH1D*> iServicesDetailsTrackingVolume;
    if (name == "outer") {
      iServicesDetailsTrackingVolume = a.getHistoServicesDetailsOuterTrackingVolumeI();
      iServicesDetailsPixelTrackingVolume_ = a.getHistoServicesDetailsPixelTrackingVolumeI();
    }
    else iServicesDetailsTrackingVolume = iServicesDetailsPixelTrackingVolume_;

    servicesTrackingVolumeIndex = 1;
    for (const auto& it : iServicesDetailsTrackingVolume) {
      prof = newProfile((TH1D*)it.second, 0., a.getEtaMaxMaterial(), materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(servicesTrackingVolumeIndex));
      histo->SetFillColor(Palette::color(servicesTrackingVolumeIndex));
      histo->SetTitle(it.first.c_str());
      iServicesDetailsTrackingVolumeStack->Add(histo);
      myTable->setContent(servicesTrackingVolumeIndex++, 2, averageHistogramValues(*histo, a.getEtaMaxMaterial()), 5);
    }
    iServicesDetailsTrackingVolumeStack->Draw("hist");
    //rServicesCompStack->GetXaxis()->SetTitle("#eta"); 
    //myCanvas->Modified();
    servicesTrackingVolumeLegend->Draw();

    myContent->addItem(myTable);
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by component type in services");
    myImage->setName("matServicesDetailsTrackingVolume");
    myContent->addItem(myImage);



    // Work area re-init
    myCanvas.reset(new TCanvas(name_countourMaterial.c_str()));
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
    ir->SetContour(temperature_levels, nullptr);
    ir->SetXTitle("z");
    ir->SetYTitle("r");
    ir->Draw("COLZ");
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    // interaction length in isolines
    ii = (TH2D*)a.getHistoIsoI().Clone();
    ii->SetNameTitle("isoi", "Interaction Length Contours");
    ii->SetContour(temperature_levels, nullptr);
    ii->SetXTitle("z");
    ii->SetYTitle("r");
    ii->Draw("COLZ");
    // Write isoline plots to web page
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material 2D distributions");
    myImage->setName("matShadow");
    myContent->addItem(myImage);
#endif // MATERIAL_SHADOW

    // Radiation length plot
    myCanvas.reset(new TCanvas(name_mapMaterialRadiation.c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->cd();
    mapRad = (TH2D*)a.getHistoMapRadiation().Clone();
    mapRad->SetContour(vis_temperature_levels, nullptr);
    mapRad->GetZaxis()->SetLabelSize(0.02);  // palette font size
    //myCanvas->SetLogz();
    mapRad->SetStats(0);
    mapRad->Draw("COLZ");
    myImage = new RootWImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation length material map");
    myImage->setName("matMapR");
    myContent->addItem(myImage);

    // Interaction length plot
    myCanvas.reset(new TCanvas(name_mapMaterialInteraction.c_str()));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->cd();
    mapInt = (TH2D*)a.getHistoMapInteraction().Clone();
    mapInt->SetContour(vis_temperature_levels, nullptr);
    mapInt->GetZaxis()->SetLabelSize(0.02); // palette font size
    mapInt->SetStats(0);
    mapInt->Draw("COLZ");
    myImage = new RootWImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Interaction length material map");
    myImage->setName("matMapI");
    myContent->addItem(myImage);

    // Nuclear interactions
    myContent = new RootWContent("Nuclear interactions", true);
    myPage->addContent(myContent);

    // Number of hits
    myCanvas.reset(new TCanvas(name_hadronsHitsNumber.c_str()));
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
    std::vector<TGraph>& hadronGoodTracksFraction=a.getHadronGoodTracksFraction();
    std::vector<double>& hadronNeededHitsFraction=a.getHadronNeededHitsFraction();
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
      TGraph* myGraph = &hadronGoodTracksFraction.at(i);
      //std::cerr << "Good Hadrons fractions at (" << i <<") has " << myGraph->GetN() << " points" << std::endl;
      //double xx, yy;
      //myGraph->GetPoint(myGraph->GetN()-1, xx, yy);
      //std::cerr << "Last point (x,y) = ("<< xx <<", " << yy <<")" << std::endl;
      averages[i] = Analyzer::average(*myGraph, geom_range_eta_regions);
      closeGraph(*myGraph);
      myGraph->SetFillColor(Palette::color(i+1));
      myGraph->SetBit(1);
      myGraph->Draw("same F");
      tempSS.str("");
      if (hadronNeededHitsFraction.at(i)!=Analyzer::ZeroHitsRequired) {
        if (hadronNeededHitsFraction.at(i)==Analyzer::OneHitRequired)
          tempSS << "1 hit required";
        else
          tempSS << int(hadronNeededHitsFraction.at(i)*100)
            << "%% hits required";
        fractionTitles[i]=tempSS.str();
        myLegend->AddEntry(myGraph, fractionTitles[i].c_str(), "F");
      }
    }
    ranger->Draw("sameaxis");
    myLegend->Draw();
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
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
          j++;
        }
      } else delta--;
    }
    summaryContent.addItem(materialSummaryTable);
    //} else {
    //  delete materialSummaryTable;
    //}

  }

  /**
   * This geometry function computes the transformation matrix that describes the position of a given module
   * in space. It also sets the shape of the provided template volume to correspond to that of the module.
   * @param m A pointer to the module object that needs to be visualised
   * @param v A pointer to the template volume that will represent the module in the visualisation
   * @return A pointer to the finished transformation matrix object
   */
  TGeoCombiTrans* Vizard::modulePlacement(Module* m, TGeoVolume* v) {
    //XYZVector ex, ey, ez, b, c, d, p;
    //TGeoArb8* arb;
    //TGeoRotation* rot;
    TGeoCombiTrans* tr = nullptr;
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


  /*
   * Create all Outer Tracker cabling plots and csv files, and display them on website.
   */
  bool Vizard::outerCablingSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
    bool isPixelTracker = tracker.isPixelTracker();

    if (!isPixelTracker) {
      std::string name = "Outer";

      std::string pageTitle = "Cabling";
      pageTitle += " (" + name + ")";
      RootWPage* myPage = new RootWPage(pageTitle);

      std::string pageAddress ="cabling" + name + ".html";
      myPage->setAddress(pageAddress);

      site.addPage(myPage);
      RootWContent* myContent;

      //********************************//
      //*                              *//
      //*       Plots                  *//
      //*                              *//
      //********************************//
      RootWImage* myImage;

      // Modules to Bundles
      std::unique_ptr<TCanvas> RZBundleCanvas = nullptr;
      std::unique_ptr<TCanvas> XYBundleNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYBundleCanvas = nullptr;   
      std::vector<std::unique_ptr<TCanvas> > XYPosBundlesDisks;
      std::vector<std::unique_ptr<TCanvas> > XYPosBundlesDiskSurfaces;
      std::vector<std::unique_ptr<TCanvas> > XYNegBundlesDisks;
      std::vector<std::unique_ptr<TCanvas> > XYNegBundlesDiskSurfaces;
   
      myContent = new RootWContent("Modules to Bundles");
      myPage->addContent(myContent);

      createOuterCablingPlotsBundles(tracker, RZBundleCanvas, XYBundleCanvas, XYBundleNegCanvas, 
					    XYPosBundlesDisks, XYPosBundlesDiskSurfaces, XYNegBundlesDisks, XYNegBundlesDiskSurfaces);

      if (RZBundleCanvas) {
	const int windowWidth = RZBundleCanvas->GetWindowWidth();
	const int windowHeight = RZBundleCanvas->GetWindowHeight();
	myImage = new RootWImage(std::move(RZBundleCanvas), windowWidth, windowHeight);
	myImage->setComment("(RZ) View : Tracker modules colored by their connections to Bundles. Different colors = different bundles.");
	myContent->addItem(myImage);
      }
      if (XYBundleNegCanvas) {
	myImage = new RootWImage(std::move(XYBundleNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, Negative cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYBundleCanvas) {
	myImage = new RootWImage(std::move(XYBundleCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      // POSITIVE CABLING SIDE
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      RootWTable* positiveSideName = new RootWTable();
      positiveSideName->setContent(0, 0, "Positive cabling side:");
      myContent->addItem(positiveSideName);
      for (auto& XYPosDisk : XYPosBundlesDisks) {
	const std::string comment = XYPosDisk->GetTitle();
	myImage = new RootWImage(std::move(XYPosDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      for (auto& XYPosSurface : XYPosBundlesDiskSurfaces) {
	const std::string comment = XYPosSurface->GetTitle();
	myImage = new RootWImage(std::move(XYPosSurface), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      // NEGATIVE CABLING SIDE
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      RootWTable* negativeSideName = new RootWTable();
      negativeSideName->setContent(0, 0, "Negative cabling side:");
      myContent->addItem(negativeSideName);
      for (auto& XYNegDisk : XYNegBundlesDisks) {
	const std::string comment = XYNegDisk->GetTitle();
	myImage = new RootWImage(std::move(XYNegDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      for (auto& XYNegSurface : XYNegBundlesDiskSurfaces) {
	const std::string comment = XYNegSurface->GetTitle();
	myImage = new RootWImage(std::move(XYNegSurface), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }


      // Modules to DTCs
      std::unique_ptr<TCanvas> RZDTCCanvas = nullptr;
      std::unique_ptr<TCanvas> XYDTCNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYDTCNegFlatCanvas = nullptr;
      std::unique_ptr<TCanvas> XYDTCCanvas = nullptr; 
      std::unique_ptr<TCanvas> XYDTCFlatCanvas = nullptr; 
      std::vector<std::unique_ptr<TCanvas> > XYDTCCanvasesDisk;
       
      myContent = new RootWContent("Modules to DTCs");
      myPage->addContent(myContent);

      createOuterCablingPlotsDTCs(tracker, RZDTCCanvas, XYDTCNegCanvas, XYDTCNegFlatCanvas, XYDTCCanvas, XYDTCFlatCanvas, XYDTCCanvasesDisk);

      if (RZDTCCanvas) {
	const int windowWidth = RZDTCCanvas->GetWindowWidth();
	const int windowHeight = RZDTCCanvas->GetWindowHeight();
	myImage = new RootWImage(std::move(RZDTCCanvas), windowWidth,  windowHeight);
	myImage->setComment("(RZ) View : Tracker modules colored by their connections to DTCs. 1 color = 1 DTC.");
	myContent->addItem(myImage);
      }
      if (XYDTCNegCanvas) {
	myImage = new RootWImage(std::move(XYDTCNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Negative cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYDTCNegFlatCanvas) {
	myImage = new RootWImage(std::move(XYDTCNegFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Negative cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYDTCCanvas) {
	myImage = new RootWImage(std::move(XYDTCCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYDTCFlatCanvas) {
	myImage = new RootWImage(std::move(XYDTCFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      for (auto& XYDTCCanvasDisk : XYDTCCanvasesDisk ) {
	const std::string comment = XYDTCCanvasDisk->GetTitle();
	myImage = new RootWImage(std::move(XYDTCCanvasDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      
      
      const OuterCablingMap* myCablingMap = tracker.getOuterCablingMap();

      // CSV files
      RootWContent* filesContent = new RootWContent("Cabling files", true);
      myPage->addContent(filesContent);   
      RootWTextFile* myTextFile;
      RootWInfo* myInfo = nullptr;
      // POSITIVE CABLING SIDE
      bool isPositiveCablingSide = true;
      filesContent->addItem(positiveSideName);
      // Modules to DTCs
      myTextFile = new RootWTextFile(Form("ModulesToDTCsPos%s.csv", name.c_str()), "Modules to DTCs");
      myTextFile->addText(createModulesToDTCsCsv(tracker, isPositiveCablingSide));
      filesContent->addItem(myTextFile);
      // DTCs to modules
      myTextFile = new RootWTextFile(Form("DTCsToModulesPos%s.csv", name.c_str()), "DTCs to modules");
      myTextFile->addText(createDTCsToModulesCsv(myCablingMap, isPositiveCablingSide));
      filesContent->addItem(myTextFile);
      // Bundles to Modules: Aggregation Patterns in TEDD
      /*This is used for bundle assembly.
	For example, for a given buddle, the pattern 3-4-3-2 means that the bundle is connected to:
	- 3 modules from fanout branch index 1.
	- 4 modules from fanout branch index 2.
	- 3 modules from fanout branch index 3.
	- 2 modules from fanout branch index 4.
      */
      myTextFile = new RootWTextFile(Form("AggregationPatternsPos%s.csv", name.c_str()), "Bundles to Modules: Aggregation Patterns in TEDD");
      myTextFile->addText(createBundlesToEndcapModulesCsv(myCablingMap, isPositiveCablingSide));
      filesContent->addItem(myTextFile);
      // Power cables distribution
      myTextFile = new RootWTextFile(Form("PowerCablesDistributionPos%s.csv", name.c_str()), "Power cables distribution");
      myTextFile->addText(createPowerCablesDistributionCsv(myCablingMap, isPositiveCablingSide));
      filesContent->addItem(myTextFile);

      // NEGATIVE CABLING SIDE
      isPositiveCablingSide = false;
      RootWTable* spacer = new RootWTable();
      spacer->setContent(0, 0, " ");
      spacer->setContent(1, 0, " ");
      spacer->setContent(2, 0, " ");
      filesContent->addItem(spacer);
      filesContent->addItem(negativeSideName);
      // Modules to DTCs
      myTextFile = new RootWTextFile(Form("ModulesToDTCsNeg%s.csv", name.c_str()), "Modules to DTCs");
      myTextFile->addText(createModulesToDTCsCsv(tracker, isPositiveCablingSide));
      filesContent->addItem(myTextFile);
      // DTCs to modules
      myTextFile = new RootWTextFile(Form("DTCsToModulesNeg%s.csv", name.c_str()), "DTCs to modules");
      myTextFile->addText(createDTCsToModulesCsv(myCablingMap, isPositiveCablingSide));
      filesContent->addItem(myTextFile);
      // Bundles to Modules: Aggregation Patterns in TEDD
      myTextFile = new RootWTextFile(Form("AggregationPatternsNeg%s.csv", name.c_str()), "Bundles to Modules: Aggregation Patterns in TEDD");
      myTextFile->addText(createBundlesToEndcapModulesCsv(myCablingMap, isPositiveCablingSide));
      filesContent->addItem(myTextFile);

      // BOTH SIDES, SUMMARY
      filesContent->addItem(spacer);
      RootWTable* bothSidesName = new RootWTable();
      bothSidesName->setContent(0, 0, "Both cabling sides, summary:");
      filesContent->addItem(bothSidesName);
      // CMSSW MODULES DETIDS TO DTC IDS
      myTextFile = new RootWTextFile(Form("CMSSWCablingMap%s.csv", name.c_str()), "CMMSW: Modules DetIds to DTCs Ids");
      myTextFile->addText(createCMSSWOuterTrackerCablingMapCsv(tracker));
      filesContent->addItem(myTextFile);


      // Cabling efficiency
      RootWContent* efficiencyContent = new RootWContent("Cabling efficiency (one side)", true);
      myPage->addContent(efficiencyContent);
      // Links
      myInfo = new RootWInfo("Total number of modules (one side)");
      int numLinks = tracker.modules().size() / 2;
      myInfo->setValue(numLinks);
      efficiencyContent->addItem(myInfo);
      // Bundles
      myInfo = new RootWInfo("Total number of fiber bundles (one side)");
      int numBundles = myCablingMap->getBundles().size();
      myInfo->setValue(numBundles);
      efficiencyContent->addItem(myInfo);
      // Bundles efficiency
      myInfo = new RootWInfo("Fiber bundles efficiency (%)");
      double bundleEfficiency = numLinks / (double)(numBundles * outer_cabling_maxNumModulesPerBundle);
      myInfo->setValue(bundleEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);
      // Cables
      myInfo = new RootWInfo("Total number of fiber cables (one side)");
      int numCables = myCablingMap->getCables().size();
      myInfo->setValue(numCables);
      efficiencyContent->addItem(myInfo);
      // Cables efficiency
      myInfo = new RootWInfo("Fiber cables efficiency (%)");
      double cableEfficiency = numBundles / (double)(numCables * outer_cabling_maxNumBundlesPerCable);
      myInfo->setValue(cableEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);
      // Overall efficiency
      myInfo = new RootWInfo("Overall cabling efficiency (%)");
      double overallEfficiency = cableEfficiency * bundleEfficiency;
      myInfo->setValue(overallEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);


      // Distinct DTCs 2D map
      RootWContent* dtcMapContent = new RootWContent("DTCs per track", false);
      myPage->addContent(dtcMapContent);
      
      std::unique_ptr<TCanvas> hitMapDTCCanvas(new TCanvas("hitmapDTCcanvas", "Hit Map DTC", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
      hitMapDTCCanvas->cd();
      hitMapDTCCanvas->SetFillColor(color_plot_background);
      hitMapDTCCanvas->SetBorderMode(0);
      hitMapDTCCanvas->SetBorderSize(0);
      analyzer.getMapPhiEtaDTC().Draw("colz");
      analyzer.getMapPhiEtaDTC().SetStats(0);
      hitMapDTCCanvas->Modified();
      myImage = new RootWImage(std::move(hitMapDTCCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Number of distinct DTCs per track");
      dtcMapContent->addItem(myImage);


      // Modules to Services Channels (optical)
      std::unique_ptr<TCanvas> XYChannelOpticalNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYChannelOpticalNegFlatCanvas = nullptr;
      std::unique_ptr<TCanvas> XYChannelOpticalCanvas = nullptr; 
      std::unique_ptr<TCanvas> XYChannelOpticalFlatCanvas = nullptr; 
      std::vector<std::unique_ptr<TCanvas> > XYChannelOpticalCanvasesDisk;
       
      myContent = new RootWContent("Modules to Services Channels (optical)");
      myPage->addContent(myContent);

      createOuterCablingPlotsServicesChannelsOptical(tracker, myCablingMap, XYChannelOpticalNegCanvas, XYChannelOpticalNegFlatCanvas, XYChannelOpticalCanvas, XYChannelOpticalFlatCanvas, XYChannelOpticalCanvasesDisk);

      if (XYChannelOpticalNegCanvas) {
	myImage = new RootWImage(std::move(XYChannelOpticalNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Negative cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYChannelOpticalNegFlatCanvas) {
	myImage = new RootWImage(std::move(XYChannelOpticalNegFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Negative cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYChannelOpticalCanvas) {
	myImage = new RootWImage(std::move(XYChannelOpticalCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYChannelOpticalFlatCanvas) {
	myImage = new RootWImage(std::move(XYChannelOpticalFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      for (auto& XYChannelOpticalCanvasDisk : XYChannelOpticalCanvasesDisk ) {
	const std::string comment = XYChannelOpticalCanvasDisk->GetTitle();
	myImage = new RootWImage(std::move(XYChannelOpticalCanvasDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }


      // Modules to Services Channels (powering)
      std::unique_ptr<TCanvas> XYChannelPowerNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYChannelPowerNegFlatCanvas = nullptr;
      std::unique_ptr<TCanvas> XYChannelPowerCanvas = nullptr; 
      std::unique_ptr<TCanvas> XYChannelPowerFlatCanvas = nullptr; 
      std::vector<std::unique_ptr<TCanvas> > XYChannelPowerCanvasesDisk;
      std::vector<std::unique_ptr<TCanvas> > XYNegChannelPowerCanvasesDisk;
       
      myContent = new RootWContent("Modules to Services Channels (powering)");
      myPage->addContent(myContent);

      createOuterCablingPlotsServicesChannelsPower(tracker, myCablingMap, XYChannelPowerNegCanvas, XYChannelPowerNegFlatCanvas, XYChannelPowerCanvas, XYChannelPowerFlatCanvas, XYChannelPowerCanvasesDisk, XYNegChannelPowerCanvasesDisk);

      // POSITIVE CABLING SIDE
      myContent->addItem(positiveSideName);
      if (XYChannelPowerCanvas) {
	myImage = new RootWImage(std::move(XYChannelPowerCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYChannelPowerFlatCanvas) {
	myImage = new RootWImage(std::move(XYChannelPowerFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Positive cabling side. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      for (auto& XYChannelPowerCanvasDisk : XYChannelPowerCanvasesDisk ) {
	const std::string comment = XYChannelPowerCanvasDisk->GetTitle();
	myImage = new RootWImage(std::move(XYChannelPowerCanvasDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
     
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      // NEGATIVE CABLING SIDE
      myContent->addItem(negativeSideName);
      if (XYChannelPowerNegCanvas) {
	myImage = new RootWImage(std::move(XYChannelPowerNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel. Negative cabling side. (CMS +Z points towards the depth of the screen)");
	myContent->addItem(myImage);
      }
      if (XYChannelPowerNegFlatCanvas) {
	myImage = new RootWImage(std::move(XYChannelPowerNegFlatCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : Tracker barrel, untilted modules. Negative cabling side. (CMS +Z points towards the depth of the screen)");
	myContent->addItem(myImage);
      }
      for (auto& XYNegChannelPowerCanvasDisk : XYNegChannelPowerCanvasesDisk ) {
	const std::string comment = XYNegChannelPowerCanvasDisk->GetTitle();
	myImage = new RootWImage(std::move(XYNegChannelPowerCanvasDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }


      // SERVICES CHANNELS TABLES
      RootWContent* channelsContent = new RootWContent("Services per PP1 and channel", false);
      myPage->addContent(channelsContent);

      RootWTable* opticalName = new RootWTable();
      opticalName->setContent(0, 0, "Optical:");
      RootWTable* poweringName = new RootWTable();
      poweringName->setContent(0, 0, "Powering:");

      // POSITIVE CABLING SIDE
      isPositiveCablingSide = true;
      channelsContent->addItem(positiveSideName);
      // SLOT B
      channelsContent->addItem(opticalName);
      ChannelSlot requestedSlot = ChannelSlot::B;
      RootWTable* channelsTablePlusB = opticalServicesChannels(myCablingMap, isPositiveCablingSide, requestedSlot);
      channelsContent->addItem(channelsTablePlusB);
      // SLOTS A AND C
      channelsContent->addItem(spacer);
      channelsContent->addItem(poweringName);
      std::vector<ChannelSlot> slots;
      slots.push_back(ChannelSlot::A);
      slots.push_back(ChannelSlot::C);
      RootWTable* channelsTablePlusAC = powerServicesChannels(myCablingMap, isPositiveCablingSide, slots);
      channelsContent->addItem(channelsTablePlusAC);

      // NEGATIVE CABLING SIDE
      isPositiveCablingSide = false;
      channelsContent->addItem(spacer);
      channelsContent->addItem(spacer);
      channelsContent->addItem(negativeSideName);
      // SLOT B
      channelsContent->addItem(opticalName);
      requestedSlot = ChannelSlot::B;
      RootWTable* channelsTableMinusB = opticalServicesChannels(myCablingMap, isPositiveCablingSide, requestedSlot);
      channelsContent->addItem(channelsTableMinusB);
      // SLOTS A AND C
      channelsContent->addItem(spacer);
      channelsContent->addItem(poweringName);
      RootWTable* channelsTableMinusAC = powerServicesChannels(myCablingMap, isPositiveCablingSide, slots);
      channelsContent->addItem(channelsTableMinusAC);
    }
    return true;
  }


  /*
   * Create all Inner Tracker cabling plots and csv files, and display them on website.
   */
  bool Vizard::innerCablingSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
    bool isPixelTracker = tracker.isPixelTracker();

    if (isPixelTracker) {
      std::string name = "Inner";

      std::string pageTitle = "Cabling";
      pageTitle += " (" + name + ")";
      RootWPage* myPage = new RootWPage(pageTitle);

      std::string pageAddress ="cabling" + name + ".html";
      myPage->setAddress(pageAddress);

      site.addPage(myPage);
      RootWContent* myContent;

      //********************************//
      //*                              *//
      //*       Plots                  *//
      //*                              *//
      //********************************//
      RootWImage* myImage;

      RootWTable* barrelName = new RootWTable();
      barrelName->setContent(0, 0, "BPIX");
      RootWTable* forwardName = new RootWTable();
      forwardName->setContent(0, 0, "FPIX and EPIX, (+Z) End");


      // MODULES TO POWER CHAINS
      std::vector<std::unique_ptr<TCanvas> > ZPhiPowerChainLayerPlots;
      std::unique_ptr<TCanvas> XYPowerChainNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYPowerChainCentralCanvas = nullptr;
      std::unique_ptr<TCanvas> XYPowerChainCanvas = nullptr;
      std::vector<std::unique_ptr<TCanvas> > XYPosPowerChainsDiskSurfaces;
   
      myContent = new RootWContent("Modules to Serial Power Chains");
      myPage->addContent(myContent);  
      
      createInnerCablingPlotsPowerChains(tracker,
					 ZPhiPowerChainLayerPlots,
					 XYPowerChainNegCanvas, XYPowerChainCentralCanvas, XYPowerChainCanvas, 
					 XYPosPowerChainsDiskSurfaces);

      // bpix
      myContent->addItem(barrelName);
      for (auto& ZPhiPlot : ZPhiPowerChainLayerPlots) {
	const std::string comment = ZPhiPlot->GetTitle();
	myImage = new RootWImage(std::move(ZPhiPlot), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      if (XYPowerChainNegCanvas) {
	myImage = new RootWImage(std::move(XYPowerChainNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, (-Z) end. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYPowerChainCentralCanvas) {
	myImage = new RootWImage(std::move(XYPowerChainCentralCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, sensors at Z = 0 only. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYPowerChainCanvas) {
	myImage = new RootWImage(std::move(XYPowerChainCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, (+Z) end. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      // fpix and epix, (+z) end
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      myContent->addItem(forwardName);
      for (auto& XYPosSurface : XYPosPowerChainsDiskSurfaces) {
	const std::string comment = XYPosSurface->GetTitle();
	myImage = new RootWImage(std::move(XYPosSurface), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
 

      // MODULES TO GBTS
      std::vector<std::unique_ptr<TCanvas> > ZPhiGBTLayerPlots;
      std::vector<std::unique_ptr<TCanvas> > XYPosGBTsDiskSurfaces;
   
      myContent = new RootWContent("Modules to LP GBTs");
      myPage->addContent(myContent);
      
      createInnerCablingPlotsGBTs(tracker,
				  ZPhiGBTLayerPlots,
				  XYPosGBTsDiskSurfaces);
      // bpix
      myContent->addItem(barrelName);
      for (auto& ZPhiPlot : ZPhiGBTLayerPlots) {
	const std::string comment = ZPhiPlot->GetTitle();
	myImage = new RootWImage(std::move(ZPhiPlot), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }
      // fpix and epix, (+z) end
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      myContent->addItem(forwardName);
      for (auto& XYPosSurface : XYPosGBTsDiskSurfaces) {
	const std::string comment = XYPosSurface->GetTitle();
	myImage = new RootWImage(std::move(XYPosSurface), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }


      // MODULES TO BUNDLES
      std::unique_ptr<TCanvas> XYBundleNegCanvas = nullptr;
      std::unique_ptr<TCanvas> XYBundlePosCanvas = nullptr;   
      std::vector<std::unique_ptr<TCanvas> > XYPosBundlesDisks;
   
      myContent = new RootWContent("Modules to Fiber Bundles");
      myPage->addContent(myContent);  
      
      createInnerCablingPlotsBundles(tracker,
				     XYBundleNegCanvas, XYBundlePosCanvas, 
				     XYPosBundlesDisks);
      // bpix
      myContent->addItem(barrelName);
      if (XYBundleNegCanvas) {
	myImage = new RootWImage(std::move(XYBundleNegCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, (-Z) end. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      if (XYBundlePosCanvas) {
	myImage = new RootWImage(std::move(XYBundlePosCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, (+Z) end. (CMS +Z points towards you)");
	myContent->addItem(myImage);
      }
      // fpix and epix, (+z) end
      myContent = new RootWContent("");
      myPage->addContent(myContent);
      myContent->addItem(forwardName);
      for (auto& XYPosDisk : XYPosBundlesDisks) {
	const std::string comment = XYPosDisk->GetTitle();
	myImage = new RootWImage(std::move(XYPosDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }


      // MODULES TO DTCs
      std::unique_ptr<TCanvas> RZDTCCanvas = nullptr;
      std::unique_ptr<TCanvas> XYDTCPosCanvas = nullptr;   
      std::vector<std::unique_ptr<TCanvas> > XYPosDTCsDisks;
       
      myContent = new RootWContent("Modules to DTCs");
      myPage->addContent(myContent);

      createInnerCablingPlotsDTCs(tracker, RZDTCCanvas, XYDTCPosCanvas, XYPosDTCsDisks);

      if (RZDTCCanvas) {
	const int windowWidth = RZDTCCanvas->GetWindowWidth();
	const int windowHeight = RZDTCCanvas->GetWindowHeight();
	myImage = new RootWImage(std::move(RZDTCCanvas), windowWidth, windowHeight);
	myImage->setComment("(RZ) View : Inner Tracker modules colored by their connections to DTCs. 1 color <=> 1 DTC.");
	myContent->addItem(myImage);
      }
      if (XYDTCPosCanvas) {
	myImage = new RootWImage(std::move(XYDTCPosCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment("(XY) Section : BPIX, (+Z) end. (CMS +Z points towards you). 1 color <=> 1 DTC.");
	myContent->addItem(myImage);
      }
      for (auto& XYPosDisk : XYPosDTCsDisks) {
	const std::string comment = XYPosDisk->GetTitle();
	myImage = new RootWImage(std::move(XYPosDisk), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
	myImage->setComment(comment);
	myContent->addItem(myImage);
      }

      
      const InnerCablingMap* myInnerCablingMap = tracker.getInnerCablingMap();

      // CSV files
      RootWContent* filesContent = new RootWContent("Cabling files", true);
      myPage->addContent(filesContent);   
      RootWTextFile* myTextFile;
      // Modules to DTCs
      myTextFile = new RootWTextFile(Form("%sTrackerModulesToDTCs.csv", name.c_str()), "Modules to DTCs");
      myTextFile->addText(createInnerTrackerModulesToDTCsCsv(tracker));
      filesContent->addItem(myTextFile);
      // DTCs to modules
      myTextFile = new RootWTextFile(Form("%sTrackerDTCsToModules.csv", name.c_str()), "DTCs to modules");
      myTextFile->addText(createInnerTrackerDTCsToModulesCsv(myInnerCablingMap));
      filesContent->addItem(myTextFile);


      // CABLING COUNT
      int numSensorsOneXSide = 0;
      int numSensorsPlusXSidePlusZEnd = 0;
      int numSensorsPlusXSideMinusZEnd = 0;
      int numPowerChainsOneXSide = 0;
      int numPowerChainsPlusXSidePlusZEnd = 0;
      int numPowerChainsPlusXSideMinusZEnd = 0;
      int numELinksOneXSide = 0;
      int numELinksPlusXSidePlusZEnd = 0;
      int numELinksPlusXSideMinusZEnd = 0;
      int numGBTsOneXSide = 0;
      int numGBTsPlusXSidePlusZEnd = 0;
      int numGBTsPlusXSideMinusZEnd = 0;
      int numBundlesOneXSide = 0;
      int numBundlesPlusXSidePlusZEnd = 0;
      int numBundlesPlusXSideMinusZEnd = 0;
      int numDTCsOneXSide = 0;
      int numDTCsPlusXSidePlusZEnd = 0;
      int numDTCsPlusXSideMinusZEnd = 0;
      computeInnerCablingCount(myInnerCablingMap,
			       numSensorsOneXSide, numSensorsPlusXSidePlusZEnd, numSensorsPlusXSideMinusZEnd,
			       numPowerChainsOneXSide, numPowerChainsPlusXSidePlusZEnd, numPowerChainsPlusXSideMinusZEnd,
			       numELinksOneXSide, numELinksPlusXSidePlusZEnd, numELinksPlusXSideMinusZEnd,
			       numBundlesOneXSide, numBundlesPlusXSidePlusZEnd, numBundlesPlusXSideMinusZEnd,
			       numGBTsOneXSide, numGBTsPlusXSidePlusZEnd, numGBTsPlusXSideMinusZEnd,
			       numDTCsOneXSide, numDTCsPlusXSidePlusZEnd, numDTCsPlusXSideMinusZEnd);

      // CABLING COUNT PER (X) SIDE
      RootWContent* countContent = new RootWContent("Cabling count (one X side, per Z end)", true);
      myPage->addContent(countContent);
      RootWInfo* myInfo = nullptr;
      // (+Z) end, (+X) side
      RootWTable* plusXSidePlusZEndName = new RootWTable();
      plusXSidePlusZEndName->setContent(0, 0, "(+X) side, (+Z) end.");
      countContent->addItem(plusXSidePlusZEndName);
      // Sensors
      myInfo = new RootWInfo("Total number of Sensors");
      myInfo->setValue(numSensorsPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // PowerChains
      myInfo = new RootWInfo("Total number of Serial Power Chains");
      myInfo->setValue(numPowerChainsPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // ELinks
      myInfo = new RootWInfo("Total number of ELinks");
      myInfo->setValue(numELinksPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // GBTs
      myInfo = new RootWInfo("Total number of GBTs");
      myInfo->setValue(numGBTsPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // Bundles
      myInfo = new RootWInfo("Total number of Fiber Bundles");
      myInfo->setValue(numBundlesPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // DTCs
      myInfo = new RootWInfo("Total number of DTCs");
      myInfo->setValue(numDTCsPlusXSidePlusZEnd);
      countContent->addItem(myInfo);
      // (+Z) end, (-X) side
      RootWTable* spacer = new RootWTable();
      spacer->setContent(0, 0, " ");
      spacer->setContent(1, 0, " ");
      spacer->setContent(2, 0, " ");
      countContent->addItem(spacer);
      RootWTable* plusXSideMinusZEndName = new RootWTable();
      plusXSideMinusZEndName->setContent(0, 0, "(+X) side, (-Z) end.");
      countContent->addItem(plusXSideMinusZEndName);
      // Sensors
      myInfo = new RootWInfo("Total number of Sensors");
      myInfo->setValue(numSensorsPlusXSideMinusZEnd);
      countContent->addItem(myInfo);
      // PowerChains
      myInfo = new RootWInfo("Total number of Serial Power Chains");
      myInfo->setValue(numPowerChainsPlusXSideMinusZEnd);
      countContent->addItem(myInfo);
      // ELinks
      myInfo = new RootWInfo("Total number of ELinks");
      myInfo->setValue(numELinksPlusXSideMinusZEnd);
      countContent->addItem(myInfo);
      // GBTs
      myInfo = new RootWInfo("Total number of GBTs");
      myInfo->setValue(numGBTsPlusXSideMinusZEnd);
      countContent->addItem(myInfo);
      // Bundles
      myInfo = new RootWInfo("Total number of Fiber Bundles");
      myInfo->setValue(numBundlesPlusXSideMinusZEnd);
      countContent->addItem(myInfo);
      // DTCs
      myInfo = new RootWInfo("Total number of DTCs");
      myInfo->setValue(numDTCsPlusXSideMinusZEnd);
      countContent->addItem(myInfo);

      // CABLING EFFICIENCY
      RootWContent* efficiencyContent = new RootWContent("Cabling efficiency (one X side)", true);
      myPage->addContent(efficiencyContent);
      // Sensors
      myInfo = new RootWInfo("Total number of Sensors (one X side)");
      myInfo->setValue(numSensorsOneXSide);
      efficiencyContent->addItem(myInfo);
      // PowerChains
      myInfo = new RootWInfo("Total number of Serial Power Chains (one X side)");
      myInfo->setValue(numPowerChainsOneXSide);
      efficiencyContent->addItem(spacer);
      efficiencyContent->addItem(myInfo);
      // PowerChains efficiency
      myInfo = new RootWInfo("Sensors <-> Serial Power Chains efficiency (%)");
      const double powerChainEfficiency = (double)numSensorsOneXSide / (numPowerChainsOneXSide * inner_cabling_maxNumModulesPerPowerChain);
      myInfo->setValue(powerChainEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);
      // ELinks
      myInfo = new RootWInfo("Total number of ELinks (one X side)");
      myInfo->setValue(numELinksOneXSide);
      efficiencyContent->addItem(spacer);
      efficiencyContent->addItem(myInfo);
      // GBTs
      myInfo = new RootWInfo("Total number of LP GBTs (one X side)");
      myInfo->setValue(numGBTsOneXSide);
      efficiencyContent->addItem(spacer);
      efficiencyContent->addItem(myInfo);
      // GBTs efficiency
      myInfo = new RootWInfo("ELinks <-> LP GBTs efficiency (%)");
      const double GBTEfficiency = (double)numELinksOneXSide / (numGBTsOneXSide * inner_cabling_maxNumELinksPerGBT);
      myInfo->setValue(GBTEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);
      // Bundles
      myInfo = new RootWInfo("Total number of Fiber Bundles (one X side)");
      myInfo->setValue(numBundlesOneXSide);
      efficiencyContent->addItem(spacer);
      efficiencyContent->addItem(myInfo);
      // Bundles efficiency
      myInfo = new RootWInfo("LP GBTs <-> Fiber Bundles efficiency (%)");
      const double bundleEfficiency = (double)numGBTsOneXSide / (numBundlesOneXSide * inner_cabling_maxNumGBTsPerBundle);
      myInfo->setValue(bundleEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);
      // DTCs
      myInfo = new RootWInfo("Total number of DTCs (one X side)");
      myInfo->setValue(numDTCsOneXSide);
      efficiencyContent->addItem(spacer);
      efficiencyContent->addItem(myInfo);
      // DTCs efficiency
      myInfo = new RootWInfo("Fiber Bundles <-> DTCs efficiency (%)");
      const double dtcEfficiency = (double)numBundlesOneXSide / (numDTCsOneXSide * inner_cabling_maxNumBundlesPerCable);
      myInfo->setValue(dtcEfficiency * 100, 0);
      efficiencyContent->addItem(myInfo);

    } // end of isPixelTracker
    return true;
  }


  /**
   * This function draws the profile of hits obtained by the analysis of the geometry
   * together with the summaries in tables with the rootweb library. It also actually does a couple of
   * calculations to count modules and such, to put the results in the tables.
   * @param analyzer A reference to the analysing class that examined the material budget and filled the histograms
   * @param site the RootWSite object for the output
   */
  bool Vizard::geometrySummary(Analyzer& analyzer, Tracker& tracker, InactiveSurfaces* inactive, RootWSite& site, bool& debugResolution, std::string name) {
    trackers_.push_back(&tracker);

    std::map<std::string, double>& tagMapWeight = analyzer.getTagWeigth();


    std::string pageTitle = "Geometry";
    if (name!="") pageTitle+=" (" +name+")";
    const bool isPixelTracker = (name == "pixel");

    RootWPage* myPage = new RootWPage(pageTitle);

    std::string pageAddress="layout"+name+".html";
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

    LayerDiskSummaryVisitor v;
    v.preVisit();
    tracker.accept(v);
    v.postVisit();

    // TO DO : replace this with a better spacer
    RootWTable* spacer = new RootWTable();
    spacer->setContent(0, 0, " ");
    spacer->setContent(1, 0, " ");
    spacer->setContent(2, 0, " ");
    spacer->setContent(3, 0, " ");

    myContent->addItem(v.layerTable);
    myContent->addItem(spacer);
    myContent->addItem(v.diskTable);


    //***************************************//
    //*                                     *//
    //*              Barrel :               *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//
    myContent = new RootWContent("Barrel : additional info", false);

    SkewedLayersVisitor sv;
    tracker.accept(sv);

    if (sv.numSkewedLayers > 0) {
      myPage->addContent(myContent);

      for (int i = 0; i < sv.numSkewedLayers; i++) {
	// layer
	myContent->addItem(sv.tables.at(i));
	// spacer
	if (i < sv.numSkewedLayers - 1) { myContent->addItem(spacer); }
      }
    }


    //***************************************//
    //*                                     *//
    //* Automatic-placement tilted layers : *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//
    myContent = new RootWContent("Tilted layers with automatic placement : additional info", false);

    TiltedLayersVisitor tv;
    tracker.accept(tv);

    // If the layer is tilted, print flat and tilted parts tables.
    if (tv.numTiltedLayers > 0) {
      myPage->addContent(myContent);

      for (int i = 0; i < tv.numTiltedLayers; i++) {
	// layer name
	myContent->addItem(tv.tiltedLayerNames.at(i));
	// flat part
	myContent->addItem(tv.flatPartNames.at(i));
	myContent->addItem(tv.flatPartTables.at(i));
	// tilted part
	myContent->addItem(tv.tiltedPartNames.at(i));
	myContent->addItem(tv.tiltedPartTables.at(i));
	// spacer
	if (i < tv.numTiltedLayers - 1) { myContent->addItem(spacer); }
      }
    }


    //***************************************//
    //*                                     *//
    //*              Endcaps :              *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//
    myContent = new RootWContent("Endcaps : additional info", false);
    myPage->addContent(myContent);
    for (int i = 0; i < v.nEndcaps; i++) {
      if (i > 0) myContent->addItem(spacer);
      myContent->addItem(v.endcapNames.at(i));
      myContent->addItem(v.endcapTables.at(i));
    }
    for (int i = 0; i < v.nDisks; i++) {
      myContent->addItem(spacer);
      myContent->addItem(v.diskNames.at(i));
      myContent->addItem(v.zErrorTables.at(i));
    }


    //********************************//
    //*                              *//
    //*       Modules                *//
    //*                              *//
    //********************************//
    double totalPower=0;
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
    std::ostringstream aWeight;
    int barrelCount=0;
    int endcapCount=0;

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
    static const int moduleWeightRow = 26;
    static const int inactiveWeightRow = 27;
    static const int totalWeightRow = 28;

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
    moduleTable->setContent(moduleWeightRow, 0, "Weight (av, g)");
    moduleTable->setContent(inactiveWeightRow, 0, "Service Weight");
    moduleTable->setContent(totalWeightRow, 0, "Total Weight");

    int loPitch;
    int hiPitch;


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
      moduleTable->setContent(moduleWeightRow, iType, aWeight.str());

      moduleTable->setContent(thicknessRow, iType, aThickness.str());
      moduleTable->setContent(moduleAreaRow, iType, aModuleArea.str());
      moduleTable->setContent(totalAreaRow, iType, aTotalArea.str());
      moduleTable->setContent(channelRow, iType, aChannel.str());
      // moduleTable->setContent(areaRow, iType, anArea.str());

    }

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
    aPower << std::fixed << std::setprecision(totalPowerPrecision) << totalPower * 1e-6;
    if (v.totalSensorPower > 1e-6) { // non-zero check for double
      aSensorPower << std::fixed << std::setprecision(totalPowerPrecision) << v.totalSensorPower * 1e-3;
    } else {
      aSensorPower << "n/a";
    }
    moduleTable->setContent(powerRow, iType, aPower.str());
    moduleTable->setContent(powerPerModuleRow, iType, aPowerPerModule.str());
    moduleTable->setContent(sensorPowerRow, iType, aSensorPower.str());
    moduleTable->setContent(sensorPowerPerModuleAvgRow, iType, aSensorPowerPerModuleAvg.str());
    moduleTable->setContent(sensorPowerPerModuleMaxRow, iType, aSensorPowerPerModuleMax.str());
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
    //*   Module types breakdown     *//
    //*                              *//
    //********************************//

    ReportModuleCount rmc;
    rmc.analyze(tracker);
    rmc.visualizeTo(*myContent);

    //********************************//
    //*                              *//
    //*    Modules parametrized      *//
    //*     spatial resolution       *//
    //*                              *//
    //********************************//
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

      // Modules' incident angles distributions (view from modules)
      std::map<std::string, TH1D>& incidentAngleLocalXBarrelDistribution = analyzer.getIncidentAngleLocalXBarrelDistribution();
      std::map<std::string, TH1D>& incidentAngleLocalYBarrelDistribution = analyzer.getIncidentAngleLocalYBarrelDistribution();
      std::map<std::string, TH1D>& incidentAngleLocalXEndcapsDistribution = analyzer.getIncidentAngleLocalXEndcapsDistribution();
      std::map<std::string, TH1D>& incidentAngleLocalYEndcapsDistribution = analyzer.getIncidentAngleLocalYEndcapsDistribution();

      // Tracks angles distributions (global coordinates)
      std::map<std::string, TH1D>& trackPhiBarrelDistribution = analyzer.getTrackPhiBarrelDistribution();
      std::map<std::string, TH1D>& trackEtaBarrelDistribution = analyzer.getTrackEtaBarrelDistribution();
      std::map<std::string, TH1D>& trackPhiEndcapsDistribution = analyzer.getTrackPhiEndcapsDistribution();
      std::map<std::string, TH1D>& trackEtaEndcapsDistribution = analyzer.getTrackEtaEndcapsDistribution();

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

	gStyle->SetOptStat("emr");
	if (parametrizedResolutionLocalXBarrelMap[tag].GetEntries() != 0) {
	  std::unique_ptr<TCanvas> resoXBarCanvas(new TCanvas());
	  resoXBarCanvas->SetFillColor(color_plot_background);
	  resoXBarCanvas->Divide(2,2);
	  TVirtualPad* myPad;
	  myPad = resoXBarCanvas->GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoXBarCanvas->GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalXBarrelMap[tag].Draw();
	  myPad = resoXBarCanvas->GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalXBarrelDistribution[tag].SetStats(1);
	  const double normA = 1. / parametrizedResolutionLocalXBarrelDistribution[tag].Integral();
	  parametrizedResolutionLocalXBarrelDistribution[tag].Scale(normA, "width");
	  parametrizedResolutionLocalXBarrelDistribution[tag].Draw();
	  myPad = resoXBarCanvas->GetPad(3);
	  myPad->cd();
	  incidentAngleLocalXBarrelDistribution[tag].SetStats(1);
	  const double normB = 1. / incidentAngleLocalXBarrelDistribution[tag].Integral();
	  incidentAngleLocalXBarrelDistribution[tag].Scale(normB, "width");
	  incidentAngleLocalXBarrelDistribution[tag].Draw();
	  myPad = resoXBarCanvas->GetPad(4);
	  myPad->cd();
	  trackPhiBarrelDistribution[tag].SetStats(1);
	  const double normC = 1. / trackPhiBarrelDistribution[tag].Integral();
	  trackPhiBarrelDistribution[tag].Scale(normC, "width");
	  trackPhiBarrelDistribution[tag].Draw();
	  RootWImage& resoXBarImage = parametrizedResolutionContent.addImage(std::move(resoXBarCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoXBarImage.setComment(Form("Resolution on local X coordinate for %s barrel modules", tag.c_str()));
	  resoXBarImage.setName(Form("Resolution on local X coordinate for %s barrel modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalYBarrelMap[tag].GetEntries() != 0) {
	  std::unique_ptr<TCanvas> resoYBarCanvas(new TCanvas());
	  resoYBarCanvas->SetFillColor(color_plot_background);
	  resoYBarCanvas->Divide(2,2);
	  TVirtualPad* myPad;
	  myPad = resoYBarCanvas->GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoYBarCanvas->GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalYBarrelMap[tag].Draw();
	  myPad = resoYBarCanvas->GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalYBarrelDistribution[tag].SetStats(1);
	  const double normA = 1. / parametrizedResolutionLocalYBarrelDistribution[tag].Integral();
	  parametrizedResolutionLocalYBarrelDistribution[tag].Scale(normA, "width");
	  parametrizedResolutionLocalYBarrelDistribution[tag].Draw();
	  myPad = resoYBarCanvas->GetPad(3);
	  myPad->cd();
	  incidentAngleLocalYBarrelDistribution[tag].SetStats(1);
	  const double normB = 1. / incidentAngleLocalYBarrelDistribution[tag].Integral();
	  incidentAngleLocalYBarrelDistribution[tag].Scale(normB, "width");
	  incidentAngleLocalYBarrelDistribution[tag].Draw();
	  myPad = resoYBarCanvas->GetPad(4);
	  myPad->cd();
	  trackEtaBarrelDistribution[tag].SetStats(1);
	  const double normC = 1. / trackEtaBarrelDistribution[tag].Integral();
	  trackEtaBarrelDistribution[tag].Scale(normC, "width");
	  trackEtaBarrelDistribution[tag].Draw();
	  RootWImage& resoYBarImage = parametrizedResolutionContent.addImage(std::move(resoYBarCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoYBarImage.setComment(Form("Resolution on local Y coordinate for %s barrel modules", tag.c_str()));
	  resoYBarImage.setName(Form("Resolution on local Y coordinate for %s barrel modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalXEndcapsMap[tag].GetEntries() != 0) {
	  std::unique_ptr<TCanvas> resoXEndCanvas(new TCanvas());
	  resoXEndCanvas->SetFillColor(color_plot_background);
	  resoXEndCanvas->Divide(2,2);
	  TVirtualPad* myPad;
	  myPad = resoXEndCanvas->GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoXEndCanvas->GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalXEndcapsMap[tag].Draw();
	  myPad = resoXEndCanvas->GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalXEndcapsDistribution[tag].SetStats(1);
	  const double normA = 1. / parametrizedResolutionLocalXEndcapsDistribution[tag].Integral();
	  parametrizedResolutionLocalXEndcapsDistribution[tag].Scale(normA, "width");
	  parametrizedResolutionLocalXEndcapsDistribution[tag].Draw();
	  myPad = resoXEndCanvas->GetPad(3);
	  myPad->cd();
	  incidentAngleLocalXEndcapsDistribution[tag].SetStats(1);
	  const double normB = 1. / incidentAngleLocalXEndcapsDistribution[tag].Integral();
	  incidentAngleLocalXEndcapsDistribution[tag].Scale(normB, "width");
	  incidentAngleLocalXEndcapsDistribution[tag].Draw();
	  myPad = resoXEndCanvas->GetPad(4);
	  myPad->cd();
	  trackPhiEndcapsDistribution[tag].SetStats(1);
	  const double normC = 1. / trackPhiEndcapsDistribution[tag].Integral();
	  trackPhiEndcapsDistribution[tag].Scale(normC, "width");
	  trackPhiEndcapsDistribution[tag].Draw();
	  RootWImage& resoXEndImage = parametrizedResolutionContent.addImage(std::move(resoXEndCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoXEndImage.setComment(Form("Resolution on local X coordinate for %s endcaps modules", tag.c_str()));
	  resoXEndImage.setName(Form("Resolution on local X coordinate for %s endcaps modules", tag.c_str()));
	}
	if (parametrizedResolutionLocalYEndcapsMap[tag].GetEntries() != 0) {
	  std::unique_ptr<TCanvas> resoYEndCanvas(new TCanvas());
	  resoYEndCanvas->SetFillColor(color_plot_background);
	  resoYEndCanvas->Divide(2,2);
	  TVirtualPad* myPad;
	  myPad = resoYEndCanvas->GetPad(0);
	  myPad->SetFillColor(color_pad_background);
	  myPad = resoYEndCanvas->GetPad(1);
	  myPad->cd();
	  parametrizedResolutionLocalYEndcapsMap[tag].Draw();
	  myPad = resoYEndCanvas->GetPad(2);
	  myPad->cd();
	  parametrizedResolutionLocalYEndcapsDistribution[tag].SetStats(1);
	  const double normA = 1. / parametrizedResolutionLocalYEndcapsDistribution[tag].Integral();
	  parametrizedResolutionLocalYEndcapsDistribution[tag].Scale(normA, "width");
	  parametrizedResolutionLocalYEndcapsDistribution[tag].Draw();
	  myPad = resoYEndCanvas->GetPad(3);
	  myPad->cd();
	  incidentAngleLocalYEndcapsDistribution[tag].SetStats(1);
	  const double normB = 1. / incidentAngleLocalYEndcapsDistribution[tag].Integral();
	  incidentAngleLocalYEndcapsDistribution[tag].Scale(normB, "width");
	  incidentAngleLocalYEndcapsDistribution[tag].Draw();
	  myPad = resoYEndCanvas->GetPad(4);
	  myPad->cd();
	  trackEtaEndcapsDistribution[tag].SetStats(1);
	  const double normC = 1. / trackEtaEndcapsDistribution[tag].Integral();
	  trackEtaEndcapsDistribution[tag].Scale(normC, "width");
	  trackEtaEndcapsDistribution[tag].Draw();
	  RootWImage& resoYEndImage = parametrizedResolutionContent.addImage(std::move(resoYEndCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
	  resoYEndImage.setComment(Form("Resolution on local Y coordinate for %s endcaps modules", tag.c_str()));
	  resoYEndImage.setName(Form("Resolution on local Y coordinate for %s endcaps modules", tag.c_str()));
	}
      }
    } // debugResolution


    //********************************//
    //*                              *//
    //*       Plots                  *//
    //*                              *//
    //********************************//
    RootWImage* myImage;
    std::unique_ptr<TCanvas> summaryCanvas = nullptr;
    std::unique_ptr<TCanvas> RZCanvas = nullptr;
    std::unique_ptr<TCanvas> RZCanvasBarrel = nullptr;
    std::unique_ptr<TCanvas> XYCanvas = nullptr;
    std::vector<std::unique_ptr<TCanvas>> XYCanvasesEC;
    std::unique_ptr<TCanvas> myCanvas = nullptr;
    createSummaryCanvasNicer(tracker, RZCanvas, RZCanvasBarrel, XYCanvas, XYCanvasesEC);
    if (isPixelTracker) {
      logINFO("PIXEL HACK for beam pipe");
      TPolyLine* beampipe  = new TPolyLine();
      beampipe->SetPoint(0, 0, 45/2.);
      beampipe->SetPoint(1, 2915/2., 45/2.);
      beampipe->SetPoint(2, 3804/2., 56.6/2.);
      beampipe->SetPoint(3, 3804/2.+1164, 91/2.);
      for (auto& XYCanvasEC : XYCanvasesEC) {
	XYCanvasEC->cd();
	drawCircle(22.5, true, 18); // "grey18"
      }
      XYCanvas->cd();
      drawCircle(22.5, true, 18); // "grey18"
      RZCanvas->cd();
      beampipe->Draw("same");
      RZCanvasBarrel->cd();
      beampipe->Draw("same");

      RZCanvas->cd();
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
      myImage = new RootWImage(std::move(summaryCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("Tracker summary: modules position in XY (endcap and barrel), YZ and number of hits vs. eta.");
      myContent->addItem(myImage);
    }
 
    if (RZCanvas) {
      const int windowWidth = RZCanvas->GetWindowWidth();
      const int windowHeight = RZCanvas->GetWindowHeight();
      myImage = new RootWImage(std::move(RZCanvas), windowWidth, windowHeight);
      myImage->setComment("RZ positions of the modules.");
      myContent->addItem(myImage);
    }
    if ((RZCanvasBarrel) && isPixelTracker) {
      myImage = new RootWImage(std::move(RZCanvasBarrel), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("RZ positions of the barrel modules.");
      myContent->addItem(myImage);
    }
    if (XYCanvas) {
      myImage = new RootWImage(std::move(XYCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment("XY Section of the tracker barrel.");
      myContent->addItem(myImage);
    }
    for (auto& XYCanvasEC : XYCanvasesEC ) {
      const std::string comment = XYCanvasEC->GetTitle();
      myImage = new RootWImage(std::move(XYCanvasEC), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      myImage->setComment(comment);
      myContent->addItem(myImage);
    }

    // Eta profile big plot
    myCanvas.reset(new TCanvas("EtaProfileHits", "Eta profile (Hit Modules)", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    drawEtaProfiles(*myCanvas.get(), analyzer);
    myImage = new RootWImage(std::move(myCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit modules across eta.");
    myContent->addItem(myImage);

    myCanvas.reset(new TCanvas("EtaProfileSensors", "Eta profile (Hits)", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    drawEtaProfilesSensors(*myCanvas.get(), analyzer);
    myImage = new RootWImage(std::move(myCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit coverage across eta.");
    myContent->addItem(myImage);

    myCanvas.reset(new TCanvas("EtaProfileStubs", "Eta profile (Stubs)", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    drawEtaProfilesStubs(*myCanvas.get(), analyzer);
    myImage = new RootWImage(std::move(myCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Stub coverage across eta.");
    myContent->addItem(myImage);

    myCanvas.reset(new TCanvas("EtaProfileNumberOfStubsRatios", "Eta profile (Stubs)", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    if (myCanvas) drawTracksDistributionPerNumberOfStubs(*myCanvas.get(), analyzer, isPixelTracker);
    myImage = new RootWImage(std::move(myCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Stub coverage across eta.");
    myContent->addItem(myImage);

    myCanvas.reset(new TCanvas("EtaProfileLayers", "Eta profile (Layers)", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    drawEtaProfilesLayers(*myCanvas.get(), analyzer);
    myImage = new RootWImage(std::move(myCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Layer coverage across eta.");
    myContent->addItem(myImage);

    if (!isPixelTracker) {
      totalEtaProfileSensors_ = &analyzer.getTotalEtaProfileSensors();
      totalEtaProfileLayers_ = &analyzer.getTotalEtaProfileLayers();
    }
    else {
      totalEtaProfileSensorsPixel_ = &analyzer.getTotalEtaProfileSensors();
      totalEtaProfileLayersPixel_ = &analyzer.getTotalEtaProfileLayers();
    }

    std::unique_ptr<TCanvas> hitMapCanvas(new TCanvas("hitmapcanvas", "Hit Map", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    hitMapCanvas->cd();
    //gStyle->SetPalette(1);
    hitMapCanvas->SetFillColor(color_plot_background);
    hitMapCanvas->SetBorderMode(0);
    hitMapCanvas->SetBorderSize(0);
    analyzer.getMapPhiEta().Draw("colz");
    analyzer.getMapPhiEta().SetStats(0);
    hitMapCanvas->Modified();
    myImage = new RootWImage(std::move(hitMapCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Hit coverage in eta, phi");
    myContent->addItem(myImage);

    drawHitCoveragePerLayer(*myPage, analyzer, isPixelTracker);
    drawStubCoveragePerLayer(*myPage, analyzer, isPixelTracker);
    if (isPixelTracker) { drawStubWith3HitsCoveragePerLayer(*myPage, analyzer); }

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
    bool withHeader = true;
    myTextFile->addText(createAllModulesCsv(tracker, withHeader));
    filesContent->addItem(myTextFile);

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

  /*
   * Draw the distribution of tracks which have exactly a given number of stubs.
   * If the stubs number is higher than plotMaxNumberOfStubs, all the tracks are gathered into the same category.
   */
  bool Vizard::drawTracksDistributionPerNumberOfStubs(TVirtualPad& myPad, Analyzer& analyzer, const bool isPixelTracker) {
    myPad.cd();
    myPad.SetFillColor(color_plot_background);
    std::map<int, TProfile>& tracksDistributionPerNumberOfStubs = analyzer.getTracksDistributionPerNumberOfStubs();
    
    const int plotMaxNumberOfStubs = (!isPixelTracker ? plotMaxNumberOfOuterTrackerStubs :  plotMaxNumberOfInnerTrackerStubs);

    TLegend* layerLegend = new TLegend(0.905, 0.5, 1., 0.9);
    int colorIndex = 1;
    for (auto& detailIt : tracksDistributionPerNumberOfStubs) {
      const int numberOfStubs = detailIt.first;
      ostringstream titleStream;
      if (numberOfStubs < plotMaxNumberOfStubs) {
	titleStream << numberOfStubs << " stub";
      }
      else { titleStream << ">=" << numberOfStubs << " stub"; }
      if (numberOfStubs >= 2) titleStream << "s";

      TProfile& detailProfile = detailIt.second;
      detailProfile.SetMinimum(0);
      detailProfile.SetMaximum(1.05);
      detailProfile.SetMarkerColor(Palette::color(colorIndex));
      detailProfile.SetLineColor(Palette::color(colorIndex));
      detailProfile.SetFillColor(Palette::color(colorIndex));
      detailProfile.SetMarkerStyle(8);
      detailProfile.SetMarkerSize(1.5);
      detailProfile.GetYaxis()->SetTitleOffset(1.3);
      detailProfile.SetStats(0);
      detailProfile.Draw("same");
      layerLegend->AddEntry(&detailProfile, titleStream.str().c_str(), "f");
      colorIndex++;
    }
    layerLegend->Draw("same");


    return true;
  }

  bool Vizard::drawEtaProfilesLayers(TVirtualPad& myPad, Analyzer& analyzer) {
    myPad.cd();
    myPad.SetFillColor(color_plot_background);
    TProfile& totalEtaProfileLayers = analyzer.getTotalEtaProfileLayers();
    std::vector<TProfile> etaProfilesLayers; // Empty vector: we count layers, therefore there is no plot "by module type" here
    return drawEtaProfilesAny(totalEtaProfileLayers, etaProfilesLayers);
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

  bool Vizard::drawTracksDistributionPerNumberOfStubs(TCanvas& myCanvas, Analyzer& analyzer, const bool isPixelTracker) {
    TVirtualPad* myVirtualPad = myCanvas.GetPad(0);
    if (!myVirtualPad) return false;
    return drawTracksDistributionPerNumberOfStubs(*myVirtualPad, analyzer, isPixelTracker);
  }

  bool Vizard::drawEtaProfilesLayers(TCanvas& myCanvas, Analyzer& analyzer) {
    TVirtualPad* myVirtualPad = myCanvas.GetPad(0);
    if (!myVirtualPad) return false;
    return drawEtaProfilesLayers(*myVirtualPad, analyzer);
  }

  bool Vizard::drawHitCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer, const bool isPixelTracker) {
    std::map<std::string, std::map<int, double> > emptyStubWith3HitsCountCount;
    return drawCoveragePerlayer(myPage, isPixelTracker, "hit", analyzer.getHitCoveragePerLayer(), analyzer.getHitCoveragePerLayerDetails(), emptyStubWith3HitsCountCount);
  }

  bool Vizard::drawStubCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer, const bool isPixelTracker) {   
    const std::string type = (!isPixelTracker ? "stub" : "1 stub <-> (>= 2 hits)");
    std::map<std::string, CoveragePerNumberOfHits> emptyDetailedPlots;
    std::map<std::string, std::map<int, double> > emptyStubWith3HitsCountCount;
    return drawCoveragePerlayer(myPage, isPixelTracker, type, analyzer.getStubCoveragePerLayer(), emptyDetailedPlots, emptyStubWith3HitsCountCount);
  }

  bool Vizard::drawStubWith3HitsCoveragePerLayer(RootWPage& myPage, Analyzer& analyzer) {
    const bool isPixelTracker = true;
    std::map<std::string, CoveragePerNumberOfHits> emptyDetailedPlots;
    return drawCoveragePerlayer(myPage, isPixelTracker, "1 stub <-> (>= 3 hits)", analyzer.getStubWith3HitsCoveragePerLayer(), emptyDetailedPlots, analyzer.getStubWith3HitsCountPerDiskAndRing());
  }

  /*
   * For each layer, draw dedicated hits or stubs coverage plots.
   */
  bool Vizard::drawCoveragePerlayer(RootWPage& myPage, const bool isPixelTracker, const std::string type, 
				    std::map<std::string, TProfile>& coveragePerLayer, 
				    std::map<std::string, CoveragePerNumberOfHits>& coveragePerLayerDetails, 
				    std::map<std::string, std::map<int, double> >& stubWith3HitsCountPerDiskAndRing) {

    if (coveragePerLayer.size() == 0) return false;

    RootWContent* myContent = new RootWContent("Layer coverage (" + type + ")", false);
    myPage.addContent(myContent);

    const bool is3HitsStubStudy = (isPixelTracker && type.find("stub") != std::string::npos && type.find("3") != std::string::npos);

    // Loop on all layer(s) / disk(s)
    for (auto& layerIt : coveragePerLayer) {
      const std::string layerName = layerIt.first;

      std::unique_ptr<TCanvas> myCanvas(new TCanvas(Form("LayerCoverage%s%s", layerName.c_str(), type.c_str()), ("Layer eta coverage (" + type + ")").c_str(), vis_std_canvas_sizeX, vis_min_canvas_sizeY));
      myCanvas->cd();

      // Distribution of tracks with at least 1 hit / stub.
      TProfile& overallCoverage = layerIt.second;
      overallCoverage.SetMinimum(0);
      overallCoverage.SetMaximum(1.05);
      overallCoverage.SetMarkerColor(Palette::color(1));
      overallCoverage.SetLineColor(Palette::color(1));
      overallCoverage.SetMarkerStyle(1);    
      
      // Regarding Inner tracker 'stubs', it doesn't make sense to zoom in to [0.9 1.0].
      // The 'stub' coverage does not intend to be hermetic here.
      if (isPixelTracker && type.find("stub") != std::string::npos) {
	overallCoverage.Draw();
	TLegend* centralLegend = new TLegend(0.85, 0.8, 0.999, 0.9);
	// In IT, a stub is defined by (>=2 hits / layer) or (>=3 hits/layer).
	std::string entry = (is3HitsStubStudy ? ">= 3 hits" : ">= 2 hits"); 
	centralLegend->AddEntry(&overallCoverage, entry.c_str());
	centralLegend->Draw();
      }
      else {
	// Prepare 3 pads
	TPad* upperLeftPad = new TPad(Form("%s_upper_left", myCanvas->GetName()), "upperLeft", 0, 0.4, 0.5, 1);
	upperLeftPad->Draw();
	TPad* upperRightPad = new TPad(Form("%s_upper_right", myCanvas->GetName()), "upperRight", 0.5, 0.4, 1, 1);
	upperRightPad->Draw();
	TPad* lowerPad = new TPad(Form("%s_lower", myCanvas->GetName()), "upper", 0, 0, 1, 0.4);
	lowerPad->Draw();

	// Top left pad: overall hit coverage.
	// Also draw detailed counts, if exist.
	upperLeftPad->cd();
	overallCoverage.Draw();
	TLegend* layerLegend = new TLegend(0.85, 0.65, 0.999, 0.95);
	std::string overallCoverageTitle = ">=1 " + type + "(s)";
	layerLegend->AddEntry(&overallCoverage, overallCoverageTitle.c_str());
	// Draw detailed hit counts per layer
	if (type == "hit") {
	  const int plotMaxNumberOfHits = (!isPixelTracker ? plotMaxNumberOfOuterTrackerHitsPerLayer : plotMaxNumberOfInnerTrackerHitsPerLayer);
	  const auto& found = coveragePerLayerDetails.find(layerName);
	  if (found != coveragePerLayerDetails.end()) {
	    CoveragePerNumberOfHits& detailedInfo = found->second;
	    drawCoveragePerlayerDetails(detailedInfo, type, layerLegend, plotMaxNumberOfHits);
	  }	  
	}
	layerLegend->Draw("same");
	
	// Top right pad: statistics.
	upperRightPad->cd();
	TH1D* efficiencyHistogram = new TH1D(Form("%s_histo", overallCoverage.GetName()), overallCoverage.GetTitle(), 50, 0, .1);
	efficiencyHistogram->SetXTitle("Inefficiency");
	efficiencyHistogram->SetYTitle("Eta bins");
	efficiencyHistogram->SetFillColor(Palette::color(1));
	efficiencyHistogram->Draw();
	for (int i=1; i<=overallCoverage.GetNbinsX(); ++i) efficiencyHistogram->Fill(1-overallCoverage.GetBinContent(i));
	TPaveText* tpt;
	tpt = new TPaveText(0.65, 0.65, 0.95, 0.95, "NB NDC");
	tpt->SetBorderSize(1);
	tpt->AddText(Form("#mu = %f%%", 100*efficiencyHistogram->GetMean()));
	tpt->AddText(Form("#sigma = %f%%", 100*efficiencyHistogram->GetRMS()));	
	tpt->Draw();

	// Bottom pad: overall coverage is showed zoomed in to [0.9 1.01] window.
	lowerPad->cd();
	TProfile* zoomedProfile = (TProfile*) overallCoverage.Clone();
	zoomedProfile->SetMinimum(0.9);
	zoomedProfile->SetMaximum(1.01);
	zoomedProfile->SetTitle("");	
	zoomedProfile->Draw();
	std::ostringstream entry;
	entry << ">=1 " << type << "(s)";
	TLegend* zoomedLegend = new TLegend(0.85, 0.65, 0.999, 0.95);
	zoomedLegend->AddEntry(zoomedProfile, entry.str().c_str());
	zoomedLegend->Draw();
      }

      RootWImage* myImage = new RootWImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      const std::string title = (type.find("stub") != std::string::npos ? "stub" : type);
      myImage->setComment("Layer coverage in eta for " + title + "s.");
      myContent->addItem(myImage);
    }

    // IT only: add tables with detailed 3-hit stubs counts, per ring transition.
    if (is3HitsStubStudy && !stubWith3HitsCountPerDiskAndRing.empty()) {
      for (auto& diskIt : stubWith3HitsCountPerDiskAndRing) {
	RootWTable* stubWith3HitsCountTable = new RootWTable();
	const std::string diskName = diskIt.first;
	const std::map<int, double>& stubWith3HitsCountPerRing = diskIt.second;
	stubWith3HitsCountTable->setContent(0, 0, diskName);
	stubWith3HitsCountTable->setContent(1, 0, "Ring transition (Ring i & i+1):");
	stubWith3HitsCountTable->setContent(2, 0, "Fraction of tracks (‰)");

	for (auto& ringTransitionIt : stubWith3HitsCountPerRing) {
	  const int ringTransition = ringTransitionIt.first;
	  const double stubWith3HitsCount = ringTransitionIt.second;
	  stubWith3HitsCountTable->setContent(1, ringTransition, ringTransition);
	  stubWith3HitsCountTable->setContent(2, ringTransition, stubWith3HitsCount * 1000., 1);
	}
	myContent->addItem(stubWith3HitsCountTable);
      }
    }

    return true;
  }

  /*
   * Draw detailed hits or stubs counts per layer.
   */
  void Vizard::drawCoveragePerlayerDetails(CoveragePerNumberOfHits& detailedInfo, const std::string type, TLegend* layerLegend, const int plotMaxNumberOfHits) {
    int colorIndex = 2;

    for (auto& detailIt : detailedInfo) {
      const int numberOfHits = detailIt.first;
      TProfile& detailProfile = detailIt.second;

      if (detailProfile.GetMaximum() > insur::hits_negligible) {	  
	detailProfile.SetMinimum(0);
	detailProfile.SetMaximum(1.05);
	detailProfile.SetMarkerColor(Palette::color(colorIndex));
	detailProfile.SetLineColor(Palette::color(colorIndex));
	detailProfile.SetMarkerStyle(8);
	detailProfile.SetMarkerSize(0.3);  
	detailProfile.Draw("same");

	std::ostringstream titleStream;
	if (numberOfHits < plotMaxNumberOfHits) {
	  titleStream << "= " << numberOfHits << " " << type;
	}
	else { titleStream << ">=" << numberOfHits << " " << type; }
	if (numberOfHits >= 2) titleStream << "s";
	layerLegend->AddEntry(&detailProfile, titleStream.str().c_str());
      }
      colorIndex++;
    }

  }


  static int nLayoutCanvases = 0;
  
  std::unique_ptr<TCanvas> Vizard::drawFullLayoutRZ() {
    std::unique_ptr<TCanvas> result = nullptr;
    PlotDrawer<YZ, Type> yzDrawer;

    for (unsigned int i=0; i< trackers_.size(); ++i) {
      Tracker& tracker = *(trackers_[i]);
      yzDrawer.addModules(tracker);
    }

    int rzCanvasX = vis_max_canvas_sizeX;
    int rzCanvasY = vis_min_canvas_sizeY;
    result.reset(new TCanvas(Form("FullRZCanvas%d", nLayoutCanvases++), "RZView Canvas (full layout)", rzCanvasX, rzCanvasY ));
    result->cd();
    yzDrawer.drawFrame<SummaryFrameStyle>(*result.get());
    yzDrawer.drawModules<ContourStyle>(*result.get());

    return result;
  }

  std::unique_ptr<TCanvas> Vizard::drawFullLayoutServicesRZ() {
    std::unique_ptr<TCanvas> result = drawFullLayoutRZ();
    result->cd();
    for (auto& matBud : materialBudgets_ ) {
      std::vector<InactiveElement> allServices = matBud->getAllServices();
      for (auto& iter : allServices) {
	double z1 = iter.getZOffset();
	double z2 = iter.getZOffset()+iter.getZLength();
	double r1 = iter.getInnerRadius();
	double r2 = iter.getInnerRadius()+iter.getRWidth();
	if ((z1>=0)||(z2>0)) {
	  TBox* myBox = new TBox(z1, r1, z2, r2);
	  myBox->SetLineColor(kGray);
	  myBox->SetFillColor(kGray);
	  myBox->SetFillStyle(1);
	  myBox->SetLineStyle(1);
	  myBox->Draw("l");
	}
      }
    }
    
    return result;
  }

  
  std::unique_ptr<TCanvas> Vizard::drawFullLayoutBarrelXY() {
    std::unique_ptr<TCanvas> result = nullptr;
    PlotDrawer<XY, Type> xyDrawer;

    for (unsigned int i=0; i< trackers_.size(); ++i) {
      Tracker& tracker = *(trackers_[i]);
      xyDrawer.addModulesType(tracker, BARREL);
    }

    int xyCanvasX = vis_min_canvas_sizeX;
    int xyCanvasY = vis_min_canvas_sizeY;
    result.reset(new TCanvas("FullXYCanvas", "XYView Canvas (full layout)", xyCanvasX, xyCanvasY ));
    result->cd();
    xyDrawer.drawFrame<SummaryFrameStyle>(*result.get());
    xyDrawer.drawModules<ContourStyle>(*result.get());

    return result;
  }


  void Vizard::stackHistos(std::vector<std::pair<std::string, TH1D*>>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation) {
    TProfile* prof;
    TH1D* histo;
    for (const auto& it : histoMap) {
      prof = newProfile((TH1D*)it.second, 0., 4.0, materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(index));
      histo->SetFillColor(Palette::color(index));
      histo->SetTitle(it.first.c_str());
      if (isRadiation) legend->AddEntry(histo, it.first.c_str());
      myStack->Add(histo);
      totalStack->Add(histo);
      myTable->setContent(index, 0, it.first);
      int column = (isRadiation ? 1 : 2);
      myTable->setContent(index++, column, averageHistogramValues(*histo, 4.0), 5);
    }
  }


  void Vizard::stackHistos(std::map<std::string, TH1D*>& histoMap, RootWTable*& myTable, int& index, THStack*& totalStack, THStack*& myStack, TLegend*& legend, bool& isRadiation) {
    TProfile* prof;
    TH1D* histo;
    for (const auto& it : histoMap) {
      prof = newProfile((TH1D*)it.second, 0., 4.0, materialNBins);
      histo = prof->ProjectionX();
      histo->SetLineColor(Palette::color(index));
      histo->SetFillColor(Palette::color(index));
      histo->SetTitle(it.first.c_str());
      if (isRadiation) legend->AddEntry(histo, it.first.c_str());
      myStack->Add(histo);
      totalStack->Add(histo);
      myTable->setContent(index, 0, it.first);
      int column = (isRadiation ? 1 : 2);
      myTable->setContent(index++, column, averageHistogramValues(*histo, 4.0), 5);
    }
  }

  void Vizard::totalMaterialSummary(Analyzer& analyzer, Analyzer& pixelAnalyzer, RootWSite& site) {
    // Pointer to an image to create on the fly
    RootWImage* myImage;
    RootWPage& myPage = site.addPage("Material (total)");
    
    // Define web-page sections
    RootWContent* materialCategoriesContent = new RootWContent("Full layout Material : Categories details (tracking volume)", true);
    myPage.addContent(materialCategoriesContent);
    RootWContent* materialOverviewContent   = new RootWContent("Full layout Material : 1d overview (tracking volume)", false);
    myPage.addContent(materialOverviewContent);
    RootWContent* materialComponentsContent = new RootWContent("Full layout Material : Components details (tracking volume)", false);
    myPage.addContent(materialComponentsContent);

    // COMPONENTS DETAILS (TRACKING VOLUME)
    RootWTable* myTable = new RootWTable();
    char titleString[256];
    sprintf(titleString, "Average (eta = [0, %.1f])", analyzer.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");

    THStack* rCompBeamPipeStack = new THStack("rcompbeampipestack", "Radiation Length by Component in beam pipe");
    THStack* iCompBeamPipeStack = new THStack("icompbeampipestack", "Interaction Length by Component in beam pipe");
    THStack* rCompPixelIntersticeStack = new THStack("rcomppixelintersticestack", "Radiation Length by Component in pixel interstice");
    THStack* iCompPixelIntersticeStack = new THStack("icomppixelintersticestack", "Interaction Length by Component in pixel interstice");
    THStack* rCompPixelTrackingVolumeStack = new THStack("rcomppixeltrackingvolumestack", "Radiation Length by Component in pixel tracking volume");
    THStack* iCompPixelTrackingVolumeStack = new THStack("icomppixeltrackingvolumestack", "Interaction Length by Component in pixel tracking volume");
    THStack* rCompIntersticeStack = new THStack("rcompintersticestack", "Radiation Length by Component in interstice");
    THStack* iCompIntersticeStack = new THStack("icompintersticestack", "Interaction Length by Component in interstice");
    THStack* rCompOuterTrackingVolumeStack = new THStack("rcompoutertrackingvolumestack", "Radiation Length by Component in outer tracking volume");
    THStack* iCompOuterTrackingVolumeStack = new THStack("icompoutertrackingvolumestack", "Interaction Length by Component in outer tracking volume");
    THStack* rCompTotalTrackingVolumeStack = new THStack("rcomptotaltrackingvolumestack", "Radiation Length by Component in tracking volume");
    THStack* iCompTotalTrackingVolumeStack = new THStack("icomptotaltrackingvolumestack", "Interaction Length by Component in tracking volume");

    TLegend* compLegend = new TLegend(0.1,0.6,0.35,0.9);

    std::unique_ptr<TCanvas> myCanvas(new TCanvas("FullLayoutMaterialComponentsTrackingVolumeRI"));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    TVirtualPad* myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    bool isRadiation = true;
    int index = 1;
    stackHistos(analyzer.getHistoBeamPipeR(), myTable, index, rCompTotalTrackingVolumeStack, rCompBeamPipeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoPixelIntersticeR(), myTable, index, rCompTotalTrackingVolumeStack, rCompPixelIntersticeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoPixelTrackingVolumeR(), myTable, index, rCompTotalTrackingVolumeStack, rCompPixelTrackingVolumeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoIntersticeR(), myTable, index, rCompTotalTrackingVolumeStack, rCompIntersticeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoOuterTrackingVolumeR(), myTable, index, rCompTotalTrackingVolumeStack, rCompOuterTrackingVolumeStack, compLegend, isRadiation);
    rCompTotalTrackingVolumeStack->Draw("hist");
    compLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    isRadiation = false;
    index = 1;
    stackHistos(analyzer.getHistoBeamPipeI(), myTable, index, iCompTotalTrackingVolumeStack, iCompBeamPipeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoPixelIntersticeI(), myTable, index, iCompTotalTrackingVolumeStack, iCompPixelIntersticeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoPixelTrackingVolumeI(), myTable, index, iCompTotalTrackingVolumeStack, iCompPixelTrackingVolumeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoIntersticeI(), myTable, index, iCompTotalTrackingVolumeStack, iCompIntersticeStack, compLegend, isRadiation);
    stackHistos(analyzer.getHistoOuterTrackingVolumeI(), myTable, index, iCompTotalTrackingVolumeStack, iCompOuterTrackingVolumeStack, compLegend, isRadiation);
    iCompTotalTrackingVolumeStack->Draw("hist");
    compLegend->Draw();

    materialComponentsContent->addItem(myTable);
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by component type in total tracking volume");
    myImage->setName("fullLayoutMatComponentsTrackingVolume");
    materialComponentsContent->addItem(myImage);

    // CATEGORIES DETAILS (TRACKING VOLUME)
    myTable = new RootWTable();
    sprintf(titleString, "Average (eta = [0, %.1f])", analyzer.getEtaMaxMaterial());
    myTable->setContent(0, 0, titleString);
    myTable->setContent(0, 1, "Radiation length");
    myTable->setContent(0, 2, "Interaction length");
    compLegend = new TLegend(0.1,0.6,0.35,0.9);

    std::vector<std::pair<std::string, TH1D*>> histoPerCategoryR, histoPerCategoryI;
    if (rCompBeamPipeStack->GetHists()) histoPerCategoryR.push_back(std::make_pair("Beam pipe", (TH1D*)rCompBeamPipeStack->GetStack()->Last()));
    if (iCompBeamPipeStack->GetHists()) histoPerCategoryI.push_back(std::make_pair("Beam pipe", (TH1D*)iCompBeamPipeStack->GetStack()->Last()));
    if (rCompPixelIntersticeStack->GetHists()) histoPerCategoryR.push_back(std::make_pair("Services and Supports under Pixel Tracking Volume", (TH1D*)rCompPixelIntersticeStack->GetStack()->Last()));
    if (iCompPixelIntersticeStack->GetHists()) histoPerCategoryI.push_back(std::make_pair("Services and Supports under Pixel Tracking Volume", (TH1D*)iCompPixelIntersticeStack->GetStack()->Last()));
    if (rCompPixelTrackingVolumeStack->GetHists()) histoPerCategoryR.push_back(std::make_pair("Pixel Tracking Volume", (TH1D*)rCompPixelTrackingVolumeStack->GetStack()->Last()));
    if (iCompPixelTrackingVolumeStack->GetHists()) histoPerCategoryI.push_back(std::make_pair("Pixel Tracking Volume", (TH1D*)iCompPixelTrackingVolumeStack->GetStack()->Last()));
    if (rCompIntersticeStack->GetHists()) histoPerCategoryR.push_back(std::make_pair("Services and Supports between Pixel and Outer Tracking Volumes", (TH1D*)rCompIntersticeStack->GetStack()->Last()));
    if (iCompIntersticeStack->GetHists()) histoPerCategoryI.push_back(std::make_pair("Services and Supports between Pixel and Outer Tracking Volumes", (TH1D*)iCompIntersticeStack->GetStack()->Last()));
    if (rCompOuterTrackingVolumeStack->GetHists()) histoPerCategoryR.push_back(std::make_pair("Outer Tracking Volume", (TH1D*)rCompOuterTrackingVolumeStack->GetStack()->Last()));
    if (iCompOuterTrackingVolumeStack->GetHists()) histoPerCategoryI.push_back(std::make_pair("Outer Tracking Volume", (TH1D*)iCompOuterTrackingVolumeStack->GetStack()->Last()));
    THStack* rCompCategoryTrackingVolumeStack = new THStack("rcompcategorytrackingvolumestack", "Radiation Length by Category in tracking volume");
    THStack* iCompCategoryTrackingVolumeStack = new THStack("icompcategorytrackingvolumestack", "Interaction Length by Category in tracking volume");
    THStack* dummy = new THStack("dummy", "dummy");

    myCanvas.reset(new TCanvas("FullLayoutMaterialCategoriesTrackingVolumeRI"));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);

    myPad = myCanvas->GetPad(1);
    myPad->cd();
    isRadiation = true;
    index = 1;
    stackHistos(histoPerCategoryR, myTable, index, dummy, rCompCategoryTrackingVolumeStack, compLegend, isRadiation);
    rCompCategoryTrackingVolumeStack->Draw("hist");
    compLegend->Draw();

    myPad = myCanvas->GetPad(2);
    myPad->cd();
    isRadiation = false;
    index = 1;
    stackHistos(histoPerCategoryI, myTable, index, dummy, iCompCategoryTrackingVolumeStack, compLegend, isRadiation);
    iCompCategoryTrackingVolumeStack->Draw("hist");
    compLegend->Draw();

    materialCategoriesContent->addItem(myTable);
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Radiation and interaction length distribution in eta by category in total tracking volume");
    myImage->setName("fullLayoutMatCategoriesTrackingVolume");
    materialCategoriesContent->addItem(myImage);
   
    // 1D OVERVIEW (TRACKING VOLUME)
    TH1D *cr = nullptr, *ci = nullptr;
    TProfile *crProf, *ciProf;
    // Work area re-init
    myCanvas.reset(new TCanvas("FullLayoutMaterialOverviewTrackingVolumeRI"));
    myCanvas->SetFillColor(color_plot_background);
    myCanvas->Divide(2, 1);
    myPad = myCanvas->GetPad(0);
    myPad->SetFillColor(color_pad_background);
    myPad = myCanvas->GetPad(1);
    myPad->cd();
    // global plots in tracking volume: radiation length
    if (rCompTotalTrackingVolumeStack->GetHists()) {
      cr = (TH1D*)rCompTotalTrackingVolumeStack->GetStack()->Last()->Clone();
      crProf = newProfile(cr, 0., analyzer.getEtaMaxMaterial(), materialNBins);
      crProf->SetFillColor(kGray + 2);
      crProf->SetTitle("Radiation Length within Full Tracking Volume; #eta; x/X_{0}");
      crProf->Draw("hist");
    }
    myPad = myCanvas->GetPad(2);
    myPad->cd();
    // global plots in tracking volume: interaction length
    if (iCompTotalTrackingVolumeStack->GetHists()) {
      ci = (TH1D*)iCompTotalTrackingVolumeStack->GetStack()->Last()->Clone();
      ciProf = newProfile(ci, 0., analyzer.getEtaMaxMaterial(), materialNBins);
      ciProf->SetFillColor(kGray + 2);
      ciProf->SetTitle("Interaction Length within Full Tracking Volume; #eta; #lambda/#lambda_{0}");
      ciProf->Draw("hist");
    }
    // Write global tracking volume plots to web pag
    myImage = new RootWImage(std::move(myCanvas), 2*vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Material in total tracking volume");
    myImage->setName("fullLayoutMatOverviewTrackingVolume");
    myTable = new RootWTable();
    sprintf(titleString, "Average radiation length in tracking volume (eta = [0, %.1f])", analyzer.getEtaMaxMaterial());
    myTable->setContent(1, 1, titleString);
    sprintf(titleString, "Average interaction length in tracking volume (eta = [0, %.1f])", analyzer.getEtaMaxMaterial());
    myTable->setContent(2, 1, titleString);
    if (cr) myTable->setContent(1, 2, averageHistogramValues(*cr, analyzer.getEtaMaxMaterial()), 5);
    if (ci) myTable->setContent(2, 2, averageHistogramValues(*ci, analyzer.getEtaMaxMaterial()), 5);
    materialOverviewContent->addItem(myTable);
    materialOverviewContent->addItem(myImage);
    
  }

  bool Vizard::additionalInfoSite(const std::string& settingsfile,
                                  Analyzer& analyzer, Analyzer& pixelAnalyzer, Tracker& tracker, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Info");
    myPage->setAddress("info.html");

    site.addPage(myPage, RootWeb::most_relevant);
    RootWContent *simulationContent, *summaryContent, *fullLayoutContent, *configFilesContent;

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
    std::unique_ptr<TCanvas> aLayout = drawFullLayoutRZ();
    std::unique_ptr<TCanvas> aLayoutServices = drawFullLayoutServicesRZ();
    std::unique_ptr<TCanvas> aLayoutXY = drawFullLayoutBarrelXY();
    if (aLayout||aLayoutXY) {
      fullLayoutContent = new RootWContent("Full layout Geometry", true);
      myPage->addContent(fullLayoutContent);
      if (aLayout) {
	const int windowWidth = aLayout->GetWindowWidth();
	const int windowHeight = aLayout->GetWindowHeight();
	RootWImage* anImage = new RootWImage(std::move(aLayout), windowWidth, windowHeight);
	anImage->setComment("RZ position of the modules (full layout)");
	anImage->setName("fullLayout");
	fullLayoutContent->addItem(anImage);
      }
      if (aLayoutXY) {
	const int windowWidth = aLayoutXY->GetWindowWidth();
	const int windowHeight = aLayoutXY->GetWindowHeight();
	RootWImage* anImage = new RootWImage(std::move(aLayoutXY), windowWidth, windowHeight);
	anImage->setComment("XY position of the barrel modules (full layout)");
	anImage->setName("fullLayoutBarrelXY");
	fullLayoutContent->addItem(anImage);
      }
      if (aLayoutServices) {
	const int windowWidth = aLayoutServices->GetWindowWidth();
	const int windowHeight = aLayoutServices->GetWindowHeight();
	RootWImage* anImage = new RootWImage(std::move(aLayoutServices), windowWidth, windowHeight);
	anImage->setComment("RZ position of the modules (full layout with services)");
	anImage->setName("fullLayout");
	fullLayoutContent->addItem(anImage);
      }
    }

    simulationContent = new RootWContent("Simulation parameters");
    myPage->addContent(simulationContent);
    RootWContent* materialsTablesContent = new RootWContent("Materials tables");
    myPage->addContent(materialsTablesContent);
    summaryContent = new RootWContent("Summary");
    myPage->addContent(summaryContent);
    configFilesContent = new RootWContent("Configuration files", false);
    myPage->addContent(configFilesContent);

    THStack* totalEtaStack = new THStack();
    if (totalEtaProfileSensors_) totalEtaStack->Add(totalEtaProfileSensors_->ProjectionX());
    if (totalEtaProfileSensorsPixel_) totalEtaStack->Add(totalEtaProfileSensorsPixel_->ProjectionX());
    std::unique_ptr<TCanvas> totalEtaProfileFull(new TCanvas("TotalEtaProfileFull", "Full eta profile (Hits)", vis_std_canvas_sizeX, vis_std_canvas_sizeY));
    totalEtaProfileFull->cd();
    ((TH1D*)totalEtaStack->GetStack()->Last())->SetMarkerStyle(8);
    ((TH1D*)totalEtaStack->GetStack()->Last())->SetMarkerSize(1);
    ((TH1D*)totalEtaStack->GetStack()->Last())->SetMinimum(0.);
    ((TH1D*)totalEtaStack->GetStack()->Last())->SetStats(0.);
    totalEtaStack->GetStack()->Last()->Draw();
    // add profile for types here...#####
    drawEtaProfilesSensors(*totalEtaProfileFull, analyzer, false);
    drawEtaProfilesSensors(*totalEtaProfileFull, pixelAnalyzer, false);
    totalEtaStack->GetStack()->Last()->Draw("same");
    RootWImage* myImage = new RootWImage(std::move(totalEtaProfileFull), vis_std_canvas_sizeX, vis_std_canvas_sizeY);
    myImage->setComment("Full hit coverage across eta");
    if (fullLayoutContent) fullLayoutContent->addItem(myImage);

    // Number of layers count
    std::unique_ptr<TCanvas> totalEtaProfileLayersFull(new TCanvas("totalEtaProfileLayersFull", "Full eta profile (Layers)", vis_std_canvas_sizeX, vis_std_canvas_sizeY));
    totalEtaProfileLayersFull->cd();
    THStack* totalLayersEtaStack = new THStack();
    TH1D* totalLayersCountOuter = (totalEtaProfileLayers_ ? (TH1D*)totalEtaProfileLayers_->ProjectionX()->Clone() : nullptr);
    if (totalLayersCountOuter) {
      totalLayersCountOuter->SetBit(1);
      totalLayersEtaStack->Add(totalLayersCountOuter);
      totalLayersCountOuter->SetMarkerStyle(8);
      totalLayersCountOuter->SetMarkerSize(1);
      totalLayersCountOuter->SetMarkerColor(Palette::color(1));
      totalLayersCountOuter->SetStats(0);
    }
    TH1D* totalLayersCountInner = (totalEtaProfileLayersPixel_ ? (TH1D*)totalEtaProfileLayersPixel_->ProjectionX()->Clone() : nullptr);
    if (totalLayersCountInner) {   
      totalLayersCountInner->SetBit(1);
      totalLayersEtaStack->Add(totalLayersCountInner);
      totalLayersCountInner->SetMarkerStyle(8);
      totalLayersCountInner->SetMarkerSize(1);
      totalLayersCountInner->SetMarkerColor(Palette::color(2));
      totalLayersCountInner->SetStats(0);
    }
    ((TH1D*)totalLayersEtaStack->GetStack()->Last())->SetMarkerStyle(8);
    ((TH1D*)totalLayersEtaStack->GetStack()->Last())->SetMarkerSize(1);
    ((TH1D*)totalLayersEtaStack->GetStack()->Last())->SetMarkerColor(kBlack);
    ((TH1D*)totalLayersEtaStack->GetStack()->Last())->SetMinimum(0.);
    ((TH1D*)totalLayersEtaStack->GetStack()->Last())->SetStats(0);
    totalLayersEtaStack->GetStack()->Last()->Draw();
    // Per detector here
    if (totalLayersCountOuter) totalLayersCountOuter->Draw("same");
    if (totalLayersCountInner) totalLayersCountInner->Draw("same");
    totalLayersEtaStack->GetStack()->Last()->Draw("same"); // To overwrite where total is the same as one of the two
    RootWImage* myImageLayers = new RootWImage(std::move(totalEtaProfileLayersFull), vis_std_canvas_sizeX, vis_std_canvas_sizeY);
    myImageLayers->setComment("Full layer coverage across eta (OT = blue, pixel = red)");
    if (fullLayoutContent) fullLayoutContent->addItem(myImageLayers);
    
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
    configFilesContent->addItem(myBinaryFileList);

    RootWInfo* myInfo;
    myInfo = new RootWInfo("Minimum bias per bunch crossing");
    myInfo->setValue(SimParms::getInstance().numMinBiasEvents(), minimumBiasPrecision);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo("Integrated luminosity");
    myInfo->setValue(SimParms::getInstance().timeIntegratedLumi(), luminosityPrecision);
    myInfo->appendValue(" fb" + superStart + "-1" + superEnd);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo("Number of tracks used for material");
    myInfo->setValue(materialTracksUsed);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo("Number of tracks used for geometry");
    myInfo->setValue(geometryTracksUsed);
    simulationContent->addItem(myInfo);
    myInfo = new RootWInfo(Form("Irradiation &alpha; parameter (at reference temperature %.0f °C)", SimParms::getInstance().referenceTemp()));
    myInfo->setValueSci(SimParms::getInstance().alphaParam(), alphaParamPrecision);
    myInfo->appendValue(" A/cm");
    simulationContent->addItem(myInfo);

    //********************************//
    //*                              *//
    //*  Summary files               *//
    //*                              *//
    //********************************//

    // Summary file with all root plots
    myBinaryFile = new RootWBinaryFile("summary.root", "ROOT file with all relevant plots");
    myBinaryFile->setNoCopy(true);
    summaryContent->addItem(myBinaryFile);

    RootWTextFile* myTextFile;

    // DetId modules list with associated geometry info
    myTextFile = new RootWTextFile("DetId_modules_list.csv", "DetId modules list with associated geometry info");
    myTextFile->addText(createModulesDetIdListCsv());
    summaryContent->addItem(myTextFile);

    // DetId sensors list with associated geometry info
    myTextFile = new RootWTextFile("DetId_sensors_list.csv", "DetId sensors list with associated geometry info");
    myTextFile->addText(createSensorsDetIdListCsv());
    summaryContent->addItem(myTextFile);

    RootWGraphViz* myGv = new RootWGraphViz("include_graph.gv", "Include structure");
    myGv->addText(mainConfigHandler::instance().createGraphVizFile());
    summaryContent->addItem(myGv);


    //********************************//
    //*                              *//
    //*  Materials tables            *//
    //*                              *//
    //********************************//
 
    // Chemical elements
    RootWTextFile* chemicalElementsFile = new RootWTextFile("chemical_elements.csv", "Chemical elements");
    chemicalElementsFile->addText(createChemicalElementsCsv());
    materialsTablesContent->addItem(chemicalElementsFile);

    // Chemical compounds
    bool hasChemicalFormula = true;
    RootWTextFile* chemicalCompoundsFile = new RootWTextFile("chemical_compounds.csv", "Chemical compounds");
    chemicalCompoundsFile->addText(createChemicalMixturesCsv(hasChemicalFormula));
    materialsTablesContent->addItem(chemicalCompoundsFile);

    // Chemical mixtures
    hasChemicalFormula = false;
    RootWTextFile* chemicalMixturesFile = new RootWTextFile("chemical_mixtures.csv", "Chemical mixtures");
    chemicalMixturesFile->addText(createChemicalMixturesCsv(hasChemicalFormula));
    materialsTablesContent->addItem(chemicalMixturesFile);

    return true;
  }


  bool Vizard::bandwidthSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
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
    std::unique_ptr<TCanvas> bandWidthCanvas(new TCanvas("ModuleBandwidthC", "Modules needed bandwidthC", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    std::unique_ptr<TCanvas> moduleHitCanvas(new TCanvas("ModuleHitC", "Module hit countC", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
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
    RootWImage* myImage = new RootWImage(std::move(bandWidthCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Module bandwidth distribution in the sparsified and unsparsified model");
    myContent->addItem(myImage);

    moduleHitCanvas->cd();
    TH1D& chanHitDistribution = analyzer.getChanHitDistribution();
    chanHitDistribution.Draw();
    myImage = new RootWImage(std::move(moduleHitCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage->setComment("Distribution of number of hits per bunch crossing (each sensor is counted separately)");
    myContent->addItem(myImage);

    RootWText* myDescription = new RootWText();
    myContent->addItem(myDescription);
    myDescription->addText( "Bandwidth useage estimate:<br/>");
    myDescription->addText( "(Pt modules: ignored)<br/>");
    myDescription->addText( "Sparsified (binary) bits/event: 23 bits/chip + 9 bit/hit<br/>");
    myDescription->addText( "Unsparsified (binary) bits/event: 16 bits/chip + 1 bit/channel<br/>");
    ostringstream aStringStream; aStringStream.str("100 kHz trigger, "); aStringStream << SimParms::getInstance().numMinBiasEvents();
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

    std::unique_ptr<TCanvas> triggerDataBandwidthCanvas(new TCanvas());
    std::unique_ptr<TCanvas> triggerFrequencyPerEventCanvas(new TCanvas());

    PlotDrawer<YZ, Type, Max> yzbwDrawer(0, 0); // we take the MAX because the Analyzer only sweeps across the first quadrant (up to PI/2),
    PlotDrawer<YZ, Type, Max> yztfDrawer(0, 0); // so there's plenty modules in Phi which don't have their property set, but Max disregards all the 0's

    yzbwDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);
    yztfDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end(), BARREL | ENDCAP);

    yzbwDrawer.drawFrame<HistogramFrameStyle>(*triggerDataBandwidthCanvas.get());
    yztfDrawer.drawFrame<HistogramFrameStyle>(*triggerFrequencyPerEventCanvas.get());

    yzbwDrawer.drawModules<ContourStyle>(*triggerDataBandwidthCanvas.get());
    yztfDrawer.drawModules<ContourStyle>(*triggerFrequencyPerEventCanvas.get());

    RootWImage& triggerDataBandwidthImage = myContent->addImage(std::move(triggerDataBandwidthCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    triggerDataBandwidthImage.setComment("Map of the bandwidth for trigger data in Gbps");
    triggerDataBandwidthImage.setName("triggerDataBandwidthMap");

    RootWImage& triggerFrequencyPerEventImage = myContent->addImage(std::move(triggerFrequencyPerEventCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
      std::unique_ptr<TCanvas> graphCanvas(new TCanvas());
      graphCanvas->cd();
      totalHisto->SetLineColor(1);
      totalHisto->SetMinimum(trueHisto->GetMinimum()*.9 < 0.1 ? 0 : trueHisto->GetMinimum()*.9);
      totalHisto->Draw();
      trueHisto->SetLineColor(2);
      trueHisto->Draw("SAME");
      RootWImage& graphImage = myContent->addImage(std::move(graphCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      graphImage.setComment("Stub rate for layer " + any2str(layer.second) + " (MHz/cm^2)");
      graphImage.setName("stubRate" + layer.first + any2str(layer.second));

    }


    return true;
  }


  bool Vizard::triggerProcessorsSummary(Analyzer& analyzer, Tracker& tracker, RootWSite& site) {
    RootWPage* myPage = new RootWPage("Trigger CPUs");
    RootWTextFile* myTextFile;
    myPage->setAddress("trigger_cpus.html");
    site.addPage(myPage);

    SummaryTable& processorSummary = analyzer.getProcessorConnectionSummary();
    SummaryTable& processorCommonSummary = analyzer.getProcessorCommonConnectionSummary();
    //std::map<std::string, SummaryTable>& moduleSummaries = analyzer.getModuleConnectionSummaries();

    // Connections between modules and trigger towers
    RootWContent& summaryContent = myPage->addContent("Summary tables", false);
    createTriggerSectorMapCsv(analyzer.getTriggerSectorMap());
    myTextFile = new RootWTextFile("trigger_sector_map.csv", "Trigger Towers to Modules connections");
    myTextFile->addText(triggerSectorMapCsv_);
    summaryContent.addItem(myTextFile);
    createModuleConnectionsCsv(analyzer.getModuleConnectionMap());
    myTextFile = new RootWTextFile("module_connections.csv", "Modules to Trigger Towers connections");
    myTextFile->addText(moduleConnectionsCsv_);
    summaryContent.addItem(myTextFile);
    
    myPage->addContent("Processor inbound connections").addTable().setContent(processorSummary.getContent());
    RootWContent& sharedConnContent = myPage->addContent("Processor shared inbound connections", false);
    std::unique_ptr<TCanvas> sharedConnCanvas(new TCanvas());
    sharedConnCanvas->cd();
    TH2I& sharedConnMap = analyzer.getProcessorCommonConnectionMap();
    sharedConnMap.GetXaxis()->LabelsOption("v");
    sharedConnMap.GetXaxis()->SetLabelSize(0.03);
    sharedConnMap.GetYaxis()->SetLabelSize(0.03);
    sharedConnMap.Draw("colz");
    RootWImage& sharedConnImage = sharedConnContent.addImage(std::move(sharedConnCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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

    std::unique_ptr<TCanvas> moduleConnectionEtaCanvas(new TCanvas());
    std::unique_ptr<TCanvas> moduleConnectionPhiCanvas(new TCanvas());
    std::unique_ptr<TCanvas> moduleConnectionEndcapPhiCanvas(new TCanvas());

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

    yzDrawer.drawFrame<HistogramFrameStyle>(*moduleConnectionEtaCanvas.get());
    xyDrawer.drawFrame<HistogramFrameStyle>(*moduleConnectionPhiCanvas.get());
    std::pair<Circle, Circle> petal = analyzer.getSampleTriggerPetal();
    TArc* a1 = new TArc(petal.first.x0, petal.first.y0, petal.first.r, (XYPoint(petal.first.x0, petal.first.y0)).Phi()*180./M_PI + 180.);
    TArc* a2 = new TArc(petal.second.x0, petal.second.y0, petal.second.r, 0., (XYPoint(petal.second.x0, petal.second.y0)).Phi()*180./M_PI + 180.);
    a1->SetFillStyle(0);
    a2->SetFillStyle(0);
    moduleConnectionPhiCanvas->cd();
    a1->Draw("only");
    a2->Draw("only");

    xyecDrawer.drawFrame<HistogramFrameStyle>(*moduleConnectionEndcapPhiCanvas.get());
    moduleConnectionEndcapPhiCanvas->cd();
    a1->Draw("only");
    a2->Draw("only");
    a1->SetBit(1);
    a2->SetBit(1);

    yzDrawer.drawModules<ContourStyle>(*moduleConnectionEtaCanvas.get());
    xyDrawer.drawModules<ContourStyle>(*moduleConnectionPhiCanvas.get());
    xyecDrawer.drawModules<ContourStyle>(*moduleConnectionEndcapPhiCanvas.get());



    /*
       moduleConnectionEtaCanvas->SetFillColor(color_plot_background);
       moduleConnectionPhiCanvas->SetFillColor(color_plot_background);
    //moduleConnectionEndcapPhiCanvas->SetFillColor(color_plot_background);

    moduleConnectionEtaCanvas->cd();
    moduleConnectionEtaMap.Draw("colz");
    */
    RootWImage& moduleConnectionEtaImage = myContent.addImage(std::move(moduleConnectionEtaCanvas), vis_max_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionEtaImage.setComment("Map of the number of connections to trigger processors per module (eta section)");
    moduleConnectionEtaImage.setName("moduleConnectionEtaMap");
    /*
       moduleConnectionPhiCanvas->cd();
       moduleConnectionPhiMap.Draw("colz");
       */
    RootWImage& moduleConnectionPhiImage = myContent.addImage(std::move(moduleConnectionPhiCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionPhiImage.setComment("Map of the number of connections to trigger processors per barrel module (phi section)");
    moduleConnectionPhiImage.setName("moduleConnectionPhiMap");

    // moduleConnectionEndcapPhiCanvas->cd();
    // moduleConnectionEndcapPhiMap.Draw("colz");

    RootWImage& moduleConnectionEndcapPhiImage = myContent.addImage(std::move(moduleConnectionEndcapPhiCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    moduleConnectionEndcapPhiImage.setComment("Map of the number of connections to trigger processors per endcap module (phi section)");
    moduleConnectionEndcapPhiImage.setName("moduleConnectionEndcapPhiMap");

    //    myContent = myPage->addContent("Module Connections distribution", true);

    std::unique_ptr<TCanvas> moduleConnectionsCanvas(new TCanvas("ModuleConnectionsC", "Modules connectionsC", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    moduleConnectionsCanvas->cd();
    TH1I& moduleConnectionsDistribution = analyzer.getModuleConnectionsDistribution();
    moduleConnectionsDistribution.SetFillColor(Palette::color(2));
    moduleConnectionsDistribution.Draw();
    RootWImage& myImage = myContent.addImage(std::move(moduleConnectionsCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
    myImage.setComment("Module connections distribution");

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
      //double verticalScale=1;
      if (additionalTag!="") {
        pageTitle += " ("+additionalTag+")";
        additionalSummaryTag = "_"+additionalTag+"_";
        //verticalScale = 10;
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

        std::unique_ptr<TCanvas> linearMomentumCanvas(new TCanvas());
        std::unique_ptr<TCanvas> momentumCanvas(new TCanvas());
        std::unique_ptr<TCanvas> distanceCanvas(new TCanvas());
        std::unique_ptr<TCanvas> angleCanvas(new TCanvas());
        std::unique_ptr<TCanvas> ctgThetaCanvas(new TCanvas());
        std::unique_ptr<TCanvas> etaCanvas(new TCanvas());
        std::unique_ptr<TCanvas> z0Canvas(new TCanvas());
        std::unique_ptr<TCanvas> pCanvas(new TCanvas());

        int myColor=0;
        int nRebin = 2;
        int markerStyle = 21;
        double markerSize = 1.;
        double lineWidth = 2.;

        linearMomentumCanvas->SetGrid(1,1);
        momentumCanvas->SetGrid(1,1);
        distanceCanvas->SetGrid(1,1);
        angleCanvas->SetGrid(1,1);
        ctgThetaCanvas->SetGrid(1,1);
        etaCanvas->SetGrid(1,1);
        z0Canvas->SetGrid(1,1);
        pCanvas->SetGrid(1,1);
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
          linearMomentumCanvas->SetLogy(0);
          momentumCanvas->SetLogy(1);
          momentumProfile.SetLineColor(momentumColor(myColor));
          momentumProfile.SetMarkerColor(momentumColor(myColor));
          momentumProfile.SetLineWidth(lineWidth);
          myColor++;
          momentumProfile.SetMarkerStyle(markerStyle);
          momentumProfile.SetMarkerSize(markerSize);
          momentumCanvas->SetFillColor(color_plot_background);
          linearMomentumCanvas->SetFillColor(color_plot_background);
          if (momentumGraph.GetN()>0) {
            momentumCanvas->cd();
            momentumProfile.Draw(plotOption.c_str());
            linearMomentumCanvas->cd();
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
          distanceCanvas->SetLogy();
          distanceProfile.SetLineColor(momentumColor(myColor));
          distanceProfile.SetMarkerColor(momentumColor(myColor));
          distanceProfile.SetLineWidth(lineWidth);
          myColor++;
          distanceProfile.SetMarkerStyle(markerStyle);
          distanceProfile.SetMarkerSize(markerSize);
          distanceCanvas->SetFillColor(color_plot_background);
          if (distanceGraph.GetN()>0) {
            distanceCanvas->cd();
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
          angleCanvas->SetLogy();
          angleProfile.SetLineColor(momentumColor(myColor));
          angleProfile.SetMarkerColor(momentumColor(myColor));
          angleProfile.SetLineWidth(lineWidth);
          myColor++;
          angleProfile.SetMarkerStyle(markerStyle);
          angleProfile.SetMarkerSize(markerSize);
          angleCanvas->SetFillColor(color_plot_background);
          if (angleGraph.GetN() > 0) {
            angleCanvas->cd();
            angleProfile.Draw(plotOption.c_str());
            plotOption = "same";
          }
        }
        plotOption = "";
        myColor=0;
        // ctgTheta and eta canvas loop
        g_guard = a.getCtgThetaGraphs(idealMaterial, isTrigger).end();
        for (g_iter = a.getCtgThetaGraphs(idealMaterial, isTrigger).begin(); g_iter != g_guard; g_iter++) {
          TGraph& ctgThetaGraph = g_iter->second;
          TProfile& ctgThetaProfile = newProfile(ctgThetaGraph, 0, a.getEtaMaxTracker(), nRebin);
	  TProfile& etaProfile = newProfile_timesSin(ctgThetaGraph, 0, a.getEtaMaxTracker(), nRebin);
	  etaProfile.SetTitle("Pseudorapidity error - const P_{T} across #eta;#eta;#delta #eta");
          ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);
          ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);
          ctgThetaCanvas->SetLogy();
          ctgThetaProfile.SetLineColor(momentumColor(myColor));
          ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
          etaProfile.SetMinimum(vis_min_dCtgTheta);
          etaProfile.SetMaximum(vis_max_dCtgTheta);
          etaCanvas->SetLogy();
          etaProfile.SetLineColor(momentumColor(myColor));
          etaProfile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          ctgThetaProfile.SetMarkerStyle(markerStyle);
          ctgThetaProfile.SetMarkerSize(markerSize);
          ctgThetaCanvas->SetFillColor(color_plot_background);
          etaProfile.SetMarkerStyle(markerStyle);
          etaProfile.SetMarkerSize(markerSize);
          etaCanvas->SetFillColor(color_plot_background);
          if (ctgThetaGraph.GetN() > 0) {
            ctgThetaCanvas->cd();
            ctgThetaProfile.Draw(plotOption.c_str());
            etaCanvas->cd();
            etaProfile.Draw(plotOption.c_str());
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
          z0Canvas->SetLogy();
          z0Profile.SetLineColor(momentumColor(myColor));
          z0Profile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          z0Profile.SetMarkerStyle(markerStyle);
          z0Profile.SetMarkerSize(markerSize);
          z0Canvas->SetFillColor(color_plot_background);
          if (z0Graph.GetN() > 0) {
            z0Canvas->cd();
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
          pCanvas->SetLogy();
          pProfile.SetLineColor(momentumColor(myColor));
          pProfile.SetMarkerColor(momentumColor(myColor));
          myColor++;
          pProfile.SetMarkerStyle(markerStyle);
          pProfile.SetMarkerSize(markerSize);
          pCanvas->SetFillColor(color_plot_background);
          if (pGraph.GetN() > 0) {
            pCanvas->cd();
            pProfile.Draw(plotOption.c_str());
            plotOption = "p same";
          }
        }
        RootWImage& linearMomentumImage = myContent->addImage(std::move(linearMomentumCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        linearMomentumImage.setComment("Transverse momentum resolution vs. eta (linear scale)");
        linearMomentumImage.setName(Form("linptres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& momentumImage = myContent->addImage(std::move(momentumCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        momentumImage.setComment("Transverse momentum resolution vs. eta");
        momentumImage.setName(Form("ptres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& distanceImage = myContent->addImage(std::move(distanceCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        distanceImage.setComment("Distance of closest approach resolution vs. eta");
        distanceImage.setName(Form("dxyres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& angleImage = myContent->addImage(std::move(angleCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        angleImage.setComment("Angle resolution vs. eta");
        angleImage.setName(Form("phires_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& ctgThetaImage = myContent->addImage(std::move(ctgThetaCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        ctgThetaImage.setComment("CtgTheta resolution vs. eta");
        ctgThetaImage.setName(Form("cotThetares_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& etaImage = myContent->addImage(std::move(etaCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        etaImage.setComment("Eta resolution vs. eta");
        etaImage.setName(Form("etares_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
	RootWImage& z0Image = myContent->addImage(std::move(z0Canvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
        z0Image.setComment("z0 resolution vs. eta");
        z0Image.setName(Form("dzres_%s_%s",additionalTag.c_str(), scenarioStr.c_str()));
        RootWImage& pImage = myContent->addImage(std::move(pCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
      //double verticalScale=1;

      pageTitle              += " ("+tag+")";
      additionalSummaryTag    = "_"+tag+"_";
      //verticalScale           = 10;
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
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::LGraph_Pt       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::BetaGraph_Pt    , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::OmegaGraph_Pt   , tag).empty() &&
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

          std::unique_ptr<TCanvas> linMomCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> logMomCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> d0Canvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> phiCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> ctgThetaCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> etaCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> z0Canvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> pCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> lCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> betaCanvas_Pt(new TCanvas());
          std::unique_ptr<TCanvas> omegaCanvas_Pt(new TCanvas());

          // Default attributes
          int myColor            = 0;
          int nBins              = insur::vis_n_bins;
          int markerStyle        = 21;
          double markerSize      = 1.;
          double lineWidth       = 2.;
          std::string plotOption = "";

          linMomCanvas_Pt->SetGrid(1,1);
          logMomCanvas_Pt->SetGrid(1,1);
          d0Canvas_Pt->SetGrid(1,1);
          phiCanvas_Pt->SetGrid(1,1);
          ctgThetaCanvas_Pt->SetGrid(1,1);
          etaCanvas_Pt->SetGrid(1,1);
          z0Canvas_Pt->SetGrid(1,1);
          pCanvas_Pt->SetGrid(1,1);
          lCanvas_Pt->SetGrid(1,1);
          betaCanvas_Pt->SetGrid(1,1);
          omegaCanvas_Pt->SetGrid(1,1);

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
            linMomCanvas_Pt->SetLogy(0);
            logMomCanvas_Pt->SetLogy(1);
            linMomCanvas_Pt->SetFillColor(color_plot_background);
            logMomCanvas_Pt->SetFillColor(color_plot_background);

            momentumProfile.SetLineColor(momentumColor(myColor));
            momentumProfile.SetMarkerColor(momentumColor(myColor));
            momentumProfile.SetLineWidth(lineWidth);
            myColor++;
            momentumProfile.SetMarkerStyle(markerStyle);
            momentumProfile.SetMarkerSize(markerSize);

            if (momentumGraph.GetN()>0) {
              linMomCanvas_Pt->cd();
              momentumProfile.Draw(plotOption.c_str());
              logMomCanvas_Pt->cd();
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
            pCanvas_Pt->SetLogy();
            pCanvas_Pt->SetFillColor(color_plot_background);

            pProfile.SetLineColor(momentumColor(myColor));
            pProfile.SetMarkerColor(momentumColor(myColor));
            pProfile.SetLineWidth(lineWidth);
            myColor++;
            pProfile.SetMarkerStyle(markerStyle);
            pProfile.SetMarkerSize(markerSize);

            if (pGraph.GetN() > 0) {
              pCanvas_Pt->cd();
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
            d0Canvas_Pt->SetLogy();
            d0Canvas_Pt->SetFillColor(color_plot_background);

            distanceProfile.SetLineColor(momentumColor(myColor));
            distanceProfile.SetMarkerColor(momentumColor(myColor));
            distanceProfile.SetLineWidth(lineWidth);
            myColor++;
            distanceProfile.SetMarkerStyle(markerStyle);
            distanceProfile.SetMarkerSize(markerSize);

            if (distanceGraph.GetN()>0) {
              d0Canvas_Pt->cd();
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
            phiCanvas_Pt->SetLogy();
            phiCanvas_Pt->SetFillColor(color_plot_background);

            angleProfile.SetLineColor(momentumColor(myColor));
            angleProfile.SetMarkerColor(momentumColor(myColor));
            angleProfile.SetLineWidth(lineWidth);
            myColor++;
            angleProfile.SetMarkerStyle(markerStyle);
            angleProfile.SetMarkerSize(markerSize);

            if (angleGraph.GetN() > 0) {
              phiCanvas_Pt->cd();
              angleProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw ctgTheta and eta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::CtgthetaGraph_Pt | idealMaterial, tag)) {

            const TGraph& ctgThetaGraph = mapel.second;
            TProfile& ctgThetaProfile   = newProfile(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
	    TProfile& etaProfile        = newProfile_timesSin(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
	    etaProfile.SetTitle("Pseudorapidity error - const P_{T} across #eta;#eta;#delta #eta");

            ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);
            ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);
            ctgThetaCanvas_Pt->SetLogy();
            ctgThetaCanvas_Pt->SetFillColor(color_plot_background);
            etaProfile.SetMinimum(vis_min_dCtgTheta);
            etaProfile.SetMaximum(vis_max_dCtgTheta);
            etaCanvas_Pt->SetLogy();
            etaCanvas_Pt->SetFillColor(color_plot_background);

            ctgThetaProfile.SetLineColor(momentumColor(myColor));
            ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
            ctgThetaProfile.SetLineWidth(lineWidth);
            etaProfile.SetLineColor(momentumColor(myColor));
            etaProfile.SetMarkerColor(momentumColor(myColor));
            etaProfile.SetLineWidth(lineWidth);
            myColor++;
            ctgThetaProfile.SetMarkerStyle(markerStyle);
            ctgThetaProfile.SetMarkerSize(markerSize);
            etaProfile.SetMarkerStyle(markerStyle);
            etaProfile.SetMarkerSize(markerSize);

            if (ctgThetaGraph.GetN() > 0) {
              ctgThetaCanvas_Pt->cd();
              ctgThetaProfile.Draw(plotOption.c_str());
              etaCanvas_Pt->cd();
              etaProfile.Draw(plotOption.c_str());
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
            z0Canvas_Pt->SetLogy();
            z0Canvas_Pt->SetFillColor(color_plot_background);

            z0Profile.SetLineColor(momentumColor(myColor));
            z0Profile.SetMarkerColor(momentumColor(myColor));
            z0Profile.SetLineWidth(lineWidth);
            myColor++;
            z0Profile.SetMarkerStyle(markerStyle);
            z0Profile.SetMarkerSize(markerSize);
            z0Canvas_Pt->SetFillColor(color_plot_background);

            if (z0Graph.GetN() > 0) {
              z0Canvas_Pt->cd();
              z0Profile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          // Draw l0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::LGraph_Pt | idealMaterial, tag)) {

            const TGraph& LGraph = mapel.second;
            TProfile& lProfile   = newProfile(LGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            lProfile.SetMinimum(vis_min_dL);
            lProfile.SetMaximum(vis_max_dL);//*verticalScale);
            lCanvas_Pt->SetLogy();
            lCanvas_Pt->SetFillColor(color_plot_background);

            lProfile.SetLineColor(momentumColor(myColor));
            lProfile.SetMarkerColor(momentumColor(myColor));
            lProfile.SetLineWidth(lineWidth);
            myColor++;
            lProfile.SetMarkerStyle(markerStyle);
            lProfile.SetMarkerSize(markerSize);
            lCanvas_Pt->SetFillColor(color_plot_background);

            if (LGraph.GetN() > 0) {
              lCanvas_Pt->cd();
              lProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          // Draw beta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::BetaGraph_Pt | idealMaterial, tag)) {

            const TGraph& BetaGraph = mapel.second;
            TProfile& betaProfile   = newProfile(BetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            betaProfile.SetMinimum(vis_min_beta);
            betaProfile.SetMaximum(vis_max_beta);
            betaCanvas_Pt->SetLogy(0);
            betaCanvas_Pt->SetFillColor(color_plot_background);

            betaProfile.SetLineColor(momentumColor(myColor));
            betaProfile.SetMarkerColor(momentumColor(myColor));
            betaProfile.SetLineWidth(lineWidth);
            myColor++;
            betaProfile.SetMarkerStyle(markerStyle);
            betaProfile.SetMarkerSize(markerSize);
            betaCanvas_Pt->SetFillColor(color_plot_background);

            if (BetaGraph.GetN() > 0) {
              betaCanvas_Pt->cd();
              betaProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          // Draw omega
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::OmegaGraph_Pt | idealMaterial, tag)) {

            const TGraph& OmegaGraph = mapel.second;
            TProfile& omegaProfile   = newProfile(OmegaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            omegaProfile.SetMinimum(vis_min_omega);
            omegaProfile.SetMaximum(vis_max_omega);
            omegaCanvas_Pt->SetLogy(0);
            omegaCanvas_Pt->SetFillColor(color_plot_background);

            omegaProfile.SetLineColor(momentumColor(myColor));
            omegaProfile.SetMarkerColor(momentumColor(myColor));
            omegaProfile.SetLineWidth(lineWidth);
            myColor++;
            omegaProfile.SetMarkerStyle(markerStyle);
            omegaProfile.SetMarkerSize(markerSize);
            omegaCanvas_Pt->SetFillColor(color_plot_background);

            if (OmegaGraph.GetN() > 0) {
              omegaCanvas_Pt->cd();
              omegaProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          RootWImage& linMomImage_Pt = myContent->addImage(std::move(linMomCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          linMomImage_Pt.setComment("Transverse momentum resolution vs. "+etaLetter+" (linear scale) - const Pt across "+etaLetter);
          linMomImage_Pt.setName(Form("linptres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& logMomImage_Pt = myContent->addImage(std::move(logMomCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          logMomImage_Pt.setComment("Transverse momentum resolution vs. "+etaLetter+" (log scale) - const Pt across "+etaLetter);
          logMomImage_Pt.setName(Form("ptres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& pImage_Pt = myContent->addImage(std::move(pCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          pImage_Pt.setComment("Momentum resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          pImage_Pt.setName(Form("pres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& d0Image_Pt = myContent->addImage(std::move(d0Canvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          d0Image_Pt.setComment("d0 resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          d0Image_Pt.setName(Form("dxyres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& z0Image_Pt = myContent->addImage(std::move(z0Canvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          z0Image_Pt.setComment("z0 resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          z0Image_Pt.setName(Form("dzres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& lImage_Pt = myContent->addImage(std::move(lCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          lImage_Pt.setComment("L resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          lImage_Pt.setName(Form("lres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& betaImage_Pt = myContent->addImage(std::move(betaCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          betaImage_Pt.setComment("Relative influence of #sigma_{d0} (max for #beta=#pi/2) and #sigma_{z0} (max for #beta=0) - const Pt across "+etaLetter);
          betaImage_Pt.setName(Form("beta_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& omegaImage_Pt = myContent->addImage(std::move(omegaCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          omegaImage_Pt.setComment("Pixel aspect ratio optimzation (for fixed pixel area). #Omega=#pi/2 means you need longer pixels - const Pt across "+etaLetter);
          omegaImage_Pt.setName(Form("omega_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& phiImage_Pt = myContent->addImage(std::move(phiCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          phiImage_Pt.setComment("Angle resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          phiImage_Pt.setName(Form("phires_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& ctgThetaImage_Pt = myContent->addImage(std::move(ctgThetaCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          ctgThetaImage_Pt.setComment("Ctg("+thetaLetter+") resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          ctgThetaImage_Pt.setName(Form("cotThetares_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& etaImage_Pt = myContent->addImage(std::move(etaCanvas_Pt), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          etaImage_Pt.setComment(etaLetter+" resolution vs. "+etaLetter+" - const Pt across "+etaLetter);
          etaImage_Pt.setName(Form("etares_%s_%s", tag.c_str(), scenarioStr.c_str()));
        }

        // Draw case II with const P across eta
        if (!gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::RhoGraph_P     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::DGraph_P       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::Z0Graph_P      , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PhiGraph_P     , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::CtgthetaGraph_P, tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::PGraph_P       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::LGraph_P       , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::BetaGraph_P    , tag).empty() &&
            !gb.getTaggedGraphs(GraphBag::RealGraph | GraphBag::OmegaGraph_P   , tag).empty()) {

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

          std::unique_ptr<TCanvas> linMomCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> logMomCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> d0Canvas_P(new TCanvas());
          std::unique_ptr<TCanvas> phiCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> ctgThetaCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> etaCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> z0Canvas_P(new TCanvas());
          std::unique_ptr<TCanvas> pCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> lCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> betaCanvas_P(new TCanvas());
          std::unique_ptr<TCanvas> omegaCanvas_P(new TCanvas());

          // Default attributes
          int myColor            = 0;
          int nBins              = insur::vis_n_bins;
          int markerStyle        = 21;
          double markerSize      = 1.;
          double lineWidth       = 2.;
          std::string plotOption = "";

          linMomCanvas_P->SetGrid(1,1);
          logMomCanvas_P->SetGrid(1,1);
          d0Canvas_P->SetGrid(1,1);
          phiCanvas_P->SetGrid(1,1);
          ctgThetaCanvas_P->SetGrid(1,1);
          etaCanvas_P->SetGrid(1,1);
          z0Canvas_P->SetGrid(1,1);
          pCanvas_P->SetGrid(1,1);
          lCanvas_P->SetGrid(1,1);
          betaCanvas_P->SetGrid(1,1);
          omegaCanvas_P->SetGrid(1,1);

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
            linMomCanvas_P->SetLogy(0);
            logMomCanvas_P->SetLogy(1);
            linMomCanvas_P->SetFillColor(color_plot_background);
            logMomCanvas_P->SetFillColor(color_plot_background);

            momentumProfile.SetLineColor(momentumColor(myColor));
            momentumProfile.SetMarkerColor(momentumColor(myColor));
            momentumProfile.SetLineWidth(lineWidth);
            myColor++;
            momentumProfile.SetMarkerStyle(markerStyle);
            momentumProfile.SetMarkerSize(markerSize);

            if (momentumGraph.GetN()>0) {
              linMomCanvas_P->cd();
              momentumProfile.Draw(plotOption.c_str());
              logMomCanvas_P->cd();
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
            pCanvas_P->SetLogy();
            pCanvas_P->SetFillColor(color_plot_background);

            pProfile.SetLineColor(momentumColor(myColor));
            pProfile.SetMarkerColor(momentumColor(myColor));
            pProfile.SetLineWidth(lineWidth);
            myColor++;
            pProfile.SetMarkerStyle(markerStyle);
            pProfile.SetMarkerSize(markerSize);

            if (pGraph.GetN() > 0) {
              pCanvas_P->cd();
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
            d0Canvas_P->SetLogy();
            d0Canvas_P->SetFillColor(color_plot_background);

            distanceProfile.SetLineColor(momentumColor(myColor));
            distanceProfile.SetMarkerColor(momentumColor(myColor));
            distanceProfile.SetLineWidth(lineWidth);
            myColor++;
            distanceProfile.SetMarkerStyle(markerStyle);
            distanceProfile.SetMarkerSize(markerSize);

            if (distanceGraph.GetN()>0) {
              d0Canvas_P->cd();
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
            phiCanvas_P->SetLogy();
            phiCanvas_P->SetFillColor(color_plot_background);

            angleProfile.SetLineColor(momentumColor(myColor));
            angleProfile.SetMarkerColor(momentumColor(myColor));
            angleProfile.SetLineWidth(lineWidth);
            myColor++;
            angleProfile.SetMarkerStyle(markerStyle);
            angleProfile.SetMarkerSize(markerSize);

            if (angleGraph.GetN() > 0) {
              phiCanvas_P->cd();
              angleProfile.Draw(plotOption.c_str());
              plotOption = "same";
            }
          }
          // Draw ctgTheta and eta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::CtgthetaGraph_P | idealMaterial, tag)) {

            const TGraph& ctgThetaGraph = mapel.second;
            TProfile& ctgThetaProfile   = newProfile(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
            TProfile& etaProfile        = newProfile_timesSin(ctgThetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);
	    etaProfile.SetTitle("Pseudorapidity error - const P across #eta;#eta;#delta #eta");

            ctgThetaProfile.SetMinimum(vis_min_dCtgTheta);
            ctgThetaProfile.SetMaximum(vis_max_dCtgTheta);
            ctgThetaCanvas_P->SetLogy();
            ctgThetaCanvas_P->SetFillColor(color_plot_background);
            etaProfile.SetMinimum(vis_min_dCtgTheta);
            etaProfile.SetMaximum(vis_max_dCtgTheta);
            etaCanvas_P->SetLogy();
            etaCanvas_P->SetFillColor(color_plot_background);

            ctgThetaProfile.SetLineColor(momentumColor(myColor));
            ctgThetaProfile.SetMarkerColor(momentumColor(myColor));
            ctgThetaProfile.SetLineWidth(lineWidth);
            etaProfile.SetLineColor(momentumColor(myColor));
            etaProfile.SetMarkerColor(momentumColor(myColor));
            etaProfile.SetLineWidth(lineWidth);
            myColor++;
            ctgThetaProfile.SetMarkerStyle(markerStyle);
            ctgThetaProfile.SetMarkerSize(markerSize);
            etaProfile.SetMarkerStyle(markerStyle);
            etaProfile.SetMarkerSize(markerSize);

            if (ctgThetaGraph.GetN() > 0) {
              ctgThetaCanvas_P->cd();
              ctgThetaProfile.Draw(plotOption.c_str());
              etaCanvas_P->cd();
              etaProfile.Draw(plotOption.c_str());
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
            z0Canvas_P->SetLogy();
            z0Canvas_P->SetFillColor(color_plot_background);

            z0Profile.SetLineColor(momentumColor(myColor));
            z0Profile.SetMarkerColor(momentumColor(myColor));
            z0Profile.SetLineWidth(lineWidth);
            myColor++;
            z0Profile.SetMarkerStyle(markerStyle);
            z0Profile.SetMarkerSize(markerSize);
            z0Canvas_P->SetFillColor(color_plot_background);

            if (z0Graph.GetN() > 0) {
              z0Canvas_P->cd();
              z0Profile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }
          plotOption = "";
          myColor=0;

          // Draw l0
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::LGraph_P | idealMaterial, tag)) {

            const TGraph& LGraph = mapel.second;
            TProfile& lProfile   = newProfile(LGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            lProfile.SetMinimum(vis_min_dL);
            lProfile.SetMaximum(vis_max_dL);//*verticalScale);
            lCanvas_P->SetLogy();
            lCanvas_P->SetFillColor(color_plot_background);

            lProfile.SetLineColor(momentumColor(myColor));
            lProfile.SetMarkerColor(momentumColor(myColor));
            lProfile.SetLineWidth(lineWidth);
            myColor++;
            lProfile.SetMarkerStyle(markerStyle);
            lProfile.SetMarkerSize(markerSize);
            lCanvas_P->SetFillColor(color_plot_background);

            if (LGraph.GetN() > 0) {
              lCanvas_P->cd();
              lProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          // Draw beta
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::BetaGraph_P | idealMaterial, tag)) {

            const TGraph& BetaGraph = mapel.second;
            TProfile& betaProfile   = newProfile(BetaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            betaProfile.SetMinimum(vis_min_beta);
            betaProfile.SetMaximum(vis_max_beta);
            betaCanvas_P->SetLogy(0);
            betaCanvas_P->SetFillColor(color_plot_background);

            betaProfile.SetLineColor(momentumColor(myColor));
            betaProfile.SetMarkerColor(momentumColor(myColor));
            betaProfile.SetLineWidth(lineWidth);
            myColor++;
            betaProfile.SetMarkerStyle(markerStyle);
            betaProfile.SetMarkerSize(markerSize);
            betaCanvas_P->SetFillColor(color_plot_background);

            if (BetaGraph.GetN() > 0) {
              betaCanvas_P->cd();
              betaProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          // Draw omega
          plotOption = "";
          myColor    = 0;
          for (const auto& mapel : gb.getTaggedGraphs(GraphBag::OmegaGraph_P | idealMaterial, tag)) {

            const TGraph& OmegaGraph = mapel.second;
            TProfile& omegaProfile   = newProfile(OmegaGraph, 0, analyzer.getEtaMaxTracker(), 1, nBins);

            omegaProfile.SetMinimum(vis_min_omega);
            omegaProfile.SetMaximum(vis_max_omega);
            omegaCanvas_P->SetLogy(0);
            omegaCanvas_P->SetFillColor(color_plot_background);

            omegaProfile.SetLineColor(momentumColor(myColor));
            omegaProfile.SetMarkerColor(momentumColor(myColor));
            omegaProfile.SetLineWidth(lineWidth);
            myColor++;
            omegaProfile.SetMarkerStyle(markerStyle);
            omegaProfile.SetMarkerSize(markerSize);
            omegaCanvas_P->SetFillColor(color_plot_background);

            if (OmegaGraph.GetN() > 0) {
              omegaCanvas_P->cd();
              omegaProfile.Draw(plotOption.c_str());
              plotOption = "p same";
            }
          }

          RootWImage& linMomImage_P = myContent->addImage(std::move(linMomCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          linMomImage_P.setComment("Transverse momentum resolution vs. "+etaLetter+" (linear scale) - const P across "+etaLetter);
          linMomImage_P.setName(Form("linptres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& logMomImage_P = myContent->addImage(std::move(logMomCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          logMomImage_P.setComment("Transverse momentum resolution vs. "+etaLetter+" (log scale) - const P across "+etaLetter);
          logMomImage_P.setName(Form("ptres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& pImage_P = myContent->addImage(std::move(pCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          pImage_P.setComment("Momentum resolution vs. "+etaLetter+" - const P across "+etaLetter);
          pImage_P.setName(Form("pres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& d0Image_P = myContent->addImage(std::move(d0Canvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          d0Image_P.setComment("d0 resolution vs. "+etaLetter+" - const P across "+etaLetter);
          d0Image_P.setName(Form("dxyres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& z0Image_P = myContent->addImage(std::move(z0Canvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          z0Image_P.setComment("z0 resolution vs. "+etaLetter+" - const P across "+etaLetter);
          z0Image_P.setName(Form("dzres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& lImage_P = myContent->addImage(std::move(lCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          lImage_P.setComment("L resolution vs. "+etaLetter+" - const P across "+etaLetter);
          lImage_P.setName(Form("lres_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& betaImage_P = myContent->addImage(std::move(betaCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          betaImage_P.setComment("Relative influence of #sigma_{d0} (max for #beta=#pi/2) and #sigma_{z0} (max for #beta=0) - const P across "+etaLetter);
          betaImage_P.setName(Form("beta_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& omegaImage_P = myContent->addImage(std::move(omegaCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          omegaImage_P.setComment("Pixel aspect ratio optimzation (for fixed pixel area). #Omega=#pi/2 means you need longer pixels - const P across "+etaLetter);
          omegaImage_P.setName(Form("omega_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& phiImage_P = myContent->addImage(std::move(phiCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          phiImage_P.setComment("Angle resolution vs. "+etaLetter+" - const P across "+etaLetter);
          phiImage_P.setName(Form("phires_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& ctgThetaImage_P = myContent->addImage(std::move(ctgThetaCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          ctgThetaImage_P.setComment("Ctg("+thetaLetter+") resolution vs. "+etaLetter+" - const P across "+etaLetter);
          ctgThetaImage_P.setName(Form("cotThetares_%s_%s", tag.c_str(), scenarioStr.c_str()));

          RootWImage& etaImage_P = myContent->addImage(std::move(etaCanvas_P), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
          etaImage_P.setComment(etaLetter+") resolution vs. "+etaLetter+" - const P across "+etaLetter);
          etaImage_P.setName(Form("etares_%s_%s", tag.c_str(), scenarioStr.c_str()));
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
      plotNames.push_back(deltaLetter+"c"+tauLetter+" ["+muLetter+"m]:");

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
      fillTaggedPlotMap(gb, plotNames[6], GraphBag::LGraph_Pt, tag, myPlotMap_Pt);


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
        if (iMom<color_names.size()) color = color_names[iMom];

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
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3 || iItem==6) myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionStd);
              else                                                          myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionHigh);
              myTable->setColor(2, baseColumn+j, myColor);
            }
            if (averagesIdeal.size() > j) {

              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3 || iItem==6) myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionStd);
              else                                                          myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionHigh);
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
      fillTaggedPlotMap(gb, plotNames[6], GraphBag::LGraph_P       , tag, myPlotMap_P);


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
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3 || iItem==6) myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionStd);
              else                                                          myTable->setContent(2, baseColumn+j,averagesReal[j],tableResolutionPrecisionHigh);
              myTable->setColor(2, baseColumn+j, myColor);
            }
            if (averagesIdeal.size() > j) {

              // Check item iterated and set precision
              int iItem = plotNameIt - plotNames.begin();
              if (iItem==0 || iItem==1 || iItem==2 || iItem==3 || iItem==6) myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionStd);
              else                                                          myTable->setContent(3, baseColumn+j,averagesIdeal[j],tableResolutionPrecisionHigh);
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

  bool Vizard::patternRecoSummary(Analyzer& a, mainConfigHandler& mainConfig, RootWSite& site) {

    bool isVisOK = true;

    std::string pageTitle   = "PatternReco";
    std::string pageAddress = "patternReco.html";

    RootWPage* myPage = new RootWPage(pageTitle);
    myPage->setAddress(pageAddress);
    site.addPage(myPage);

    // Canvases
    gStyle->SetGridStyle(style_grid);
    gStyle->SetOptStat(0);

    //
    // Pt option
    RootWContent& myContentPlotsPt = myPage->addContent("Pt: Track purity & probability of contamination by bkg hits when propagating track through tracker", true);

    // Summary
    RootWInfo& myInfoPt = myContentPlotsPt.addInfo("Number of minimum bias events per bunch crossing");
    myInfoPt.setValue(SimParms::getInstance().numMinBiasEvents(), minimumBiasPrecision);

    // a) Bkg contamination probability -> In-Out approach
    std::unique_ptr<TCanvas> canvasPtBkgContInOut(new TCanvas("canvasPtBkgContInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
    canvasPtBkgContInOut->SetGrid(1,1);
    canvasPtBkgContInOut->SetLogy(0);
    canvasPtBkgContInOut->SetFillColor(color_plot_background);
    canvasPtBkgContInOut->SetObjectStat(false);

    // For each momentum/transverse momentum compute
    int iMomentum = 0;
    TLegend* legendPtInOut = new TLegend(0.11,0.66,0.36,0.89,"p_{T} options (InOut+IP, full TRK):");
    legendPtInOut->SetTextSize(0.025);

    for (const auto& pIter : mainConfig.getMomenta()) {

      a.hisPatternRecoInOutPt[iMomentum]->SetNameTitle(std::string("hisPatternRecoInOutPt"+any2str(pIter/Units::GeV)+"GeV").c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
      a.hisPatternRecoInOutPt[iMomentum]->SetLineWidth(2.);
      a.hisPatternRecoInOutPt[iMomentum]->SetMarkerStyle(21);
      a.hisPatternRecoInOutPt[iMomentum]->SetMarkerSize(1.);
      a.hisPatternRecoInOutPt[iMomentum]->SetLineColor(momentumColor(iMomentum));
      a.hisPatternRecoInOutPt[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

      if (iMomentum==0) a.hisPatternRecoInOutPt[iMomentum]->Draw("PE1");
      else              a.hisPatternRecoInOutPt[iMomentum]->Draw("SAME PE1");

      a.hisPatternRecoInOutPt[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
      legendPtInOut->AddEntry(a.hisPatternRecoInOutPt[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

      iMomentum++;
    }
    legendPtInOut->Draw("SAME");

    RootWImage& myImagePtBkgContInOut = myContentPlotsPt.addImage(std::move(canvasPtBkgContInOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtBkgContInOut.setComment("In-Out approach for pT: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
    myImagePtBkgContInOut.setName("bkg_pt_pContam_inout");

    // b) Bkg contamination probability -> OutIn approach
    std::unique_ptr<TCanvas> canvasPtBkgContOutIn(new TCanvas("canvasPtBkgContOutIn","",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
    canvasPtBkgContOutIn->SetGrid(1,1);
    canvasPtBkgContOutIn->SetLogy(0);
    canvasPtBkgContOutIn->SetFillColor(color_plot_background);
    canvasPtBkgContOutIn->SetObjectStat(false);

    // For each momentum/transverse momentum compute
    iMomentum = 0;
    TLegend* legendPtOutIn = new TLegend(0.11,0.66,0.36,0.89,"p_{T} options (OutIn, full TRK):");
    legendPtOutIn->SetTextSize(0.025);

    for (const auto& pIter : mainConfig.getMomenta()) {

      a.hisPatternRecoOutInPt[iMomentum]->SetNameTitle(std::string("hisPatternRecoOutInPt"+any2str(pIter/Units::GeV)+"GeV").c_str(),"Out-In: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
      a.hisPatternRecoOutInPt[iMomentum]->SetLineWidth(2.);
      a.hisPatternRecoOutInPt[iMomentum]->SetMarkerStyle(21);
      a.hisPatternRecoOutInPt[iMomentum]->SetMarkerSize(1.);
      a.hisPatternRecoOutInPt[iMomentum]->SetLineColor(momentumColor(iMomentum));
      a.hisPatternRecoOutInPt[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

      if (iMomentum==0) a.hisPatternRecoOutInPt[iMomentum]->Draw("PE1");
      else              a.hisPatternRecoOutInPt[iMomentum]->Draw("SAME PE1");

      a.hisPatternRecoOutInPt[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
      legendPtOutIn->AddEntry(a.hisPatternRecoOutInPt[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

      iMomentum++;
    }
    legendPtOutIn->Draw("SAME");

    RootWImage& myImagePtBkgContOutIn = myContentPlotsPt.addImage(std::move(canvasPtBkgContOutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePtBkgContOutIn.setComment("Out-In approach for pT: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
    myImagePtBkgContOutIn.setName("bkg_pt_pContam_outin");

    //
    // Detail on pt in-out studies

    // a) D0
    RootWContent& myContentPlotsPtD0InOut = myPage->addContent("Pt in-out - D: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitDProjInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtD0InOut(new TCanvas(std::string("canvasPtD0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtD0InOut->SetGrid(1,1);
      canvasPtD0InOut->SetLogy(0);
      canvasPtD0InOut->SetFillColor(color_plot_background);
      canvasPtD0InOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtD0InOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPtInOut->SetX1(0.64);
      legendPtInOut->SetX2(0.89);
      legendPtInOut->SetEntryLabel("p_{T} options (InOut+IP, full TRK):");
      legendPtInOut->Draw("SAME");

      RootWImage& myImagePtD0InOut = myContentPlotsPtD0InOut.addImage(std::move(canvasPtD0InOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtD0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in in-out approach.").c_str());
      myImagePtD0InOut.setName(std::string("bkg_pt_d0_inout_"+name).c_str());
    }

    // b) Z0
    RootWContent& myContentPlotsPtZ0InOut = myPage->addContent("Pt in-out - Z: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitZProjInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtZ0InOut(new TCanvas(std::string("canvasPtZ0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtZ0InOut->SetGrid(1,1);
      canvasPtZ0InOut->SetLogy(0);
      canvasPtZ0InOut->SetFillColor(color_plot_background);
      canvasPtZ0InOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtZ0InOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPtInOut->Draw("SAME");

      RootWImage& myImagePtZ0InOut = myContentPlotsPtZ0InOut.addImage(std::move(canvasPtZ0InOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtZ0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in in-out approach.").c_str());
      myImagePtZ0InOut.setName(std::string("bkg_pt_z0_inout_"+name).c_str());
    }

    // c) pContamination
    RootWContent& myContentPlotsPtProbContamInOut = myPage->addContent("Pt in-out - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitProbContamInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtProbContamInOut(new TCanvas(std::string("canvasPtProbContamInOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtProbContamInOut->SetGrid(1,1);
      canvasPtProbContamInOut->SetLogy(0);
      canvasPtProbContamInOut->SetFillColor(color_plot_background);
      canvasPtProbContamInOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtProbContamInOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,0.1);

        iMomentum++;
      }
      legendPtInOut->SetX1(0.11);
      legendPtInOut->SetX2(0.36);
      legendPtInOut->Draw("SAME");

      RootWImage& myImagePtProbContamInOut = myContentPlotsPtProbContamInOut.addImage(std::move(canvasPtProbContamInOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtProbContamInOut.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in in-out approach.").c_str());
      myImagePtProbContamInOut.setName(std::string("bkg_pt_pContam_inout_"+name).c_str());
    }

    //
    // Detail on pt out-in studies

    // a) D0
    RootWContent& myContentPlotsPtD0OutIn = myPage->addContent("Pt out-in - D: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitDProjOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtD0OutIn(new TCanvas(std::string("canvasPtD0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtD0OutIn->SetGrid(1,1);
      canvasPtD0OutIn->SetLogy(0);
      canvasPtD0OutIn->SetFillColor(color_plot_background);
      canvasPtD0OutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtD0OutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPtOutIn->SetX1(0.64);
      legendPtOutIn->SetX2(0.89);
      legendPtOutIn->SetEntryLabel("p_{T} options (OutIn, full TRK):");
      legendPtOutIn->Draw("SAME");

      RootWImage& myImagePtD0OutIn = myContentPlotsPtD0OutIn.addImage(std::move(canvasPtD0OutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtD0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in out-in approach.").c_str());
      myImagePtD0OutIn.setName(std::string("bkg_pt_d0_outin_"+name).c_str());
    }

    // b) Z0
    RootWContent& myContentPlotsPtZ0OutIn = myPage->addContent("Pt out-in - Z: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitZProjOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtZ0OutIn(new TCanvas(std::string("canvasPtZ0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtZ0OutIn->SetGrid(1,1);
      canvasPtZ0OutIn->SetLogy(0);
      canvasPtZ0OutIn->SetFillColor(color_plot_background);
      canvasPtZ0OutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtZ0OutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPtOutIn->Draw("SAME");

      RootWImage& myImagePtZ0OutIn = myContentPlotsPtZ0OutIn.addImage(std::move(canvasPtZ0OutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtZ0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in out-in approach.").c_str());
      myImagePtZ0OutIn.setName(std::string("bkg_pt_z0_outin_"+name).c_str());
    }

    // c) pContamination
    RootWContent& myContentPlotsPtProbContamOutIn = myPage->addContent("Pt out-in - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPtHitProbContamOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPtProbContamOutIn(new TCanvas(std::string("canvasPtProbContamOutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPtProbContamOutIn->SetGrid(1,1);
      canvasPtProbContamOutIn->SetLogy(0);
      canvasPtProbContamOutIn->SetFillColor(color_plot_background);
      canvasPtProbContamOutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPtProbContamOutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,0.1);

        iMomentum++;
      }
      legendPtOutIn->SetX1(0.11);
      legendPtOutIn->SetX2(0.36);
      legendPtOutIn->Draw("SAME");

      RootWImage& myImagePtProbContamOutIn = myContentPlotsPtProbContamOutIn.addImage(std::move(canvasPtProbContamOutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePtProbContamOutIn.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in out-in approach.").c_str());
      myImagePtProbContamOutIn.setName(std::string("bkg_pt_pContam_outin_"+name).c_str());
    }

    //
    // P option
    RootWContent& myContentPlotsP = myPage->addContent("P: Track purity & probability of contamination by bkg hits when propagating track through tracker", true);

    // Summary
    RootWInfo& myInfoP = myContentPlotsP.addInfo("Number of minimum bias events per bunch crossing");
    myInfoP.setValue(SimParms::getInstance().numMinBiasEvents(), minimumBiasPrecision);

    // a) Bkg contamination probability - In-Out approach
    std::unique_ptr<TCanvas> canvasPBkgContInOut(new TCanvas("canvasPBkgContInOut","",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
    canvasPBkgContInOut->SetGrid(1,1);
    canvasPBkgContInOut->SetLogy(0);
    canvasPBkgContInOut->SetFillColor(color_plot_background);
    canvasPBkgContInOut->SetObjectStat(false);

    // For each momentum/transverse momentum compute
    iMomentum = 0;
    TLegend* legendPInOut = new TLegend(0.11,0.66,0.36,0.89,"p options (InOut+IP, full TRK):");
    legendPInOut->SetTextSize(0.025);

    for (const auto& pIter : mainConfig.getMomenta()) {

      a.hisPatternRecoInOutP[iMomentum]->SetNameTitle(std::string("hisPBkgContInOut"+any2str(pIter/Units::GeV)+"GeV").c_str(),"In-Out: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
      a.hisPatternRecoInOutP[iMomentum]->SetLineWidth(2.);
      a.hisPatternRecoInOutP[iMomentum]->SetMarkerStyle(21);
      a.hisPatternRecoInOutP[iMomentum]->SetMarkerSize(1.);
      a.hisPatternRecoInOutP[iMomentum]->SetLineColor(momentumColor(iMomentum));
      a.hisPatternRecoInOutP[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

      if (iMomentum==0) a.hisPatternRecoInOutP[iMomentum]->Draw("PE1");
      else              a.hisPatternRecoInOutP[iMomentum]->Draw("SAME PE1");

      a.hisPatternRecoInOutP[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
      legendPInOut->AddEntry(a.hisPatternRecoInOutP[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

      iMomentum++;
    }
    legendPInOut->Draw("SAME");

    RootWImage& myImagePBkgContInOut = myContentPlotsP.addImage(std::move(canvasPBkgContInOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePBkgContInOut.setComment("In-Out approach for p: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
    myImagePBkgContInOut.setName("bkg_p_pContam_inout");

    // b) Bkg contamination probability - Out-In approach
    std::unique_ptr<TCanvas> canvasPBkgContOutIn(new TCanvas("canvasPBkgContOutIn","",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
    canvasPBkgContOutIn->SetGrid(1,1);
    canvasPBkgContOutIn->SetLogy(0);
    canvasPBkgContOutIn->SetFillColor(color_plot_background);
    canvasPBkgContOutIn->SetObjectStat(false);

    // For each momentum/transverse momentum compute
    iMomentum = 0;
    TLegend* legendPOutIn = new TLegend(0.11,0.66,0.36,0.89,"p options (OutIn, full TRK):");
    legendPOutIn->SetTextSize(0.025);

    for (const auto& pIter : mainConfig.getMomenta()) {

      a.hisPatternRecoOutInP[iMomentum]->SetNameTitle(std::string("hisPBkgContOutIn"+any2str(pIter/Units::GeV)+"GeV").c_str(),"Out-In: Bkg contamination prob. in 95% area of 2D error ellipse accumulated accross N layers;#eta;1 - #Pi_{i=1}^{N} (1-p^{i}_{bkg_95%})");
      a.hisPatternRecoOutInP[iMomentum]->SetLineWidth(2.);
      a.hisPatternRecoOutInP[iMomentum]->SetMarkerStyle(21);
      a.hisPatternRecoOutInP[iMomentum]->SetMarkerSize(1.);
      a.hisPatternRecoOutInP[iMomentum]->SetLineColor(momentumColor(iMomentum));
      a.hisPatternRecoOutInP[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

      if (iMomentum==0) a.hisPatternRecoOutInP[iMomentum]->Draw("PE1");
      else              a.hisPatternRecoOutInP[iMomentum]->Draw("SAME PE1");

      a.hisPatternRecoOutInP[iMomentum]->GetYaxis()->SetRangeUser(0,1.3);
      legendPOutIn->AddEntry(a.hisPatternRecoOutInP[iMomentum],std::string(any2str(pIter/Units::GeV)+"GeV").c_str(),"lp");

      iMomentum++;
    }
    legendPOutIn->Draw("SAME");

    RootWImage& myImagePBkgContOutIn = myContentPlotsP.addImage(std::move(canvasPBkgContOutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    myImagePBkgContOutIn.setComment("Out-In approach for p: Bkg contamination prob. in 95% area of 2D error ellipse accumulated across N layers");
    myImagePBkgContOutIn.setName("bkg_p_pContam_outin");

    //
    // Detail on p in-out studies

    // a) D0
    RootWContent& myContentPlotsPD0InOut = myPage->addContent("P in-out - D: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitDProjInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPD0InOut(new TCanvas(std::string("canvasPD0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPD0InOut->SetGrid(1,1);
      canvasPD0InOut->SetLogy(0);
      canvasPD0InOut->SetFillColor(color_plot_background);
      canvasPD0InOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPD0InOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPInOut->SetX1(0.64);
      legendPInOut->SetX2(0.89);
      legendPInOut->SetEntryLabel("p options (InOut+IP, full TRK):");
      legendPInOut->Draw("SAME");

      RootWImage& myImagePD0InOut = myContentPlotsPD0InOut.addImage(std::move(canvasPD0InOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePD0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in in-out approach.").c_str());
      myImagePD0InOut.setName(std::string("bkg_p_d0_inout_"+name).c_str());
    }

    // b) Z0
    RootWContent& myContentPlotsPZ0InOut = myPage->addContent("P in-out - Z: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitZProjInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPZ0InOut(new TCanvas(std::string("canvasPZ0InOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPZ0InOut->SetGrid(1,1);
      canvasPZ0InOut->SetLogy(0);
      canvasPZ0InOut->SetFillColor(color_plot_background);
      canvasPZ0InOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPZ0InOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPInOut->SetX1(0.64);
      legendPInOut->SetX2(0.89);
      legendPInOut->Draw("SAME");

      RootWImage& myImagePZ0InOut = myContentPlotsPZ0InOut.addImage(std::move(canvasPZ0InOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePZ0InOut.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in in-out approach.").c_str());
      myImagePZ0InOut.setName(std::string("bkg_p_z0_inout_"+name).c_str());
    }

    // c) pContamination
    RootWContent& myContentPlotsPProbContamInOut = myPage->addContent("P in-out - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitProbContamInOut) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPProbContamInOut(new TCanvas(std::string("canvasPProbContamInOut"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPProbContamInOut->SetGrid(1,1);
      canvasPProbContamInOut->SetLogy(0);
      canvasPProbContamInOut->SetFillColor(color_plot_background);
      canvasPProbContamInOut->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPProbContamInOut"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1.1);

        iMomentum++;
      }
      legendPInOut->SetX1(0.11);
      legendPInOut->SetX2(0.36);
      legendPInOut->Draw("SAME");

      RootWImage& myImagePProbContamInOut = myContentPlotsPProbContamInOut.addImage(std::move(canvasPProbContamInOut), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePProbContamInOut.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in in-out approach.").c_str());
      myImagePProbContamInOut.setName(std::string("bkg_p_pContam_inout_"+name).c_str());
    }
    //
    // Detail on p out-in studies

    // a) D0
    RootWContent& myContentPlotsPD0OutIn = myPage->addContent("P out-in - D: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitDProjOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPD0OutIn(new TCanvas(std::string("canvasPD0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPD0OutIn->SetGrid(1,1);
      canvasPD0OutIn->SetLogy(0);
      canvasPD0OutIn->SetFillColor(color_plot_background);
      canvasPD0OutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPD0OutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{R-#Phi} from previous layers/discs;#eta; #sigma_{R-#Phi} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPOutIn->SetX1(0.64);
      legendPOutIn->SetX2(0.89);
      legendPOutIn->SetEntryLabel("p options (OutIn+IP, full TRK):");
      legendPOutIn->Draw("SAME");

      RootWImage& myImagePD0OutIn = myContentPlotsPD0OutIn.addImage(std::move(canvasPD0OutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePD0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in R-Phi from previous layers/discs in out-in approach.").c_str());
      myImagePD0OutIn.setName(std::string("bkg_p_d0_outin_"+name).c_str());
    }

    // b) Z0
    RootWContent& myContentPlotsPZ0OutIn = myPage->addContent("P out-in - Z: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitZProjOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPZ0OutIn(new TCanvas(std::string("canvasPZ0OutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPZ0OutIn->SetGrid(1,1);
      canvasPZ0OutIn->SetLogy(0);
      canvasPZ0OutIn->SetFillColor(color_plot_background);
      canvasPZ0OutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPZ0OutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: an extrapolated #sigma_{Z} from previous layers/discs;#eta; #sigma_{Z} [#mum]").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1000);

        iMomentum++;
      }
      legendPOutIn->SetX1(0.64);
      legendPOutIn->SetX2(0.89);
      legendPOutIn->Draw("SAME");

      RootWImage& myImagePZ0OutIn = myContentPlotsPZ0OutIn.addImage(std::move(canvasPZ0OutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePZ0OutIn.setComment(std::string("Detector: "+name+" - an extrapolated sigma in Z from previous layers/discs in out-in approach.").c_str());
      myImagePZ0OutIn.setName(std::string("bkg_p_z0_outin_"+name).c_str());
    }

    // c) pContamination
    RootWContent& myContentPlotsPProbContamOutIn = myPage->addContent("P out-in - Bkg contamination prob.: Details of track purity & probability of contamination by bkg", false);
    for (auto const & iterMap : a.hisPHitProbContamOutIn) {

      // Get name -> remove first 2 artificial characters used to sort correctly map
      std::string name = iterMap.first;
      name = name.erase(0,2);

      int         iMomentum = 0;

      std::unique_ptr<TCanvas> canvasPProbContamOutIn(new TCanvas(std::string("canvasPProbContamOutIn"+name).c_str(),"",vis_std_canvas_sizeY,vis_min_canvas_sizeY));
      canvasPProbContamOutIn->SetGrid(1,1);
      canvasPProbContamOutIn->SetLogy(0);
      canvasPProbContamOutIn->SetFillColor(color_plot_background);
      canvasPProbContamOutIn->SetObjectStat(false);

      for (const auto& pIter : mainConfig.getMomenta()) {

        auto & profHis = iterMap.second;
        profHis[iMomentum]->SetNameTitle(std::string("canvasPProbContamOutIn"+name+any2str(pIter/Units::GeV)+"GeV").c_str(),std::string(name+" In-Out approach: bkg contamination prob. as error ellipse extrap. from previous layers/discs;#eta; probability").c_str());
        profHis[iMomentum]->SetLineWidth(2.);
        profHis[iMomentum]->SetMarkerStyle(21);
        profHis[iMomentum]->SetMarkerSize(1.);
        profHis[iMomentum]->SetLineColor(momentumColor(iMomentum));
        profHis[iMomentum]->SetMarkerColor(momentumColor(iMomentum));

        if (iMomentum==0) profHis[iMomentum]->Draw("PE1");
        else              profHis[iMomentum]->Draw("SAME PE1");

        profHis[iMomentum]->GetYaxis()->SetRangeUser(0,1.1);

        iMomentum++;
      }
      legendPOutIn->SetX1(0.11);
      legendPOutIn->SetX2(0.36);
      legendPOutIn->Draw("SAME");

      RootWImage& myImagePProbContamOutIn = myContentPlotsPProbContamOutIn.addImage(std::move(canvasPProbContamOutIn), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      myImagePProbContamOutIn.setComment(std::string("Detector: "+name+" - bkg contamination prob. as error ellipse extrap. from previous layers/discs in out-in approach.").c_str());
      myImagePProbContamOutIn.setName(std::string("bkg_p_pContam_outin_"+name).c_str());
    }

    return isVisOK;
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

    profileBag& aProfileBag = a.getProfileBag();
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
      std::unique_ptr<TCanvas> pointsCanvas(new TCanvas());
      pointsCanvas->SetGrid(1,1);
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
          tempSS << tempString.c_str() << any2str(myPt/Units::GeV,2); tempString = ", ";
        }

        npointsProfile.SetMinimum(1E-2);
        //npointsProfile.GetXaxis()->SetLimits(0, 2.4);
        npointsProfile.SetLineColor(Palette::color(myColor));
        npointsProfile.SetMarkerColor(Palette::color(myColor));
        npointsProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        npointsProfile.SetMarkerStyle(8);
        pointsCanvas->SetFillColor(color_plot_background);

        pointsCanvas->cd();
        // std::cerr << "About to draw plot " << myPt << std::endl; // debug
        npointsProfile.Draw(plotOption.c_str());
        //plotOption = "E6 same";
        plotOption = "E1 same";
        //plotOption = "same";
      }
      tempSS << " GeV";

      std::unique_ptr<TCanvas> pointsLogCanvas((TCanvas*)pointsCanvas->DrawClone());
      pointsLogCanvas->SetLogy();
           
      RootWImage& npointsImage = myContent.addImage(std::move(pointsCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      npointsImage.setComment(tempSS.str().c_str());
      npointsImage.setName("ntrigpoints");

      RootWImage& npointsLogImage = myContent.addImage(std::move(pointsLogCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS << " (log scale)";
      npointsLogImage.setComment(tempSS.str().c_str());
      npointsLogImage.setName("ntrigpointsLog");

      std::unique_ptr<TCanvas> fractionCanvas(new TCanvas());
      fractionCanvas->SetGrid(1,1);
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
          tempSS << tempString.c_str() << any2str(myPt/Units::GeV,2) ; tempString = ", ";
        }

        fractionProfile.SetMinimum(1E-2);
        fractionProfile.SetMaximum(100);
        fractionProfile.SetLineColor(Palette::color(myColor));
        fractionProfile.SetMarkerColor(Palette::color(myColor));
        fractionProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        fractionProfile.SetMarkerStyle(8);
        fractionCanvas->SetFillColor(color_plot_background);

        fractionCanvas->cd();
        // std::cerr << "About to draw fraction plot " << myPt << std::endl; // debug
        fractionProfile.Draw(plotOption.c_str());
        plotOption = "E1 same";
        //aValue = fractionProfile.GetBinContent(fractionProfile.GetMaximumBin());
        //if (aValue>maxy) maxy=aValue;
        //aValue = fractionProfile.GetBinContent(fractionProfile.GetMinimumBin());
        //if (aValue<miny) miny=aValue;
        //std::cerr << "Fraction plots between " << miny << " and " << maxy << std::endl;
      }
      tempSS << " GeV";

      std::unique_ptr<TCanvas> fractionLogCanvas((TCanvas*)fractionCanvas->DrawClone());
      fractionLogCanvas->SetLogy();

      RootWImage& fractionImage = myContent.addImage(std::move(fractionCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      fractionImage.setComment(tempSS.str().c_str());
      fractionImage.setName("fractiontrigpoints");
      
      RootWImage& fractionLogImage = myContent.addImage(std::move(fractionLogCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS << " (log scale)";
      fractionLogImage.setComment(tempSS.str().c_str());
      fractionLogImage.setName("fractiontrigpointsLog");


      std::unique_ptr<TCanvas> purityCanvas(new TCanvas());
      purityCanvas->SetGrid(1,1);
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
          tempSS << tempString.c_str() << any2str(myPt/Units::GeV,2) ; tempString = ", ";
        }

        purityProfile.SetMinimum(1E-2);
        purityProfile.SetMaximum(100);
        purityProfile.SetLineColor(Palette::color(myColor));
        purityProfile.SetMarkerColor(Palette::color(myColor));
        purityProfile.SetFillColor(Palette::color(myColor));
        myColor++;
        purityProfile.SetMarkerStyle(8);
        purityCanvas->SetFillColor(color_plot_background);

        purityCanvas->cd();
        // std::cerr << "About to draw purity plot " << myPt << std::endl; // debug
        purityProfile.Draw(plotOption.c_str());
        plotOption = "E1 same";
        //aValue = purityProfile.GetBinContent(purityProfile.GetMaximumBin());
        //if (aValue>maxy) maxy=aValue;
        //aValue = purityProfile.GetBinContent(purityProfile.GetMinimumBin());
        //if (aValue<miny) miny=aValue;
        //std::cerr << "Purity plots between " << miny << " and " << maxy << std::endl;
      }
      tempSS << " GeV";

      std::unique_ptr<TCanvas> purityLogCanvas((TCanvas*)purityCanvas->DrawClone());
      purityLogCanvas->SetLogy();

      RootWImage& purityImage = myContent.addImage(std::move(purityCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      purityImage.setComment(tempSS.str().c_str());
      purityImage.setName("puritytrigpoints");
      
      RootWImage& purityLogImage = myContent.addImage(std::move(purityLogCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
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
    mapBag& myMapBag = a.getMapBag();
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
        std::unique_ptr<TCanvas> myCanvas(new TCanvas());
        double myPt = it->first;
        if (myPt>maxPt) maxPt=myPt;
        TH2D& myMap = it->second;
        myCanvas->SetFillColor(color_plot_background);
        myCanvas->cd();

        // Actually plot the map
        myMap.SetMinimum(0);
        if (myPt<1.5) { // TODO: make this 1.5 a global constant (also in Analyzer)
          myMap.SetMaximum(0.025);
        } else {
          myMap.SetMaximum(1.0);
        }
        myMap.Draw("colz");

        // Create the image object
        RootWImage& myImage = myContent.addImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
        std::unique_ptr<TCanvas> myCanvas(new TCanvas(Form("StubEfficiencyCoverageCanvas%s", lmel.first.c_str()), "Stub efficiency eta coverage", vis_std_canvas_sizeX, vis_min_canvas_sizeY));
        //myCanvas->SetFillColor(color_plot_background);
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
        RootWImage* myImage = new RootWImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
        std::unique_ptr<TCanvas> myCanvas(new TCanvas());
        double myEfficiency = it->first;
        TH2D& myMap = it->second;
        myCanvas->SetFillColor(color_plot_background);
        myCanvas->cd();

        // Actually plot the map
        myMap.SetMinimum(0);
        if (maxPt>0) myMap.SetMaximum(maxPt);
        else myMap.SetMaximum(10);
        myMap.Draw("colz");

        // Create the image object
        RootWImage& myImage = myContent.addImage(std::move(myCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
    std::unique_ptr<TCanvas> thickCanvas(new TCanvas());
    std::unique_ptr<TCanvas> windowCanvas(new TCanvas());
    std::unique_ptr<TCanvas> suggestedSpacingCanvas(new TCanvas());
    std::unique_ptr<TCanvas> suggestedSpacingAWCanvas(new TCanvas());
    std::unique_ptr<TCanvas> nominalCutCanvas(new TCanvas());
    thickCanvas->SetFillColor(color_plot_background);
    windowCanvas->SetFillColor(color_plot_background);
    suggestedSpacingCanvas->SetFillColor(color_plot_background);
    suggestedSpacingAWCanvas->SetFillColor(color_plot_background);
    nominalCutCanvas->SetFillColor(color_plot_background);

    struct Spacing { double operator()(const Module& m) { return m.dsDistance(); } };
    PlotDrawer<YZ, Spacing> thicknessDrawer;
    thicknessDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    thicknessDrawer.drawFrame<HistogramFrameStyle>(*thickCanvas.get());
    thicknessDrawer.drawModules<ContourStyle>(*thickCanvas.get());

    struct TriggerWindow { double operator()(const Module& m) { return m.triggerWindow(); } };
    PlotDrawer<YZ, TriggerWindow> windowDrawer;
    windowDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
    windowDrawer.drawFrame<HistogramFrameStyle>(*windowCanvas.get());
    windowDrawer.drawModules<ContourStyle>(*windowCanvas.get());


    // Actually plot the maps
    //thickCanvas->cd();
    //thicknessMap.Draw("colz");
    //windowCanvas->cd();
    //windowMap.Draw("colz");
    if (extended) {
      suggestedSpacingCanvas->cd();
      suggestedSpacingMap.Draw("colz");
      suggestedSpacingAWCanvas->cd();
      suggestedSpacingMapAW.Draw("colz");
      nominalCutCanvas->cd();
      //struct PtCut { double operator()(const Module& m) { return PtErrorAdapter(m).getPtCut(); } };
      //PlotDrawer<YZ, PtCut> cutDrawer;
      //cutDrawer.addModulesType(tracker.modules().begin(), tracker.modules().end());
      //cutDrawer.drawFrame<HistogramFrameStyle>(*nominalCutCanvas.get());
      //cutDrawer.drawModules<ContourStyle>(*nominalCutCanvas.get());
      nominalCutCanvas->SetLogz();
      nominalCutMap.Draw("colz");
    }

    // Create the image objects
    RootWImage& thicknessImage = myContent.addImage(std::move(thickCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    thicknessImage.setComment("Map of sensor distances");
    thicknessImage.setName("ThicknessMap");
    RootWImage& windowImage = myContent.addImage(std::move(windowCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
    windowImage.setComment("Map of selection windows");
    windowImage.setName("WindowMap");
    if (extended) {
      RootWImage& suggestedSpacingImage = myContent.addImage(std::move(suggestedSpacingCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      suggestedSpacingImage.setComment("Map of selection suggestedSpacings [default window]");
      suggestedSpacingImage.setName("SuggestedSpacingMap");
      RootWImage& suggestedSpacingAWImage = myContent.addImage(std::move(suggestedSpacingAWCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      suggestedSpacingAWImage.setComment("Map of selection suggestedSpacings [selected windows]");
      suggestedSpacingAWImage.setName("SuggestedSpacingMapAW");
      RootWImage& nominalCutImage = myContent.addImage(std::move(nominalCutCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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

      TH1D* myFrame = new TH1D("myFrame", "", nBins, 0, nBins);
      myFrame->SetYTitle("Optimal distance range [mm]");
      myFrame->SetMinimum(0);
      myFrame->SetMaximum(6);
      TAxis* xAxis = myFrame->GetXaxis();

      TGraphErrors* rangeGraphBad = new TGraphErrors();
      TGraphErrors* rangeGraph = new TGraphErrors();
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
	  rangeGraphPoints=rangeGraph->GetN();
	  rangeGraph->SetPoint(rangeGraphPoints, i+0.5, (min+max)/2.);
	  rangeGraph->SetPointError(rangeGraphPoints, 0.25, (max-min)/2.);
	} else {
	  rangeGraphPoints=rangeGraph->GetN();
	  rangeGraphBad->SetPoint(rangeGraphPoints, i+0.5, (min+max)/2.);
	  rangeGraphBad->SetPointError(rangeGraphPoints, 0.25, (min-max)/2.);
	}
      }

      std::unique_ptr<TCanvas> rangeCanvas(new TCanvas());
      rangeCanvas->cd();
      rangeCanvas->SetFillColor(color_plot_background);
      rangeCanvas->SetGrid(0,1);
      myFrame->SetBit(1); // rangeCanvas takes ownership of it
      myFrame->Draw();
      rangeGraph->SetFillColor(Palette::color(1));
      rangeGraph->SetBit(1);  // rangeCanvas takes ownership of it
      rangeGraph->Draw("same 2");
      rangeGraphBad->SetFillColor(Palette::color(2));
      rangeGraphBad->SetBit(1);  // rangeCanvas takes ownership of it
      rangeGraphBad->Draw("same 2");
      availableSpacings.SetMarkerStyle(0);
      availableSpacings.SetBit(0);
      availableSpacings.Draw("p same");

      RootWImage& RangeImage = spacingSummaryContent.addImage(std::move(rangeCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS.str(""); tempSS << "Sensor distance range tuning";
      RangeImage.setComment(tempSS.str());
      tempSS.str(""); tempSS << "TriggerRangeTuning";
      RangeImage.setName(tempSS.str());

      std::unique_ptr<TCanvas> tuningCanvas(new TCanvas());
      tuningCanvas->cd();
      tuningCanvas->SetFillColor(color_plot_background);
      tuningCanvas->SetGrid(0,1);
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

      RootWImage& tuningImage = spacingSummaryContent.addImage(std::move(tuningCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
      tempSS.str(""); tempSS << "Sensor distance range tuning for different windows";
      tuningImage.setComment(tempSS.str());
      tempSS.str(""); tempSS << "TriggerRangeTuningWindows";
      tuningImage.setName(tempSS.str());

      TH1D& spacingDistribution = a.getHistoOptimalSpacing(false);
      std::unique_ptr<TCanvas> spacingCanvas(new TCanvas());
      spacingCanvas->SetFillColor(color_plot_background);
      spacingCanvas->cd();
      spacingDistribution.SetFillColor(Palette::color(1));
      spacingDistribution.Draw();
      RootWImage& spacingImage = spacingSummaryContent.addImage(std::move(spacingCanvas), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
      spacingImage.setComment("Distribution of minimal spacing for low pT rejection @ standard window");
      spacingImage.setName("SpacingDistribution");
      TH1D& spacingDistributionAW = a.getHistoOptimalSpacing(true);
      std::unique_ptr<TCanvas> spacingCanvasAW(new TCanvas());
      spacingCanvasAW->SetFillColor(color_plot_background);
      spacingCanvasAW->cd();
      spacingDistributionAW.SetFillColor(Palette::color(1));
      spacingDistributionAW.Draw();
      RootWImage& spacingImageAW = spacingSummaryContent.addImage(std::move(spacingCanvasAW), vis_min_canvas_sizeX, vis_min_canvas_sizeY);
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
        std::unique_ptr<TCanvas> tuningCanvas(new TCanvas());
        tuningCanvas->SetFillColor(color_plot_background);
        tuningCanvas->cd();

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

        RootWImage& tuningImage = spacingDetailedContent.addImage(std::move(tuningCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
        std::unique_ptr<TCanvas> turnonCanvas(new TCanvas());
        turnonCanvas->SetFillColor(color_plot_background);
        turnonCanvas->cd();

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
        RootWImage& turnonImage = turnOnDetailedContent.addImage(std::move(turnonCanvas), vis_std_canvas_sizeX, vis_min_canvas_sizeY);
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
                                   std::unique_ptr<TCanvas> &YZCanvas, std::unique_ptr<TCanvas> &XYCanvas,
                                   std::unique_ptr<TCanvas> &XYCanvasEC) {
    Int_t irep;
    TVirtualPad* myPad;

    YZCanvas.reset(new TCanvas("YZCanvas", "YZView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvasEC.reset(new TCanvas("XYCanvasEC", "XYView Canvas (Endcap)", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));

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
                                        std::unique_ptr<TCanvas> &RZCanvas, std::unique_ptr<TCanvas> &RZCanvasBarrel, std::unique_ptr<TCanvas> &XYCanvas,
                                        std::vector<std::unique_ptr<TCanvas> > &XYCanvasesEC) {

    RZCanvas.reset(new TCanvas("RZCanvas", "RZView Canvas",  insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY));
    RZCanvas->cd();
    PlotDrawer<YZ, Type> yzDrawer;
    yzDrawer.addModules(tracker);
    yzDrawer.drawFrame<SummaryFrameStyle>(*RZCanvas.get());
    yzDrawer.drawModules<ContourStyle>(*RZCanvas.get());

    double viewPortMax = MAX(tracker.barrels().at(0).maxR() * 1.1, tracker.barrels().at(0).maxZ() * 1.1); // Style to improve. Calculate (with margin) the barrel geometric extremum
    RZCanvasBarrel.reset(new TCanvas("RZCanvasBarrel", "RZView CanvasBarrel", vis_min_canvas_sizeX, vis_min_canvas_sizeY));
    RZCanvasBarrel->cd();
    PlotDrawer<YZ, Type> yzDrawerBarrel(viewPortMax, viewPortMax);
    yzDrawerBarrel.addModulesType(tracker, BARREL);
    yzDrawerBarrel.drawFrame<SummaryFrameStyle>(*RZCanvasBarrel.get());
    yzDrawerBarrel.drawModules<ContourStyle>(*RZCanvasBarrel.get());

    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, Type> xyBarrelDrawer;
    xyBarrelDrawer.addModulesType(tracker, BARREL);
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());

    for (auto& anEndcap : tracker.endcaps() ) {
      std::unique_ptr<TCanvas> XYCanvasEC(new TCanvas(Form("XYCanvasEC_%s", anEndcap.myid().c_str()),
					Form("XY projection of Endcap %s", anEndcap.myid().c_str()),
						      vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
      XYCanvasEC->cd();
      PlotDrawer<XY, Type> xyEndcapDrawer;
      xyEndcapDrawer.addModules(anEndcap);
      xyEndcapDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasEC.get());
      xyEndcapDrawer.drawModules<ContourStyle>(*XYCanvasEC.get());
      XYCanvasesEC.push_back(std::move(XYCanvasEC));
    }

    // And now one per disk surface
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size()>0) {
	const Disk& lastDisk = anEndcap.disks().back();
     
	const std::map<int, std::vector<const Module*> >& allSurfaceModules = lastDisk.getSurfaceModules();
	for (int surfaceIndex = 1; surfaceIndex <= 4; surfaceIndex++) {
	  auto found = allSurfaceModules.find(surfaceIndex);
	  if (found != allSurfaceModules.end()) {
	    const std::vector<const Module*>& surfaceModules = found->second;
	    std::unique_ptr<TCanvas> XYCanvasEC(new TCanvas(Form("XYCanvasEC_%s_%d", anEndcap.myid().c_str(), surfaceIndex),
					      Form("XY section of Endcap %s -- surface %d", anEndcap.myid().c_str(), surfaceIndex),
					      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	    XYCanvasEC->cd();
	    PlotDrawer<XY, Type> xyEndcapDrawer;

	    xyEndcapDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	    xyEndcapDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasEC.get());
	    xyEndcapDrawer.drawModules<ContourStyle>(*XYCanvasEC.get());
	    xyEndcapDrawer.drawModuleContours<ContourStyle>(*XYCanvasEC.get());
	    XYCanvasesEC.push_back(std::move(XYCanvasEC));
	  }
	  else logERROR("Tried to access modules belonging to one of the 4 disk surfaces, but empty container.");
	}
      }
    }

  }


  /* OT CABLING */

  /* 
   * Bundle cabling plots.
   */
  void Vizard::createOuterCablingPlotsBundles(const Tracker& tracker,
					      std::unique_ptr<TCanvas> &RZCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYNegCanvas,
					      std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDisks, std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDiskSurfaces,
					      std::vector<std::unique_ptr<TCanvas> > &XYNegBundlesDisks, std::vector<std::unique_ptr<TCanvas> > &XYNegBundlesDiskSurfaces) {
    
    const std::set<Module*>& trackerModules = tracker.modules();
    RZCanvas.reset(new TCanvas("RZCanvas", "RZView Canvas", insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY));
    RZCanvas->cd();
    PlotDrawer<YZFull, TypeBundleTransparentColor> yzDrawer;
    // TBPS flat part: only draw 2 consecutives MFBs in phi, otherwise confusing drawing effect happens on bundle coloring.
    yzDrawer.addModules(trackerModules.begin(), trackerModules.end(), [] (const Module& m ) { 
	return ( (!m.getBundle()->isBarrelPSFlatPart()) 
		 || (m.getBundle()->isBarrelPSFlatPart() && (m.getBundle()->phiPosition().phiSegmentRef() == 0 
							     || m.getBundle()->phiPosition().phiSegmentRef() == 1)) 
		 ); 
      } );
    yzDrawer.drawFrame<SummaryFrameStyle>(*RZCanvas.get());
    yzDrawer.drawModules<ContourStyle>(*RZCanvas.get());
   
    // NEGATIVE CABLING SIDE. BARREL.
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNeg, TypeBundleColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveCablingSide() < 0); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines

    // POSITIVE CABLING SIDE. BARREL.
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, TypeBundleColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveCablingSide() > 0); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines

    // POSITIVE CABLING SIDE.
    // ENDCAPS DISK.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();
	std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosBundleEndcap_%sAnyDisk", anEndcap.myid().c_str()),
					    Form("(XY) Projection : Endcap %s, any Disk. (CMS +Z points towards you)", anEndcap.myid().c_str()),
					    vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	XYCanvasDisk->cd();
	PlotDrawer<XY, TypeBundleTransparentColor> xyDiskDrawer;
	xyDiskDrawer.addModules(lastDisk);
	xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
	XYPosBundlesDisks.push_back(std::move(XYCanvasDisk));
      }
    }

    // ENDCAPS DISK SURFACE.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();	
	const std::map<int, std::vector<const Module*> >& allSurfaceModules = lastDisk.getSurfaceModules();
	for (int surfaceIndex = 1; surfaceIndex <= 4; surfaceIndex++) {
	  auto found = allSurfaceModules.find(surfaceIndex);
	  if (found != allSurfaceModules.end()) {  
	    // Surface seen rotated: (+Z) towards the depth of the screen
	    if ((surfaceIndex % 2) == 1) {
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosRotateY180BundleEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : Endcap %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      PlotDrawer<XYRotateY180, TypeFanoutBranchTransparentColor> xyDiskDrawer;
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      const bool isRotatedY180 = true;
	      drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
	      XYPosBundlesDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	    // (+Z) towards you
	    else {
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosBundleEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : Endcap %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      PlotDrawer<XY, TypeFanoutBranchTransparentColor> xyDiskDrawer;
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
	      XYPosBundlesDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	  }
	  else logERROR("Tried to access modules belonging to one of the 4 disk surfaces, but empty container.");
	}
      }
    }

    // NEGATIVE CABLING SIDE.
    // ENDCAPS DISK.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& firstDisk = anEndcap.disks().front();
	std::unique_ptr<TCanvas> XYNegCanvasDisk(new TCanvas(Form("XYNegBundleEndcap_%sAnyDisk", anEndcap.myid().c_str()),
					       Form("(XY) Projection : Endcap %s, any Disk. (CMS +Z points towards the depth of the screen)", 
						    anEndcap.myid().c_str()),
					       vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	XYNegCanvasDisk->cd();
	PlotDrawer<XYNegRotateY180, TypeBundleTransparentColor> xyDiskDrawer;
	xyDiskDrawer.addModules(firstDisk);
	xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvasDisk.get());
	xyDiskDrawer.drawModules<ContourStyle>(*XYNegCanvasDisk.get());
	const bool isRotatedY180 = true;
	drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
	XYNegBundlesDisks.push_back(std::move(XYNegCanvasDisk));
      }
    }

    // ENDCAPS DISK SURFACE.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& firstDisk = anEndcap.disks().front();	
	const std::map<int, std::vector<const Module*> >& allSurfaceModules = firstDisk.getSurfaceModules();
	for (int surfaceIndex = 1; surfaceIndex <= 4; surfaceIndex++) {
	  auto found = allSurfaceModules.find(surfaceIndex);
	  if (found != allSurfaceModules.end()) {
	    // (+Z) towards you
	    if ((surfaceIndex % 2) == 1) {
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYNegSurfaceDisk(new TCanvas(Form("XYNegBundleEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						      Form("(XY) Section : Endcap %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", 
							   anEndcap.myid().c_str(), surfaceIndex),
						      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYNegSurfaceDisk->cd();
	      PlotDrawer<XYNeg, TypeFanoutBranchTransparentColor> xyDiskDrawer;
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYNegSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYNegSurfaceDisk.get());
	      drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
	      XYNegBundlesDiskSurfaces.push_back(std::move(XYNegSurfaceDisk));
	    }
	    // Surface seen rotated: (+Z) towards the depth of the screen
	    else {
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYNegSurfaceDisk(new TCanvas(Form("XYNegBundleEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						      Form("(XY) Section : Endcap %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", 
							   anEndcap.myid().c_str(), surfaceIndex),
						      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYNegSurfaceDisk->cd();
	      PlotDrawer<XYNegRotateY180, TypeFanoutBranchTransparentColor> xyDiskDrawer;
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYNegSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYNegSurfaceDisk.get());
	      const bool isRotatedY180 = true;
	      drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
	      XYNegBundlesDiskSurfaces.push_back(std::move(XYNegSurfaceDisk));
	    }
	  }
	  else logERROR("Tried to access modules belonging to one of the 4 disk surfaces, but empty container.");
	}
      }
    }

  }


  /* 
   * DTC cabling plots.
   */
  void Vizard::createOuterCablingPlotsDTCs(Tracker& tracker,
					   std::unique_ptr<TCanvas> &RZCanvas, 
					   std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, 
					   std::vector<std::unique_ptr<TCanvas> > &XYCanvasesDisk) {

    const std::set<Module*>& trackerModules = tracker.modules();
    RZCanvas.reset(new TCanvas("RZCanvas", "RZView Canvas", insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY));
    RZCanvas->cd();
    PlotDrawer<YZFull, TypeDTCTransparentColor> yzDrawer;
    yzDrawer.addModules(trackerModules.begin(), trackerModules.end(), [] (const Module& m ) { 
	return ( (m.isPositiveCablingSide() > 0 && m.dtcPhiSectorRef() == 1) || (m.isPositiveCablingSide() < 0 && m.dtcPhiSectorRef() == 2) ); 
      } );
    yzDrawer.drawFrame<SummaryFrameStyle>(*RZCanvas.get());
    yzDrawer.drawModules<ContourStyle>(*RZCanvas.get());

    // NEGATIVE CABLING SIDE. BARREL.
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNeg, TypeDTCColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0)); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines

    // NEGATIVE CABLING SIDE. BARREL FLAT PART.
    XYNegFlatCanvas.reset(new TCanvas("XYNegFlatCanvas", "XYNegFlatView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegFlatCanvas->cd();
    PlotDrawer<XYNeg, TypeDTCColor> xyNegFlatBarrelDrawer;
    xyNegFlatBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0) && !m.isTilted()); } );
    xyNegFlatBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegFlatCanvas.get());
    xyNegFlatBarrelDrawer.drawModules<ContourStyle>(*XYNegFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines

    // POSITIVE CABLING SIDE. BARREL.
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, TypeDTCColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0)); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines

    // POSITIVE CABLING SIDE. BARREL FLAT PART.
    XYFlatCanvas.reset(new TCanvas("XYFlatCanvas", "XYView FlatCanvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYFlatCanvas->cd();
    PlotDrawer<XY, TypeDTCColor> xyBarrelFlatDrawer;
    xyBarrelFlatDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0) && !m.isTilted()); } );
    xyBarrelFlatDrawer.drawFrame<SummaryFrameStyle>(*XYFlatCanvas.get());
    xyBarrelFlatDrawer.drawModules<ContourStyle>(*XYFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    
    // ENDCAPS DISK.
    for (auto& anEndcap : tracker.endcaps() ) {
      for (auto& aDisk : anEndcap.disks() ) {
	if (aDisk.side()) {
	  std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosDTCEndcap_%sDisk_%d", anEndcap.myid().c_str(), aDisk.myid()),
					      Form("(XY) Projection : Endcap %s Disk %d. (CMS +Z points towards you)", anEndcap.myid().c_str(), aDisk.myid()),
					      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	  XYCanvasDisk->cd();
	  PlotDrawer<XY, TypeDTCColor> xyDiskDrawer;
	  xyDiskDrawer.addModules(aDisk);
	  xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	  xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	  XYCanvasesDisk.push_back(std::move(XYCanvasDisk));
	  drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
	}
      }
    }
  }


  /*
   * Optical cables channels assignments plots.
   */
  void Vizard::createOuterCablingPlotsServicesChannelsOptical(Tracker& tracker, const OuterCablingMap* myCablingMap,
							   std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, 
							   std::vector<std::unique_ptr<TCanvas> > &XYCanvasesDisk) {
    bool isPowerCabling = false;
    bool isPositiveCablingSide = true;
    TLegend* channelsLegendPos = new TLegend(0.905,0.3,1.0,0.8);
    computeServicesChannelsLegend(channelsLegendPos, myCablingMap, isPositiveCablingSide, isPowerCabling);
    isPositiveCablingSide = false;
    TLegend* channelsLegendNeg = new TLegend(0.905,0.3,1.0,0.8);
    computeServicesChannelsLegend(channelsLegendNeg, myCablingMap, isPositiveCablingSide, isPowerCabling);

    // NEGATIVE CABLING SIDE. BARREL.
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNeg, TypeOpticalChannelColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0)); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendNeg->Draw("same"); 

    // NEGATIVE CABLING SIDE. BARREL FLAT PART.
    XYNegFlatCanvas.reset(new TCanvas("XYNegFlatCanvas", "XYNegFlatView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegFlatCanvas->cd();
    PlotDrawer<XYNeg, TypeOpticalChannelColor> xyNegFlatBarrelDrawer;
    xyNegFlatBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0) && !m.isTilted()); } );
    xyNegFlatBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegFlatCanvas.get());
    xyNegFlatBarrelDrawer.drawModules<ContourStyle>(*XYNegFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendNeg->Draw("same");

    // POSITIVE CABLING SIDE. BARREL.
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, TypeOpticalChannelColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0)); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendPos->Draw("same");

    // POSITIVE CABLING SIDE. BARREL FLAT PART.
    XYFlatCanvas.reset(new TCanvas("XYFlatCanvas", "XYView FlatCanvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYFlatCanvas->cd();
    PlotDrawer<XY, TypeOpticalChannelColor> xyBarrelFlatDrawer;
    xyBarrelFlatDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0) && !m.isTilted()); } );
    xyBarrelFlatDrawer.drawFrame<SummaryFrameStyle>(*XYFlatCanvas.get());
    xyBarrelFlatDrawer.drawModules<ContourStyle>(*XYFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendPos->Draw("same");
    
    // ENDCAPS DISK.
    for (auto& anEndcap : tracker.endcaps() ) {
      for (auto& aDisk : anEndcap.disks() ) {
	if (aDisk.side()) {
	  std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosOpticalChannelsEndcap_%sDisk_%d", anEndcap.myid().c_str(), aDisk.myid()),
					      Form("(XY) Projection : Endcap %s Disk %d. (CMS +Z points towards you)", anEndcap.myid().c_str(), aDisk.myid()),
					      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	  XYCanvasDisk->cd();
	  PlotDrawer<XY, TypeOpticalChannelColor> xyDiskDrawer;
	  xyDiskDrawer.addModules(aDisk);
	  xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	  xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	  XYCanvasesDisk.push_back(std::move(XYCanvasDisk));
	  drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
	  channelsLegendPos->Draw("same");
	}
      }
    }
  }


  /*
   * Power cables channels assignments plots.
   */
  void Vizard::createOuterCablingPlotsServicesChannelsPower(Tracker& tracker, const OuterCablingMap* myCablingMap,
							    std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYNegFlatCanvas, std::unique_ptr<TCanvas> &XYCanvas, std::unique_ptr<TCanvas> &XYFlatCanvas, 
							    std::vector<std::unique_ptr<TCanvas> > &XYCanvasesDisk, std::vector<std::unique_ptr<TCanvas> > &XYNegCanvasesDisk) {
    bool isPowerCabling = true;

    bool isPositiveCablingSide = true;
    TLegend* channelsLegendPos = new TLegend(0.905, 0., 1., 1.);
    computeServicesChannelsLegend(channelsLegendPos, myCablingMap, isPositiveCablingSide, isPowerCabling);
    isPositiveCablingSide = false;
    TLegend* channelsLegendNeg = new TLegend(0.905, 0., 1., 1.);
    computeServicesChannelsLegend(channelsLegendNeg, myCablingMap, isPositiveCablingSide, isPowerCabling);

    // NEGATIVE CABLING SIDE. BARREL.
    bool isRotatedY180 = true;
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNegRotateY180, TypePowerChannelColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0)); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
    channelsLegendNeg->Draw("same");

    // NEGATIVE CABLING SIDE. BARREL FLAT PART.
    isRotatedY180 = true;
    XYNegFlatCanvas.reset(new TCanvas("XYNegFlatCanvas", "XYNegFlatView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegFlatCanvas->cd();
    PlotDrawer<XYNegRotateY180, TypePowerChannelColor> xyNegFlatBarrelDrawer;
    xyNegFlatBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() < 0) && !m.isTilted()); } );
    xyNegFlatBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegFlatCanvas.get());
    xyNegFlatBarrelDrawer.drawModules<ContourStyle>(*XYNegFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
    channelsLegendNeg->Draw("same");

    // POSITIVE CABLING SIDE. BARREL.
    isRotatedY180 = false;
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, TypePowerChannelColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0)); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendPos->Draw("same");

    // POSITIVE CABLING SIDE. BARREL FLAT PART.
    XYFlatCanvas.reset(new TCanvas("XYFlatCanvas", "XYView FlatCanvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYFlatCanvas->cd();
    PlotDrawer<XY, TypePowerChannelColor> xyBarrelFlatDrawer;
    xyBarrelFlatDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return ((m.subdet() == BARREL) && (m.isPositiveCablingSide() > 0) && !m.isTilted()); } );
    xyBarrelFlatDrawer.drawFrame<SummaryFrameStyle>(*XYFlatCanvas.get());
    xyBarrelFlatDrawer.drawModules<ContourStyle>(*XYFlatCanvas.get());
    drawPhiSectorsBoundaries(outer_cabling_nonantWidth);  // Spider lines
    channelsLegendPos->Draw("same");
    
    for (auto& anEndcap : tracker.endcaps() ) {
      for (auto& aDisk : anEndcap.disks() ) {
	// POSITIVE CABLING SIDE. ENDCAPS DISK.
	if (aDisk.side()) {
	  isRotatedY180 = false;
	  std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosPowerChannelsEndcap_%sDisk_%d", anEndcap.myid().c_str(), aDisk.myid()),
					      Form("(XY) Projection : Endcap %s Disk %d. (CMS +Z points towards you)", anEndcap.myid().c_str(), aDisk.myid()),
					      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	  XYCanvasDisk->cd();
	  PlotDrawer<XY, TypePowerChannelColor> xyDiskDrawer;
	  xyDiskDrawer.addModules(aDisk);
	  xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	  xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	  XYCanvasesDisk.push_back(std::move(XYCanvasDisk));
	  drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
	  channelsLegendPos->Draw("same");
	}
	// NEGATIVE CABLING SIDE. ENDCAPS DISK.
	else {
	  isRotatedY180 = true;
	  std::unique_ptr<TCanvas> XYNegCanvasDisk(new TCanvas(Form("XYNegPowerChannelsEndcap_%sDisk_%d", anEndcap.myid().c_str(), aDisk.myid()),
					      Form("(XY) Projection : Endcap %s Disk %d. (CMS +Z points towards the depth of the screen)", anEndcap.myid().c_str(), aDisk.myid()),
					      vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	  XYNegCanvasDisk->cd();
	  PlotDrawer<XYNegRotateY180, TypePowerChannelColor> xyDiskDrawer;
	  xyDiskDrawer.addModules(aDisk);
	  xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvasDisk.get());
	  xyDiskDrawer.drawModules<ContourStyle>(*XYNegCanvasDisk.get());
	  XYNegCanvasesDisk.push_back(std::move(XYNegCanvasDisk));
	  drawPhiSectorsBoundaries(outer_cabling_nonantWidth, isRotatedY180);  // Spider lines
	  channelsLegendNeg->Draw("same");
	}
      }
    }
  }


  /* Interface to gather information on services channels, and create a table storing it.
   */
  RootWTable* Vizard::opticalServicesChannels(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const ChannelSlot requestedSlot) {
    std::map<int, std::vector<int> > cablesPerChannel;
    std::map<int, int> tbpsBundlesPerChannel;
    std::map<int, int> tbssBundlesPerChannel;
    std::map<int, int> teddpsBundlesPerChannel;
    std::map<int, int> teddssBundlesPerChannel;

    // Fill services channels maps.
    analyzeOpticalServicesChannels(myCablingMap, cablesPerChannel, tbpsBundlesPerChannel, tbssBundlesPerChannel, teddpsBundlesPerChannel, teddssBundlesPerChannel, isPositiveCablingSide, requestedSlot); 

    // Create table.
    RootWTable* channelsTable = createOpticalServicesChannelTable(cablesPerChannel, tbpsBundlesPerChannel, tbssBundlesPerChannel, teddpsBundlesPerChannel, teddssBundlesPerChannel, isPositiveCablingSide, requestedSlot); 

    return channelsTable;
  }


  /* Get the requested Services Channels info from the cabling map.
   */
  void Vizard::analyzeOpticalServicesChannels(const OuterCablingMap* myCablingMap, std::map<int, std::vector<int> > &cablesPerChannel, 
					      std::map<int, int> &tbpsBundlesPerChannel, std::map<int, int> &tbssBundlesPerChannel, std::map<int, int> &teddpsBundlesPerChannel, std::map<int, int> &teddssBundlesPerChannel,
					      const bool isPositiveCablingSide, const ChannelSlot requestedSlot) {

    // MFCs
    const std::map<const int, std::unique_ptr<OuterCable> >& cables = (isPositiveCablingSide ? 
								       myCablingMap->getCables() 
								       : myCablingMap->getNegCables());
    // loop on all MFCs
    for (const auto& myCableIt : cables) {
      const OuterCable* myCable = myCableIt.second.get();
      const ChannelSection* mySection = myCable->opticalChannelSection();
      const ChannelSlot& myChannelSlot = mySection->channelSlot();

      // If a specific slot is requested, select that slot only!
      if ( requestedSlot == ChannelSlot::UNKNOWN 
	   || (requestedSlot != ChannelSlot::UNKNOWN && myChannelSlot == requestedSlot)
	   ) {
	// Channel number
	const int channelNumber = mySection->channelNumber();

	// MFC
	const int cableId = myCableIt.first;
	cablesPerChannel[channelNumber].push_back(cableId);
      }
    }

    // MFBs
    const std::map<const int, std::unique_ptr<OuterBundle> >& bundles = (isPositiveCablingSide ? myCablingMap->getBundles() : myCablingMap->getNegBundles());
    // loop on all MFBs
    for (const auto& myBundle : bundles) {
      const ChannelSection* mySection = myBundle.second->opticalChannelSection();
      const ChannelSlot& myChannelSlot = mySection->channelSlot();

      // If a specific slot is requested, select that slot only!
      if ( requestedSlot == ChannelSlot::UNKNOWN 
	   || (requestedSlot != ChannelSlot::UNKNOWN && myChannelSlot == requestedSlot)
	   ) {

	// Channel number
	const int channelNumber = mySection->channelNumber();
	const std::string subDetectorName = myBundle.second->subDetectorName();
	const Category bundleType = myBundle.second->type();   
	// TBPS   
	if (subDetectorName == outer_cabling_tbps) tbpsBundlesPerChannel[channelNumber] += 1;
	// TB2S
	else if (subDetectorName == outer_cabling_tb2s) tbssBundlesPerChannel[channelNumber] += 1;
	// TEDD
	else {
	  if (bundleType == Category::PS10G 
	      || bundleType == Category::PS10GA 
	      || bundleType == Category::PS10GB 
	      || bundleType == Category::PS5G) { teddpsBundlesPerChannel[channelNumber] += 1; }
	  else if (bundleType == Category::SS) { teddssBundlesPerChannel[channelNumber] += 1; }
	  else { logERROR("analyzeServicesChannels : subdetector " + any2str(subDetectorName) 
			  + " has undetected bundle type" + any2str(bundleType)); }
	}
      }
    }

  }


  /* Create the table with Services Channel information.
   */
  RootWTable* Vizard::createOpticalServicesChannelTable(const std::map<int, std::vector<int> > &cablesPerChannel, 
							std::map<int, int> &tbpsBundlesPerChannel, std::map<int, int> &tbssBundlesPerChannel, std::map<int, int> &teddpsBundlesPerChannel, std::map<int, int> &teddssBundlesPerChannel,
							const bool isPositiveCablingSide, const ChannelSlot requestedSlot) {

    RootWTable* channelsTable = new RootWTable();

    // Header table
    channelsTable->setContent(0, 1, any2str(requestedSlot));
    channelsTable->setContent(0, 2, "# MFC");
    channelsTable->setContent(0, 3, "# MFB TBPS");
    channelsTable->setContent(0, 4, "# MFB TB2S");
    channelsTable->setContent(0, 5, "# MFB TEDD PS");
    channelsTable->setContent(0, 6, "# MFB TEDD 2S");
    channelsTable->setContent(0, 7, "# MFB Total");

    int totalCables = 0;
    int totalTbpsBundles = 0;
    int totalTbssBundles = 0;
    int totalTeddpsBundles = 0;
    int totalTeddssBundles = 0;
    int totalBundles = 0;

    // Fill table
    for (int i = 1; i <= 12; i++) {
      const int channelNumber = (isPositiveCablingSide ? i : -i);
      int numCablesPerChannel = (cablesPerChannel.count(channelNumber) != 0 ? cablesPerChannel.at(channelNumber).size() : 0);
      int numTbpsBundlesPerChannel = (tbpsBundlesPerChannel.count(channelNumber) != 0 ? tbpsBundlesPerChannel.at(channelNumber) : 0);
      int numTbssBundlesPerChannel = (tbssBundlesPerChannel.count(channelNumber) != 0 ? tbssBundlesPerChannel.at(channelNumber) : 0);
      int numTeddpsBundlesPerChannel = (teddpsBundlesPerChannel.count(channelNumber) != 0 ? teddpsBundlesPerChannel.at(channelNumber) : 0);
      int numTeddssBundlesPerChannel = (teddssBundlesPerChannel.count(channelNumber) != 0 ? teddssBundlesPerChannel.at(channelNumber) : 0);
      int numBundlesPerChannel = numTbpsBundlesPerChannel + numTbssBundlesPerChannel + numTeddpsBundlesPerChannel + numTeddssBundlesPerChannel;     

      // PP1 name
      const int pp1 = channelNumber 
	+ (channelNumber >= 0 ? (fabs(channelNumber) <= 6 ? 2 : 5) 
	   : -(fabs(channelNumber) <= 6 ? 2 : 5) );
      std::stringstream pp1Name;
      std::string sign = (pp1 >= 0 ? "+" : "");
      pp1Name << "PP1" << sign << pp1;
      if (requestedSlot != ChannelSlot::UNKNOWN) pp1Name << " " << any2str(requestedSlot);
      channelsTable->setContent(i, 0, pp1Name.str());

      // Channel name
      std::stringstream channelName;
      channelName << "OT" << channelNumber;
      if (requestedSlot != ChannelSlot::UNKNOWN) channelName << " " << any2str(requestedSlot);
      channelsTable->setContent(i, 1, channelName.str());

      channelsTable->setContent(i, 2, numCablesPerChannel);
      channelsTable->setContent(i, 3, numTbpsBundlesPerChannel);
      channelsTable->setContent(i, 4, numTbssBundlesPerChannel);
      channelsTable->setContent(i, 5, numTeddpsBundlesPerChannel);
      channelsTable->setContent(i, 6, numTeddssBundlesPerChannel);
      channelsTable->setContent(i, 7, numBundlesPerChannel);

      // Gather totals
      totalCables += numCablesPerChannel;
      totalTbpsBundles += numTbpsBundlesPerChannel;
      totalTbssBundles += numTbssBundlesPerChannel;
      totalTeddpsBundles += numTeddpsBundlesPerChannel;
      totalTeddssBundles += numTeddssBundlesPerChannel;
      totalBundles += numBundlesPerChannel;
    }
    channelsTable->setContent(13, 1, "Total");
    channelsTable->setContent(13, 2, totalCables);
    channelsTable->setContent(13, 3, totalTbpsBundles);
    channelsTable->setContent(13, 4, totalTbssBundles);
    channelsTable->setContent(13, 5, totalTeddpsBundles);
    channelsTable->setContent(13, 6, totalTeddssBundles);
    channelsTable->setContent(13, 7, totalBundles);

    return channelsTable;
  }


  /* Interface to gather information on powerServices channels, and create a table storing it.
   */
  RootWTable* Vizard::powerServicesChannels(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const std::vector<ChannelSlot>& slots) {

    RootWTable* channelsTable = new RootWTable();

    for (const auto& requestedSlot : slots) {
      std::map<int, int> tbpsBundlesPerChannel;
      std::map<int, int> tbssBundlesPerChannel;
      std::map<int, int> teddpsBundlesPerChannel;
      std::map<int, int> teddssBundlesPerChannel;

      // Fill powerServices channels maps.
      analyzePowerServicesChannels(myCablingMap, tbpsBundlesPerChannel, tbssBundlesPerChannel, teddpsBundlesPerChannel, teddssBundlesPerChannel, isPositiveCablingSide, requestedSlot); 

      // Create table.
      createPowerServicesChannelTable(channelsTable, tbpsBundlesPerChannel, tbssBundlesPerChannel, teddpsBundlesPerChannel, teddssBundlesPerChannel, isPositiveCablingSide, requestedSlot); 
    }

    return channelsTable;
  }


  /* Get the requested PowerServices Channels info from the cabling map.
   */
  void Vizard::analyzePowerServicesChannels(const OuterCablingMap* myCablingMap, 
					    std::map<int, int> &tbpsBundlesPerChannel, std::map<int, int> &tbssBundlesPerChannel, std::map<int, int> &teddpsBundlesPerChannel, std::map<int, int> &teddssBundlesPerChannel,
					    const bool isPositiveCablingSide, const ChannelSlot requestedSlot) {

    // Power cables = MFBs scheme
    const std::map<const int, std::unique_ptr<OuterBundle> >& bundles = (isPositiveCablingSide ? 
									 myCablingMap->getBundles() 
									 : myCablingMap->getNegBundles());

    for (const auto& myBundleIt : bundles) {
      const OuterBundle* myBundle = myBundleIt.second.get();
      const ChannelSection* mySection = myBundle->powerChannelSection();
      const ChannelSlot& myChannelSlot = mySection->channelSlot();

      // If necessary, can select the PowerServices Channels corresponding to the requested slot.
      if ( requestedSlot == ChannelSlot::UNKNOWN 
	   || (requestedSlot != ChannelSlot::UNKNOWN && myChannelSlot == requestedSlot)
	   ) {

	// Channel number
	const int channelNumber = mySection->channelNumber();

	const std::string subDetectorName = myBundle->subDetectorName();
	const Category bundleType = myBundle->type();    

	// TBPS
	if (subDetectorName == outer_cabling_tbps) { tbpsBundlesPerChannel[channelNumber] += 1; }
	// TB2S
	else if (subDetectorName == outer_cabling_tb2s) { tbssBundlesPerChannel[channelNumber] += 1; }
	// TEDD
	else {
	  if (bundleType == Category::PS10G 
	      || bundleType == Category::PS10GA 
	      || bundleType == Category::PS10GB 
	      || bundleType == Category::PS5G) { teddpsBundlesPerChannel[channelNumber] += 1; }
	  else if (bundleType == Category::SS) { teddssBundlesPerChannel[channelNumber] += 1; }
	  else { logERROR("analyzePowerServicesChannels : subdetector " + any2str(subDetectorName) 
			  + " has undetected bundle type" + any2str(bundleType)); }
	}
      }
    }
  }


  /* Create the table with PowerServices Channel information.
   */
  void Vizard::createPowerServicesChannelTable(RootWTable* channelsTable, 
					       std::map<int, int> &tbpsBundlesPerChannel, std::map<int, int> &tbssBundlesPerChannel, std::map<int, int> &teddpsBundlesPerChannel, std::map<int, int> &teddssBundlesPerChannel,
					       const bool isPositiveCablingSide, const ChannelSlot requestedSlot) {

    const int maxCol = channelsTable->maxCol();
    const int startCol = (maxCol == 0 ? 0 : maxCol + 1);

    // Header table
    channelsTable->setContent(0, startCol + 1, any2str(requestedSlot));
    channelsTable->setContent(0, startCol + 2, "# PWR TBPS");
    channelsTable->setContent(0, startCol + 3, "# PWR TB2S");
    channelsTable->setContent(0, startCol + 4, "# PWR TEDD PS");
    channelsTable->setContent(0, startCol + 5, "# PWR TEDD 2S");
    channelsTable->setContent(0, startCol + 6, "# PWR Total");
    
    int totalTbpsBundles = 0;
    int totalTbssBundles = 0;
    int totalTeddpsBundles = 0;
    int totalTeddssBundles = 0;
    int totalBundles = 0;

    // Fill table
    for (int i = 1; i <= 12; i++) {
      const int channelNumber = (isPositiveCablingSide ? i : -i);
      int numTbpsBundlesPerChannel = (tbpsBundlesPerChannel.count(channelNumber) != 0 ? tbpsBundlesPerChannel.at(channelNumber) : 0);
      int numTbssBundlesPerChannel = (tbssBundlesPerChannel.count(channelNumber) != 0 ? tbssBundlesPerChannel.at(channelNumber) : 0);
      int numTeddpsBundlesPerChannel = (teddpsBundlesPerChannel.count(channelNumber) != 0 ? teddpsBundlesPerChannel.at(channelNumber) : 0);
      int numTeddssBundlesPerChannel = (teddssBundlesPerChannel.count(channelNumber) != 0 ? teddssBundlesPerChannel.at(channelNumber) : 0);
      int numBundlesPerChannel = numTbpsBundlesPerChannel + numTbssBundlesPerChannel + numTeddpsBundlesPerChannel + numTeddssBundlesPerChannel;

      // PP1 name
      const int pp1 = channelNumber 
	+ (channelNumber >= 0 ? (fabs(channelNumber) <= 6 ? 2 : 5) 
	   : -(fabs(channelNumber) <= 6 ? 2 : 5) );
      std::stringstream pp1Name;
      std::string sign = (pp1 >= 0 ? "+" : "");
      pp1Name << "PP1" << sign << pp1;
      if (requestedSlot != ChannelSlot::UNKNOWN) pp1Name << " " << any2str(requestedSlot);
      channelsTable->setContent(i, startCol, pp1Name.str());

      // Channel name
      std::stringstream channelName;
      channelName << "OT" << channelNumber;
      if (requestedSlot != ChannelSlot::UNKNOWN) channelName << " " << any2str(requestedSlot);
      channelsTable->setContent(i, startCol + 1, channelName.str());

      channelsTable->setContent(i, startCol + 2, numTbpsBundlesPerChannel);
      channelsTable->setContent(i, startCol + 3, numTbssBundlesPerChannel);
      channelsTable->setContent(i, startCol + 4, numTeddpsBundlesPerChannel);
      channelsTable->setContent(i, startCol + 5, numTeddssBundlesPerChannel);
      channelsTable->setContent(i, startCol + 6, numBundlesPerChannel);

      totalTbpsBundles += numTbpsBundlesPerChannel;
      totalTbssBundles += numTbssBundlesPerChannel;
      totalTeddpsBundles += numTeddpsBundlesPerChannel;
      totalTeddssBundles += numTeddssBundlesPerChannel;
      totalBundles += numBundlesPerChannel;
    }
    channelsTable->setContent(13, startCol + 1, "Total");
    channelsTable->setContent(13, startCol + 2, totalTbpsBundles);
    channelsTable->setContent(13, startCol + 3, totalTbssBundles);
    channelsTable->setContent(13, startCol + 4, totalTeddpsBundles);
    channelsTable->setContent(13, startCol + 5, totalTeddssBundles);
    channelsTable->setContent(13, startCol + 6, totalBundles);
  }


  /* IT CABLING */

  /*
   * Power chains plots
   */
  void Vizard::createInnerCablingPlotsPowerChains(const Tracker& tracker,
						  std::vector<std::unique_ptr<TCanvas> > &ZPhiLayerPlots,
						  std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYCentralCanvas, std::unique_ptr<TCanvas> &XYCanvas,
						  std::vector<std::unique_ptr<TCanvas> > &XYPosPowerChainsDiskSurfaces) {

    const std::pair<double, double> maxRadii = computeInnerCablingPlotsMaxRadii(tracker);
    const double forwardViewPort = maxRadii.second;

    const std::pair<double, double> scalingFactors = computeInnerCablingPlotsScalingFactors(tracker);
    const double barrelScalingFactor = scalingFactors.first;
    const double forwardScalingFactor = scalingFactors.second;

    // BARREL LAYERS, (ZPhi).
    const int numLayers = tracker.barrels().at(0).layers().size(); // TO DO : ugly
    for (int layerNumber = 1; layerNumber <= numLayers; layerNumber++) {
      // POSITIVE X SIDE
      std::unique_ptr<TCanvas> ZPhiCanvasPos(new TCanvas(Form("ZPhiPowerChainBarrelLayer%d_positiveXSide", layerNumber),
					Form("(ZPhi), Barrel Layer %d. (+X) side.", layerNumber), vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
      ZPhiCanvasPos->cd();
      PlotDrawer<ZPhi, TypePowerChainTransparentColor> zphiBarrelDrawerPos;
      zphiBarrelDrawerPos.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && m.isPositiveXSide()
		  ); 
	} );
      zphiBarrelDrawerPos.drawFrame<SummaryFrameStyle>(*ZPhiCanvasPos.get());
      zphiBarrelDrawerPos.drawModules<ContourStyle>(*ZPhiCanvasPos.get());
      ZPhiLayerPlots.push_back(std::move(ZPhiCanvasPos));
      // NEGATIVE X SIDE
      std::unique_ptr<TCanvas> ZPhiCanvasNeg(new TCanvas(Form("ZPhiPowerChainBarrelLayer%d_negativeXSide", layerNumber),
					Form("(ZPhi), Barrel Layer %d. (-X) side.", layerNumber), vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
      ZPhiCanvasNeg->cd();
      PlotDrawer<ZPhi, TypePowerChainTransparentColor> zphiBarrelDrawerNeg;
      zphiBarrelDrawerNeg.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && !m.isPositiveXSide()
		  ); 
	} );
      zphiBarrelDrawerNeg.drawFrame<SummaryFrameStyle>(*ZPhiCanvasNeg.get());
      zphiBarrelDrawerNeg.drawModules<ContourStyle>(*ZPhiCanvasNeg.get());
      ZPhiLayerPlots.push_back(std::move(ZPhiCanvasNeg));
    }
       
    // NEGATIVE CABLING SIDE. BARREL (XY).
    bool isRotatedY180 = false;
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNeg, TypePowerChainTransparentColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveZEnd() < 0); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // POSITIVE CABLING SIDE. BARREL CENTRAL MODULES (XY).
    isRotatedY180 = false;
    XYCentralCanvas.reset(new TCanvas("XYCentralCanvas", "XYCentralView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCentralCanvas->cd();
    PlotDrawer<XY, TypePowerChainTransparentColor> xyCentralBarrelDrawer;
    xyCentralBarrelDrawer.addModules( tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.uniRef().ring == 1); } );
    xyCentralBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCentralCanvas.get());
    xyCentralBarrelDrawer.drawModules<ContourStyle>(*XYCentralCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // POSITIVE CABLING SIDE. BARREL (XY).
    isRotatedY180 = false;
    XYCanvas.reset(new TCanvas("XYCanvas", "XYView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYCanvas->cd();
    PlotDrawer<XY, TypePowerChainTransparentColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveZEnd() > 0); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // ENDCAPS DISK SURFACE.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();	
	const std::map<int, std::vector<const Module*> >& allSurfaceModules = lastDisk.getSurfaceModules();
	for (int surfaceIndex = 1; surfaceIndex <= 4; surfaceIndex++) {
	  auto found = allSurfaceModules.find(surfaceIndex);
	  if (found != allSurfaceModules.end()) {  
	    // Surface seen rotated: (+Z) towards the depth of the screen
	    if ((surfaceIndex % 2) == 1) {
	      isRotatedY180 = true;;
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosRotateY180PowerChainEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      PlotDrawer<XYRotateY180, TypePowerChainTransparentColor> xyDiskDrawer(forwardViewPort, forwardViewPort);
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	      XYPosPowerChainsDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	    // (+Z) towards you
	    else {
	      isRotatedY180 = false;
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosPowerChainEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      PlotDrawer<XY, TypePowerChainTransparentColor> xyDiskDrawer(forwardViewPort, forwardViewPort);
	      xyDiskDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { return (m.subdet() == ENDCAP); } );
	      xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	      XYPosPowerChainsDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	  }
	  else logERROR("Tried to access modules belonging to one of the 4 disk surfaces, but empty container.");
	}
      }
    }

  }


  /*
   * GBT plots.
   */
  void Vizard::createInnerCablingPlotsGBTs(const Tracker& tracker,
					   std::vector<std::unique_ptr<TCanvas> > &ZPhiLayerPlots,
					   std::vector<std::unique_ptr<TCanvas> > &XYPosGBTsDiskSurfaces) {

    const std::pair<double, double> maxRadii = computeInnerCablingPlotsMaxRadii(tracker);
    const double forwardViewPort = maxRadii.second;

    const std::pair<double, double> scalingFactors = computeInnerCablingPlotsScalingFactors(tracker);
    const double forwardScalingFactor = scalingFactors.second;

    // BARREL LAYERS, (ZPhi).
    const int numLayers = tracker.barrels().at(0).layers().size(); // TO DO : ugly
    for (int layerNumber = 1; layerNumber <= numLayers; layerNumber++) {
      // POSITIVE X SIDE
      std::unique_ptr<TCanvas> ZPhiCanvasPos(new TCanvas(Form("ZPhiGBTBarrelLayer%d_positiveXSide", layerNumber),
					   Form("(ZPhi), Barrel Layer %d. (+X) side. (≠ colors) => (≠ power chains). Alternance of groups of filled/contoured module(s) is used to show the alternance of GBTs.", layerNumber), vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
      ZPhiCanvasPos->cd();
      // Contour modules
      PlotDrawer<ZPhi, TypeGBTTransparentColor> zphiBarrelContourDrawerPos;
      zphiBarrelContourDrawerPos.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && m.isPositiveXSide()
		  && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 0)
		  ); 
	} );
      zphiBarrelContourDrawerPos.drawFrame<SummaryFrameStyle>(*ZPhiCanvasPos.get());
      zphiBarrelContourDrawerPos.drawModules<ContourStyle>(*ZPhiCanvasPos.get());
      // Filled modules
      PlotDrawer<ZPhi, TypeGBTTransparentColor> zphiBarrelFillDrawerPos;
      zphiBarrelFillDrawerPos.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && m.isPositiveXSide()
		  && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 1)
		  ); 
	} );
      zphiBarrelFillDrawerPos.drawModules<FillStyle>(*ZPhiCanvasPos.get());
      ZPhiLayerPlots.push_back(std::move(ZPhiCanvasPos));
      // NEGATIVE X SIDE
      std::unique_ptr<TCanvas> ZPhiCanvasNeg(new TCanvas(Form("ZPhiGBTBarrelLayer%d_negativeXSide", layerNumber),
					   Form("(ZPhi), Barrel Layer %d. (+X) side. (≠ colors) => (≠ power chains). Alternance of groups of filled/contoured module(s) is used to show the alternance of GBTs.", layerNumber), vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
      ZPhiCanvasNeg->cd();
      // Contour modules
      PlotDrawer<ZPhi, TypeGBTTransparentColor> zphiBarrelContourDrawerNeg;
      zphiBarrelContourDrawerNeg.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && !m.isPositiveXSide()
		  && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 0)
		  ); 
	} );
      zphiBarrelContourDrawerNeg.drawFrame<SummaryFrameStyle>(*ZPhiCanvasNeg.get());
      zphiBarrelContourDrawerNeg.drawModules<ContourStyle>(*ZPhiCanvasNeg.get());
      // Filled modules
      PlotDrawer<ZPhi, TypeGBTTransparentColor> zphiBarrelFillDrawerNeg;
      zphiBarrelFillDrawerNeg.addModules(tracker.modules().begin(), tracker.modules().end(), [layerNumber] (const Module& m ) { 
	  return (m.subdet() == BARREL 
		  && m.uniRef().layer == layerNumber
		  && !m.isPositiveXSide()
		  && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 1)
		  ); 
	} );
      zphiBarrelFillDrawerNeg.drawModules<FillStyle>(*ZPhiCanvasNeg.get());
      ZPhiLayerPlots.push_back(std::move(ZPhiCanvasNeg));
    }


    // POSITIVE CABLING SIDE.
    // ENDCAPS DISK SURFACE.
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();	
	const std::map<int, std::vector<const Module*> >& allSurfaceModules = lastDisk.getSurfaceModules();
	for (int surfaceIndex = 1; surfaceIndex <= 4; surfaceIndex++) {
	  auto found = allSurfaceModules.find(surfaceIndex);
	  if (found != allSurfaceModules.end()) {  
	    // Surface seen rotated: (+Z) towards the depth of the screen
	    if ((surfaceIndex % 2) == 1) {
	      bool isRotatedY180 = true;;
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosRotateY180GBTEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      // Filled modules
	      PlotDrawer<XYRotateY180, TypeGBTTransparentColor> xyDiskFillDrawer(forwardViewPort, forwardViewPort);
	      xyDiskFillDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { 
		  return ( (m.subdet() == ENDCAP)
			   && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 0)
			   );
		} );
	      xyDiskFillDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskFillDrawer.drawModules<FillStyle>(*XYSurfaceDisk.get());
	      // Contour modules
	      PlotDrawer<XYRotateY180, TypeGBTTransparentColor> xyDiskContourDrawer(forwardViewPort, forwardViewPort);
	      xyDiskContourDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { 
		  return ( (m.subdet() == ENDCAP)
			   && ((m.getGBT() ? m.getGBT()->indexColor() : 0) == 1)
			   );
		} );
	      xyDiskContourDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	      XYPosGBTsDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	    // (+Z) towards you
	    else {
	      bool isRotatedY180 = false;
	      const std::vector<const Module*>& surfaceModules = found->second;
	      std::unique_ptr<TCanvas> XYSurfaceDisk(new TCanvas(Form("XYPosGBTEndcap_%sAnyDiskSurface_%d", anEndcap.myid().c_str(), surfaceIndex),
						   Form("(XY) Section : %s, any Disk, Surface %d. (The 4 surfaces of a disk are indexed such that |zSurface1| < |zSurface2| < |zSurface3| < |zSurface4|)", anEndcap.myid().c_str(), surfaceIndex),
						   vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	      XYSurfaceDisk->cd();
	      // Filled modules
	      PlotDrawer<XY, TypeGBTTransparentColor> xyDiskFillDrawer(forwardViewPort, forwardViewPort);
	      xyDiskFillDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { 
		  return ( (m.subdet() == ENDCAP)
			   && (femod((m.getGBT() ? m.getGBT()->GBTPhiIndex() : 0), 2) == 0)
			   );
		} );
	      xyDiskFillDrawer.drawFrame<SummaryFrameStyle>(*XYSurfaceDisk.get());
	      xyDiskFillDrawer.drawModules<FillStyle>(*XYSurfaceDisk.get());
	      // Contour modules
	      PlotDrawer<XY, TypeGBTTransparentColor> xyDiskContourDrawer(forwardViewPort, forwardViewPort);
	      xyDiskContourDrawer.addModules(surfaceModules.begin(), surfaceModules.end(), [] (const Module& m ) { 
		  return ( (m.subdet() == ENDCAP)
			   && (femod((m.getGBT() ? m.getGBT()->GBTPhiIndex() : 0), 2) == 1)
			   );
		} );
	      xyDiskContourDrawer.drawModules<ContourStyle>(*XYSurfaceDisk.get());
	      drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	      XYPosGBTsDiskSurfaces.push_back(std::move(XYSurfaceDisk));
	    }
	  }
	  else logERROR("Tried to access modules belonging to one of the 4 disk surfaces, but empty container.");
	}
      }
    }

  }


  /*
   * Bundles plots.
   */
  void Vizard::createInnerCablingPlotsBundles(const Tracker& tracker,
					      std::unique_ptr<TCanvas> &XYNegCanvas, std::unique_ptr<TCanvas> &XYPosCanvas,
					      std::vector<std::unique_ptr<TCanvas> > &XYPosBundlesDisks) {

    const std::pair<double, double> maxRadii = computeInnerCablingPlotsMaxRadii(tracker);
    const double forwardViewPort = maxRadii.second;

    const std::pair<double, double> scalingFactors = computeInnerCablingPlotsScalingFactors(tracker);
    const double barrelScalingFactor = scalingFactors.first;
    const double forwardScalingFactor = scalingFactors.second;
       
    // NEGATIVE CABLING SIDE. BARREL.
    bool isRotatedY180 = false;
    XYNegCanvas.reset(new TCanvas("XYNegCanvas", "XYNegView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYNegCanvas->cd();
    PlotDrawer<XYNeg, TypeInnerBundleTransparentColor> xyNegBarrelDrawer;
    xyNegBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveZEnd() < 0); } );
    xyNegBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYNegCanvas.get());
    xyNegBarrelDrawer.drawModules<ContourStyle>(*XYNegCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // POSITIVE CABLING SIDE. BARREL.
    isRotatedY180 = false;
    XYPosCanvas.reset(new TCanvas("XYPosCanvas", "XYPosView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYPosCanvas->cd();
    PlotDrawer<XY, TypeInnerBundleTransparentColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(tracker.modules().begin(), tracker.modules().end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveZEnd() > 0); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYPosCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYPosCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // POSITIVE CABLING SIDE.
    // ENDCAPS DISK.
    isRotatedY180 = false;
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();
	std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosBundleEndcap_%sAnyDisk", anEndcap.myid().c_str()),
					    Form("(XY) Projection : %s, any Disk. (CMS +Z points towards you). 1 color <=> 1 Fibre Bundle.", 
						 anEndcap.myid().c_str()),
					    vis_min_canvas_sizeX, vis_min_canvas_sizeY) );
	XYCanvasDisk->cd();
	PlotDrawer<XY, TypeInnerBundleTransparentColor> xyDiskDrawer(forwardViewPort, forwardViewPort);
	xyDiskDrawer.addModules(lastDisk);
	xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	XYPosBundlesDisks.push_back(std::move(XYCanvasDisk));
      }
    }

  }


  /*
   * DTC plots.
   */
  void Vizard::createInnerCablingPlotsDTCs(const Tracker& tracker,
					   std::unique_ptr<TCanvas> &RZCanvas,
					   std::unique_ptr<TCanvas> &XYPosCanvas,
					   std::vector<std::unique_ptr<TCanvas> > &XYPosDTCsDisks) {

    const std::set<Module*>& trackerModules = tracker.modules();
    RZCanvas.reset(new TCanvas("RZCanvas", "RZView Canvas", insur::vis_max_canvas_sizeX, insur::vis_min_canvas_sizeY));
    RZCanvas->cd();
    PlotDrawer<YZFull, TypeInnerDTCTransparentColor> yzDrawer;
    yzDrawer.addModules(tracker);
    yzDrawer.drawFrame<SummaryFrameStyle>(*RZCanvas.get());
    yzDrawer.drawModules<ContourStyle>(*RZCanvas.get());

    const std::pair<double, double> maxRadii = computeInnerCablingPlotsMaxRadii(tracker);
    const double forwardViewPort = maxRadii.second;

    const std::pair<double, double> scalingFactors = computeInnerCablingPlotsScalingFactors(tracker);
    const double barrelScalingFactor = scalingFactors.first;
    const double forwardScalingFactor = scalingFactors.second;

    // POSITIVE CABLING SIDE. BARREL.
    bool isRotatedY180 = false;
    XYPosCanvas.reset(new TCanvas("XYPosCanvas", "XYPosView Canvas", vis_min_canvas_sizeX, vis_min_canvas_sizeY ));
    XYPosCanvas->cd();
    PlotDrawer<XY, TypeInnerDTCTransparentColor> xyBarrelDrawer;
    xyBarrelDrawer.addModules(trackerModules.begin(), trackerModules.end(), [] (const Module& m ) { return (m.subdet() == BARREL && m.isPositiveZEnd() > 0); } );
    xyBarrelDrawer.drawFrame<SummaryFrameStyle>(*XYPosCanvas.get());
    xyBarrelDrawer.drawModules<ContourStyle>(*XYPosCanvas.get());
    drawFrameOfReference(isRotatedY180, barrelScalingFactor);

    // POSITIVE CABLING SIDE.
    // ENDCAPS DISK.
    isRotatedY180 = false;
    for (auto& anEndcap : tracker.endcaps() ) {
      if (anEndcap.disks().size() > 0) {
	const Disk& lastDisk = anEndcap.disks().back();
	std::unique_ptr<TCanvas> XYCanvasDisk(new TCanvas(Form("XYPosDTCEndcap_%sAnyDisk", anEndcap.myid().c_str()),
					    Form("(XY) Projection : %s, one Disk. (CMS +Z points towards you). 1 color <=> 1 DTC.", 
						 anEndcap.myid().c_str()),
					    vis_min_canvas_sizeX, vis_min_canvas_sizeY));
	XYCanvasDisk->cd();
	PlotDrawer<XY, TypeInnerDTCTransparentColor> xyDiskDrawer(forwardViewPort, forwardViewPort);
	xyDiskDrawer.addModules(lastDisk);
	xyDiskDrawer.drawFrame<SummaryFrameStyle>(*XYCanvasDisk.get());
	xyDiskDrawer.drawModules<ContourStyle>(*XYCanvasDisk.get());
	drawFrameOfReference(isRotatedY180, forwardScalingFactor);
	XYPosDTCsDisks.push_back(std::move(XYCanvasDisk));
      }
    }
  }


  /*
   * Count all Power chains, GBTs, Bundles, and DTCs in the Inner Tracker cabling map.
   */
  void Vizard::computeInnerCablingCount(const InnerCablingMap* myInnerCablingMap,
					int& numSensorsOneXSide, int& numSensorsPlusXSidePlusZEnd, int& numSensorsPlusXSideMinusZEnd,
					int& numPowerChainsOneXSide, int& numPowerChainsPlusXSidePlusZEnd, int& numPowerChainsPlusXSideMinusZEnd,
					int& numELinksOneXSide, int& numELinksPlusXSidePlusZEnd, int& numELinksPlusXSideMinusZEnd,
					int& numBundlesOneXSide, int& numBundlesPlusXSidePlusZEnd, int& numBundlesPlusXSideMinusZEnd,
					int& numGBTsOneXSide, int& numGBTsPlusXSidePlusZEnd, int& numGBTsPlusXSideMinusZEnd,
					int& numDTCsOneXSide, int& numDTCsPlusXSidePlusZEnd, int& numDTCsPlusXSideMinusZEnd) const {

    // PowerChains
    const std::map<int, std::unique_ptr<PowerChain> >& powerChains = myInnerCablingMap->getPowerChains();
    for (const auto& it : powerChains) {
      const PowerChain* myPowerChain = it.second.get();
      if (myPowerChain->isPositiveXSide()) {
	numPowerChainsOneXSide++;
	if (myPowerChain->isPositiveZEnd()) numPowerChainsPlusXSidePlusZEnd++; else numPowerChainsPlusXSideMinusZEnd++;

	const int numSensors = myPowerChain->numModules();
	numSensorsOneXSide += numSensors;
	if (myPowerChain->isPositiveZEnd()) numSensorsPlusXSidePlusZEnd += numSensors; else numSensorsPlusXSideMinusZEnd += numSensors;
      }
    }

    // GBTs
    const std::map<std::string, std::unique_ptr<GBT> >& gbts = myInnerCablingMap->getGBTs();
    for (const auto& it : gbts) {
      const GBT* myGBT = it.second.get();
      if (myGBT->isPositiveXSide()) {
	numGBTsOneXSide++;
	if (myGBT->isPositiveZEnd()) numGBTsPlusXSidePlusZEnd++; else numGBTsPlusXSideMinusZEnd++;

	const int numELinks = myGBT->numELinks();
	numELinksOneXSide += numELinks;
	if (myGBT->isPositiveZEnd()) numELinksPlusXSidePlusZEnd += numELinks; else numELinksPlusXSideMinusZEnd += numELinks;
      }
    }

    // Bundles
    const std::map<int, std::unique_ptr<InnerBundle> >& bundles = myInnerCablingMap->getBundles();
    for (const auto& it : bundles) {
      const InnerBundle* myBundle = it.second.get();
      if (myBundle->isPositiveXSide()) {
	numBundlesOneXSide++;
	if (myBundle->isPositiveZEnd()) numBundlesPlusXSidePlusZEnd++; else numBundlesPlusXSideMinusZEnd++;
      }
    }

    // DTCs
    const std::map<int, std::unique_ptr<InnerDTC> >& dtcs = myInnerCablingMap->getDTCs();
    for (const auto& it : dtcs) {
      const InnerDTC* myDTC = it.second.get();
      if (myDTC->isPositiveXSide()) {
	numDTCsOneXSide++;
	if (myDTC->isPositiveZEnd()) numDTCsPlusXSidePlusZEnd++; else numDTCsPlusXSideMinusZEnd++;
      }
    }
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
  TProfile* Vizard::newProfile(TH1D* sourceHistogram, double xlow, double xup, int desiredNBins /* = 0 */) {
    if (!sourceHistogram) return nullptr;
    int nBins = sourceHistogram->GetNbinsX();
    if (desiredNBins != 0 && desiredNBins < nBins) nBins = desiredNBins;
    TProfile* resultProfile = new TProfile(Form("%s_profile",sourceHistogram->GetName()),
                                 sourceHistogram->GetTitle(),
                                 nBins,
                                 xlow,
                                 xup);
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

  TProfile& Vizard::newProfile_timesSin(const TGraph& sourceGraph, double xlow, double xup, int rebin /* = 1 */, int nBins) {
    TProfile* resultProfile;
    int nPoints = sourceGraph.GetN();
    // Rebin by factor 1 or user defined factor
    if (nBins==0) nPoints /= rebin;
    // Or set new number of bins
    else if (nBins <= nPoints) nPoints = nBins;
    resultProfile = new TProfile(Form("%s_timesSin_profile", sourceGraph.GetName()), sourceGraph.GetTitle(), nPoints, xlow, xup);
    double x, y;
    //double sintheta;
    for (int i=0; i<sourceGraph.GetN(); ++i) {
      sourceGraph.GetPoint(i, x, y);
      resultProfile->Fill(x, y/cosh(x));
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
    ss << "subdetectorId, z, rho, phi, detId, tt_list" << csv_eol;
    for (const auto& mapel : moduleConnections) {
      auto pos = mapel.first->posRef();
      ss << pos.subdetectorId << csv_separator << pos.z << csv_separator << pos.rho << csv_separator << pos.phi << csv_separator << mapel.second.detId();
      for (const auto& conn : mapel.second.connectedProcessors) {
        ss << csv_separator << 't' << conn.first << '_' << conn.second;
      }
      ss << csv_eol;
    }
    moduleConnectionsCsv_ = ss.str();
  }

  std::string Vizard::createAllModulesCsv(const Tracker& t, bool& withHeader) {
    TrackerVisitor v;
    if (withHeader) v.preVisit();
    t.accept(v);
    return v.output();
  }

  std::string Vizard::createBarrelModulesCsv(const Tracker& t) {
    BarrelVisitor v;
    v.preVisit();
    t.accept(v);
    return v.output();
  }

  std::string Vizard::createEndcapModulesCsv(const Tracker& t) {
    EndcapVisitor v;
    v.preVisit();
    t.accept(v);
    return v.output();
  }

  std::string Vizard::createModulesDetIdListCsv() {
    std::stringstream header;
    header << "DetId/U, BinaryDetId/B, Section/C, Layer/I, Ring/I, r_mm/D, z_mm/D, tiltAngle_deg/D, skewAngle_deg/D, phi_deg/D, meanWidth_mm/D, length_mm/D, sensorSpacing_mm/D, sensorThickness_mm/D" << std::endl;
    std::string detIdsListCsv = header.str();

    for (unsigned int i=0; i< trackers_.size(); ++i) {
      Tracker& tracker = *(trackers_.at(i));
      bool withHeader = false;
      detIdsListCsv += createAllModulesCsv(tracker, withHeader);
    }
    return detIdsListCsv;
  }

  std::string Vizard::createSensorsDetIdListCsv() {
    std::stringstream header;
    header << "DetId/U, BinaryDetId/B, Section/C, Layer/I, Ring/I, r_mm/D, z_mm/D, phi_deg/D" << std::endl;
    std::string detIdsListCsv = header.str();

    for (unsigned int i=0; i< trackers_.size(); ++i) {
      Tracker& tracker = *(trackers_.at(i));
      TrackerSensorVisitor v;
      tracker.accept(v);
      detIdsListCsv += v.output();
    }
    return detIdsListCsv;
  }


  std::string Vizard::createChemicalElementsCsv() {

    std::stringstream myCsv;
    myCsv << "Atomic Symbol /C, Density (g/cm^3) /D, Atomic Number /I, Standard Atomic Weight (u) /D, Radiation length (g/cm^2) /D, Interaction length (g/cm^2) /D" << std::endl;

    const MaterialsTable& myTable = MaterialsTable::instance();
    /*
    double density = myTable.density("CO2");
    double rl = myTable.radiationLength("CO2");
    double il = myTable.interactionLength("CO2");
    std::cout << "CO2" << " density = " << density << " rl = " << rl << " il = " << il << std::endl;
    */
    const ChemicalElementMap& allChemicalElements = myTable.getAllChemicalElements();

    for (const auto& elemIt : allChemicalElements) {
      const std::string elementName = elemIt.first;
      const ChemicalElement& elem = elemIt.second;
      myCsv << elementName << ","
	    << (elem.getDensity() * 1000.) << ","
	    << elem.getAtomicNumber() << ","
	    << elem.getAtomicWeight() << ","
	    << elem.getRadiationLength() << ","
	    << elem.getInteractionLength() 
	    << std::endl;
      //myCsv << std::endl;
    }

    return myCsv.str();
  }


  std::string Vizard::createChemicalMixturesCsv(const bool hasChemicalFormula) {

    std::stringstream myCsv;
    if (hasChemicalFormula) myCsv << "Compound Name /C, Density (g/cm^3) /D, Radiation length (g/cm^2) /D, Interaction length (g/cm^2) /D" << std::endl;
    else { myCsv << "Mixture Name /C,  Substance name /C, Substance mass fraction /D, Mixture Density (g/cm^3) /D, Mixture Radiation length (g/cm^2) /D, Mixture Interaction length (g/cm^2) /D" << std::endl; }

    const MaterialsTable& myTable = MaterialsTable::instance();
    const ChemicalMixtureMap& allChemicalMixtures = myTable.getAllChemicalMixtures();

    for (const auto& mixIt : allChemicalMixtures) {
      const std::string mixtureName = mixIt.first;
      const ChemicalMixture& mix = mixIt.second;

      if (mix.hasChemicalFormula() == hasChemicalFormula) {
	myCsv << mixtureName << ",";
	if (!hasChemicalFormula) myCsv << "," << ",";	      
	myCsv << (mix.getDensity() * 1000.) << ","     // g/ cm3
	      << mix.getRadiationLength() << ","
	      << mix.getInteractionLength() 
	      << std::endl;

	if (!hasChemicalFormula) {
	  const MassComposition& fractions = mix.getMassComposition();
	  for (const auto& fractionIt : fractions) {
	    myCsv << ","
		  << fractionIt.first << ","
		  << fractionIt.second
		  << std::endl;
	  }
	  myCsv << std::endl;
	  myCsv << std::endl;
	}
	
      }
    }

    return myCsv.str();
  }


  /* Create csv file (Outer Tracker), navigating from Module hierarchy level to DTC hierarchy level.
   */
  std::string Vizard::createModulesToDTCsCsv(const Tracker& tracker, const bool isPositiveCablingSide) {
    ModulesToDTCsVisitor v(isPositiveCablingSide);
    v.preVisit();
    tracker.accept(v);
    return v.output();
  }


  /* Create csv file (Outer Tracker), navigating from DTC hierarchy level to Module hierarchy level.
   */
  std::string Vizard::createDTCsToModulesCsv(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide) {

    std::stringstream modulesToDTCsCsv;
    modulesToDTCsCsv << "DTC name/C, DTC CMSSW Id/U, DTC Phi Sector Ref/I, type /C, DTC Slot/I, DTC Phi Sector Width_deg/D, MFC #/I, MFC type/C, MFB #/I, OPT Services Channel/I, PWR Services Channel/I, Module DetId/U, Module Section/C, Module Layer/I, Module Ring/I, Module phi_deg/D" << std::endl;

    const std::map<const std::string, std::unique_ptr<const OuterDTC> >& myDTCs = (isPositiveCablingSide ? 
										   myCablingMap->getDTCs() 
										   : myCablingMap->getNegDTCs());
    for (const auto& dtcIt : myDTCs) {
      const OuterDTC* myDTC = dtcIt.second.get();
      if (myDTC) {
	std::stringstream DTCInfo;
	DTCInfo << myDTC->name() << ", "
		<< myDTC->getCMSSWId() << ", "
		<< myDTC->phiSectorRef() << ", "
		<< any2str(myDTC->type()) << ", "
		<< myDTC->slot() << ", "
		<< std::fixed << std::setprecision(6)
		<< myDTC->phiSectorWidth() * 180. / M_PI << ", ";

	const std::vector<OuterCable*>& myCables = myDTC->cable();
	for (const auto& cable : myCables) {
	  std::stringstream cableInfo;
	  cableInfo << cable->myid() << ", "
		    << any2str(cable->type()) << ", ";
	  const ChannelSection* myOpticalSection = cable->opticalChannelSection();
	  const int opticalChannelNumber = myOpticalSection->channelNumber();
	  const ChannelSlot& opticalChannelSlot = myOpticalSection->channelSlot();

	  const std::vector<OuterBundle*>& myBundles = cable->bundles();
	  for (const auto& bundle : myBundles) {
	    std::stringstream bundleInfo;
	    bundleInfo << bundle->myid() << ", "
		       << opticalChannelNumber << " " 
		       << any2str(opticalChannelSlot) << ", "
		       << bundle->powerChannelSection()->channelNumber() << " " 
		       << any2str(bundle->powerChannelSection()->channelSlot()) << ", ";

	    const std::vector<Module*>& myModules = bundle->modules();
	    for (const auto& module : myModules) {
	      std::stringstream moduleInfo;
	      moduleInfo << module->myDetId() << ", "
			 << module->uniRef().subdetectorName << ", "
			 << module->uniRef().layer << ", "
			 << module->moduleRing() << ", "
			 << module->center().Phi() * 180. / M_PI;
	      modulesToDTCsCsv << DTCInfo.str() << cableInfo.str() << bundleInfo.str() << moduleInfo.str() << std::endl;
	    }
	    if (myModules.size() == 0) modulesToDTCsCsv << DTCInfo.str() << cableInfo.str() << bundleInfo.str() << std::endl;
	  }
	  if (myBundles.size() == 0) modulesToDTCsCsv << DTCInfo.str() << cableInfo.str() << std::endl;
	}
	if (myCables.size() == 0) modulesToDTCsCsv << DTCInfo.str() << std::endl;
      }
    }
    if (myDTCs.size() == 0) modulesToDTCsCsv << std::endl;

    return modulesToDTCsCsv.str();
  }


  /* Create csv file, navigating, in TEDD, from Bundle hierarchy level to Module hierarchy level.
     This also provides modules aggregation patterns. 
     A pattern provides the number of modules connected to each MFB fanout branch.
     The MFB fanout branches are sorted by increasing index.
     For example, for a given MFB, the pattern 3-4-3-2 means that the MFB is connected to:
     - 3 modules from fanout branch index 1.
     - 4 modules from fanout branch index 2.
     - 3 modules from fanout branch index 3.
     - 2 modules from fanout branch index 4.
   */
  std::string Vizard::createBundlesToEndcapModulesCsv(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide) {
    std::stringstream bundlesToEndcapModulesCsv;

    const std::string& summaryText = countBundlesToEndcapModulesCombinations(myCablingMap, isPositiveCablingSide);
    bundlesToEndcapModulesCsv << summaryText << std::endl;

    bundlesToEndcapModulesCsv << "MFB #/I, # Modules per fanout branch index, Module DetId/U, Module Section/C, Module Disk/I, Module Ring/I, Module phi_deg/D, Module Z (mm)/D, Fanout branch index/U" << std::endl;

    const std::map<const std::string, std::unique_ptr<const OuterDTC> >& myDTCs = (isPositiveCablingSide ? 
										   myCablingMap->getDTCs() 
										   : myCablingMap->getNegDTCs());
    for (const auto& dtcIt : myDTCs) {
      const OuterDTC* myDTC = dtcIt.second.get();
      if (myDTC) {

	const std::vector<OuterCable*>& myCables = myDTC->cable();
	for (const auto& cable : myCables) {

	  const std::vector<OuterBundle*>& myBundles = cable->bundles();
	  for (const auto& bundle : myBundles) {

	    std::string subDetectorName = bundle->subDetectorName();
	    // Only in TEDD.
	    if (subDetectorName == outer_cabling_tedd1 || subDetectorName == outer_cabling_tedd2) {
	      // Bundle related info.
	      std::stringstream bundleInfo;
	      bundleInfo << bundle->myid() << ",";

	      // Create pattern related to the bundle.
	      std::map<int, int> pattern;
	      std::vector<std::string> modulesInBundleInfo;
	      const std::vector<Module*>& myModules = bundle->modules();
	      for (const auto& module : myModules) {
		// Get which MFB fanout branch the module is connected to.
		const int fanoutBranchIndex = module->getEndcapFiberFanoutBranch();

		// Module related info.
		std::stringstream moduleInfo;
		moduleInfo << module->myDetId() << ", "
			   << module->uniRef().subdetectorName << ", "
			   << module->uniRef().layer << ", "
			   << module->moduleRing() << ", "
			   << module->center().Phi() * 180. / M_PI << ", "
			   << module->center().Z() << ", "
			   << fanoutBranchIndex;
		modulesInBundleInfo.push_back(moduleInfo.str());
		
		// Count the number of modules per fanout branch index.
		pattern[fanoutBranchIndex] += 1; 
	      }

	      // Checks whether pattern makes sense, and put it in a-b-c-d format.
	      std::stringstream patternInfo;
	      for (int fanoutBranchIndex = 1; fanoutBranchIndex <= outer_cabling_maxNumFanoutBranchesPerEndcapBundle; fanoutBranchIndex++) {
		auto found = pattern.find(fanoutBranchIndex);
		if (found != pattern.end()) {
		  if (fanoutBranchIndex != 1) patternInfo << "-";
		  const int numModulesPerDiskSurface = found->second;
		  patternInfo << numModulesPerDiskSurface;
		}
		else { logERROR("In TEDD, bundle " + any2str(bundle->myid()) 
				+ "has 0 module belonging to fanout branch " + any2str(fanoutBranchIndex)); }
	      }
	      patternInfo << ", ";
  
	      // Print info in csv file: bundle info + pattern info + associated modules info.
	      bundlesToEndcapModulesCsv << bundleInfo.str() << patternInfo.str();
	      const int numModulesInBundle = modulesInBundleInfo.size();
	      for (int i = 0; i < numModulesInBundle; i++) {
		// Set empty the first 2 columns, since they are the same for all modules belonging to a given bundle.
		if (i != 0) bundlesToEndcapModulesCsv << ", " << ", ";
		// List info from all modules belonging to the same bundle.
		bundlesToEndcapModulesCsv << modulesInBundleInfo.at(i) << std::endl;
	      }
	      if (myModules.size() == 0) bundlesToEndcapModulesCsv << std::endl;
	    }
	  }
	  if (myBundles.size() == 0) bundlesToEndcapModulesCsv << std::endl;
	}
	if (myCables.size() == 0) bundlesToEndcapModulesCsv << std::endl;
      }
    }
    if (myDTCs.size() == 0) bundlesToEndcapModulesCsv << std::endl;

    return bundlesToEndcapModulesCsv.str();
  }


  /* Provide a summary text file, with the distribution of modules aggregation patterns which are encountered in TEDD.
     A pattern is, for a given bundle, the number of connected modules per fanout branch.
     Here, patterns are irrespective of the fanout branch ordering.
     For example, 1-2-3-4 or 3-4-1-2 are both considered to be combination 1-2-3-4.
     All this is because Electronics/Mechanics will need, in TEDD, custom aggregation patch cords, 
     to group the fibers from each fanout branch into one bundle.
     One need to know how many customs aggregation patch cords are needed!
  */
  std::string Vizard::countBundlesToEndcapModulesCombinations(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide) {
    std::stringstream summaryText;
    summaryText << "Fanouts: # Modules per branch (Irrespective of branches ordering)" << std::endl;

    std::map<std::multiset<int>, int> combinationsDistribution;

    const std::map<const std::string, std::unique_ptr<const OuterDTC> >& myDTCs = (isPositiveCablingSide ? 
										   myCablingMap->getDTCs() 
										   : myCablingMap->getNegDTCs());
    for (const auto& dtcIt : myDTCs) {
      const OuterDTC* myDTC = dtcIt.second.get();
      if (myDTC) {

	const std::vector<OuterCable*>& myCables = myDTC->cable();
	for (const auto& cable : myCables) {

	  const std::vector<OuterBundle*>& myBundles = cable->bundles();
	  for (const auto& bundle : myBundles) {

	    std::string subDetectorName = bundle->subDetectorName();
	    // Only in TEDD.
	    if (subDetectorName == outer_cabling_tedd1 || subDetectorName == outer_cabling_tedd2) {
	      // Create pattern related to the bundle.
	      std::map<int, int> pattern;

	      const std::vector<Module*>& myModules = bundle->modules();
	      for (const auto& module : myModules) {
		// Get which MFB fanout branch the module belongs to.
		const int fanoutBranchIndex = module->getEndcapFiberFanoutBranch();
		// Count the number of modules per MFB fanout branch.
		pattern[fanoutBranchIndex] += 1; 
	      }

	      // Checks pattern makes sense, and create the corresponding combination.
	      // A combination is the number of modules per fanout branch, irrespective of branches ordering.
	      // One wants 1-2-3-4 and 3-4-1-2 to end up in the same combination: 1-2-3-4.
	      // Duplicates are allowed: combination 1-2-3-3 can happen!
	      std::multiset<int> combination;  
	      for (int fanoutBranchIndex = 1; fanoutBranchIndex <= outer_cabling_maxNumFanoutBranchesPerEndcapBundle; fanoutBranchIndex++) {
		auto found = pattern.find(fanoutBranchIndex);
		if (found != pattern.end()) {
		  const int numModulesPerDiskSurface = found->second;
		  // Create combination
		  combination.insert(numModulesPerDiskSurface);
		}
		else { logERROR("In TEDD, bundle " + any2str(bundle->myid()) 
				+ "has 0 module belonging to fanout branch " + any2str(fanoutBranchIndex)); }
	      }
	      // Count the occurences of each combination.
	      combinationsDistribution[combination] += 1;
	    }
	  }	 
	}
      }
    }

    // Print the different encountered combinations, and their occurences.
    for (const auto& comb : combinationsDistribution) {
      summaryText << "Combination";
      std::multiset<int> combination = comb.first;
      for (const auto& numModules : combination) {
	summaryText << " " << numModules;
      }
      summaryText << " appears " << comb.second << " times." << std::endl;
    } 
    summaryText << std::endl;

    return summaryText.str();
  }


  /**
   * Create Power Cables distribution csv file.
   * This provides a count of power cables based on their connections (how many modules per cable, of which type).
   */
  std::string Vizard::createPowerCablesDistributionCsv(const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide) {
    std::stringstream powerCablesDistributionCsv;
    powerCablesDistributionCsv << " , # Modules per power cable /U, # Power Cables /U" << std::endl;

    std::map<int, int> countBundlesPerPS10GConnections;
    std::map<int, int> countBundlesPerPS5GConnections;
    std::map<int, int> countBundlesPerSSConnections;

    // IN THE OT CABLING MAP, 1 OPTICAL BUNDLE <=> 1 POWER CABLE!!!!
    const std::map<const int, std::unique_ptr<OuterBundle> >& myBundles = (isPositiveCablingSide ? 
									   myCablingMap->getBundles() 
									   : myCablingMap->getNegBundles());
    // Fill counters
    for (const auto& myBundleIt : myBundles) {
      const OuterBundle* myBundle = myBundleIt.second.get();

      const Category& bundleType = myBundle->type();      
      const int numModulesPerBundle = myBundle->numModules();	  

      // PS 10G
      if (bundleType == Category::PS10G || bundleType == Category::PS10GA || bundleType == Category::PS10GB) {
	countBundlesPerPS10GConnections[numModulesPerBundle] += 1;
      }
      // PS 5G
      else if (bundleType == Category::PS5G) { 
	countBundlesPerPS5GConnections[numModulesPerBundle] += 1;
      }
      // 2S
      else if (bundleType == Category::SS) {
	countBundlesPerSSConnections[numModulesPerBundle] += 1;
      }
      else { 
	logERROR("Unknown bundle type: " + any2str(bundleType));
      }
    }

    int totalNumBundles = 0;
    // Print info
    // PS 10G
    powerCablesDistributionCsv << std::endl;
    for (const auto& countBundlesIt : countBundlesPerPS10GConnections) {
      powerCablesDistributionCsv << any2str(Category::PS10G) << ", " << countBundlesIt.first << ", " << countBundlesIt.second << std::endl;
      totalNumBundles += countBundlesIt.second; 
    }
    // PS 5G
    powerCablesDistributionCsv << std::endl;
    for (const auto& countBundlesIt : countBundlesPerPS5GConnections) {
      powerCablesDistributionCsv << any2str(Category::PS5G) << ", " << countBundlesIt.first << ", " << countBundlesIt.second << std::endl;
      totalNumBundles += countBundlesIt.second; 
    }
    // 2S
    powerCablesDistributionCsv << std::endl;
    for (const auto& countBundlesIt : countBundlesPerSSConnections) {
      powerCablesDistributionCsv << any2str(Category::SS) << ", " << countBundlesIt.first << ", " << countBundlesIt.second << std::endl;
      totalNumBundles += countBundlesIt.second; 
    }

    powerCablesDistributionCsv << std::endl;
    powerCablesDistributionCsv << " , Total, " << totalNumBundles << std::endl;

    return powerCablesDistributionCsv.str();
  }


  /* Create csv file (Outer Tracker), summary on both cabling sides. Info needed by CMSSW: Modules DetIds to DTCIds.
   */
  std::string Vizard::createCMSSWOuterTrackerCablingMapCsv(const Tracker& tracker) {
    CMSSWOuterTrackerCablingMapVisitor v;
    v.preVisit();
    tracker.accept(v);
    return v.output();
  }


  /* Create csv file (Inner Tracker), navigating from Module hierarchy level to DTC hierarchy level.
   */
  std::string Vizard::createInnerTrackerModulesToDTCsCsv(const Tracker& tracker) {
    InnerTrackerModulesToDTCsVisitor v;
    v.preVisit();
    tracker.accept(v);
    return v.output();
  }


  /* Create csv file (Inner Tracker), navigating from DTC hierarchy level to Module hierarchy level.
   */
  std::string Vizard::createInnerTrackerDTCsToModulesCsv(const InnerCablingMap* myInnerCablingMap) {

    std::stringstream dtcsToModulesCsv;
    dtcsToModulesCsv << "(+Z) End ?/Boolean, (+X) Side?/Boolean, DTC #/I, MFB #/I, LP GBT #/C, # ELinks Per Module/I, Power Chain #/I, Power Chain Type/C, Long Barrel ?/Boolean, Module DetId/U, Module Section/C, Module Layer/I, Module Ring/I, Module phi_deg/D" << std::endl;

    const std::map<int, std::unique_ptr<InnerDTC> >& myDTCs = myInnerCablingMap->getDTCs();
    for (const auto& itDTC : myDTCs) {
      InnerDTC* myDTC = itDTC.second.get();
      if (myDTC) {
	std::stringstream DTCInfo;
	DTCInfo << myDTC->isPositiveZEnd() << ","
		<< myDTC->isPositiveXSide() << ","
		<< myDTC->myid() << ",";

	const std::vector<InnerBundle*>& myBundles = myDTC->bundles();
	for (const auto& myBundle : myBundles) {
	  std::stringstream bundleInfo;
	  bundleInfo << myBundle->myid() << ",";

	  const std::vector<GBT*>& myGBTs = myBundle->GBTs();
	  for (const auto& myGBT : myGBTs) {
	    std::stringstream GBTInfo;
	    GBTInfo << any2str(myGBT->GBTId()) << ","
		    << myGBT->numELinksPerModule() << ",";

	    const PowerChain* myPowerChain = myGBT->getPowerChain();
	    if (myPowerChain != nullptr) {
	      std::stringstream powerChainInfo;
	      powerChainInfo << myPowerChain->myid() << ","
			     << any2str(myPowerChain->powerChainType()) << ","
			     << any2str(myPowerChain->isBarrelLong()) << ",";

	      const std::vector<Module*>& myModules = myGBT->modules();
	      for (const auto& module : myModules) {
		std::stringstream moduleInfo;
		moduleInfo << module->myDetId() << ", "
			   << module->uniRef().subdetectorName << ", "
			   << module->uniRef().layer << ", "
			   << module->moduleRing() << ", "
			   << module->center().Phi() * 180. / M_PI;
		dtcsToModulesCsv << DTCInfo.str() << bundleInfo.str() << GBTInfo.str() << powerChainInfo.str() << moduleInfo.str() << std::endl;
	      }
	      if (myModules.size() == 0) dtcsToModulesCsv << DTCInfo.str() << bundleInfo.str() << GBTInfo.str() << powerChainInfo.str() << std::endl;
	    }
	    else dtcsToModulesCsv << DTCInfo.str() << bundleInfo.str() << GBTInfo.str() << std::endl;
	  }
	  if (myGBTs.size() == 0) dtcsToModulesCsv << DTCInfo.str() << bundleInfo.str() << std::endl;
	}
	if (myBundles.size() == 0) dtcsToModulesCsv << DTCInfo.str() << std::endl;
      }
    }
    if (myDTCs.size() == 0) dtcsToModulesCsv << std::endl;

    return dtcsToModulesCsv.str();
  }


  /*
   * Draw circle of a given radius, fill, and color.
   */
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


  /*
  * Draw spider net to delimit the Phi Sectors.
  * bool isRotatedY180 : 
  - false: draws in CMS global frame of reference.
  - true: draws in CMS global frame of reference rotated of 180° around CMS_Y.
  */
  void Vizard::drawPhiSectorsBoundaries(const double phiSectorWidth, const bool isRotatedY180) {
    int numPhiSectors = round(2. * M_PI / phiSectorWidth);
    double phiSectorBoundaryRadius = 2 * vis_min_canvas_sizeX;

    for (int i = 0; i < numPhiSectors; i++) {
      const double angle = i * phiSectorWidth;
      const double rotatedAngle = (isRotatedY180 ?  M_PI - angle : angle);
      TLine* line = new TLine(0., 0., phiSectorBoundaryRadius * cos(rotatedAngle), phiSectorBoundaryRadius * sin(rotatedAngle)); 
      line->SetLineWidth(2); 
      line->Draw("same");     
    }
    drawFrameOfReference(isRotatedY180);
  }


  /*
   *  Draw frame of reference reminder.
   *  It is drawn on the top right corner of the plot.
   */
  void Vizard::drawFrameOfReference(const bool isRotatedY180, const double scalingFactor) {

    const double arrowMin = 900 * scalingFactor;
    const double arrowMax = 1100 * scalingFactor;
    const double circleZRadius = 50 * scalingFactor;
    const double arrowWidth = 0.02;
    const double textSize = 0.025;
    
    // CMS reference frame of reference.
    if (!isRotatedY180) {
      TArrow* arrowX = new TArrow(arrowMin, arrowMin, arrowMax, arrowMin, arrowWidth, "|>");
      arrowX->Draw();
      const double textXAbs = 1000 * scalingFactor;
      const double textXOrd = 820 * scalingFactor;
      TLatex* textX = new TLatex(textXAbs, textXOrd, "X");
      textX->SetTextSize(textSize);
      textX->Draw("same");

      TArrow* arrowY = new TArrow(arrowMin, arrowMin, arrowMin, arrowMax, arrowWidth, "|>");
      arrowY->Draw();
      const double textYAbs = 820 * scalingFactor;
      const double textYOrd = 1000 * scalingFactor;
      TLatex* textY = new TLatex(textYAbs, textYOrd, "Y");
      textY->SetTextSize(textSize);
      textY->Draw("same");

      const double circleZCentre = 1050 * scalingFactor;
      const double pointZRadius = 10 * scalingFactor;
      TEllipse* circleZ = new TEllipse(circleZCentre, circleZCentre, circleZRadius, circleZRadius);
      circleZ->SetLineWidth(2);
      circleZ->Draw("same");
      TEllipse* pointZ = new TEllipse(circleZCentre, circleZCentre, pointZRadius, pointZRadius);
      pointZ->SetFillColor(kBlack);
      pointZ->Draw("same");
      const double textZAbs = 1030 * scalingFactor;
      const double textZOrd = 1110 * scalingFactor;
      TLatex* textZ = new TLatex(textZAbs, textZOrd, "Z");
      textZ->SetTextSize(textSize);
      textZ->Draw("same");
    }

    // CMS frame of reference rotated by 180 degrees around CMS_Y.
    else {
      TArrow* arrowX = new TArrow(arrowMax, arrowMin, arrowMin, arrowMin, arrowWidth, "|>");
      arrowX->Draw();
      const double textXAbs = 950 * scalingFactor;
      const double textXOrd = 820 * scalingFactor;
      TLatex* textX = new TLatex(textXAbs, textXOrd, "X");
      textX->SetTextSize(textSize);
      textX->Draw("same");

      TArrow* arrowY = new TArrow(arrowMax, arrowMin, arrowMax, arrowMax, arrowWidth, "|>");
      arrowY->Draw();
      const double textYAbs = 1020 * scalingFactor;
      const double textYOrd = 1000 * scalingFactor;
      TLatex* textY = new TLatex(textYAbs, textYOrd, "Y");
      textY->SetTextSize(textSize);
      textY->Draw("same");

      const double circleZCentreAbs = 850 * scalingFactor;
      const double circleZCentreOrd = 1050 * scalingFactor;
      TEllipse* circleZ = new TEllipse(circleZCentreAbs, circleZCentreOrd, circleZRadius, circleZRadius);
      circleZ->SetLineWidth(2);
      circleZ->Draw("same");
      const double crossMinAbs = 815 * scalingFactor;
      const double crossMaxAbs = 885 * scalingFactor;
      const double crossMinOrd = 1015 * scalingFactor;
      const double crossMaxOrd = 1085 * scalingFactor;
      TLine* lineU = new TLine(crossMinAbs, crossMinOrd, crossMaxAbs, crossMaxOrd);
      lineU->SetLineWidth(2);
      lineU->Draw("same");
      TLine* lineD = new TLine(crossMinAbs, crossMaxOrd, crossMaxAbs, crossMinOrd);
      lineD->SetLineWidth(2);
      lineD->Draw("same");
      const double textZAbs = 830 * scalingFactor;
      const double textZOrd = 1110 * scalingFactor;
      TLatex* textZ = new TLatex(textZAbs, textZOrd, "Z");
      textZ->SetTextSize(textSize);
      textZ->Draw("same");
    }
  }


  /*
   * Compute scaling factor to place the frame of reference reminder at the right position (top right corner of the plot).
   */
  const std::pair<double, double> Vizard::computeInnerCablingPlotsScalingFactors(const Tracker& tracker) {
    const std::pair<double, double> maxRadii = computeInnerCablingPlotsMaxRadii(tracker);
    const double barrelViewPort = maxRadii.first;
    const double forwardViewPort = maxRadii.second;

    const double barrelScalingFactor = barrelViewPort / geom_outer_strip_radius;  
    const double forwardScalingFactor = forwardViewPort / geom_outer_strip_radius; 

    return std::make_pair(barrelScalingFactor, forwardScalingFactor);
  }


  /*
   * Compute the max radii to be drawn on Barrel / Endcaps plots.
   */
  const std::pair<double, double> Vizard::computeInnerCablingPlotsMaxRadii(const Tracker& tracker) {
    double barrelViewPort, forwardViewPort;
    if (tracker.isPixelTracker()) {
      if (tracker.barrels().size() > 0 && tracker.endcaps().size() >= 2) {
	barrelViewPort = tracker.barrels().at(0).maxR() * 1.1;
	forwardViewPort = tracker.endcaps().at(1).maxR() * 1.1;
      }
      else { logERROR("Vizard::computeInnerCablingPlotsMaxRadii : Unexpected number of subdetectors in Inner Tracker."); }
    }
    else { logERROR("Vizard::computeInnerCablingPlotsMaxRadii : Should be used on Inner Tracker only."); }

    return std::make_pair(barrelViewPort, forwardViewPort);
  }


  /*
   * Compute colored legend for services channels.
   * This just lists all the cables associated to a given cabling side, and sum up all the encountered colors.
   * isPowerCabling is added because in case of non-optical cabling, one wants also the possibility of transparent colors.
   * Transparent colors are used to distinguish channels sections A and C, which are specific to power cabling.
   * TO DO: Would be nicer to have this drawn on the fly while the plots are created.
   */
  void Vizard::computeServicesChannelsLegend(TLegend* legend, const OuterCablingMap* myCablingMap, const bool isPositiveCablingSide, const bool isPowerCabling) {
    std::map<std::string, int > channelsColors;

    if (!isPowerCabling) {
      // Only consider the relevant cables: cables from (+Z) side or (-Z) side.
      const std::map<const int, std::unique_ptr<OuterCable> >& cables = (isPositiveCablingSide ? 
									 myCablingMap->getCables() 
									 : myCablingMap->getNegCables());

      // Loop on all the encountered cables
      for (const auto& myCableIt : cables) {
	const OuterCable* myCable = myCableIt.second.get();
	const ChannelSection* mySection = myCable->opticalChannelSection();
	const int& myChannelNumber = mySection->channelNumber();
	const int& myPlotColor = mySection->plotColor();

	// This is simply to add 0 in front of single-digit numbers, so that the sorting directly makes sense.
	std::stringstream channelNameStream;
	channelNameStream << "OT";
	// Find single-digit numbers
	if (fabs(myChannelNumber) <= 9) {	
	  if (myChannelNumber >= 0) channelNameStream << "0" << myChannelNumber; // Add 0 in front of positive digit	
	  else channelNameStream << "-0" << fabs(myChannelNumber); // Add -0 in front of negative digit
	}
	else channelNameStream << myChannelNumber;

	// If the legend is for power cabling, one need to distinguish sections A and C.
	const ChannelSlot& mySlot = mySection->channelSlot();
	channelNameStream << any2str(mySlot);
	channelNameStream << std::endl;
	const std::string channelName = channelNameStream.str();

	// ADD CHANNEL COLOR TO THE MAP
	if (channelsColors.find(channelName) == channelsColors.end()) {
	  channelsColors[channelName] = myPlotColor;
	}
      }
    }

    else {
      // Only consider the relevant bundles: bundles from (+Z) side or (-Z) side.
      const std::map<const int, std::unique_ptr<OuterBundle> >& bundles = (isPositiveCablingSide ? 
									   myCablingMap->getBundles() 
									   : myCablingMap->getNegBundles());

      // Loop on all the encountered bundles
      for (const auto& myBundleIt : bundles) {
	const OuterBundle* myBundle = myBundleIt.second.get();
	const ChannelSection* mySection = myBundle->powerChannelSection();
	const int& myChannelNumber = mySection->channelNumber();
	const int& myPlotColor = mySection->plotColor();

	// This is simply to add 0 in front of single-digit numbers, so that the sorting directly makes sense.
	std::stringstream channelNameStream;
	channelNameStream << "OT";
	// Find single-digit numbers
	if (fabs(myChannelNumber) <= 9) {	
	  if (myChannelNumber >= 0) channelNameStream << "0" << myChannelNumber; // Add 0 in front of positive digit	
	  else channelNameStream << "-0" << fabs(myChannelNumber); // Add -0 in front of negative digit
	}
	else channelNameStream << myChannelNumber;

	// If the legend is for power cabling, one need to distinguish sections A and C.
	const ChannelSlot& mySlot = mySection->channelSlot();
	channelNameStream << any2str(mySlot);
	channelNameStream << std::endl;
	const std::string channelName = channelNameStream.str();

	// ADD CHANNEL COLOR TO THE MAP
	if (channelsColors.find(channelName) == channelsColors.end()) {
	  channelsColors[channelName] = myPlotColor;
	}
      }
    }


    // Create legenda
    for (const auto& it : channelsColors) {
      const std::string& channelName = it.first;
      const int& myPlotColor = it.second;

      // Just fakely used to add an entry, not drawn!
      Double_t x[1] = {0.};
      Double_t y[1] = {0.};
      TPolyLine* line = new TPolyLine(1, x, y);
      // Obtain the channel color.
      const bool isTransparentActivated = isPowerCabling;
      line->SetLineColor(Palette::colorChannel(myPlotColor, isTransparentActivated));
      line->SetFillColor(Palette::colorChannel(myPlotColor, isTransparentActivated));
      legend->AddEntry(line, channelName.c_str(), "f");
    }
  }


  WeightsPerSubdetector Vizard::computeDetailedWeights(MaterialBudget& materialBudget, RootWPage& myPage) {
    Tracker& myTracker = materialBudget.getTracker();
    //const bool isIT = myTracker.isPixelTracker();
    std::string myTrackerName = myTracker.myid();
    RootWContent& myContent = myPage.addContent("All volumes details");
 

    // Graphic representation of the materials volumes in the (RZ) plane
    double maxR = myTracker.maxR()*1.2;
    double maxZ = myTracker.maxZ()*1.2;
    std::unique_ptr<TCanvas> allVolumesCanvas(new TCanvas("allVolumesCanvas", "allVolumesCanvas")); // TODO Factory for canvases?!
    allVolumesCanvas->cd();
    TH2D* allVolumesPlot = new TH2D("allVolumesPlot", ";z [mm];r [mm]", 200, -maxZ, maxZ, 100, 0, maxR);
    allVolumesPlot->Draw();
  

    std::stringstream allVolumesStream;
    allVolumesStream << "volume_ID/I,sim_category/C,z1 (mm)/D,z2 (mm)/D,r1 (mm)/D,r2 (mm)/D,volume_RL/D,volume_IL/D,elementID/I,subdetector/C,mechanical_category/C,Component/C,Element/C,mass (g)/D,mass_per_length (g/mm)/D" << std::endl;

    std::stringstream modulesStream;
    modulesStream << "subdetector/C,layer_or_disk/I,ring/I,sim_category/C,volume_RL/D,volume_IL/D, mechanical_category/C,Component/C,Element/C,mass (g)/D" << std::endl;


    WeightsPerSubdetector totalWeights;
    //std::vector<MaterialProperties*> totalMaterials;
    //Vector<unique_ptr<Base> >.  vec.emplace_back(new Derived())

    std::vector<int> allColors;
    allColors.push_back(kOrange);
    allColors.push_back(kCyan);
    allColors.push_back(kRed);
    allColors.push_back(kGreen);
    allColors.push_back(kAzure + 1);
    int colorIndex = 0;
    std::map<std::string, int> subdetectorColors;

    // SERVICES
    bool isModule = false;
    std::vector<InactiveElement> allServices = materialBudget.getAllServices(); 
   
    // Counting services with an ad-hoc index
    int serviceId = 0;
    maxZ = 0.; maxR = 0.;
    for (auto& serviceIt : allServices) {
      double z1 = serviceIt.getZOffset();
      double z2 = serviceIt.getZOffset()+serviceIt.getZLength();
      double r1 = serviceIt.getInnerRadius();
      double r2 = serviceIt.getInnerRadius()+serviceIt.getRWidth();
      double serviceLength = serviceIt.getLength();
      double rl = serviceIt.getRadiationLength();
      double il = serviceIt.getInteractionLength();

      // Update the maxZ and maxR with respect to the inactive surfaces
      if (fabs(z1)>maxZ) maxZ=fabs(z1);
      if (fabs(z2)>maxZ) maxZ=fabs(z2);
      if (fabs(r1)>maxR) maxR=fabs(r1);
      if (fabs(r2)>maxR) maxR=fabs(r2);
  
     
      const std::map<LocalElement, double, ElementNameCompare>& allMasses = serviceIt.getLocalElementsDetails();
      plotAndPrintVolumeMaterials(totalWeights, allVolumesStream, modulesStream, 
				  allMasses, z1, z2, r1, r2, rl, il,
				  subdetectorColors, allColors, colorIndex,
				  isModule, serviceId, serviceLength
				  );
    
      serviceId++;
    }


    // MODULE CAPS
    isModule = true;
    const std::vector<std::vector<ModuleCap> >& barrelModules = materialBudget.getBarrelModuleCaps();
    const std::vector<std::vector<ModuleCap> >& endcapModules = materialBudget.getEndcapModuleCaps(); 
    std::vector<ModuleCap> allModules;
    for (const auto& barrelIt : barrelModules) allModules.insert(allModules.end(), barrelIt.begin(),  barrelIt.end());
    for (const auto& endcapIt : endcapModules) allModules.insert(allModules.end(), endcapIt.begin(),  endcapIt.end());

    std::set<std::tuple<std::string, int, int> > allModuleMaterialsRefs;

    for (auto& moduleIt : allModules) {
      const Module* detectorModule = &(moduleIt.getModule());
      double z1 = detectorModule->minZ();
      double z2 = detectorModule->maxZ();
      double r1 = detectorModule->minR();
      double r2 = detectorModule->maxR();
      //length = iter.getLength();
      double rl = moduleIt.getRadiationLength();
      double il = moduleIt.getInteractionLength();

      // Update the maxZ and maxR with respect to the inactive surfaces
      if (fabs(z1)>maxZ) maxZ=fabs(z1);
      if (fabs(z2)>maxZ) maxZ=fabs(z2);
      if (fabs(r1)>maxR) maxR=fabs(r1);
      if (fabs(r2)>maxR) maxR=fabs(r2);

      // Find out whether the module info will be printed in the modules dedicated csv file.
      // Indeed, only one module per subdetector + layer + ring is printed.
      bool printModulesCsv = false;
      if (detectorModule != nullptr) {
	const std::string subdetectorName = detectorModule->uniRef().subdetectorName;
	const int layerOrDiskIndex = detectorModule->uniRef().layer;
	const int ringIndex = detectorModule->uniRef().ring;

	const std::tuple<std::string, int, int> moduleMaterialsRef = std::make_tuple(subdetectorName, layerOrDiskIndex, ringIndex);
	if (allModuleMaterialsRefs.find(moduleMaterialsRef) == allModuleMaterialsRefs.end()) {
	  printModulesCsv = true;
	  allModuleMaterialsRefs.insert(moduleMaterialsRef);
	}
      }

      const std::map<LocalElement, double, ElementNameCompare>& allMasses = moduleIt.getLocalElementsDetails();
      plotAndPrintVolumeMaterials(totalWeights, allVolumesStream, modulesStream, 
				  allMasses, z1, z2, r1, r2, rl, il,
				  subdetectorColors, allColors, colorIndex,
				  isModule, 0, 0., detectorModule, printModulesCsv
				  );
    }


    // ADD PLOT AND CSV FILE TO WEBSITE
    allVolumesPlot->GetXaxis()->SetRangeUser(-maxZ, maxZ);
    allVolumesPlot->GetYaxis()->SetRangeUser(0, maxR);

    RootWImage& servicesImage = myContent.addImage(std::move(allVolumesCanvas), vis_max_canvas_sizeX, vis_min_canvas_sizeY);
    servicesImage.setComment("All material volumes, (RZ) view.");
    servicesImage.setName("AllMaterialVolumesRZ");

    RootWTextFile* allVolumesFile = new RootWTextFile(Form("allVolumesMaterials_%s.csv", myTrackerName.c_str()), "All volumes weights");
    allVolumesFile->addText(allVolumesStream.str());
    myContent.addItem(allVolumesFile);

    RootWTextFile* modulesFile = new RootWTextFile(Form("modulesMaterials_%s.csv", myTrackerName.c_str()), "Modules weights");
    modulesFile->addText(modulesStream.str());
    myContent.addItem(modulesFile);


    return totalWeights;
  }




  void Vizard::plotAndPrintVolumeMaterials(WeightsPerSubdetector& totalWeights, std::stringstream& allVolumesStream, std::stringstream& modulesStream, 
					   const std::map<LocalElement, double, ElementNameCompare>& allMasses, 
					   const double z1, const double z2, const double r1, const double r2, const double rl, const double il,
					   std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
					   const bool isModule, const int serviceId, const double serviceLength, 
					   const Module* detectorModule, const bool printModulesCsv) {

    bool isEmpty = true;
    int elementId = 0;
    std::set<std::string> volumeSubdetectorNames;
    const uint32_t detId = ( (isModule && detectorModule != nullptr) ? detectorModule->myDetId() : 0);
    const int layerOrDiskIndex = ( (isModule && detectorModule != nullptr) ? detectorModule->uniRef().layer : 0);
    const int ringIndex = ( (isModule && detectorModule != nullptr) ? detectorModule->uniRef().ring : 0);
    double moduleMass = 0.;

    // LOOP ON ALL LOCAL ELEMENTS
    for (const auto& massIt : allMasses) {

      const LocalElement& myElement = massIt.first;

      const std::string subdetectorName = myElement.subdetectorName();
      volumeSubdetectorNames.insert(subdetectorName);

      const MechanicalCategory& mechanicalCategory = myElement.mechanicalCategory();
      std::string tableMechanicalCategory = any2str(mechanicalCategory);
      if (mechanicalCategory == MechanicalCategory::COOLING || mechanicalCategory == MechanicalCategory::SUPPORT) tableMechanicalCategory = "SUPPORTS & COOLING";

      const std::string componentName = myElement.componentName();
      const std::string elementName = myElement.elementName();

      const double mass = massIt.second;
	    
      if (mass > 0.) isEmpty = false;


      // ALL VOLUMES: COMPUTE WEIGHTS TOTALS
      totalWeights[subdetectorName][tableMechanicalCategory][componentName] += mass;
	

      // ALL VOLUMES: DETAILED CSV OUTPUT
      allVolumesStream << (isModule ? detId : serviceId) << ","
		       << (isModule ? "Module" : "Service") << ","
		       << z1 << ","
		       << z2 << ","
		       << r1 << ","
		       << r2 << ","
		       << rl << ","
		       << il << ","
		       << elementId++ << ","
		       << subdetectorName << ","
		       << any2str(mechanicalCategory) << ","
		       << componentName << ","	       
		       << elementName << ","
		       << mass;
      if (!isModule) allVolumesStream << "," << (mass / serviceLength);
      allVolumesStream << std::endl;


      if (printModulesCsv) {
	// MODULE VOLUMES: DETAILED CSV OUTPUT
	modulesStream << subdetectorName << ","
		      << layerOrDiskIndex << ","
		      << ringIndex << ","
		      << "Module" << ","
		      << rl << ","
		      << il << ","
		      << any2str(mechanicalCategory) << ","
		      << componentName << ","	       
		      << elementName << ","
		      << mass
		      << std::endl;
	moduleMass += mass;
      }
    }
    allVolumesStream << std::endl;
    if (printModulesCsv) {
      modulesStream <<  "," << "," << "," << "," << "," << "," << "," << "," << "Total " << "," << moduleMass << std::endl;
      modulesStream << std::endl;
    }


    // ALL VOLUMES: PLOT
    if (volumeSubdetectorNames.size() > 4) { 
      std::cout << "!!! More than 4 subdetectors assigned to a materials volume." << std::endl; 
    }
      
    else if (volumeSubdetectorNames.size() == 0) { plotVolumeBox("", subdetectorColors, allColors, colorIndex, isEmpty, z1, z2, r1, r2); }

    else {
      double plotR1 = r1;
      double plotR2 = r2;
      for (const auto& subdetectorName : volumeSubdetectorNames) {  

	if (isModule && detectorModule != nullptr) {
	  if (detectorModule->isTilted() 
	      && detectorModule->subdet() == ModuleSubdetector::BARREL
	      ) {
	    const double rhoAtMinZ = (z1 > 0. ? r2 : r1); 
	    const double rhoAtMaxZ = (z1 > 0. ? r1 : r2);
	    const bool isFilled = false;
	    plotVolumeBox(subdetectorName, subdetectorColors, allColors, colorIndex, isEmpty, z1, z2, rhoAtMinZ, rhoAtMaxZ, isFilled);
	  }
	  else plotVolumeBox(subdetectorName, subdetectorColors, allColors, colorIndex, isEmpty, z1, z2, r1, r2);
	}

	else {
	  plotVolumeBox(subdetectorName, subdetectorColors, allColors, colorIndex, isEmpty, z1, z2, plotR1, plotR2);
	  if (!isModule) {
	    plotR1 += 2;
	    plotR2 += 2;
	  }
	}
      }
    }
    

  }

  void Vizard::plotVolumeBox(const std::string subdetectorName, 
			     std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
			     const bool isEmpty, 
			     const double z1, const double z2, const double r1, const double r2, const bool isFilled) {

    const int color = computeSubdetectorColor(subdetectorName, subdetectorColors, allColors, colorIndex, isEmpty);
   
    if (isFilled) {
      TBox* myBox = new TBox(z1, r1, z2, r2);
      myBox->SetLineColor(color);
      myBox->SetFillStyle(3003);
      myBox->SetFillColor(color);
      myBox->Draw("l");
    }
    else {
      TLine* myLine = new TLine(z1, r1, z2, r2);
      myLine->SetLineColor(color);
      myLine->Draw("l");
    }
  }



  const int Vizard::computeSubdetectorColor(const std::string subdetectorName,
					    std::map<std::string, int>& subdetectorColors, const std::vector<int>& allColors, int& colorIndex,
					    const bool isEmpty) {
    int color;

    if (!isEmpty) {
      if (subdetectorName == "") color = kGray;
      else {
	const auto& found = subdetectorColors.find(subdetectorName);
	if (found != subdetectorColors.end()) {
	  color = found->second;
	}
	else {
	  const int numColors = allColors.size();
	  if (colorIndex < numColors) {
	    color = allColors.at(colorIndex);
	    colorIndex++;
	    subdetectorColors.insert(std::make_pair(subdetectorName, color));
	  }
	  else { logERROR("Not enough colors are defined with respect to the total number of subdetectors"); }
	}
	//if (subdetectorName == "TBPS" || subdetectorName == "PXB") color = kAzure + 1;
	//else if (subdetectorName == "TB2S") color = kCyan;
	//else if (subdetectorName == "TEDD_1" || subdetectorName == "FPIX_1") color = kRed;
	// else if (subdetectorName == "TEDD_2" || subdetectorName == "FPIX_2") color = kOrange;
	//else if (subdetectorName == "OTST" || subdetectorName == "ITST") color = kGreen;
	//else color = kMagenta;
      }
    }
    else { color = kBlack; }

    return color;
  }





  // Create an extra tab for XML files linking
  void Vizard::createXmlSite(RootWSite& site, std::string xmlDir, std::string layoutDir) {
    RootWPage* myPage = new RootWPage("XML");
    myPage->setAddress("xml.html");
    site.addPage(myPage);

    std::vector<std::string> origMetadataXmlFiles, origPixelXmlFiles, origOuterTrackerXmlFiles;
    std::vector<std::string> metadataXmlFileNames, pixelXmlFileNames, outerTrackerXmlFileNames;

    try {
      boost::filesystem::directory_iterator end_iter;
      for (boost::filesystem::directory_iterator dir_iter(xmlDir); dir_iter != end_iter; dir_iter++) {
	if (boost::filesystem::is_regular_file(dir_iter->path())) {
	  std::string origFile = dir_iter->path().string();
	  std::string fileName = dir_iter->path().filename().string();

	  if (fileName.find(".cfg") != std::string::npos ) {
	    origMetadataXmlFiles.push_back(origFile);
	    metadataXmlFileNames.push_back(fileName);
	  }
	  else if (fileName.find(".xml") != std::string::npos ) {
	    if (fileName.find("pixel") != std::string::npos ) {
	      origPixelXmlFiles.push_back(origFile);
	      pixelXmlFileNames.push_back(fileName);
	    }
	    else {
	      origOuterTrackerXmlFiles.push_back(origFile);
	      outerTrackerXmlFileNames.push_back(fileName);
	    }
	  }
	}
      }
    }
    catch (boost::filesystem::filesystem_error e) {
      cerr << e.what() << " when trying to copy XML files from XML directory to website directory." << endl;
    }

    RootWContent* content = new RootWContent("XML files");
    if (!metadataXmlFileNames.empty()) {
      RootWBinaryFileList* metadataXmlFileList = new RootWBinaryFileList(metadataXmlFileNames.begin(), metadataXmlFileNames.end(), "XML generation Metadata", origMetadataXmlFiles.begin(), origMetadataXmlFiles.end());
      content->addItem(metadataXmlFileList);
    }
    if (!pixelXmlFileNames.empty()) {
      RootWBinaryFileList* pixelXmlFileList = new RootWBinaryFileList(pixelXmlFileNames.begin(), pixelXmlFileNames.end(), "XML for Pixel", origPixelXmlFiles.begin(), origPixelXmlFiles.end());
      content->addItem(pixelXmlFileList);
    }
    if (!outerTrackerXmlFileNames.empty()) {
      RootWBinaryFileList* outerTrackerXmlFileList = new RootWBinaryFileList(outerTrackerXmlFileNames.begin(), outerTrackerXmlFileNames.end(), "XML for Outer Tracker", origOuterTrackerXmlFiles.begin(), origOuterTrackerXmlFiles.end());
      content->addItem(outerTrackerXmlFileList);
    }
    myPage->addContent(content);
  }
  // NB : XML files are copied from the XML directory to the www/layout directory with RootWBinaryFileList::dump.

}
