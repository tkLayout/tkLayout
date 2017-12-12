// Standard c++ libraries
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <string>
#include <cstdio>

// BOOST
#include <boost/program_options.hpp>

// Private libraries
#include <PtError.hh>
#include <module.hh>

using namespace std;
namespace po = boost::program_options;

void syntax(char* programName) {
  //void syntax(char* po::options_description& visible, char* programName) {
  std::cout << "Syntax: " << programName << endl
            << " [Barrel|Endcap]" << std::endl
            << " mod_d_start mod_d_stop mod_d_step" << std::endl
            << " pitch strip_l mod_z mod_r window efficient_pt"<< std::endl;
  std::cout << "Pt units : GeV/c" << std::endl;
  std::cout << "Pitch    : micrometers" << std::endl;
  std::cout << "mod_*    : millimeters" << std::endl;
  std::cout << "window   : number of strips" << std::endl;
  //std::cout << "And then the show: " << visible << std::endl;
}

enum { Momentum, Distance };

int main(int argc, char*argv[]) {
  cout << "# Tuning Thinkness" << endl;
  if (argc!=11) {
    syntax(argv[0]);
    return -1;
  }

  int ModuleType;
  if (strcmp(argv[1], "Barrel")==0) {
    ModuleType = Module::Barrel;
    std::cout << "# Module type = Barrel" << std::endl;
  } else if (strcmp(argv[1], "Endcap")==0) {
    ModuleType = Module::Endcap;
    std::cout << "# Module type = Endcap" << std::endl;
  } else {
    syntax(argv[0]);
    return 0;
  }

  double start = atof(argv[2]);
  double stop = atof(argv[3]);
  double step = atof(argv[4]);

  double pitch = atof(argv[5])/1000.;
  double strip_l = atof(argv[6]);
  double mod_z = atof(argv[7]);
  double mod_r = atof(argv[8]);
  double window = atof(argv[9]);
  double efficient_pt = atof(argv[10]);

  ptError myError;
  myError.setModuleType(ModuleType);
  myError.setPitch(pitch);
  myError.setStripLength(strip_l);
  myError.setZ(mod_z);
  myError.setR(mod_r);

  double p;
  double curvature, cur_error;
  double p_cut; //, cur_cut;  
 
  double mod_d; 

  double theta = atan(mod_r/mod_z);
  cout << "# eta = " << -1 * log(tan(theta/2)) << endl;
  cout << "mod_d p_cut eff_1 inef_"<< efficient_pt << " pt(1%) pt(90%) geom_ineff" << endl;
  for (mod_d=start; mod_d<stop; mod_d+=step) {
    myError.setDistance(mod_d);
    
    p_cut = myError.stripsToP(window/2.);
    //cur_cut = 1/p_cut;

    //cout << mod_d << " "
    //     << p_cut << " ";
    printf("%.2f %.2f ", mod_d, p_cut);

    p = 1;
    curvature = 1 / p;
    cur_error = myError.computeError(p) * curvature;
    //cout << 100 * myError.probabilityInside(cur_cut, curvature, cur_error) * myError.geometricEfficiency() << " ";
    printf("%.2f ", 100 * myError.probabilityInside(1/p_cut, curvature, cur_error) * myError.geometricEfficiency());

    p = efficient_pt;
    curvature = 1 / p;
    cur_error = myError.computeError(p) * curvature;
    //cout << 100-100 * myError.probabilityInside(cur_cut, curvature, cur_error) * myError.geometricEfficiency() << " ";
    printf("%.2f ",100-100 * myError.probabilityInside(1/p_cut, curvature, cur_error) * myError.geometricEfficiency());

    double myValue;
    myValue = myError.find_probability(0.01, p_cut);
    //std::cout << 1/myValue << " ";
    printf("%.2f ", myValue);
    myValue = myError.find_probability(0.90, p_cut);
    //std::cout << 1/myValue;
    printf("%.2f ", myValue);


    myValue = 100*(1-myError.geometricEfficiency());
    printf("%.2f ", myValue);
    std::cout << std::endl;
  }

  return(0);
}

/*

int main(int ac, char* av[]) {
  try {
    // Declare the variables later used in the parameter scan
    std::cout << "Syntax: " << programName << endl
	      << " [Barrel|Endcap]" << std::endl
	      << " mod_d_start mod_d_stop mod_d_step" << std::endl
	      << " pitch strip_l mod_z mod_r window"<< std::endl;
    std::cout << "Pt units : GeV/c" << std::endl;
    std::cout << "Pitch    : micrometers" << std::endl;
    std::cout << "mod_*    : millimeters" << std::endl;
    std::cout << "window   : number of strips" << std::endl;
    std::cout << "And then the show: " << visible << std::endl;
    
    string ModuleType = "Barrel";
    double pitch = 90;     // micrometers
    double strip_l = 23.2; // mm
    double mod_z = 954;    // mm
    double mod_r = 348;    // mm
    double mod_d = 1.3;    // mm
    double window = 3;     // number of strips
    double ptCut;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic");
    generic.add_options()
      ("help", "produce help message")
      ;

    // Declare optional options
    po::options_description optional("Material plot");
    optional.add_options()
      ("trackerOnly,t", "takes only into account the outer tracker")
      ("maxRange,r", po::value<double>(&range), "sets the vertical range for the material plot")
      ("output,o", po::value<string>(&outputName), "output file name (.png will be added)")
      ;

    // Create the command line option set
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(optional);

    // Create the list to be shown
    po::options_description visible("Options");
    visible.add(generic).add(optional);

    // Actually read the options from command line
    po::variables_map vm;
    store(po::command_line_parser(ac, av).
	  options(cmdline_options).run(), vm);
    notify(vm);

    // List the options
    if (vm.count("help")) {
      printUsage(visible, av[0]);
      return 0;
    }
    
    bool usePixel = (vm.count("trackerOnly")==0);
    int myColor;

    // Look into the input-file list
    if (vm.count("input-file")) {
      string dirName;
      prepareCanvas();
      Palette::skipColors(100);
      Palette::prepare(fileV.size(), 210, 0.6,.6);
      for (unsigned int iMat=0 ; iMat<fileV.size(); ++iMat) {
	dirName = (mainConfiguration.getLayoutDirectory()
		   + "/" + fileV.at(iMat));
	cerr << dirName << endl; // debug
	myColor = Palette::color(iMat);
	plotMaterial(dirName, myColor, usePixel, fileV.at(iMat));
      }
      if (lastPixel) {
	myLegend->AddEntry(lastPixel, "pixel", "F");
      }
      //Palette::prepare(10, 210, 0.9, 0.4);
      printCanvas(outputName);
    } else {
      cout << "Error: at least a layout name should be given" << endl << endl;

      printUsage(visible, av[0]);
      return -1;
    }

  } catch(exception& e) {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
  
*/
