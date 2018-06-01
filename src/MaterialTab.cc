#include <fstream>
#include <sstream>
#include "MaterialTab.hh"
#include "global_constants.hh"
#include "MainConfigHandler.hh"
#include <MessageLogger.hh>
#include <stdexcept>

namespace material {


  ChemicalBase::ChemicalBase(const double density) : 
    density_(density) 
  { }




  ChemicalElement::ChemicalElement(const double density, const int atomicNumber, const double atomicMass) : 
    ChemicalBase(density),
    atomicNumber_(atomicNumber), 
    atomicMass_(atomicMass) 
  {
    radiationLength_ = computeRadiationLength(atomicNumber, atomicMass);
    interactionLength_ = computeInteractionLength(atomicMass);
  }


  /**
   * Calculate the estimated radiation length associated to an estimated atomic number Z and estimated atomic weight A.
   */
  const double ChemicalElement::computeRadiationLength(const int atomicNumber, const int atomicMass) {
    const double alpha = 1/137.035999139;
    const double a = alpha * atomicNumber;
    const double a_2 = pow(a, 2.);
    const double a_3 = pow(a, 3.);
    const double a_4 = pow(a, 4.);
    const double a_6 = pow(a, 6.);
    const double f_Z = a_2 * ( 1/(1 + a_2) + 0.20206 - 0.0369 * a_3 + 0.0083 * a_4 - 0.002 * a_6);
    double L_Z, L_prime_Z;
    switch (atomicNumber) {
    case 1:
      // H
      L_Z       = 5.31;
      L_prime_Z = 6.144;
      break;
    case 2:
      // He
      L_Z       = 4.79;
      L_prime_Z = 5.621;
      break;
    case 3:
      // Li
      L_Z       = 4.74;
      L_prime_Z = 5.805;
      break;
    case 4:
      // Be:
      L_Z       = 4.71;
      L_prime_Z = 5.924;
      break;
    default:
      // Z >4
      L_Z = log(184.15) - 1/3. * log(atomicNumber);
      L_prime_Z = log(1194) - 2/3. * log(atomicNumber);
    }
    
    // radiationLength
    const double radiationLength = 716.408 * atomicMass / (pow(atomicNumber, 2.) * (L_Z - f_Z) + atomicNumber * L_prime_Z);
    return radiationLength;
  }



  const double ChemicalElement::computeInteractionLength(const int atomicMass) {
    const double A_a = 31.645;
    const double A_b = 11.57238;

    const double interactionLength = pow(atomicMass, 1./3.) * A_a + A_b;
    return interactionLength;
  }



  ChemicalMixture::ChemicalMixture(const double density, const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements) :
    ChemicalBase(density),
    formula_(formula)
  {
    ratios_ = computeMassicComposition(formula, allChemicalElements);

    //ChemicalBaseMap alreadyDefinedMaterials(allChemicalElements.begin(), allChemicalElements.end());
    ChemicalBaseMap alreadyDefinedMaterials;
    alreadyDefinedMaterials.insert(allChemicalElements.begin(), allChemicalElements.end());

    const std::pair<double, double>& radiationAndInteractionLengths = computeRadiationAndInteractionLengths(ratios_, alreadyDefinedMaterials);
    radiationLength_ = radiationAndInteractionLengths.first;
    interactionLength_ = radiationAndInteractionLengths.second;
  }


  ChemicalMixture::ChemicalMixture(const double density, const MassicComposition& ratios, const ChemicalBaseMap& alreadyDefinedMaterials) :
    ChemicalBase(density),
    ratios_(ratios)
  {
    const std::pair<double, double>& radiationAndInteractionLengths = computeRadiationAndInteractionLengths(ratios_, alreadyDefinedMaterials);
    radiationLength_ = radiationAndInteractionLengths.first;
    interactionLength_ = radiationAndInteractionLengths.second;
  }



  const MassicComposition ChemicalMixture::computeMassicComposition(const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements) const {
    MassicComposition ratios;
    double totalMoleculeMass = 0.;

    for (const auto& elementIt : formula) {
      const std::string chemicalElementName = elementIt.first;
      const int chemicalElementNumber = elementIt.second;

      const auto found = allChemicalElements.find(chemicalElementName);
      if (found == allChemicalElements.end()) {
	std::cout << "Tried to create molecule made of unknow atom(s) name." << std::endl;
      }
      else {
	const ChemicalElement& element = found->second;
	const double elementAtomicMass = element.getAtomicMass();
	const double mass = chemicalElementNumber * elementAtomicMass;
	ratios.push_back(std::make_pair(chemicalElementName, mass));
	totalMoleculeMass += mass;

	std::cout << "ChemicalMixture::computeMassicComposition " 
		  << "chemicalElementName = " << chemicalElementName 
		  << "chemicalElementNumber = " << chemicalElementNumber
		  << "chemicalElementName = " << chemicalElementName
		  << "elementAtomicMass = " << elementAtomicMass
		  << "mass = " << mass
		  << std::endl;


      }
    }

    for (auto& ratioIt : ratios) {
      ratioIt.second /= totalMoleculeMass;
    }

    return ratios;
  }



  const std::pair<double, double> ChemicalMixture::computeRadiationAndInteractionLengths(const MassicComposition& ratios, const ChemicalBaseMap& alreadyDefinedMaterials) const {
    double invertedRadiationLength = 0.;
    double invertedInteractionLength = 0.;

    for (const auto& ratioIt : ratios) {
      const std::string chemicalBaseName = ratioIt.first;
      const double chemicalBaseMassicWeight = ratioIt.second;

      const auto found = alreadyDefinedMaterials.find(chemicalBaseName);
      if (found == alreadyDefinedMaterials.end()) {
	std::cout << "Tried to create mixture made of unknow chemical element / chemical compound / mixture:" << chemicalBaseName << std::endl;
      }
      else {
	const ChemicalBase& base = found->second;
	const double radiationLength = base.getRadiationLength();

	if (fabs(radiationLength) < insur::mat_negligible) std::cout << "Found a null radiation length for " << chemicalBaseName << std::endl;
	else {
	  invertedRadiationLength += chemicalBaseMassicWeight / radiationLength;
	}

	const double interactionLength = base.getInteractionLength();
	if (fabs(interactionLength) <  insur::mat_negligible) std::cout << "Found a null interaction length for " << chemicalBaseName << std::endl;
	else {
	  invertedInteractionLength += chemicalBaseMassicWeight / interactionLength;
	}

	std::cout << "ChemicalMixture::computeRadiationAndInteractionLengths " 
		  << "chemicalBaseName = " << chemicalBaseName 
		  << "base density = " << base.getDensity()
		  << "base radiationLength = " << radiationLength
		  << "base interactionLength = " << interactionLength
		  << std::endl;

      }
    }

    double radiationLength = 0.;
    if (fabs(invertedRadiationLength) < insur::mat_negligible) std::cout << "Mixture with infinite radiation length!!" << std::endl;
    else { radiationLength = 1. / invertedRadiationLength; }

    double interactionLength = 0.;
    if (fabs(invertedInteractionLength) < insur::mat_negligible) std::cout << "Mixture with infinite interaction length!!" << std::endl;
    else { interactionLength = 1. / invertedInteractionLength; }


    return std::make_pair(radiationLength, interactionLength);
  }



  MaterialsTable::MaterialsTable() {
    //std::string mattabFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_mattabfile);

    // CHEMICAL ELEMENTS
    std::ifstream chemicalElementsFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_chemicalElementsFile);
 
    ChemicalElementMap allChemicalElements; 

    if (chemicalElementsFile.good()) {
      std::string lineString;

      while (std::getline(chemicalElementsFile, lineString)) {
	std::istringstream myLine;
        myLine.str(lineString);
	
	std::string elementName;
        myLine >> elementName;

        //check if is a comment
        if (elementName != "" && elementName[0] != '#') {
	  // TO DO: check whether input is of expected type and 3 values
	  double elementDensity;
	  int atomicNumber;
	  double atomicMass;
          myLine >> elementDensity >> atomicNumber >> atomicMass;
          elementDensity /= 1000.;   // convert g/cm3 in g/mm3

	  std::cout << "MaterialsTable::MaterialsTable() create eleemtray table " 
		    << " elementName = " <<  elementName
		    << " elementDensity = " << elementDensity 
		    << " atomicNumber = " << atomicNumber
		    << " atomicMass  =" << atomicMass 
		    << std::endl;


	  ChemicalElement element = ChemicalElement(elementDensity, atomicNumber, atomicMass);
          allChemicalElements.insert(std::make_pair(elementName, element));
        }
        myLine.clear();
      }
    } else {
      logERROR("Could not open chemical elements file.");
    }


    for (const auto& elemIt : allChemicalElements) {
      std::cout << "MaterialsTable::MaterialsTable finsihed computing all pure elem. load Elementary mat = " << elemIt.first;
      const ChemicalElement& elem = elemIt.second;
      std::cout << " elem.getDensity() = " << elem.getDensity()
		<< " elem.getRadiationLength() = " << elem.getRadiationLength()
		<< " elem.getInteractionLength() = " << elem.getInteractionLength() 
		<< " elem.getAtomicNumber() = " << elem.getAtomicNumber()
		<< " elem.getAtomicMass() = " << elem.getAtomicMass()
		<< " elem.isChemicalElement() = " << elem.isChemicalElement()
		<< std::endl;
    }



    // CHEMICAL COMPOUNDS
    std::ifstream chemicalCompoundsFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_chemicalCompoundsFile); 
   
    ChemicalMixtureMap allChemicalMixtures;

    if (chemicalCompoundsFile.good()) {
      std::string lineString;
      while (std::getline(chemicalCompoundsFile, lineString)) {
	std::istringstream myLine;
        myLine.str(lineString);

	std::string compoundName;
        myLine >> compoundName;

        //check if is a comment
        if (compoundName != "" && compoundName[0] != '#') {
	  double compoundDensity;
          myLine >> compoundDensity;
          compoundDensity /= 1000.;   // convert g/cm3 in g/mm3

	  ChemicalFormula compoundFormula;
	  std::string element;

	  std::cout << "MaterialsTable::MaterialsTable() create compound table " 
		    << " compoundName = " <<  compoundName
		    << " compoundDensity = " << compoundDensity;


	  while (myLine >> element) {

	    const auto delimiterPosition = element.find(insur::default_composition_delimiter);
	    if (delimiterPosition != std::string::npos) {
	      const std::string elementName = element.substr(0, delimiterPosition);
	      const std::string elementNumberString = element.substr(delimiterPosition + insur::default_composition_delimiter.length());
	      const int elementNumber = atoi(elementNumberString.c_str());
	      std::cout << "elementName = " << elementName << " elementNumber = " << elementNumber;
	      compoundFormula.push_back(std::make_pair(elementName, elementNumber));
	    }
	    else { std::cout << "Chemical compound: could not find the : delimiter." << std::endl; }

	    element.clear();
	  }

	  std::cout << "." << std::endl;

	  ChemicalMixture coumpound = ChemicalMixture(compoundDensity, compoundFormula, allChemicalElements);
	  allChemicalMixtures.insert(std::make_pair(compoundName, coumpound));	  
        }
        myLine.clear();
      }
    } 
    else {
      logERROR("Could not open chemical compounds file.");
    }


    for (const auto& mixIt : allChemicalMixtures) {
      std::cout << "MaterialsTable::MaterialsTable finsihed computing all chemical composites. load Composite = " << mixIt.first;
      const ChemicalMixture& mix = mixIt.second;
      std::cout << " mix.getDensity() = " << mix.getDensity()
		<< " mix.getRadiationLength() = " << mix.getRadiationLength()
		<< " mix.getInteractionLength() = " << mix.getInteractionLength() 
		<< " mix.hasChemicalFormula() = " << mix.hasChemicalFormula()
		<< " mix.isChemicalElement() = " << mix.isChemicalElement();
      const ChemicalFormula& formula = mix.getChemicalFormula();
      for (const auto& formulaIt : formula) {
	std::cout << " formulaIt.first = " << formulaIt.first
		  << " formulaIt.second = " << formulaIt.second;
      }
      const MassicComposition& ratio = mix.getMassicComposition();
      for (const auto& ratioIt : ratio) {
	std::cout << " ratioIt.first = " << ratioIt.first
		  << " ratioIt.second = " << ratioIt.second;
      }
      std::cout << "." << std::endl;
    }



    // CHEMICAL MIXTURES
    std::ifstream chemicalMixturesFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_chemicalMixturesFile); 

    ChemicalBaseMap alreadyDefinedMaterials;
    alreadyDefinedMaterials.insert(allChemicalElements.begin(), allChemicalElements.end());
    alreadyDefinedMaterials.insert(allChemicalMixtures.begin(), allChemicalMixtures.end());

    if (chemicalMixturesFile.good()) {
      std::string lineString;
      while (std::getline(chemicalMixturesFile, lineString)) {
	std::istringstream myLine;
        myLine.str(lineString);

	std::string mixtureName;
        myLine >> mixtureName;

        //check if is a comment
        if (mixtureName != "" && mixtureName[0] != '#') {
	  double mixtureDensity;
          myLine >> mixtureDensity;
          mixtureDensity /= 1000.;   // convert g/cm3 in g/mm3


	  std::cout << "MaterialsTable::MaterialsTable() create mixture table " 
		    << " mixtureName = " <<  mixtureName
		    << " mixtureDensity = " << mixtureDensity;


	  MassicComposition mixtureComposition;
	  std::string constituant;
	  while (myLine >> constituant) { // TO DO: cross-check this

	    const auto delimiterPosition = constituant.find(insur::default_composition_delimiter);
	    if (delimiterPosition != std::string::npos) {
	      const std::string constituantName = constituant.substr(0, delimiterPosition);
	      const std::string constituantMassicWeightString = constituant.substr(delimiterPosition + insur::default_composition_delimiter.length());
	      const double constituantMassicWeight = std::stod(constituantMassicWeightString);
	      std::cout << "constituantName = " << constituantName << " constituantMassicWeight = " << constituantMassicWeight;
	      mixtureComposition.push_back(std::make_pair(constituantName, constituantMassicWeight));
	    }
	    else { std::cout << "Chemical mixture: could not find the : delimiter." << std::endl; }

	    constituant.clear();
	  }

	  std::cout << "." << std::endl;

	  ChemicalMixture mixture = ChemicalMixture(mixtureDensity, mixtureComposition, alreadyDefinedMaterials);
	  allChemicalMixtures.insert(std::make_pair(mixtureName, mixture));
	  alreadyDefinedMaterials.insert(std::make_pair(mixtureName, mixture));	  
        }
        myLine.clear();
      }
    } 
    else {
      logERROR("Could not open chemical mixtures file.");
    }


    for (const auto& mixIt : allChemicalMixtures) {
      std::cout << "MaterialsTable::MaterialsTable finsihed computing all mixtures. load mixture = " << mixIt.first;
      const ChemicalMixture& mix = mixIt.second;
      std::cout << " mix.getDensity() = " << mix.getDensity()
		<< " mix.getRadiationLength() = " << mix.getRadiationLength()
		<< " mix.getInteractionLength() = " << mix.getInteractionLength() 
		<< " mix.hasChemicalFormula() = " << mix.hasChemicalFormula()
		<< " mix.isChemicalElement() = " << mix.isChemicalElement();
      const ChemicalFormula& formula = mix.getChemicalFormula();
      for (const auto& formulaIt : formula) {
	std::cout << " formulaIt.first = " << formulaIt.first
		  << " formulaIt.second = " << formulaIt.second;
      }
      const MassicComposition& ratio = mix.getMassicComposition();
      for (const auto& ratioIt : ratio) {
	std::cout << " ratioIt.first = " << ratioIt.first
		  << " ratioIt.second = " << ratioIt.second;
      }
      std::cout << "." << std::endl;
    }




    this->first = allChemicalElements;
    this->second = allChemicalMixtures;
  }

  const MaterialsTable& MaterialsTable::instance() {
    static MaterialsTable instance_;
    return instance_;
  }

  double MaterialsTable::density(const std::string materialName) const {
    ChemicalElementMap allChemicalElements = this->first;
    ChemicalMixtureMap allChemicalMixtures = this->second;

    double density = 0.;

    const auto& found = allChemicalElements.find(materialName);
    if (found != allChemicalElements.end()) { density = found->second.getDensity(); }
    else { std::cout << "MaterialsTable::density: material " << found->first << " could not be found in MaterialsTable." << std::endl; }
    return density;
  }

  double MaterialsTable::radiationLength(const std::string materialName) const {
    ChemicalElementMap allChemicalElements = this->first;
    ChemicalMixtureMap allChemicalMixtures = this->second;

    double radiationLength = 0.;

    const auto& found = allChemicalElements.find(materialName);
    if (found != allChemicalElements.end()) { radiationLength = found->second.getRadiationLength(); }
    else { std::cout << "MaterialsTable::radiationLength: material " << found->first << " could not be found in MaterialsTable." << std::endl; }
    return radiationLength;
  }

  double MaterialsTable::interactionLength(const std::string materialName) const {
    ChemicalElementMap allChemicalElements = this->first;
    ChemicalMixtureMap allChemicalMixtures = this->second;

    double interactionLength = 0.;

    const auto& found = allChemicalElements.find(materialName);
    if (found != allChemicalElements.end()) { interactionLength = found->second.getInteractionLength(); }
    else { std::cout << "MaterialsTable::interactionLength: material " << found->first << " could not be found in MaterialsTable." << std::endl; }
    return interactionLength;
  }










  const std::string MaterialTab::msg_no_mat_file = "Material tab file does not exist.";
  const std::string MaterialTab::msg_no_mat_file_entry1 = "Material '";
  const std::string MaterialTab::msg_no_mat_file_entry2 = "' not found in Material tab file.";


  MaterialTab::MaterialTab() {
    std::string mattabFile(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_mattabfile);
    std::ifstream mattabStream(mainConfigHandler::instance().getMattabDirectory() + "/" + insur::default_mattabfile);
    std::string line;
    std::string material;
    std::istringstream lineStream;
    double density, radiationLength, interactionLength;

    if (mattabStream.good()) {
      while (!mattabStream.eof()) {
        std::getline(mattabStream, line);
        lineStream.str(line);
        lineStream >> material;

        //check if is a comment
        if (material[0] != '#') {
          lineStream >> density >> radiationLength >> interactionLength;
          density /= 1000; // convert g/cm3 in g/mm3
          insert(make_pair(material, make_tuple(density, radiationLength, interactionLength)));
        }

        lineStream.clear();
      }
    } else {
      logERROR(msg_no_mat_file);
    }
  }

  const MaterialTab& MaterialTab::instance() {
    static MaterialTab instance_;
    return instance_;
  }

  double MaterialTab::density(std::string material) const {
    double val = 0;
    try {
      val = std::get<0>(at(material));
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
      val = -1;
    }
    return val;
  }

  double MaterialTab::radiationLength(std::string material) const {
    double val = 0;
    try {
      val = std::get<1>(at(material));
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
      val = -1;
    }
    return val;
  }

  double MaterialTab::interactionLength(std::string material) const {
    double val = 0;
    try {
      val = std::get<2>(at(material));
    } catch (const std::out_of_range& ex) {
      logERROR(msg_no_mat_file_entry1 + material + msg_no_mat_file_entry2);
      val = -1;
    }
    return val;
  }
} /* namespace material */
