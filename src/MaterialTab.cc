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




  ChemicalElement::ChemicalElement(const double density, const int atomicNumber, const double atomicWeight) : 
    ChemicalBase(density),
    atomicNumber_(atomicNumber), 
    atomicWeight_(atomicWeight) 
  {
    radiationLength_ = computeRadiationLength(atomicNumber, atomicWeight);
    interactionLength_ = computeInteractionLength(atomicWeight);
  }


  /**
   * Calculate the estimated radiation length associated to an estimated atomic number Z and estimated atomic weight A.
   */
  const double ChemicalElement::computeRadiationLength(const int atomicNumber, const double atomicWeight) {
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
    const double radiationLength = 716.408 * atomicWeight / (pow(atomicNumber, 2.) * (L_Z - f_Z) + atomicNumber * L_prime_Z);
    return radiationLength;
  }



  const double ChemicalElement::computeInteractionLength(const double atomicWeight) {
    const double A_a = 31.645;
    const double A_b = 11.57238;

    const double interactionLength = pow(atomicWeight, 1./3.) * A_a + A_b;
    //const double invertedInteractionLength = 1. / 35. * pow(atomicWeight, -1./3.);
    //const double interactionLength = 1. / invertedInteractionLength;

    return interactionLength;
  }



  ChemicalMixture::ChemicalMixture(const double density, const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements) :
    ChemicalBase(density),
    formula_(formula)
  {
    fractions_ = computeMassComposition(formula, allChemicalElements);

    //ChemicalBaseMap alreadyDefinedMaterials(allChemicalElements.begin(), allChemicalElements.end());
    ChemicalBaseMap alreadyDefinedMaterials;
    alreadyDefinedMaterials.insert(allChemicalElements.begin(), allChemicalElements.end());

    const std::pair<double, double>& radiationAndInteractionLengths = computeRadiationAndInteractionLengths(fractions_, alreadyDefinedMaterials);
    radiationLength_ = radiationAndInteractionLengths.first;
    interactionLength_ = radiationAndInteractionLengths.second;
  }


  ChemicalMixture::ChemicalMixture(const double density, const MassComposition& fractions, const ChemicalBaseMap& alreadyDefinedMaterials) :
    ChemicalBase(density),
    fractions_(fractions)
  {
    checkMassFractionsSum(fractions);

    const std::pair<double, double>& radiationAndInteractionLengths = computeRadiationAndInteractionLengths(fractions_, alreadyDefinedMaterials);
    radiationLength_ = radiationAndInteractionLengths.first;
    interactionLength_ = radiationAndInteractionLengths.second;
  }


  void ChemicalMixture::checkMassFractionsSum(const MassComposition& fractions) const {
    double fractionSum = 0.;
    for (const auto& fractionIt : fractions) {
      const double chemicalBaseMassFraction = fractionIt.second;
      fractionSum += chemicalBaseMassFraction;
    }

    if (fabs(fractionSum - 1.) > insur::mat_negligible) { 
      std::cout << "Error defining Chemical mixture: sum of massic weights is not equal to 1." << std::endl;
      std::cout << "Inadequate mixture composition is:";
      for (const auto& fractionIt : fractions) {
	const std::string chemicalBaseName = fractionIt.first;
	const double chemicalBaseMassFraction = fractionIt.second;
	std::cout << " " << chemicalBaseName << ":" << chemicalBaseMassFraction;
      }
      std::cout << "." << std::endl;
    }

  }



  const MassComposition ChemicalMixture::computeMassComposition(const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements) const {
    MassComposition fractions;
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
	const double elementAtomicWeight = element.getAtomicWeight();
	const double mass = chemicalElementNumber * elementAtomicWeight;
	fractions.push_back(std::make_pair(chemicalElementName, mass));
	totalMoleculeMass += mass;

	/*std::cout << "ChemicalMixture::computeMassComposition " 
		  << "chemicalElementName = " << chemicalElementName 
		  << "chemicalElementNumber = " << chemicalElementNumber
		  << "chemicalElementName = " << chemicalElementName
		  << "elementAtomicWeight = " << elementAtomicWeight
		  << "mass = " << mass
		  << std::endl;*/


      }
    }

    for (auto& fractionIt : fractions) {
      fractionIt.second /= totalMoleculeMass;
    }
    checkMassFractionsSum(fractions);

    return fractions;
  }



  const std::pair<double, double> ChemicalMixture::computeRadiationAndInteractionLengths(const MassComposition& fractions, const ChemicalBaseMap& alreadyDefinedMaterials) const {
    double invertedRadiationLength = 0.;
    double invertedInteractionLength = 0.;

    for (const auto& fractionIt : fractions) {
      const std::string chemicalBaseName = fractionIt.first;
      const double chemicalBaseMassFraction = fractionIt.second;

      const auto found = alreadyDefinedMaterials.find(chemicalBaseName);
      if (found == alreadyDefinedMaterials.end()) {
	std::cout << "Tried to create mixture made of unknow chemical element / chemical compound / mixture:" << chemicalBaseName << std::endl;
      }
      else {
	const ChemicalBase& base = found->second;
	const double radiationLength = base.getRadiationLength();

	if (fabs(radiationLength) < insur::mat_negligible) std::cout << "Found a null radiation length for " << chemicalBaseName << std::endl;
	else {
	  invertedRadiationLength += chemicalBaseMassFraction / radiationLength;
	}

	const double interactionLength = base.getInteractionLength();
	if (fabs(interactionLength) <  insur::mat_negligible) std::cout << "Found a null interaction length for " << chemicalBaseName << std::endl;
	else {
	  invertedInteractionLength += chemicalBaseMassFraction / interactionLength;
	}

	/*std::cout << "ChemicalMixture::computeRadiationAndInteractionLengths " 
		  << "chemicalBaseName = " << chemicalBaseName 
		  << "base density = " << base.getDensity()
		  << "base radiationLength = " << radiationLength
		  << "base interactionLength = " << interactionLength
		  << std::endl;*/

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
	  double atomicWeight;
          myLine >> elementDensity >> atomicNumber >> atomicWeight;
          elementDensity /= 1000.;   // convert g/cm3 in g/mm3

	  /*std::cout << "MaterialsTable::MaterialsTable() create eleemtray table " 
		    << " elementName = " <<  elementName
		    << " elementDensity = " << elementDensity 
		    << " atomicNumber = " << atomicNumber
		    << " atomicWeight  =" << atomicWeight 
		    << std::endl;*/


	  ChemicalElement element = ChemicalElement(elementDensity, atomicNumber, atomicWeight);
          allChemicalElements.insert(std::make_pair(elementName, element));
        }
        myLine.clear();
      }
    } else {
      logERROR("Could not open chemical elements file.");
    }

   
    /* 
    for (const auto& elemIt : allChemicalElements) {
      const MaterialTab& oldTable = MaterialTab::instance();
      const double oldRad = oldTable.radiationLength(elemIt.first);
      const double oldInt = oldTable.interactionLength(elemIt.first);
      const ChemicalElement& elem = elemIt.second;
      const double radRatio = (elem.getRadiationLength() - oldRad) / oldRad * 100.;
      const double intRatio = (elem.getInteractionLength() - oldInt) / oldInt * 100.;
      std::cout << elemIt.first << " radRatio = " << radRatio << "intRatio = " << intRatio << std::endl;

      std::cout << "MaterialsTable::MaterialsTable finsihed computing all pure elem. load Elementary mat = " << elemIt.first;
      const ChemicalElement& elem = elemIt.second;
      std::cout << " elem.getDensity() = " << elem.getDensity()
		<< " elem.getRadiationLength() = " << elem.getRadiationLength()
		<< " elem.getInteractionLength() = " << elem.getInteractionLength() 
		<< " elem.getAtomicNumber() = " << elem.getAtomicNumber()
		<< " elem.getAtomicWeight() = " << elem.getAtomicWeight()
		<< " elem.isChemicalElement() = " << elem.isChemicalElement()
		<< std::endl;
    }*/
    
    



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

	  /*std::cout << "MaterialsTable::MaterialsTable() create compound table " 
		    << " compoundName = " <<  compoundName
		    << " compoundDensity = " << compoundDensity;*/


	  while (myLine >> element) {

	    const auto delimiterPosition = element.find(insur::default_composition_delimiter);
	    if (delimiterPosition != std::string::npos) {
	      const std::string elementName = element.substr(0, delimiterPosition);
	      const std::string elementNumberString = element.substr(delimiterPosition + insur::default_composition_delimiter.length());
	      const int elementNumber = atoi(elementNumberString.c_str());
	      //std::cout << "elementName = " << elementName << " elementNumber = " << elementNumber;
	      compoundFormula.push_back(std::make_pair(elementName, elementNumber));
	    }
	    else { std::cout << "Chemical compound: could not find the : delimiter." << std::endl; }

	    element.clear();
	  }

	  //std::cout << "." << std::endl;

	  ChemicalMixture coumpound = ChemicalMixture(compoundDensity, compoundFormula, allChemicalElements);
	  allChemicalMixtures.insert(std::make_pair(compoundName, coumpound));	  
        }
        myLine.clear();
      }
    } 
    else {
      logERROR("Could not open chemical compounds file.");
    }


    /*
    for (const auto& mixIt : allChemicalMixtures) {
      const MaterialTab& oldTable = MaterialTab::instance();
      const double oldRad = oldTable.radiationLength(mixIt.first);
      const double oldInt = oldTable.interactionLength(mixIt.first);
      const ChemicalMixture& mix = mixIt.second;
      const double radRatio = (mix.getRadiationLength() - oldRad) / oldRad * 100.;
      const double intRatio = (mix.getInteractionLength() - oldInt) / oldInt * 100.;
      std::cout << mixIt.first << " radRatio = " << radRatio << "intRatio = " << intRatio << std::endl;

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
      const MassComposition& fraction = mix.getMassComposition();
      for (const auto& fractionIt : fraction) {
	std::cout << " fractionIt.first = " << fractionIt.first
		  << " fractionIt.second = " << fractionIt.second;
      }
      std::cout << "." << std::endl;
      }*/
    
    


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


	  /*std::cout << "MaterialsTable::MaterialsTable() create mixture table " 
		    << " mixtureName = " <<  mixtureName
		    << " mixtureDensity = " << mixtureDensity;*/


	  MassComposition mixtureComposition;
	  std::string constituant;
	  while (myLine >> constituant) { // TO DO: !!!!! should check for (myLine >> constituant) error

	    const auto delimiterPosition = constituant.find(insur::default_composition_delimiter);
	    if (delimiterPosition != std::string::npos) {
	      const std::string constituantName = constituant.substr(0, delimiterPosition);
	      const std::string constituantMassFractionString = constituant.substr(delimiterPosition + insur::default_composition_delimiter.length());
	      const double constituantMassFraction = std::stod(constituantMassFractionString);
	      //std::cout << "constituantName = " << constituantName << " constituantMassFraction = " << constituantMassFraction;
	      mixtureComposition.push_back(std::make_pair(constituantName, constituantMassFraction));
	    }
	    else { 
	      std::cout << "Chemical mixture: could not find the : delimiter. Please set name:value, with no space between name and value." << std::endl; 
	    }

	    constituant.clear();
	  }

	  //std::cout << "." << std::endl;

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
      if (!mixIt.second.hasChemicalFormula()) {
	const MaterialTab& oldTable = MaterialTab::instance();
	const double oldRad = oldTable.radiationLength(mixIt.first);
	const double oldInt = oldTable.interactionLength(mixIt.first);

	std::string closestElementName = "none";
	double closestElementRad = -1.;
	double closestElementInt = -1.;
	double distance = 10000.;
	for (const auto& elemIt : allChemicalElements) {
	  const std::string elementName = elemIt.first;
	  const ChemicalElement& elem = elemIt.second;
	  const double elementRad = elem.getRadiationLength();
	  const double elementInt = elem.getInteractionLength();

	  double elementRadDistance = fabs(elementRad - oldRad) / oldRad;
	  double elementIntDistance = fabs(elementInt - oldInt) / oldInt;
	  if ( (elementRadDistance + elementIntDistance) < distance) {
	    closestElementName = elementName;
	    closestElementRad = elementRad;
	    closestElementInt = elementInt;
	    distance = elementRadDistance + elementIntDistance;
	  }
	}

	const double radRatio = (closestElementRad - oldRad) / oldRad * 100.;
	const double intRatio = (closestElementInt - oldInt) / oldInt * 100.;
	//std::cout << mixIt.first << " radRatio = " << radRatio << "intRatio = " << intRatio << std::endl;
	std::cout << mixIt.first << " " << closestElementName << std::endl;

      
	/*
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
	  const MassComposition& fraction = mix.getMassComposition();
	  for (const auto& fractionIt : fraction) {
	  std::cout << " fractionIt.first = " << fractionIt.first
	  << " fractionIt.second = " << fractionIt.second;
	  }
	  std::cout << "." << std::endl;*/
      }
    }
    
    



    this->first = allChemicalElements;
    this->second = allChemicalMixtures;
  }


  const MaterialsTable& MaterialsTable::instance() {
    static MaterialsTable instance_;
    return instance_;
  }

  double MaterialsTable::getDensity(const std::string materialName) const {
    const ChemicalElementMap& allChemicalElements = this->first;
    const ChemicalMixtureMap& allChemicalMixtures = this->second;

    double density = 0.;

    const auto foundElem = allChemicalElements.find(materialName);
    if (foundElem != allChemicalElements.end()) { density = foundElem->second.getDensity(); }
    else { 
      const auto foundMix = allChemicalMixtures.find(materialName);
      if (foundMix != allChemicalMixtures.end()) { density = foundMix->second.getDensity(); }
      else { std::cout << "MaterialsTable::density: material " << materialName << " could not be found in MaterialsTable." << std::endl; }
    }
    return density;
  }

  double MaterialsTable::getRadiationLength(const std::string materialName) const {
    const ChemicalElementMap& allChemicalElements = this->first;
    const ChemicalMixtureMap& allChemicalMixtures = this->second;

    double radiationLength = 0.;

    const auto foundElem = allChemicalElements.find(materialName);
    if (foundElem != allChemicalElements.end()) { radiationLength = foundElem->second.getRadiationLength(); }
    else { 
      const auto foundMix = allChemicalMixtures.find(materialName);
      if (foundMix != allChemicalMixtures.end()) { radiationLength = foundMix->second.getRadiationLength(); }
      else { std::cout << "MaterialsTable::radiationLength: material " << materialName << " could not be found in MaterialsTable." << std::endl; }
    }
    return radiationLength;
  }

  double MaterialsTable::getInteractionLength(const std::string materialName) const {
    const ChemicalElementMap& allChemicalElements = this->first;
    const ChemicalMixtureMap& allChemicalMixtures = this->second;

    double interactionLength = 0.;

    const auto foundElem = allChemicalElements.find(materialName);
    if (foundElem != allChemicalElements.end()) { interactionLength = foundElem->second.getInteractionLength(); }
    else { 
      const auto foundMix = allChemicalMixtures.find(materialName);
      if (foundMix != allChemicalMixtures.end()) { interactionLength = foundMix->second.getInteractionLength(); }
      else { std::cout << "MaterialsTable::interactionLength: material " << materialName << " could not be found in MaterialsTable." << std::endl; }
    }
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
