#ifndef MATERIALTAB_HH
#define MATERIALTAB_HH

//#include <map>
#include <tuple>

#include <global_constants.hh>
#include <global_funcs.hh>

namespace material {

  //class TableChemicalElement : public TableMaterial {


  class ChemicalBase {
  public:
    ChemicalBase(const int materialsTreeHierarchyLevel, const double density);
    const double getDensity() const { return density_; }
    const double getRadiationLength() const { return radiationLength_; }
    const double getInteractionLength() const { return interactionLength_; }
    const int getMaterialsTreeHierarchyLevel() const { return materialsTreeHierarchyLevel_; }

    //virtual const bool isChemicalElement() const = 0;

  protected:
    int materialsTreeHierarchyLevel_;
    double density_;
    double radiationLength_;
    double interactionLength_;    
  };
  

  class ChemicalElement : public ChemicalBase {
  public:
    ChemicalElement(const double density, const int atomicNumber, const double atomicWeight);
    const int getAtomicNumber() const { return atomicNumber_; } 
    const double getAtomicWeight() const { return atomicWeight_; }   // standard atomic weight (u)
    //const bool isChemicalElement() const override { return true; }
    const bool isChemicalElement() const { return true; }

  private:
    // TO DO: should probably set more accurate RL and IL values directly in the cfg file?
    const double computeRadiationLength(const int atomicNumber, const double atomicWeight); 
    const double computeInteractionLength(const double atomicWeight);

    int atomicNumber_;       // atome's Z
    double atomicWeight_;    // atome's A    
  };


  typedef std::map<std::string, ChemicalBase> ChemicalBaseMap; // to do: define ChemicalBase directly inside MaterialsTable to avoid this?
  typedef std::map<std::string, ChemicalElement> ChemicalElementMap;

  typedef std::vector< std::pair<std::string, int> > ChemicalFormula;
  typedef std::vector< std::pair<std::string, double> > MassComposition;

  class ChemicalMixture : public ChemicalBase {
  public:
    ChemicalMixture(const double density, const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements);
    ChemicalMixture(const double density, const MassComposition& fractions, const ChemicalBaseMap& alreadyDefinedMaterials);

    //const bool isChemicalElement() const override { return false; }
    const bool isChemicalElement() const { return false; }
    const bool hasChemicalFormula() const { return (formula_.size() != 0); }
    const ChemicalFormula getChemicalFormula() const { return formula_; }
    const MassComposition getMassComposition() const { return fractions_; }

  private:
    // TO DO: should probably place more accurate RL and IL values directly in the file?
    const MassComposition computeMassComposition(const ChemicalFormula& formula, const ChemicalElementMap& allChemicalElements) const;
    const MassComposition computeMassComposition(const MassComposition& fractionsFromCfg) const;
    void normalizeMassComposition(MassComposition& fractions) const;
    void checkMassCompositionSum(const MassComposition& fractions) const;

    const int computeMaterialsTreeHierarchyLevel(const MassComposition& fractions, const ChemicalBaseMap& alreadyDefinedMaterials) const;

    const std::pair<double, double> computeRadiationAndInteractionLengths(const MassComposition& fractions, const ChemicalBaseMap& alreadyDefinedMaterials) const;

    ChemicalFormula formula_;
    MassComposition fractions_;
  };

  typedef std::map<std::string, ChemicalMixture> ChemicalMixtureMap;



  typedef std::pair<ChemicalElementMap, ChemicalMixtureMap > MaterialsTableType;
  class MaterialsTable : public MaterialsTableType {
  private:
    MaterialsTable();
 
  public:
    static const MaterialsTable& instance();

    int getMaxMaterialsTreeHierarchyLevel() const;

    double getDensity(const std::string materialName) const;  // TO DO: should return const double!!!
    double getRadiationLength(const std::string materialName) const;
    double getInteractionLength(const std::string materialName) const;

    const ChemicalElementMap getAllChemicalElements() const { return this->first; }
    const ChemicalMixtureMap getAllChemicalMixtures() const { return this->second; }
  };












  typedef std::map<std::string, std::tuple<double, double, double> > MaterialTabType;
  class MaterialTab : public MaterialTabType {
  private:
    MaterialTab();
    static const std::string msg_no_mat_file;
    static const std::string msg_no_mat_file_entry1;
    static const std::string msg_no_mat_file_entry2;

  public:
    static const MaterialTab& instance();

    double density(std::string material) const;
    double radiationLength(std::string material) const;
    double interactionLength(std::string material) const;
  };
} /* namespace material */

#endif /* MATERIALTAB_HH */
