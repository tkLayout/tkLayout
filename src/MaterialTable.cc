/**
 * @file MaterialTable.cc
 * @brief This is the implementation of the internal data structure for material properties.
 */

#include <MaterialTable.hh>
namespace insur {
    
  void MaterialTable2::parseMaterial(std::string str) {
    str = trim(str);
    BaseMaterial* mat = 0;
    if (str.find("Elementary") == 0) { 
      mat = new ElementaryMaterial();
      mat->fromCfg(str.substr(10));
    } else if (str.find("Composite") == 0) {
      mat = new CompositeMaterial(this);
      mat->fromCfg(str.substr(9));
    }
    materials_[mat->getTag()] = mat;
  }


  void ElementaryMaterial::fromCfg(std::string str) {
    //Example string: "Carbon; C; C_twp" 0.4  1  2
    
    str = ltrim(str);

    int tagEnd; 
    if (str[0] == '\"') {  // multiple tags or a tag with spaces inside
      tagEnd = str.find_first_of("\"", 1);
      std::vector<std::string> tags = split(str.substr(1, tagEnd-1), ";");
      std::transform(tags.begin(), tags.end(), tags.begin(), &trim);
      tag_ = tags.at(0);
      aliases_.insert(tags.begin()+1, tags.end());
    } else { // only 1 tag with no spaces inside
      tagEnd = str.find_first_of(" \t", 0);
      tag_ = str.substr(0, tagEnd);
    }
    str.erase(0, tagEnd+1);

    std::stringstream ss(str);

    ss >> density_;
    ss >> z_;
    ss >> a_;
  }


  void CompositeMaterial::fromCfg(std::string str) {
    //Example string: Composite "Kapton" 0.6  "Pyroid":0.  "C":0.3
   
    str = ltrim(str);

    int tagEnd; 
    if (str[0] == '\"') {  // multiple tags or a tag with spaces inside
      tagEnd = str.find_first_of("\"", 1);
      std::vector<std::string> tags = split(str.substr(1, tagEnd-1), ";");
      std::transform(tags.begin(), tags.end(), tags.begin(), &trim);
      tag_ = tags.at(0);
      aliases_.insert(tags.begin()+1, tags.end());
    } else { // only 1 tag with no spaces inside
      tagEnd = str.find_first_of(" \t", 0);
      tag_ = str.substr(0, tagEnd);
    }
    str.erase(0, tagEnd+1);
    
    std::stringstream ss(str);
    ss >> density_;

    std::string tag;
    double fraction;
    while (!str.empty()) {
      str = ltrim(str);
      if (str[0] == '\"') {
        tagEnd = str.find_first_of("\"", 1);
        tag = trim(str.substr(1, tagEnd-1));
        tagEnd = str.find_first_of(":", tagEnd+1);
      } else {
        tagEnd = str.find_first_of(":", 0);
        tag = trim(str.substr(0, tagEnd));
      }
      BaseMaterial* mat = mTable_->getMaterial(tag);
      
      str.erase(0, tagEnd+1);

      std::stringstream ss2(str);
      ss2 >> fraction;

      materialFractions_.push_back((MaterialFraction){mat, fraction});
    }

  }


  double ElementaryMaterial::getRadiationLength() {
    if (rlength_ <= 0) {
      rlength_ = 716.4*a_/(z_*(z_+1)*log(287/sqrt(z_)));
    }
    return rlength_;
  }

  double ElementaryMaterial::getInteractionLength() {
    if (ilength_ <= 0) {
      ilength_ = 35*pow(a_,1/3)/density_; 
    }
    return ilength_;
  }

  double CompositeMaterial::getRadiationLength() {
    if (rlength_ <= 0) {
      for (std::vector<MaterialFraction>::const_iterator it = materialFractions_.begin(); it != materialFractions_.end(); ++it) {
        rlength_ += it->fraction/it->material->getRadiationLength();
      }
      rlength_ = 1/rlength_;
    }
    return rlength_;
  }

  double CompositeMaterial::getInteractionLength() {
    if (ilength_ <= 0) {
      for (std::vector<MaterialFraction>::const_iterator it = materialFractions_.begin(); it != materialFractions_.end(); ++it) {
        ilength_ += it->fraction/it->material->getInteractionLength();
      }
      ilength_ = 1/ilength_;
    }
    return ilength_;
  }



  /**
   * Adds a material already bundled into a <i>MaterialRow</i> struct to the container class.
   * @param mat The struct that will be copied and appended to the internal vector container
   */
  void MaterialTable::addMaterial(MaterialRow mat) {
    materials.push_back(mat);
  }

  /**
   * Add a material to the container class by bundling the individual parts into a struct.
   * @param tag The name of the material
   * @param density The material density
   * @param rlength The radiation length of the material
   * @param ilength The interaction length of the material
   */
  void MaterialTable::addMaterial(std::string tag, double density, double rlength, double ilength) {
    MaterialRow row;
    row.tag = tag;
    row.density = density;
    row.rlength = rlength;
    row.ilength = ilength;
    materials.push_back(row);
  }

  /**
   * Find a material by its name. If no material is found, the function throws an exception.
   * @param tag The name of the requested material
   * @return A reference to the requested <i>MaterialRow</i> struct
   */
  MaterialRow& MaterialTable::getMaterial(std::string tag) { // throws std::runtime_error
    int index = findIndex(tag);
    if (index < 0) throw std::runtime_error("MaterialTable::getMaterialByTag(): " + errEntryNotFound);
    return getMaterial(index);
  }

  /**
   * Find a material by its index within the collection. If the index is out of range, the function throws an exception.
   * @param index The material index
   * @return A reference to the requested <i>MaterialRow</i> struct
   */
  MaterialRow& MaterialTable::getMaterial(int index) { // throws std::runtime_error
    if (index >= (int)materials.size()) throw std::runtime_error("MaterialTable::getMaterialByIndex(): " + errIndexOutOfRange);
    return materials.at(index);
  }

  /**
   * Replace a material as identified by its tag, if it exists.
   * @param oldtag The name of the material that is to be replaced
   * @param newmat The <i>MaterialRow</i> struct containing the replacement
   * @return True if the operation was successful, false if the material that was to be replaced does not exist
   */
  bool MaterialTable::replaceMaterial(std::string oldtag, MaterialRow newmat) {
    int index = findIndex(oldtag);
    return replaceMaterial(index, newmat);
  }

  /**
   * Replace a material as identified by its index, if it exists.
   * @param index The material index in the internal vector container
   * @param newmat The <i>MaterialRow</i> struct containing the replacement
   * @return True if the operation was successful, false if the index was out of range
   */
  bool MaterialTable::replaceMaterial(int index, MaterialRow newmat) {
    if (index < 0 || index >= (int)materials.size()) return false;
    materials.at(index) = newmat;
    return false;
  }

  /**
   * This is a bookkeeping function that wraps around a basic function of <i>std::vector</i>.
   * @return The number of entries in the table.
   */
  unsigned int MaterialTable::rowCount() {
    return materials.size();
  }

  /**
   * This is a bookkeeping function that wraps around a basic function of <i>std::vector</i>.
   * @return True if there are no entries in the table, false otherwise
   */
  bool MaterialTable::empty() {
    return materials.empty();
  }

  /**
   * Print the contents of the material table to <i>cout</i>.
   */
  void MaterialTable::print() {
    std::cout << "The material table contains " << materials.size() << " elements:" << std::endl;
    for (unsigned int i = 0; i < materials.size(); i++) {
      std::cout << "Material: " << materials.at(i).tag << ", density: " << materials.at(i).density << ", radiation length: ";
      std::cout << materials.at(i).rlength << ", interaction length: " << materials.at(i).ilength << std::endl;
    }
    std::cout << std::endl;
  }

  /**
   * Find the internal index of a material as identified by its tag
   * @param tag The name of the material
   * @return The index of the requested material; -1 if no such material exists in the table
   */
  int MaterialTable::findIndex(std::string tag) {
    bool found = false;
    int index = 0;
    while (index < (int)materials.size() && !found) {
      if (tag.compare(materials.at(index).tag) == 0) found = true;
      else index++;
    }
    if (!found) return -1;
    return index;
  }
}
