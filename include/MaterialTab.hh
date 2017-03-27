/**
 * @file MaterialTab.h
 *
 * @date Jul 29, 2014
 * @author Stefano Martina
 */

#ifndef MATERIALTAB_H_
#define MATERIALTAB_H_

#include <map>
#include <tuple>

namespace material {
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

#endif /* MATERIALTAB_H_ */
