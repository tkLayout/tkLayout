#include "include/SvnRevision.h"
#include <iostream>

using namespace std;
int main() {
  cout << "Revision is __" << SvnRevision::revisionNumber << "__" << endl;
  return 0;
}
