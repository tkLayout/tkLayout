#include <mainConfigHandler.h>

int main() {
  mainConfigHandler m;
  bool result = m.getConfiguration();
  if (result) return 0;
  else return -1;
}
