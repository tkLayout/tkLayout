#include <iostream>
#include <string>
#include <memory>

int main(int argc, char* argv[]) {
  std::string someThing="Hello, world!";
  auto xx = std::make_unique<std::string>(someThing);
  std::cout << *xx << std::endl;
  return 0;
}
