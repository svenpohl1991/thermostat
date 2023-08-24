#include "thermostat/properties.hpp"
int main(int argc, char* argv[]) {
    std::string path;
    if (argc > 1) {
        path = argv[1];
        auto i = create_statistics(path);
    }
    else
    {
        throw std::invalid_argument("No input json was given");
    }
    
}
