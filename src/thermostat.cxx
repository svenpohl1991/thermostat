#include "thermostat/properties.hpp"
#include <filesystem>
int main(int argc, char* argv[]) {
    std::string path;
    auto cwd = std::filesystem::current_path(); //getting path
    std::filesystem::current_path(cwd); //setting path
    if (argc > 1) {
        path = argv[1];
        auto i = create_statistics(path);
    }
    else
    {
        throw std::invalid_argument("No input json was given");
    }
    
}
