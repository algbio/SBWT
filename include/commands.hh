#pragma once

#include <vector>
#include <string>

std::vector<std::string> get_available_variants();

int build_main(int argc, char** argv);
int search_main(int argc, char** argv);
int build_from_plain_main(int argc, char** argv);
