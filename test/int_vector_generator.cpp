#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

int main(int argc, char * argv[])
{
    if (argc < 5)
    {
        cout << "Usage: " << argv[0] << " FILE SIZE WIDTH DEFAULT_VALUE [PERM_SEED]" << endl;
        cout << " (1) Generates an int_vector<>(SIZE, DEFAULT_VALUE, WIDTH)" << endl;
        cout << "     Vector will be initialized with random bits, if " << endl;
        cout << "     DEFAULT_VALUE=r. If DEFAULT_VALUE=i, v will be set to" << endl;
        cout << "     the identity." << endl;
        cout << " (2) If PERM_SEED is specified, a random_shuffle seeded with" << endl;
        cout << "     PERM_SEED will be performed." << endl;
        cout << " (3) Stores the vector to FILE." << endl;
        cout << "Only " << argc << " elements provided at current program call." << endl;
        return 1;
    }
    uint64_t size = stoull(argv[2]);
    uint64_t width = stoull(argv[3]);
    int_vector<> v(size, 0, width);
    if ('r' == argv[4][0])
    {
        util::set_random_bits(v);
    }
    else if ('i' == argv[4][0])
    {
        util::set_to_id(v);
    }
    else
    {
        uint64_t default_value = stoull(argv[4]);
        util::set_to_value(v, default_value);
    }
    if (argc > 5)
    {
        std::random_device rd;
        std::mt19937 g(rd());
        unsigned long seed = stoul(argv[5]);
        srand(seed);
        shuffle(v.begin(), v.end(), g);
    }
    store_to_file(v, argv[1]);
}
