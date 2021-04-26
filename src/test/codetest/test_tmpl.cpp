#include "test_tmpl.h"

int main() {
    std::cout << "start" << std::endl;

    // This prints standard base class
    test<int> t;

    t.print(42);
    t.print2();

    //try template specialisation
    test<bool> t2;

    t2.print(true);
    t2.print2();

    std::cout << "done" << std::endl;

    return 0;
}
