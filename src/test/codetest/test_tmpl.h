// Simple test for template specialisation, with just one specialised function
// build with  g++ test_tmpl.cpp -o test_tmpl

#include <iostream>


class test_interface
{
public:
    virtual void print2() = 0;
};


template<typename T>
class test: public test_interface
{
public:
    void print(T in) {
        std::cout << "Test: " << in << std::endl;
    };

    void print2() {
        std::cout << "Test2 " << std::endl;
    };
};

template<>
void test<bool>::print2() {
    std::cout << "Test2 specialised " << std::endl;
};
