



class phy_module_base
{
public:
    phy_module_base();

    ~phy_module_base();


    bool iniitialise_memory() = 0;
    bool initial_conditions() = 0;

    // TBD, how does it get data? friend of ESP ? grid ?
    bool loop(data) = 0; 

    // TBD
    bool store(); //? should be "get_data()"? "store_data(h5file)"?
    
    bool configure();

    bool free_memory();

    
    

}
    
