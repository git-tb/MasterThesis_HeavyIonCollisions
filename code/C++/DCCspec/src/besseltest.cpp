 #include <boost/math/special_functions/bessel.hpp> 
 #include <cmath>
 #include <chrono>
 #include <iostream>



 int main() {

    int Nx = 3;
    double x[Nx] = {0.0098436,2.231657233,213.7320129};

    std::cout << ">> J0\n\n";
    {
        std::cout << j0(x[0]) << std::endl;
        std::cout << jn(0,x[0]) << std::endl;
        std::cout << std::cyl_bessel_j(0,x[0]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(0,x[0]) << std::endl;
        std::cout << "\n";
        std::cout << j0(x[1]) << std::endl;
        std::cout << jn(0,x[1]) << std::endl;
        std::cout << std::cyl_bessel_j(0,x[1]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(0,x[1]) << std::endl;
        std::cout << "\n";
        std::cout << j0(x[2]) << std::endl;
        std::cout << jn(0,x[2]) << std::endl;
        std::cout << std::cyl_bessel_j(0,x[2]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(0,x[2]) << std::endl;   
        std::cout << "\n";
    }

    {
        std::cout << " > j0()\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = j0(x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > jn(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = jn(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > std::cyl_bessel_j(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = std::cyl_bessel_j(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > boost::math::cyl_bessel_j(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = boost::math::cyl_bessel_j(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    std::cout << "\n";

    std::cout << ">> J1\n\n";
    {
        std::cout << j1(x[0]) << std::endl;
        std::cout << jn(1,x[0]) << std::endl;
        std::cout << std::cyl_bessel_j(1,x[0]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(1,x[0]) << std::endl;
        std::cout << "\n";
        std::cout << j1(x[1]) << std::endl;
        std::cout << jn(1,x[1]) << std::endl;
        std::cout << std::cyl_bessel_j(1,x[1]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(1,x[1]) << std::endl;
        std::cout << "\n";
        std::cout << j1(x[2]) << std::endl;
        std::cout << jn(1,x[2]) << std::endl;
        std::cout << std::cyl_bessel_j(1,x[2]) << std::endl;
        std::cout << boost::math::cyl_bessel_j(1,x[2]) << std::endl;   
        std::cout << "\n";
    }

    {
        std::cout << " > j1()\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = j1(x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > jn(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = jn(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > std::cyl_bessel_j(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = std::cyl_bessel_j(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > boost::math::cyl_bessel_j(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = boost::math::cyl_bessel_j(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    std::cout << "\n";

    std::cout << ">> Y0\n\n";
    {
        std::cout << y0(x[0]) << std::endl;
        std::cout << yn(0,x[0]) << std::endl;
        std::cout << std::cyl_neumann(0,x[0]) << std::endl;
        std::cout << boost::math::cyl_neumann(0,x[0]) << std::endl;
        std::cout << "\n";
        std::cout << y0(x[1]) << std::endl;
        std::cout << yn(0,x[1]) << std::endl;
        std::cout << std::cyl_neumann(0,x[1]) << std::endl;
        std::cout << boost::math::cyl_neumann(0,x[1]) << std::endl;
        std::cout << "\n";
        std::cout << y0(x[2]) << std::endl;
        std::cout << yn(0,x[2]) << std::endl;
        std::cout << std::cyl_neumann(0,x[2]) << std::endl;
        std::cout << boost::math::cyl_neumann(0,x[2]) << std::endl;   
        std::cout << "\n";
    }

    {
        std::cout << " > y0()\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = y0(x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > yn(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = yn(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > std::cyl_neumann(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = std::cyl_neumann(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > boost::math::cyl_neumann(0,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = boost::math::cyl_neumann(0,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    std::cout << "\n";

    std::cout << ">> Y1\n\n";
    {
        std::cout << y1(x[0]) << std::endl;
        std::cout << yn(1,x[0]) << std::endl;
        std::cout << std::cyl_neumann(1,x[0]) << std::endl;
        std::cout << boost::math::cyl_neumann(1,x[0]) << std::endl;
        std::cout << "\n";
        std::cout << y1(x[1]) << std::endl;
        std::cout << yn(1,x[1]) << std::endl;
        std::cout << std::cyl_neumann(1,x[1]) << std::endl;
        std::cout << boost::math::cyl_neumann(1,x[1]) << std::endl;
        std::cout << "\n";
        std::cout << y1(x[2]) << std::endl;
        std::cout << yn(1,x[2]) << std::endl;
        std::cout << std::cyl_neumann(1,x[2]) << std::endl;
        std::cout << boost::math::cyl_neumann(1,x[2]) << std::endl;   
        std::cout << "\n";
    }

    {
        std::cout << " > y1()\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = y1(x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > yn(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = yn(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > std::cyl_neumann(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = std::cyl_neumann(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    {
        std::cout << " > boost::math::cyl_neumann(1,)\n";
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < (int)1e7; i++) {
            double y = boost::math::cyl_neumann(1,x[i%Nx]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double>(end-start);
        std::cout << duration << std::endl;
    }
    std::cout << "\n";

    return 0;
 }
